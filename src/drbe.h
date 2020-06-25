#include <assert.h>
#include <getopt.h>
#include <signal.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "hw.h"
#include "util.h"

class Band {
  public:
  float coef_per_ppu(path_proc_unit& ppu, drbe_wafer& wafer, float& ppus_per_link, bool verbose = false) {
    if(ppu._num_clusters == 0) return 0;

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500;  
    int ppus_in_full_range = 1 + total_buffer_needed / ppu._input_buffer_length;

    int objects_contributed_per_ppu 
      = ceil( (_n_obj - _n_full_range_obj)/(float)ppus_in_full_range) +
        _n_full_range_obj;

    // ciel: b/c coef_per_clust may be too large
    int clusters_required_per_ppu_per_link = objects_contributed_per_ppu * 
                        ceil((float)_avg_coef_per_object / (float)ppu._coef_per_cluster);  

    // How many replicas do we need to make
    int input_replication_factor = ceil(clusters_required_per_ppu_per_link 
                                        / (float)ppu._num_clusters);

    assert(input_replication_factor > 0);

    //if(clusters_required_per_ppu_per_link > ppu._num_clusters) return 0;

    int comp_constrainted_links_per_ppu = 
      input_replication_factor * ppu._num_clusters / clusters_required_per_ppu_per_link;

    int num_links_sharing_ppu = min(ppu._output_router._in_degree,comp_constrainted_links_per_ppu);

    int total_ppus_considered = ppus_in_full_range * input_replication_factor;

    float effective_coef_per_ppu = num_links_sharing_ppu * _n_obj * _avg_coef_per_object 
      / total_ppus_considered;


    // We now need to take into account Wafer-level effects for coefficient updates
    // Note the following:
    // Fast moving objects get charged the high-update rate (others, the slow). 
    // Dynamically, we only need to send updates for actual coefficients.
   
    float metadata_size = 32;  //metatdata bits per coef
    float coef_size = 32; //bits per coefficient
    float update_bits_per_coef = coef_size + metadata_size / (float)ppu._coef_per_cluster;
    //coef_bandwith in bit/cyc

    // If Either the source or destination is a fast object, then
    // we will need to update coefficients often!
    //  Source:        Fixed      Slow    Fast
    //                ------------------------
    //  Dest:  Fixed | no-coef    slow    fast
    //         Slow  | slow       slow    fast
    //         Fast  | fast       fast    fast

    int n_platform = _n_fixed + _n_obj;
    float frac_fixed = (float)(_n_fixed)/((float)n_platform);
    float frac_slow  = (float)(_n_slow) /((float)n_platform);
    float frac_fast  = (float)( _n_fast)/((float)n_platform);

    float frac_coef_slow  = frac_slow * (frac_slow + 2 * frac_fixed);
    float frac_coef_fast  = frac_fast * (1 + frac_slow + frac_fixed);


    float coef_per_ppu = num_links_sharing_ppu * _n_obj * _avg_coef_per_object
                                    / total_ppus_considered;

    float avg_coef_bw=0; //average bandwidth in bits/ns
    // First add in the slow bits contribution
    //printf("fixed, fast, slow: %f, %f, %f\n", frac_fixed * frac_fixed, frac_coef_slow, frac_coef_fast);
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_coef_slow / _low_update_period;
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_coef_fast / _high_update_period;

    // with the above coef bandwidth, and the max input bandwidth, compute the maximum
    // depth we can send coefficients without losing bandwidth
    float max_ppu_coef_depth = ppu._coef_router._in_degree * ppu._coef_router._bitwidth / 
                                avg_coef_bw;

    float fraction_ppus_active_due_to_coeff = min(max_ppu_coef_depth/wafer.half_depth(),1.0f);   

    if(verbose) {
      printf("Half Depth %d, max_ppu_coef_depth %f\n", wafer.half_depth(), max_ppu_coef_depth);
    }

    //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
    //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
    //    comp_constrainted_links_per_ppu);
    
    //return effective_coef_per_ppu;
    ppus_per_link = (float) total_ppus_considered / (float) num_links_sharing_ppu * fraction_ppus_active_due_to_coeff;

    return effective_coef_per_ppu * fraction_ppus_active_due_to_coeff;
  }

  int num_links() {
    int links_per_band = _n_tx * _n_rx;
    int total_links = links_per_band * _n_bands;
    return total_links;
  }

  int platforms() {
    return _n_slow + _n_fast + _n_fixed;
  }

  int reflectors() {
    return _n_obj;
  }

  bool could_be_harder_than(Band& other) {
    if(platforms() > other.platforms()) return true;
    if(reflectors() > other.reflectors()) return true;
    if(_n_bands < other._n_bands) return true;
    if(_avg_coef_per_object > other._avg_coef_per_object) return true;
    if(_n_fast > other._n_fast) return true;
    if(_range > other._range) return true;
    return false;
  }

  void increase_difficulty();
  float normalized_distance_to(Band& other);
  std::vector<float>& normalized_vec();



  void print_csv() {
    printf("%d, %d, %d, %d, %d, %d", platforms(), reflectors(),
        _n_bands, _avg_coef_per_object, _n_fast, _range);
  }

  void recalculate_txrx() {
      // n_tx and n_rx indicate the number of 
      // transmitters and receivers PER BAND
      _n_tx = ceil((float) platforms() / _n_bands);
      _n_rx = ceil((float) platforms() / _n_bands); 
  
      // objects in the scene that move are objects to model
      // Other objects get lost in clutter
      _n_obj = _n_slow + _n_fast;
  }

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_bands;

  //Objects by speed
  int _n_slow=0;
  int _n_fast=0;
  int _n_fixed=0;

  int _n_full_range_obj;
  int _avg_coef_per_object;
  int _range;

  float _low_update_period=1000000; //
  float _high_update_period=10000; //clock cycles?

  std::vector<float> _norm_features;
};


class ScenarioGen {

  public:

  static void gen_scenarios(std::vector<Band>& scene_vec) {
    // Scenario Generation
    for(int i_scene = 0; i_scene < 10000; i_scene++) {
      Band b;
      int platforms = rand_rng(min_platforms,max_platforms);
  
      int cutoff1 = rand_rng(0,platforms-1);
      int cutoff2 = rand_rng(0,platforms);
      if(cutoff1>cutoff2) {
        std::swap(cutoff1,cutoff2);
      }
     
      b._n_fixed = cutoff1;
      b._n_slow  = cutoff2 - cutoff1;
      b._n_fast = platforms - cutoff2;
      //int half_fast = b._n_fast/2;
      //b._n_fast -= half_fast;
      //b._n_slow += half_fast;
  
      b._n_bands = rand_rng(min_bands,max_bands);

      b.recalculate_txrx();
  
      b._avg_coef_per_object = rand_rng(min_coef_per_obj,max_coef_per_obj);
      b._n_full_range_obj = b._n_fast;
      b._range = rand_rng(min_range,max_range);
      
      b._high_update_period = rand_rng(10000 /*10us*/,100000 /*100us*/);
  
      printf("Stationary obj: %d, Slow Obj: %d, Fast Obj: %d, bands %d,  \
          coef_per_obj: %d, range: %d, update period: %f\n", 
          b._n_fixed, b._n_slow, b._n_fast, b._n_bands, 
          b._avg_coef_per_object, b._range, b._high_update_period);
  
      scene_vec.emplace_back(b);
    }
  }

  const static int min_platforms=40;
  const static int max_platforms=200;
  const static int min_bands=1;
  const static int max_bands=4;
  const static int min_coef_per_obj=20;
  const static int max_coef_per_obj=100;
  const static int min_range=50;
  const static int max_range=500;

};

