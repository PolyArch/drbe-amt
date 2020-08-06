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

struct ge_stats {
  // Computation Area breakdown
  float avg_relative_location_comp_area = 0.0;
  float avg_affine_transform_comp_area = 0.0;
  float avg_relative_motion_comp_area = 0.0;
  float avg_path_gain_comp_area = 0.0;
  float avg_path_delay_comp_area = 0.0;
  float avg_ge_comp_area = 0.0;
  // Total Computation Area breakdown
  float total_relative_location_comp_area = 0.0;
  float total_affine_transform_comp_area = 0.0;
  float total_relative_motion_comp_area = 0.0;
  float total_path_gain_comp_area = 0.0;
  float total_path_delay_comp_area = 0.0;
  float total_ge_comp_area = 0.0;
  // Count
  int num_ge_core = -1;
  // Total Memory Capacity
  float total_mem_bytes = 0.0;
  float avg_mem_bytes = 0.0;
  // histogram
  int chiplet_histo[100000];
  int ge_core_histo[100000];
  void print_chiplet_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",chiplet_histo[i]);
    }
  }
  void print_ge_core_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",ge_core_histo[i]);
    }
  }
};


//Note for later: Sumeet says that if either tx or rx is fast, then we need variable doppler
//yikes

class Band {
  public:

  float ppus_per_band(path_proc_unit& ppu, drbe_wafer& w, int& failures, bool verbose=false) {
    // Table for reasoning about number of slow/fast tx/rx pairs
    //          tx:   Fixed       Slow       Fast
    //                -----------------------------
    //  rx:    Fixed | fixed      slow      fast
    //         Slow  | slow       slow      fast
    //         Fast  | fast       fast      fast

    int n_platform = _n_fixed + _n_obj;
    float frac_fixed = (float)(_n_fixed)/((float)n_platform);
    float frac_slow  = (float)(_n_slow) /((float)n_platform);
    //float frac_fast  = (float)(_n_fast)/((float)n_platform);

    float frac_txrx_fixed_or_slow  = (frac_slow+frac_fixed) * (frac_slow+frac_fixed);
    float frac_txrx_fast = 1 - frac_txrx_fixed_or_slow;


    float total_ppus=0;

    //iterate over combinations of fast and slow ppus
    for(int fast = 0; fast < 2; ++fast) {
      for(int clutter = 0; clutter < 2; ++clutter) {
        float links_per_ppu = calc_links_per_ppu(ppu,w,fast,clutter,verbose);
        if(links_per_ppu==0) {
          failures+=1;
          continue;
        }
        float frac_speed = fast==0 ? frac_txrx_fixed_or_slow : frac_txrx_fast;
        float frac_clutter = clutter==0 ? (1-_frac_clutter) : _frac_clutter;
        float frac= frac_speed * frac_clutter;

        float ppus = ceil(num_links() * frac / links_per_ppu);

        //printf("%f %f %f %f; ", frac_speed, frac_clutter,frac,ppus);
        total_ppus +=ppus;
      }
    }
    //printf(" PPUs\n");
    return total_ppus;
  }
              


  float calc_links_per_ppu(path_proc_unit& ppu, drbe_wafer& wafer, 
                     bool fast_txrx, bool clutter,
                     bool verbose = false) {
    if(ppu._num_clusters == 0) {
      return 0;
    }

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500.0;  

    if(_is_direct_path || _is_scatter_sep) {
      total_buffer_needed /=2;  //Only need half the range for direct path
    }

    if(verbose) {
      printf("input buffer length: %d\n",ppu._input_buffer_length);
      printf("total_buffer_needed: %d\n",total_buffer_needed);
    }

    int ppus_in_full_range = 1 + total_buffer_needed / ppu._input_buffer_length;

    if(verbose) {
      printf("ppus_in_full_range: %d, mem_ratio %f\n",ppus_in_full_range,ppu._mem_ratio);
    }

    int objects_contributed_per_ppu;
    if(_is_direct_path) {
      // also need to expand ppus_in_full_range due to i/o limitations
      // we're going to get _n_obj inputs
      int io_ppus_per_side = ceil((float)_n_obj / (float)ppu._output_router._in_degree);
      ppus_in_full_range = max(ppus_in_full_range, io_ppus_per_side * io_ppus_per_side);

      objects_contributed_per_ppu = ceil(_n_obj/(float)ppus_in_full_range);
    } else if(_is_scatter_sep) {
      objects_contributed_per_ppu = 2;  // Funny right?
    } else { // Tapped delay case
      if(ppu._is_dyn_reconfig) {
       objects_contributed_per_ppu = 
        ceil( (_n_obj)/(float)ppus_in_full_range);
      } else {
       objects_contributed_per_ppu = 
        ceil( (_n_obj - _n_full_range_obj)/(float)ppus_in_full_range) +
        _n_full_range_obj;

      }
    } 


    // ciel: b/c coef_per_clust may be too large
    int clusters_required_per_ppu_per_link = objects_contributed_per_ppu * 
                        ceil((float)_avg_coef_per_object / (float)ppu._coef_per_cluster);  

    // How many replicas do we need to make
    int input_replication_factor = ceil(clusters_required_per_ppu_per_link 
                                        / (float)ppu._num_clusters);

    if(verbose) {
      printf("obj/ppu/link %d, clusters/ppu/link : %d\n, input repl: %d\n",
          objects_contributed_per_ppu, clusters_required_per_ppu_per_link, input_replication_factor);
    }


    assert(input_replication_factor > 0);

    //if(clusters_required_per_ppu_per_link > ppu._num_clusters) return 0;

    int comp_constrainted_links_per_ppu = 
      input_replication_factor * ppu._num_clusters / clusters_required_per_ppu_per_link;

    //TODO: for now i don't think its worth it to consider multiple links sharing the same
    //PPU for a direct-path model, b/c the I/O will almost always be the constraint.
    //similar reasoning for scatter sep, except the limitation is just memory
    int max_sharing_factor = (_is_direct_path || _is_scatter_sep) ? 1 : ppu._output_router._in_degree;

    int num_links_sharing_ppu = min(max_sharing_factor,comp_constrainted_links_per_ppu);

    int total_ppus_considered = ppus_in_full_range * input_replication_factor;

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
    //  tx/rx:        Fixed/Slow  Fast
    //                -----------------------------
    //         Slow  | slow       fast
    //         Fast  | fast       fast


    float frac_slow, frac_fast;
    if(fast_txrx) {
      frac_slow = 0;
      frac_fast = 1;
    } else {
      frac_slow  = (float)(_n_slow) / ((float)_n_obj);
      frac_fast  = (float)(_n_fast) / ((float)_n_obj);
    }

    float coef_per_ppu = num_links_sharing_ppu * objects_contributed_per_ppu * _avg_coef_per_object
                                    / total_ppus_considered;

    float avg_coef_bw=0; //average bandwidth in bits/ns
    // First add in the slow bits contribution
    //printf("fixed, fast, slow: %f, %f, %f\n", frac_fixed * frac_fixed, frac_coef_slow, frac_coef_fast);
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_slow / _low_update_period;
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_fast / _high_update_period;

    // with the above coef bandwidth, and the max input bandwidth, compute the maximum
    // depth we can send coefficients without losing bandwidth
    float max_ppu_coef_depth = ppu._coef_router._in_degree * ppu._coef_router._bitwidth / 
                                avg_coef_bw;

    float fraction_ppus_active_due_to_coeff = min(max_ppu_coef_depth/wafer.half_depth(),1.0f);   

    if(verbose) {
      printf("Half Depth %d, max_ppu_coef_depth %f, frac_active: %f\n", 
          wafer.half_depth(), max_ppu_coef_depth, fraction_ppus_active_due_to_coeff);
    }

    //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
    //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
    //    comp_constrainted_links_per_ppu);
    
    float ppus_per_link = (float) total_ppus_considered / (float) num_links_sharing_ppu / fraction_ppus_active_due_to_coeff;

    if(verbose) {
      printf("num_links_sharing_ppu %d, ppus per link %f\n\n\n",num_links_sharing_ppu,ppus_per_link);
    }

    return 1.0f/ppus_per_link;
  }

  int algorithmic_num_links() {
    int links_per_band = _n_tx * _n_rx;
    int total_links = links_per_band * _n_bands;
    return total_links;
  }


  int num_links() {
    int links_per_band;
    if(_is_direct_path || _is_scatter_sep) {
      // No need to include self links for direct path
      links_per_band = (_n_tx-1) * (_n_rx-1);
    } else {
      links_per_band = _n_tx * _n_rx;
    }

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

  bool increase_difficulty(int kind);
  float normalized_distance_to(Band& other);
  std::vector<float>& normalized_vec();
  float ge_mem_byte(ge_core & ge, ge_stats & stats);
  float ge_comp_area(ge_core & ge, ge_stats & stats);

  void print_csv() {
    printf("%d, %d, %d, %d, %d, %d", platforms(), reflectors(),
        _n_bands, _avg_coef_per_object, _n_fast, _range);
  }

  void print_norm_csv();

  void recalculate_txrx() {
      // n_tx and n_rx indicate the number of 
      // transmitters and receivers PER BAND
      _n_tx = ceil((float) platforms() / _n_bands);
      _n_rx = ceil((float) platforms() / _n_bands); 
  
      // objects in the scene that move are objects to model
      // Other objects get lost in clutter
      _n_obj = _n_slow + _n_fast;
      _n_full_range_obj = _n_fast;
  }

  bool _is_direct_path=false;
  bool _is_scatter_sep=false;

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_bands;

  //Objects by speed
  int _n_platform=0;
  int _n_slow=0;
  int _n_fast=0;
  int _n_fixed=0; 

  int _n_full_range_obj;
  int _avg_coef_per_object;
  int _range;

  float _low_update_period=1000000; //
  float _high_update_period=10000; //clock cycles?
  // Compute Area needed
  int num_chiplet = 0;
  int num_ppu_chiplet = 0;
  int num_ge_chiplet = 0;
  float ge_area = 0.0;
  float _frac_clutter = 0.0;

  // Geometry Engine Resource 
  float antenna_0th_comp = 0.0; // in 
  float antenna_4th_comp = 0.0;
  float antenna_5th_comp = 0.0;

  float antenna_0th_mem = 0.0;
  float antenna_4th_mem = 0.0;
  float antenna_5th_mem = 0.0;

  // Memory
  float mem_byte = 0.0; //in mm2

  std::vector<float> _norm_features;
};


class ScenarioGen {

  public:

  static void gen_scenarios(std::vector<Band>& scene_vec,bool direct_path, bool scatter_sep) {
    // Lets always add the most difficult scenario -- comment this if you don't want it : )
    //Band b;
    //b._n_fixed = max_platforms -  max_reflectors;
    //b._n_slow = 0;
    //b._n_fast = max_reflectors;
    //b._n_bands = 1;
    //b.recalculate_txrx();
    //b._avg_coef_per_object = 60;
    //b._n_full_range_obj = b._n_fast;
    //b._range=max_range;
    //b._high_update_period = 10000; //10us
    //b._is_direct_path=direct_path;
    //b._is_scatter_sep=scatter_sep;
    //scene_vec.emplace_back(b); 

    // Scenario Generation
    for(int i_scene = 0; i_scene < 10000; i_scene++) {
      Band b;
      int platforms = rand_rng(min_platforms,max_platforms);
      int max_rand_reflect = min(max_reflectors,max_platforms);
      int num_reflectors = rand_rng(1,max_rand_reflect);

      int slow = rand_rng(0,num_reflectors);

      b._is_direct_path=direct_path;
      b._is_scatter_sep=scatter_sep;

      b._n_fixed = platforms-num_reflectors;
      b._n_slow  = slow;
      b._n_fast = num_reflectors-slow;
      //int half_fast = b._n_fast/2;
      //b._n_fast -= half_fast;
      //b._n_slow += half_fast;
  
      b._n_bands = rand_rng(min_bands,max_bands);

      b.recalculate_txrx();
  
      b._avg_coef_per_object = (rand_rng(min_coef_per_obj,max_coef_per_obj) / 4) * 4;
      b._n_full_range_obj = b._n_fast;
      b._range = (rand_rng(min_range,max_range)/10)*10;
      
      b._high_update_period = rand_rng(10000 /*10us*/,100000 /*100us*/);

      for(int i = 0; i < 100; ++i) {
        //b.increase_difficulty(1);
        //b.increase_difficulty(2); //upgrade slow to fast
      }

  /*
      printf("Stationary obj: %d, Slow Obj: %d, Fast Obj: %d, bands %d,  \
          coef_per_obj: %d, range: %d, update period: %f\n", 
          b._n_fixed, b._n_slow, b._n_fast, b._n_bands, 
          b._avg_coef_per_object, b._range, b._high_update_period);
  */
      scene_vec.emplace_back(b);
    }
  }

  const static int clutter_links=0;
  const static int min_platforms=80;
  const static int max_platforms=200;
  const static int max_reflectors=100;
  const static int min_bands=1;
  const static int max_bands=4;
  const static int min_coef_per_obj=20;
  const static int max_coef_per_obj=100;
  const static int min_range=50;
  const static int max_range=500;

};

