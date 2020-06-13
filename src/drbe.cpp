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

using namespace std;

//using sec = chrono::seconds;
//using get_time = chrono::steady_clock;
//
//static struct option long_options[] = {
//    {"algorithm",      required_argument, nullptr, 'a',},
//    {"verbose",        no_argument,       nullptr, 'v',},
//    {0, 0, 0, 0,},
//};


class band {
public:
  float coef_per_ppu(path_proc_unit& ppu, drbe_wafer& wafer, float& ppus_per_link) {
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

    ppus_per_link = (float) total_ppus_considered / (float) num_links_sharing_ppu;

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

    float frac_fixed = (float)(_n_fixed)/((float)_n_obj);
    float frac_slow  = (float)(_n_slow)/((float)_n_obj);
    float frac_fast  = (float)( _n_fast)/((float)_n_obj);

    float frac_coef_slow  = frac_slow * (frac_slow + 2 * frac_fixed);
    float frac_coef_fast  = frac_fast * (1 + frac_slow + frac_fixed);


    float coef_per_ppu = num_links_sharing_ppu * _n_obj * _avg_coef_per_object
                                    / total_ppus_considered;

    float avg_coef_bw=0; //average bandwidth in bits/ns
    // First add in the slow bits contribution
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_coef_slow / _low_update_period;
    avg_coef_bw += coef_per_ppu * update_bits_per_coef * frac_coef_fast / _high_update_period;

    // with the above coef bandwidth, and the max input bandwidth, compute the maximum
    // depth we can send coefficients without losing bandwidth
    float max_ppu_coef_depth = ppu._coef_router._in_degree * ppu._coef_router._bitwidth / 
                                avg_coef_bw;

    float fraction_ppus_active_due_to_coeff = min(max_ppu_coef_depth/wafer.half_depth(),1.0f);   

    //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
    //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
    //    comp_constrainted_links_per_ppu);
    
    //return effective_coef_per_ppu;
    return effective_coef_per_ppu * fraction_ppus_active_due_to_coeff;
  }

  int num_links() {
    int links_per_band = _n_tx * _n_rx;
    int total_links = links_per_band * _n_bands;
    return total_links;
  }

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_bands;

  //Objects by speed
  int _n_slow;
  int _n_fast;
  int _n_fixed;

  int _n_full_range_obj;
  int _avg_coef_per_object;
  int _range;

  float _low_update_period=1000000; //
  float _high_update_period=10000; //clock cycles?
};



int get_best_ppu() {

  return 0;
}

struct ppu_stats {
  float avg_coef_per_ppu=0;
  float avg_coef_per_mm2=0;
  int total_failures=0;
  float avg_wafers=0;
  float avg_ppus_per_link=0;
  int wafer_histo[100000];
  int coef_per_mm2[100000]; 

  void print_wafer_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",wafer_histo[i]);
    }
  }

};

void evaluate_ppu(path_proc_unit* ppu, std::vector<band>& scene_vec, drbe_wafer& w, 
                   ppu_stats& stats) {
  stats.total_failures=0;

  float total_coef_per_ppu=0;
  float wafers=0, total_ppus_per_link=0;

  for(auto& band : scene_vec) {
    float ppus_per_link=0;
    int coef_per_ppu = band.coef_per_ppu(*ppu,w,ppus_per_link);
    if(coef_per_ppu==0) stats.total_failures+=1;
    total_ppus_per_link += ppus_per_link;

    total_coef_per_ppu += coef_per_ppu;

    float ppus = band.num_links() * ppus_per_link;
    float num_wafers = ppus / w.num_units();
    int int_num_wafers = ceil(num_wafers);
    wafers += int_num_wafers;
    assert(int_num_wafers<100000);
    stats.wafer_histo[int_num_wafers]++;
    int int_avg_coef_per_mm2 = stats.avg_coef_per_ppu / w.chiplet_area();
    assert(int_avg_coef_per_mm2 < 100000);
    stats.coef_per_mm2[int_num_wafers]++;
    //printf("Avg wafer: %d, PPUs: %f, band.num_links():%d, ppus_per_link:%f\n",
    //    int_num_wafers, ppus, band.num_links(), ppus_per_link);
  }
  stats.avg_ppus_per_link = total_ppus_per_link / scene_vec.size();
  stats.avg_wafers = wafers / scene_vec.size();
  stats.avg_coef_per_ppu = total_coef_per_ppu / scene_vec.size();
  stats.avg_coef_per_mm2 = stats.avg_coef_per_ppu / w.chiplet_area();
}

path_proc_unit* design_ppu_for_scenarios(std::vector<band>& scene_vec, drbe_wafer& w) {
  ppu_stats best_stats;
  path_proc_unit* best_ppu = NULL; //needs to be null so we don't delete a non-pointer
  float ppu_area = w.chiplet_area();
  tech_params& t = *w._t;

  // AggNet computaitons
  for(int agg_network = 1; agg_network < 20; agg_network+=1) {
  float area_per_side=200.0*sqrt((float)ppu_area);
  int total_ins_per_side = ceil(area_per_side/32.0); //divide by bits
  int remain=total_ins_per_side - (2*2 + agg_network*2);
  if(remain < 2) break;

    // Loop over coefficients per cluster
    for(int cpc = 20; cpc < 21; cpc+=1) {

      for(int mem_ratio = 1; mem_ratio < 100; ++mem_ratio) {
        path_proc_unit* ppu =new path_proc_unit(&t);

        ppu->_coef_per_cluster=cpc;
        ppu->_input_router._in_degree=2;
        ppu->_input_router._out_degree=2;
        ppu->_output_router._in_degree=agg_network;
        ppu->_output_router._out_degree=agg_network;
        ppu->_coef_router._in_degree=remain/2;
        ppu->_coef_router._out_degree=remain/2;

        ppu->set_params_by_mem_ratio(mem_ratio,ppu_area);
        if(ppu->_num_clusters==0) break;
        //ppu_vec[i].print_area_breakdown();
 
        ppu_stats stats; 
        evaluate_ppu(ppu,scene_vec,w,stats);

        if((stats.avg_coef_per_mm2 > best_stats.avg_coef_per_mm2) && stats.total_failures==0) {
          path_proc_unit* old_best = best_ppu;
          best_stats = stats;
          best_ppu=ppu;
          if(old_best) delete old_best;
        } else {
          delete ppu;
        }

        //printf("Mem Ratio: %d, Clusters/PPU: %d, Coef/PPU/mm2: %f\n",
        //    mem_ratio,ppu->_num_clusters,avg_coef_per_mm2);
      }
    }
  }

  return best_ppu;
}



int main(int argc, char* argv[]) {
  bool verbose = false;
  string str_schedType = string("sa");

  tech_params t;

  /*
  int opt;
  while ((opt = getopt_long(argc, argv, "va:", long_options, nullptr)) != -1) {
    switch (opt) {
      case 'a': str_schedType = string(optarg); break;
      case 'v': verbose = true; break;
      default: exit(1);
    }
  }

  argc -= optind;
  argv += optind;

  if (argc != 2) {
    //cerr << "Usage: ss_sched [FLAGS] test.drbe_hw test.drbe_scene \n";
    //exit(1);
    cout << "Assuming Default Params\n";
  }
  */

  std::vector<band> scene_vec;

  // Scenario Generation
  for(int i_scene = 0; i_scene < 10000; i_scene++) {
    band b;
    int platforms = rand_rng(80,200);

    int cutoff1 = rand_rng(0,platforms-1);
    int cutoff2 = rand_rng(0,platforms);
    if(cutoff1>cutoff2) {
      std::swap(cutoff1,cutoff2);
    }
   
    b._n_fixed = cutoff1;
    b._n_slow  = cutoff2 - cutoff1;
    b._n_fast = platforms - cutoff2;

    b._n_bands = rand_rng(1,4);

    // n_tx and n_rx indicate the number of 
    // transmitters and receivers PER BAND
    b._n_tx = ceil((float) platforms / b._n_bands);
    b._n_rx = ceil((float) platforms / b._n_bands); 

    // objects in the scene that move are objects to model
    // Other objects get lost in clutter
    b._n_obj = b._n_slow + b._n_fast;

    b._avg_coef_per_object = rand_rng(20,100);
    b._n_full_range_obj = b._n_fast;
    b._range = rand_rng(50,500);
    b._high_update_period = rand_rng(10000 /*10us*/,100000 /*100us*/);

    printf("Stationary obj: %d, Slow Obj: %d, Fast Obj: %d, bands %d,  \
        coef_per_obj: %d, range: %d, update period: %f\n", 
        b._n_fixed, b._n_slow, b._n_fast, b._n_bands, 
        b._avg_coef_per_object, b._range, b._high_update_period);

    scene_vec.emplace_back(b);
  }



  
  printf("\n  --- DSE Summary --- \n");

  for(int ppu_area = 20; ppu_area < 21; ppu_area+=5) {
  //int ppu_area=20;

  //for(float fast_update_period = 1000; fast_update_period < 1000000; 
  //    fast_update_period*=1.2589254117941672104239541063958) {
  //fast_update_period=10000;

     //for(auto& band : scene_vec) {
     //  band._high_update_period=fast_update_period;
     //}


    drbe_wafer w(&t,300,(float)ppu_area);

  //for(int spmm = 1; spmm < 200; ++spmm) {
    //t._sram_Mb_per_mm2=spmm;


    //for(int agg_network = 1; agg_network < 200; agg_network+=1) {
    //agg_network=11

    path_proc_unit* best_ppu = design_ppu_for_scenarios(scene_vec,w);

    ppu_stats stats;
    evaluate_ppu(best_ppu,scene_vec,w,stats);

    
    printf("%dmm^2 PPU (%0.2f), in-MB: %0.2f, clust: %d, coef_per_clust %d, "\
           "Agg nets: %d, Mem Ratio: %d, Coef per mm2: %0.2f, "\
           "avg_wafers: %0.2f, avg_ppus_per_link: %0.2f, fails: %d\n",
        ppu_area, 
        best_ppu->area(), 
        best_ppu->input_buf_area()*t._sram_Mb_per_mm2/8,
        best_ppu->_num_clusters, best_ppu->_coef_per_cluster,
        best_ppu->_output_router._in_degree, 
        (int)best_ppu->_mem_ratio, 
        stats.avg_coef_per_mm2,
        stats.avg_wafers,
        stats.avg_ppus_per_link,
        stats.total_failures);
    printf("Wafers Needed Histogram: ");
    stats.print_wafer_histo();
    printf("\n");
    //best_ppu->print_params();
    //best_ppu->print_area_breakdown();
  }


  return 0;
}

