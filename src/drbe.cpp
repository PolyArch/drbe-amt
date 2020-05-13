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
  float coef_per_ppu(path_proc_unit& ppu, drbe_wafer& wafer) {
    if(ppu._num_clusters == 0) return 0;

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500;  
    int ppus_in_full_range = 1 + total_buffer_needed / ppu._input_buffer_length;

    int objects_contributed_per_ppu 
      = ceil( (_n_obj - _n_full_range_obj)/(float)ppus_in_full_range) +
        _n_full_range_obj;

    // ciel: b/c coef_per_clust may be too large
    int clusters_required_per_ppu_per_link = objects_contributed_per_ppu * 
                        ceil((float)_coef_per_object / (float)ppu._coef_per_cluster);  

    // How many replicas do we need to make
    int input_replication_factor = ceil(clusters_required_per_ppu_per_link 
                                        / (float)ppu._num_clusters);

    assert(input_replication_factor > 0);

    //if(clusters_required_per_ppu_per_link > ppu._num_clusters) return 0;

    int comp_constrainted_links_per_ppu = 
      input_replication_factor * ppu._num_clusters / clusters_required_per_ppu_per_link;

    int num_links_sharing_ppu = min(ppu._output_router._in_degree,comp_constrainted_links_per_ppu);

    int total_ppus_considered = ppus_in_full_range * input_replication_factor;

    float effective_coef_per_ppu = num_links_sharing_ppu * _n_obj * _coef_per_object 
      / total_ppus_considered;

    // We now need to take into account Wafer-level effects for coefficient updates
    // Note the following:
    // Fast moving objects get charged the high-update rate (others, the slow). 
    // Dynamically, we only need to send updates for actual coefficients.
    
    float slow_coef_per_ppu = num_links_sharing_ppu * (_n_obj- _n_full_range_obj) * _coef_per_object 
                                    / total_ppus_considered;

    float fast_coef_per_ppu = num_links_sharing_ppu * _n_full_range_obj * _coef_per_object 
                                    / total_ppus_considered;

    float metadata_size = 32;  //metatdata bits per coef
    float coef_size = 32; //bits per coefficient
    float update_bits_per_coef = coef_size + metadata_size / (float)ppu._coef_per_cluster;
    //coef_bandwith in bit/cyc
    float coef_bandwidth = (fast_coef_per_ppu  / _high_update_period + 
                            slow_coef_per_ppu / _low_update_period) * update_bits_per_coef; 


    float coef_bandwidth_fast_obj = effective_coef_per_ppu / _high_update_period * 
                                     update_bits_per_coef;


    float fraction_slow = (float)(_n_obj- _n_full_range_obj)/((float)_n_obj);
    float fraction_fast = (float)( _n_full_range_obj)/((float)_n_obj);

    // have to take weighted average of not-moving and moving coef bandwidth
    float avg_coef_bandwidth = coef_bandwidth * fraction_slow + 
                               coef_bandwidth_fast_obj * fraction_fast;

    // with the above coef bandwidth, and the max input bandwidth, compute the maximum
    // depth we can send coefficients without losing bandwidth
    float max_ppu_coef_depth = ppu._coef_router._in_degree * ppu._coef_router._bitwidth / 
                                avg_coef_bandwidth;


    float fraction_ppus_active_due_to_coeff = min(max_ppu_coef_depth/wafer.half_depth(),1.0f);   

    //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
    //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
    //    comp_constrainted_links_per_ppu);
    
    //return effective_coef_per_ppu;
    return effective_coef_per_ppu * fraction_ppus_active_due_to_coeff;
  }

  band(int n_tx, int n_rx, int n_obj, int n_full_range_obj, int range, int coef_per_object) {
    _n_tx=n_tx;
    _n_rx=n_rx;
    _n_obj=n_obj;
    _range=range;
    _coef_per_object=coef_per_object;
    _n_full_range_obj=n_full_range_obj;
  }

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_full_range_obj;
  int _coef_per_object;
  int _range;

  float _low_update_period=1000000; //
  float _high_update_period=2000; //clock cycles?
};

class scene {
public:
  std::vector<band> bands;
  int _n_obj;
  int _obj_size; 
};


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

  std::vector<int>   multi_path_vec = {10,20,30,40,50,60,70,80,90,100};
  std::vector<int>   coef_per_obj = {20,30,40,50,60,70,80,90,100};
  std::vector<int>   obj_range_vec  = {50,100,150,200,250,300,350,400,450,500};
  std::vector<int>   percent_full_range_vec = {0,10,20,30,40,50,60,70,80,90,100};

  //std::vector<int> multi_path_vec = {100};
  //std::vector<int> obj_range_vec = {100};
  //std::vector<int> percent_full_range_vec = {10};
  //std::vector<int> coef_per_obj = {40};
  std::vector<band> band_vec;

  // Many scenarios!
  for(int num_obj : multi_path_vec) {
    for(int obj_range : obj_range_vec) {
      for(int per_full_range : percent_full_range_vec) {
        for(int coef_per_obj : coef_per_obj) {
          band_vec.emplace_back(100,100,num_obj,num_obj*per_full_range/100,obj_range,coef_per_obj);
        }
      }
    }
  }




  //for(int ppu_area = 5; ppu_area < 200; ppu_area+=5) {
  int ppu_area=20;

  for(float fast_update_period = 1000; fast_update_period < 1000000; 
      fast_update_period*=1.2589254117941672104239541063958) {
  //fast_update_period=10000;

     for(auto& band : band_vec) {
       band._high_update_period=fast_update_period;
     }


    drbe_wafer w(&t,300,(float)ppu_area);

  //for(int spmm = 1; spmm < 200; ++spmm) {
    //t._sram_Mb_per_mm2=spmm;


    //for(int agg_network = 1; agg_network < 200; agg_network+=1) {
    //agg_network=11

    float best_mem_ratio;
    int best_failures=0;
    float best_coef_per_mm2=0;
    path_proc_unit* best_ppu = NULL;

    // AggNet computaitons
    for(int agg_network = 1; agg_network < 20; agg_network+=1) {
    float area_per_side=200.0*sqrt((float)ppu_area);
    int total_ins_per_side = ceil(area_per_side/32.0); //divide by bits
    int remain=total_ins_per_side - (2*2 + agg_network*2);
    if(remain < 2) break;

      for(int cpc = 10; cpc < 40; cpc+=10) {

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
  
          float total_coef_per_ppu=0;
          int total_failures=0;
          for(auto& band : band_vec) {
            int coef_per_ppu = band.coef_per_ppu(*ppu,w);
            if(coef_per_ppu==0) total_failures+=1;
            total_coef_per_ppu += coef_per_ppu;
          }
          float avg_coef_per_ppu = total_coef_per_ppu / band_vec.size();
          float avg_coef_per_mm2 = avg_coef_per_ppu / ppu_area;

          if(avg_coef_per_mm2 > best_coef_per_mm2) {
            path_proc_unit* old_best = best_ppu;
            best_coef_per_mm2 = avg_coef_per_mm2;
            best_mem_ratio=mem_ratio;
            best_ppu=ppu;
            best_failures=total_failures;
            if(old_best) delete old_best;
          } else {
            delete ppu;
          }

          //printf("Mem Ratio: %d, Avg Paths/PPU: %f, Paths/PPU/mm2: %f\n",
          //    mem_ratio,avg_coef_per_ppu,avg_coef_per_mm2);
        }
      }
    }
    printf("For %dmm^2 PPU, sram_per_mm2 %f, in buffer MB %f, %fmm^2 PPU, clust: %d, coef_per_clust %d,\
            AggNet: %d, Mem Ratio: %f, Coef per mm2: %f, fails: %d, fast_period %f\n",
        ppu_area, t._sram_Mb_per_mm2,
        best_ppu->area(), 
        best_ppu->input_buf_area()*t._sram_Mb_per_mm2/8,
        best_ppu->_num_clusters, best_ppu->_coef_per_cluster,
        best_ppu->_output_router._in_degree, 
        best_mem_ratio, best_coef_per_mm2,
        best_failures, fast_update_period);
    //best_ppu->print_params();
    //best_ppu->print_area_breakdown();

  }


  scene s;

  return 0;
}

