#include <assert.h>
#include <getopt.h>
#include <signal.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

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
  float coef_per_ppu(path_proc_unit& ppu) {
    if(ppu._num_clusters == 0) return 0;

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500;  
    int ppus_in_full_range = 1 + total_buffer_needed / ppu._input_buffer_length;

    int objects_contributed_per_ppu 
      = ceil( (_n_obj - _n_full_range_obj)/(float)ppus_in_full_range) +
        _n_full_range_obj;

    int clusters_required_per_ppu_per_link = objects_contributed_per_ppu * 
                                                    _coef_per_object / ppu._coef_per_cluster;

    // How many replicas do we need to make
    int input_replication_factor = ceil(clusters_required_per_ppu_per_link 
                                        / (float)ppu._num_clusters);

    assert(input_replication_factor > 0);

    //if(clusters_required_per_ppu_per_link > ppu._num_clusters) return 0;

    int comp_constrainted_links_per_ppu = 
      input_replication_factor * ppu._num_clusters / clusters_required_per_ppu_per_link;

    int num_links_sharing_ppu = min(ppu._output_router._in_degree,comp_constrainted_links_per_ppu);

    float effective_coef_per_ppu = num_links_sharing_ppu * _n_obj * _coef_per_object 
      / ppus_in_full_range / input_replication_factor;

   //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
   //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
   //    comp_constrainted_links_per_ppu);
    return effective_coef_per_ppu;
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

  std::vector<int> multi_path_vec = {10,20,30,40,50,60,70,80,90,100};
  std::vector<int> obj_range_vec = {100,200,300,400,500};
  std::vector<int> percent_full_range_vec = {0,10,20,30,40,50,60,70,80,90,100};
  std::vector<int> coef_per_obj = {20,40,60,80,100};

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

  std::vector<path_proc_unit> ppu_vec(99,t);

  for(int ppu_area = 5; ppu_area < 150; ppu_area+=5) {
    float best_mem_ratio;
    float best_coef_per_mm2=0;

    for(int i = 0; i < ppu_vec.size(); ++i) {
      int mem_ratio = 1 + i;
      ppu_vec[i].set_params_by_mem_ratio(mem_ratio,ppu_area);
      //ppu_vec[i].print_area_breakdown();
  
      float total_coef_per_ppu=0;
      int total_failures=0;
      for(auto& band : band_vec) {
        int coef_per_ppu = band.coef_per_ppu(ppu_vec[i]);
        if(coef_per_ppu==0) total_failures+=1;
        total_coef_per_ppu += coef_per_ppu;
  
      }
      float avg_coef_per_ppu = total_coef_per_ppu / band_vec.size();
      float avg_coef_per_mm2 = avg_coef_per_ppu / ppu_area;

      if(avg_coef_per_mm2 > best_coef_per_mm2) {
        best_coef_per_mm2 = avg_coef_per_mm2;
        best_mem_ratio=mem_ratio;
      }

      //printf("Mem Ratio: %d, Avg Paths/PPU: %f, Paths/PPU/mm2: %f\n",
      //    mem_ratio,avg_coef_per_ppu,avg_coef_per_mm2);
    }

    printf("For %dmm^2 PPU, Mem Ratio: %f, Coef per mm2: %f\n",
        ppu_area, best_mem_ratio, best_coef_per_mm2);
  }


  scene s;
  drbe_wafer w(t);

  return 0;
}

