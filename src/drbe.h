#include <assert.h>
#include <getopt.h>
#include <signal.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include "dsp.h"
#include "hw.h"
#include "util.h"
#include <sstream>
#include <tuple>

#define LINK_COMPLEXITY_OVERPROV (1.25f)

struct WaferStats {
  // Wafer level
  int num_wafer = 0; // number of wafers we can provide (this is the input)
  int wafer_io_limit = 0; // the target wafer io
  float tech_scaling = 0.0;
  // chiplet level
  float chiplet_io_layer = 0.0;

  // Compute Area needed
  int num_ppu_chiplet = 0;
  int num_ge_chiplet = 0;

  // Geometry Engine Breakdown
  /*Requirement: num_ge_chiplet == num_gc_chiplet + num_gm_chiplet*/
  int num_gc_chiplet = 0;
  int num_gm_chiplet = 0;
};

// compute, memory, bandwidth, latency
struct cmbl{
  float compute = 0.0; // unit: MAC -> compute * 1e-9 to get Gig MACs
  float memory = 0.0; // unit: 64-bit -> memory * (64*1e-6/8) to get Mega Bytes
  float bandwidth = 0.0; // unit: 64-bit -> bandwidth * ( (64*1e-9) / 8 ) to get Gig Byte per second
  float latency = 0.0; // clock cycles
  void set_cmbl(float c, float m, float b, float l){
    compute = c;
    memory = m;
    bandwidth = b;
    latency = l;
  }
};

// Global Fidelity
struct st_global_fid {
  float num_path = 0.0;
  float num_obj = 0.0;
  float upd_rate = 0.0;
};

// Coordinate transformation
// Geodetic to Cartesian Tradeoff Space Analysis 
struct st_coordinate_trans : cmbl {
  float ta1_scene_upd_rate = 0.0;
};

// NR engine
struct st_nr_engine : cmbl {
  float interpolation_ord = 0.0;
  float conv_fed = 0.0;
};

// Relative Orientation
struct st_relative_orientation : cmbl{};

// Antenna Gain
struct st_antenna_gain : cmbl{
  float order = 0.0;
  float num_antenna = 0.0;
  float res_angle = 0.0;
  float dict_dim = 0.0;
  int general = 0;
};

// Path Gain and Velocity
struct st_path_velocity : cmbl {};

// RCS
struct st_rcs : cmbl{
  float order = 0.0;
  float points = 0.0;
  float angle = 0.0;
  float freq = 0.0;
  float plzn = 0.0;
  float pntAngle = 0.0;
  float plates = 0.0;
};

struct st_tu : cmbl {};

struct ge_stat_per_band {
  st_global_fid global_fid;
  st_coordinate_trans coordinate_trans;
  st_nr_engine nr_engine;
  st_relative_orientation relative_orientation;
  st_antenna_gain antenna_gain;
  st_path_velocity path_velocity;
  st_rcs rcs;
  st_tu tu;
  // Total Compute and Total Memory
  float total_compute = 0.0;
  float total_memory = 0.0;
  // Area
  float compute_area = 0.0;
  float mem_area = 0.0;
  float total_area = 0.0;
};

// This structure records the number of "MACCS" in
struct DatapathOverheads {
  typedef enum  {
    MAX=0, //actual computaiton
    FAST_OBJ,
    POINT_FULL_RATIO,
    CLUTTER,
    IMPERFECT_MAP,
    COEF_NETWORK,
    NUM_ENTRIES 
  } OVERHEAD_REASON;

  static const char* nameof(int i) {
    switch(i) {
      case MAX: return "max";
      case FAST_OBJ: return "fast";
      case POINT_FULL_RATIO: return "p/f";
      case CLUTTER: return "clut";
      case IMPERFECT_MAP: return "map";
      case COEF_NETWORK: return "net";
      default: return "util";
    }
  }

  void scan() {
    for(int i = 0; i < NUM_ENTRIES-1; ++i) {
      reason[i] = reason[i]-reason[i+1];
    }
    
  }

  void add_ov(DatapathOverheads& other, float ratio) {
    for(int i = 0; i < NUM_ENTRIES; ++i) {
      reason[i] += ratio * other.reason[i];
    }
  }

  float reason[NUM_ENTRIES] = {0};
};

struct ppu_stat_per_band {
  float num_ppu_chiplet = 0.0;
};


struct PPUStats {
  float avg_coef_per_ppu=0;
  float avg_coef_per_mm2=0;
  int total_failures=0;
  float avg_wafers=0;
  float avg_links_per_wafer=0;
  float percent_in_target_wafer=0;

  float avg_links_per_mm2=0;

  // Histogram
  int wafer_histo[100000];
  void print_wafer_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",wafer_histo[i]);
    }
  }

  float num_scenarios=0;
  DatapathOverheads overhead;

  // Per Band Statistic
  ppu_stat_per_band * ppu_stat_vec;
};



//Note for later: Sumeet says that if either tx or rx is fast, then we need variable doppler
//yikes

class Band {
  public:

  static float average_clutter(std::vector<Band>);

  //int _tdl_max_systems_per_wafer=0;
  //int _aidp_max_systems_per_wafer=0;
  //int _dp_max_systems_per_wafer=0;

  //void set_dp_max_systems_per_wafer(int wafer_in) {
  //  _tdl_max_systems_per_wafer = min(wafer_in,wafer_out); 
  //}

  //void set_dp_max_systems_per_wafer(int wafer_in) {
  //  _dp_max_systems_per_wafer = wafer
  //}
  //void set_dp_max_systems_per_wafer() {
  //  
  //}

  //This only thinks about inputs, since outputs are similar
  //float io_for_n_systems(float systems_per_wafer) { 
  //  if(_is_direct_path) {
  //    return (_n_obj-1) * systems_per_wafer;
  //  } else if (_is_aidp) {
  //    return systems_per_wafer * _k_rcs_points;
  //  } else {
  //    return systems_per_wafer;
  //  }
  //}

  //This only thinks about inputs, since outputs are similar
  float max_systems_per_wafer(drbe_wafer& w) { 
    int max_by_in = w.input_sig_io();

    int max_by_agg;
    if(_is_direct_path) {
      max_by_agg = w.agg_input_io() / (float)(max(_n_tx-1,1));                
    } else if (_is_aidp) {
      //multiply by 1.25 just to be conservative, got help us if we get this wrong
      max_by_agg = w.agg_input_io() / (1 + (_k_rcs_points-1)* _avg_frac_full_objects * 1.25);
    } else {
      max_by_agg = w.agg_input_io() / 2; //because we don't know how to get top and bottom
    }
    //printf("max_by_in %d, max_by_agg %d\n",max_by_in,max_by_agg);
    return std::min(max_by_in,max_by_agg);
  }


 
  float calc_links_per_wafer(float wafer_unconstrained_ppus, path_proc_unit* ppu, drbe_wafer& w, WaferStats & w_stats, bool verbose = false) {

    if(w_stats.num_wafer==0) w_stats.num_wafer=1; //FIXME not sure why this is zero, --tony

    //We're goint to just cut the biggest square we can out, and assume that's what we'll
    //use all over.
    float GE_chiplets_per_wafer = w_stats.num_ge_chiplet / w_stats.num_wafer;

    float num_wafers = wafer_unconstrained_ppus / 
      (w.num_units() * ppu->_ppus_per_chiplet - GE_chiplets_per_wafer);
      // We substract the number of chiplets that was used to Geometry Engine (Compute + Memory) 
    
    float links_per_wafer = num_links() / num_wafers;

    if(w.limit_wafer_io()==false) {
      return links_per_wafer;
    }

    float systems_per_wafer = sqrt(links_per_wafer); 

    systems_per_wafer = min(systems_per_wafer,max_systems_per_wafer(w));
    
    links_per_wafer = systems_per_wafer * systems_per_wafer;

    return max(links_per_wafer,1.0f);
  }

  float ppus_per_band(path_proc_unit& ppu, drbe_wafer& w, PPUStats& ppu_stats, bool verbose=false) {
    // Table for reasoning about number of slow/fast tx/rx pairs
    //          tx:   Fixed       Slow       Fast
    //                -----------------------------
    //  rx:    Fixed | fixed      slow      fast
    //         Slow  | slow       slow      fast
    //         Fast  | fast       fast      fast

    ppu_stats.num_scenarios++;

    float frac_fixed = (float)(_n_fixed)/((float)platforms());
    float frac_slow  = (float)(_n_slow) /((float)platforms());
    //float frac_fast  = (float)(_n_fast)/((float)n_platform);

    float frac_txrx_fixed=frac_fixed * frac_fixed;
    float frac_txrx_fixed_or_slow  = (frac_slow+frac_fixed)*(frac_slow+frac_fixed) - frac_txrx_fixed;
    float frac_txrx_fast = max(0.0f,1 - frac_txrx_fixed_or_slow - frac_txrx_fixed);
    float total = frac_txrx_fixed + frac_txrx_fixed_or_slow + frac_txrx_fast;

    assert(total > 0.99f && total < 1.01f);

    float frac_speed_vec[3] ={frac_txrx_fixed,frac_txrx_fixed_or_slow,frac_txrx_fast};
    float total_ppus=0;

    float ppus_affected_by_var_doppler_clutter=0;
    bool print_clutter_impact_statement=false;

    float native_frac_clutter = _frac_clutter;
    if(_is_direct_path || _is_aidp) {
      native_frac_clutter=0; //clutter handled by scond loop
    }

    total = 0;

    //iterate over combinations of fast and slow ppus
    for(int pulsed = 0; pulsed <2; ++ pulsed) {
      for(int speed = 0; speed < 3; ++speed) {
        for(int clutter = 0; clutter < 2; ++clutter) {
          float frac_pulsed = pulsed == 0 ? (1- _frac_pulsed) : _frac_pulsed;
          float frac_speed = frac_speed_vec[speed];
          float frac_clutter = clutter==0 ? (1-native_frac_clutter) : native_frac_clutter;
          float frac= frac_speed * frac_clutter * frac_pulsed;
  
          total+=frac;

          assert(frac_speed >=0 && frac_clutter >=0 && frac >=0);
          assert(frac_speed <=1 && frac_clutter <=1 && frac <=1);
          if(frac==0) continue; //save some time : )
 
          // For pulsed==0, we want full duty cycle (no pulsed tx)
          // For pulsed==1, we want pulsed duty cycle
          float duty_cycle = (pulsed==0) ? 1.0 : _pulsed_duty_cycle;

          if(verbose) {
            print_clutter_impact_statement=true;
            printf("MAPPING WITH PARAMS clutter: %d, speed %d, frac %f\n",clutter,speed, frac);
          }
  
          DatapathOverheads ov;
          float links_per_ppu = calc_links_per_ppu(ppu,ov,w,_is_direct_path, _is_aidp,                                                                   speed,clutter,false /*just clutter*/,
                                              duty_cycle, verbose);
          
          ppu_stats.overhead.add_ov(ov,frac); //accumulate total overhead statistics



          if(links_per_ppu==0) {
            ppu_stats.total_failures+=1;
            continue;
          }
  
          if(verbose) printf("\n");
  
          float ppus = ceil(num_links() * frac / links_per_ppu);


          //Need to pay for FULL tapped delay lines in direct-path models
          if(_is_direct_path || _is_aidp) {

            float clutter_links_per_ppu_for_dp = 
              calc_links_per_ppu(ppu,ov,w,false,false/*tap-delay*/,
                                 speed,true/*clutter*/,true/*just clutter*/,
                                 duty_cycle, verbose);
  
            ppu_stats.overhead.add_ov(ov,frac); //accumulate total overhead statistics

            float clutter_ppus_for_dp = ceil(num_links() * frac / clutter_links_per_ppu_for_dp);
            ppus += clutter_ppus_for_dp;
          }

          if(clutter && speed !=0) {
            ppus_affected_by_var_doppler_clutter+=ppus;
          }
          //printf("%f %f %f %f; ", frac_speed, frac_clutter,frac,ppus);
          total_ppus +=ppus;
        }
      }
    }

    ppu_stats.overhead.scan(); // 

    assert(total > 0.99f && total < 1.01f);

    if(verbose && print_clutter_impact_statement) {
      printf("PPUS: %f, Var Doppler PPUS: %f, var_doppler clutter: %f\n", 
          total_ppus, ppus_affected_by_var_doppler_clutter,
          ppus_affected_by_var_doppler_clutter / total_ppus);
    }
    if(verbose) printf("\n");
    return total_ppus;
  }
              

  int map_cluster_resources(int full_clusters_per_ppu, int point_clusters_per_ppu,
                            path_proc_unit& ppu, int input_replication_factor,
                            int num_links_sharing_ppu) {
      int extra_full_clusters_required = full_clusters_per_ppu * num_links_sharing_ppu -  
                                         ppu.num_full_clusters() * input_replication_factor;
      int remaining_ppu_full_clusters=0;
      //if we don't need any more, then we're good, don't add a penalty
      if(extra_full_clusters_required<0) {
        remaining_ppu_full_clusters=-extra_full_clusters_required; //don't add twice below
        extra_full_clusters_required=0; 
      }
      
      int extra_point_clusters_required = point_clusters_per_ppu * num_links_sharing_ppu +
                                          extra_full_clusters_required * ppu._coef_per_cluster -
                                          ppu.num_point_clusters() * input_replication_factor -
                                          remaining_ppu_full_clusters;
      return extra_point_clusters_required;
  }

  //fast==0 (fixed), fast==1
  float calc_links_per_ppu(path_proc_unit& ppu, DatapathOverheads& ov, drbe_wafer& wafer, 
                     bool is_direct_path, bool is_aidp,
                     int speed_txrx, bool clutter, bool just_clutter, float tx_duty_cycle,
                     bool verbose = false) {
    if(ppu._num_clusters == 0) {
      return 0;
    }

 
    // Can't have clutter and direct path or scatter seperable
    assert(!(clutter && (is_direct_path || is_aidp)));
    assert(!(just_clutter && !clutter));

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500.0;

    if(!_is_aidp && !_is_direct_path) {
      total_buffer_needed *= tx_duty_cycle; // multiply by duty cycle
    }

    if(is_direct_path || is_aidp) {
      total_buffer_needed /= 2;  //Only need half the range for direct path
    }

    if(verbose) {
      printf("input buffer length: %d\n",ppu._input_buffer_length);
      printf("total_buffer_needed: %d\n",total_buffer_needed);
    }

    int ppus_in_full_range = 1 + total_buffer_needed / ppu._input_buffer_length;

    if(verbose) {
      printf("ppus_in_full_range: %d, mem_ratio %f\n",ppus_in_full_range,ppu._mem_ratio);
    }

    //Calculate the number of MACCs/PPU (ideal case)
    ov.reason[DatapathOverheads::MAX] = ppu.num_full_clusters() * ppu._coef_per_cluster +
                                           ppu.num_point_clusters();


    // This is used for overhead calculation;
    int ideal_coef_per_ppu_per_link = -1;  //need to set this...


    float extra_obj_ratio_for_fast = 0;
    
    int objects_per_ppu;
    if(is_direct_path) {
      // also need to expand ppus_in_full_range due to i/o limitations
      // we're going to get _n_obj inputs
      int io_ppus_per_side = ceil((float)_n_tx / (float)ppu._output_router._in_degree);
      ppus_in_full_range = max(ppus_in_full_range, io_ppus_per_side * io_ppus_per_side);

      //objects are split up among these ppus
      objects_per_ppu = ceil(_n_obj/(float)ppus_in_full_range);
    } else if(is_aidp) {

      // Similar constraint here due to _k_rcs_points to aggregate
      int io_ppus_per_side = ceil(_k_rcs_points / (float)ppu._output_router._in_degree);
      ppus_in_full_range = max(ppus_in_full_range, io_ppus_per_side * io_ppus_per_side);

      // FOR ojects/ppu, THIS IS A BIT OF A WEIRD MODEL, 
      // SINCE COMPUTATION DOES NOT HAPPEN IN TERMS OF OBJECTS
      //instead, we will convert to effective number of object
      float coef_required = _k_rcs_points * _coef_per_rcs_point * 2;
      //compute the equivalent objects:
      objects_per_ppu = ceil(coef_required/_avg_coef_per_object);

      ov.reason[DatapathOverheads::MAX] = coef_required;

      // this will get multiplied out correctly later
      // Note also taht we don't both to split these up among ppus since we expect the
      // total clusters to be low anyways

    } else { // Tapped delay case
      float ideal_objects_per_ppu = 
        ceil( (_n_obj)/(float)ppus_in_full_range);
      
      // compute this PURELY for computing overheads
      ideal_coef_per_ppu_per_link = ideal_objects_per_ppu * _avg_frac_full_objects * 
                                    _avg_coef_per_object +
                                    ideal_objects_per_ppu * (1-_avg_frac_full_objects);

      if(ppu._is_dyn_reconfig) {
        objects_per_ppu= ideal_objects_per_ppu;  
      } else {
        objects_per_ppu = 
          ceil( (_n_obj - _n_full_range_obj)/(float)ppus_in_full_range) +
          _n_full_range_obj;
        extra_obj_ratio_for_fast = ideal_objects_per_ppu / objects_per_ppu;
      }
    } 

    ov.reason[DatapathOverheads::FAST_OBJ] = ov.reason[DatapathOverheads::MAX] * 
                                                extra_obj_ratio_for_fast;

    int full_clusters_per_ppu = 0;
    int point_clusters_per_ppu = 0;

    int full_obj = objects_per_ppu * _avg_frac_full_objects+0.999f; //round up
    int point_obj = objects_per_ppu - full_obj;

    int clusters_per_full_object = 
                ceil((float)_avg_coef_per_object / (float)ppu._coef_per_cluster);  

    // These are normal objects
    if(!just_clutter) {
       full_clusters_per_ppu += full_obj * clusters_per_full_object; 
       point_clusters_per_ppu += point_obj;
    }

    // JUST FOR OVERHEAD< Kind of a hard computation to get
    float max_links_for_full_clusters = 1000;
    float max_links_for_point_clusters = 1000;

    if(full_clusters_per_ppu) {
      max_links_for_full_clusters = (ppu.num_full_clusters()+ppu.num_point_clusters()/
                                     ppu._coef_per_cluster)/(float)full_clusters_per_ppu;
    }
    if(point_clusters_per_ppu) {
      max_links_for_point_clusters = (ppu.num_point_clusters()+ppu.num_full_clusters())
                                     /(float)point_clusters_per_ppu;
    }
    float max_links_pre_clutter=std::min(max_links_for_full_clusters,max_links_for_point_clusters);
    float max_coef_pre_clutter=std::min(max_links_pre_clutter*full_clusters_per_ppu*ppu._coef_per_cluster +
                               max_links_pre_clutter*point_clusters_per_ppu *
                               extra_obj_ratio_for_fast, ov.reason[DatapathOverheads::FAST_OBJ]);

    ov.reason[DatapathOverheads::POINT_FULL_RATIO] = max_coef_pre_clutter;



    float coef_pre_clutter = full_clusters_per_ppu * ppu._coef_per_cluster +
                             point_clusters_per_ppu;

    if(clutter) {
      full_clusters_per_ppu  += full_obj * clusters_per_full_object; 
      full_clusters_per_ppu  += point_obj;  // not a typo, we use a full cluster for one point obj
    }

    float coef_post_clutter = full_clusters_per_ppu * ppu._coef_per_cluster +
                             point_clusters_per_ppu;

    //Calculate the number of MACCs/PPU (if you had to count clutter)
    ov.reason[DatapathOverheads::CLUTTER] = ov.reason[DatapathOverheads::POINT_FULL_RATIO] * 
                                            coef_pre_clutter / coef_post_clutter;

    if(!is_direct_path && !is_aidp) {
      //finally, lets reduce the compute intensity by duty cycle
      full_clusters_per_ppu *= tx_duty_cycle;
      point_clusters_per_ppu *= tx_duty_cycle;
    }


    int input_replication_factor = ceil((float)full_clusters_per_ppu / 
                                        (float)ppu.num_full_clusters())-1;

    bool found_valid_repl_factor=false;
    int num_links_sharing_ppu=1;

    while(found_valid_repl_factor==false) {
      input_replication_factor++;

      int extra_point_clusters_required = 
        map_cluster_resources(full_clusters_per_ppu, point_clusters_per_ppu, ppu,
                              input_replication_factor, num_links_sharing_ppu);
    
      if(extra_point_clusters_required <= 0) {
        found_valid_repl_factor = true;
      }      
    } 

    if(verbose) {
      printf("obj/ppu/link %d, clusters/ppu/link : %d, obj input repl: %d\n",
          objects_per_ppu, full_clusters_per_ppu, input_replication_factor);
      printf("h/w param -- clusters %d, flex clusters %d\n",
          ppu._num_clusters, ppu._num_flexible_clusters);
    }

    assert(input_replication_factor > 0);


    if(!is_direct_path && !is_aidp) {
      //then we can potentially jam multiple links per PPU
      for(; num_links_sharing_ppu < ppu._output_router._in_degree; num_links_sharing_ppu++) {
       
        int extra_point_clusters_required = 
        map_cluster_resources(full_clusters_per_ppu, point_clusters_per_ppu, ppu,
                              input_replication_factor, num_links_sharing_ppu);

        if(extra_point_clusters_required > 0) {
          num_links_sharing_ppu--; //go back one : )
          break;
        }
      }
    }

    int total_ppus_considered = ppus_in_full_range * input_replication_factor;


    float pre_network_links_per_ppu = (float) num_links_sharing_ppu / (float) total_ppus_considered;
    ov.reason[DatapathOverheads::IMPERFECT_MAP] = 
      ideal_coef_per_ppu_per_link * pre_network_links_per_ppu;




    // We now need to take into account Wafer-level effects for coefficient updates
    // Note the following:
    // Fast moving objects get charged the high-update rate (others, the slow). 
    // Dynamically, we only need to send updates for actual coefficients.

    float metadata_size = 32;  //metatdata bits per coef
    float coef_size = 32; //bits per coefficient
    float doppler_coef_size = 32; //bits per coefficient
    //coef_bandwith in bit/cyc

    // If Either the source or destination is a fast object, then
    // we will need to update coefficients often!
    //  tx/rx:        Fixed/Slow  Fast
    //                -----------------------------
    //         Slow  | slow       fast
    //         Fast  | fast       fast
                                 

    float full_clusters_per_ppu_total = (float) full_clusters_per_ppu * num_links_sharing_ppu / 
                                        (float) input_replication_factor;

    float point_clusters_per_ppu_total = (float) point_clusters_per_ppu * num_links_sharing_ppu / 
                                         (float) input_replication_factor;

    float full_coef_bits_per_ppu = full_clusters_per_ppu_total * 
          (metadata_size + (coef_size + doppler_coef_size)  * ppu._coef_per_cluster);
                                   
    float point_coef_bits_per_ppu = point_clusters_per_ppu_total * 
          (metadata_size + coef_size + doppler_coef_size);

    float total_coef_bits_per_ppu = full_coef_bits_per_ppu + point_coef_bits_per_ppu;

    float frac_slow, frac_fast;
    if(speed_txrx==2) { //FAST sees all as fast
      frac_slow = 0;
      frac_fast = 1;
    } else {
      frac_slow  = (float)(_n_slow) / ((float)(_n_fast+_n_slow));
      frac_fast  = (float)(_n_fast) / ((float)(_n_fast+_n_slow));
    }

    float avg_coef_bw =0; //average bandwidth in bits/ns
    avg_coef_bw += total_coef_bits_per_ppu * frac_slow / _low_update_period;
    avg_coef_bw += total_coef_bits_per_ppu * frac_fast / _high_update_period;


    // with the above coef bandwidth, and the max input bandwidth, compute the maximum
    // depth we can send coefficients without losing bandwidth
    float max_ppu_coef_depth = ppu._coef_router._in_degree * ppu._coef_router._bitwidth / 
                                avg_coef_bw;


    float fraction_ppus_active_due_to_coeff = 1; //min(max_ppu_coef_depth/wafer.half_depth(),1.0f);   

    if(verbose) {
      printf("Half Depth %d, max_ppu_coef_depth %f, frac_active: %f\n", 
          wafer.half_depth(), max_ppu_coef_depth, fraction_ppus_active_due_to_coeff);
    }

    //printf("in_buf: %d, in_buf_area: %f\n", ppu._input_buffer_length, ppu.input_buf_area());
    //printf("ppus %d, repl %d, best_link/ppu %d,\n",ppus_in_full_range, input_replication_factor, 
    //    comp_constrainted_links_per_ppu);
    
    float ppus_per_link = (float) total_ppus_considered / (float) num_links_sharing_ppu / fraction_ppus_active_due_to_coeff;

    if(verbose) {
      printf("num_links_sharing_ppu %d, ppus per link %f\n",num_links_sharing_ppu,ppus_per_link);
    }

    float links_per_ppu = 1.0f/ppus_per_link;

    ov.reason[DatapathOverheads::COEF_NETWORK] = 
      ideal_coef_per_ppu_per_link * links_per_ppu;


    return links_per_ppu;
  }

  int algorithmic_num_links() {
    int links_per_band = _n_tx * _n_rx;
    int total_links = links_per_band * _n_bands;
    return total_links;
  }

  int num_links() {
    int links_per_band;
    if(_is_direct_path || _is_aidp) {
      // No need to include self links for direct path
      links_per_band = (_n_tx-1) * (_n_rx-1);
    } else {
      links_per_band = _n_tx * _n_rx;
    }

    int total_links = links_per_band * _n_bands;
    return total_links;
  }

  float frac_fast() {
    if(platforms()==0) return 0;
    return _n_fast / platforms();
  }


  float frac_max_links();
  float frac_max_link_complexity();

  float link_complexity() {
    float avg_coef_per_link = (_avg_coef_per_object * _avg_frac_full_objects +
                               1.0f            * (1 - _avg_frac_full_objects)) * _n_obj;
    return avg_coef_per_link;
  }

  float frac_full_range() {
    return _n_full_range_obj / (float)_n_obj;
  }

  int platforms() {
    return _n_slow + _n_fast + _n_fixed;
  }

  int reflectors() {
    return _n_obj;
  }

  int num_paths(){
    return _n_tx * _n_rx * reflectors() * _n_bands;
  }

  //TODO UPDATE THIS
  bool could_be_harder_than(Band& other) {
    //if(platforms() > other.platforms()) return true;
    //if(reflectors() > other.reflectors()) return true;
    //if(_n_bands < other._n_bands) return true;
    if(_n_obj > other._n_obj) return true; 
    if(frac_max_links() > other.frac_max_links()) return true;
    if(frac_max_link_complexity() > other.frac_max_link_complexity()) return true;
    if(frac_fast() > other.frac_fast()) return true;
    if(frac_full_range() > other.frac_full_range()) return true; 
    if(_frac_clutter > other._frac_clutter) return true;
    if(_range > other._range) return true;
    return false;
  }

  int increase_difficulty(int kind);
  float normalized_distance_to(Band& other);
  std::vector<float>& normalized_vec();


  void print_csv();
  //void print_csv() {
  //  printf("%d, %d, %d, %d, %d, %d", platforms(), reflectors(),
  //      _n_bands, _avg_coef_per_object, _n_fast, _range);
  //}

  void print_norm_csv();

  void recalculate_txrx() {
      // n_tx and n_rx indicate the number of 
      // transmitters and receivers PER BAND
      _n_tx = ceil((float) platforms() / _n_bands);
      _n_rx = ceil((float) platforms() / _n_bands) * 2; 
  }

  void print_detailed() {
    printf("      ");
    if(_is_direct_path) {
      printf("dp:");
    } else if(_is_aidp) {
      printf("aidp (kp:%d/%d)", _k_rcs_points, _coef_per_rcs_point);
    } else {
      printf("tdl:");
    }

    printf("tx %d, rx %d, obj %d, bands %d, slow %d, fast %d, fixed %d;"\
           "full_range %d, coef_per_obj %d, frac_full %f, range %d"\
           "frac_clutter %f, frac_pulsed %f, pulsed_duty %f, \
           update_low %f, update_high %f",
           _n_tx, _n_rx, _n_obj, _n_bands, _n_slow, _n_fast, _n_fixed, 
           _n_full_range_obj, _avg_coef_per_object, _avg_frac_full_objects, _range,
           _frac_clutter, _frac_pulsed, _pulsed_duty_cycle,
           _low_update_period, _high_update_period);
  }


  bool _is_direct_path=false;
  bool _is_aidp=false;

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_bands;

  int _k_rcs_points=20; //number of points in RCS model
  int _coef_per_rcs_point=20; //number of points in RCS model
 

  //Objects by speed
  int _n_slow=0;
  int _n_fast=0;
  int _n_fixed=0; 

  int _n_full_range_obj;
  int _avg_coef_per_object;
  float _avg_frac_full_objects;
  int _range;

  float _frac_clutter = 0.0;
  float _frac_pulsed = 0.0; // out of 1
  float _pulsed_duty_cycle = 0.5; //out of 1
  float _low_update_period=1000000; //
  float _high_update_period=5000; //clock cycles?

  std::vector<float> _norm_features;

};



struct GEStats {
// ---------------- Geometry Engine -----------------
  // Coordinate Translation
  void coordinate_trans_cmbl(ge_stat_per_band & fed_cmbl);
  // NR Engine
  void nr_engine_cmbl(ge_stat_per_band & fed_cmbl);
  // Relative Orientation
  void relative_orientation_cmbl(ge_stat_per_band & fed_cmbl);
  // Antenna Gain
  void antenna_gain_cmbl(ge_stat_per_band & fed_cmbl);
  // Path Gain and Velocity
  void path_gain_cmbl(ge_stat_per_band & fed_cmbl);
  // RCS
  void rcs_cmbl(ge_stat_per_band & fed_cmbl);
  // Update Time
  void tu_cmbl(ge_stat_per_band & fed_cmbl);
  // All
  void get_ge_cmbl(ge_stat_per_band & fed_cmbl);

  // Dump GE tradeoff
  std::string print_ge_tradeoff(Band & b, ge_stat_per_band & ge_stat, WaferStats & w_stats){
    std::stringstream ss;
    assert(!(b._is_aidp && b._is_direct_path));// can not be AIDP and DP the same time
    ss // Channel Model
      << (b._is_aidp ? "AIDP" : (b._is_direct_path ? "DP" : "TDL")) << ", "
      // Global Parameter
      << ge_stat.global_fid.upd_rate << ", "
      // Object(on sky, reflector), Transmitter, Receiver, Path and Bands
      << ge_stat.global_fid.num_obj << ", "                
      << b._n_tx << ", "
      << b._n_rx << ", "
      << ge_stat.global_fid.num_path << ", "
      << b._n_bands << ", "
      // ratio obj / tx or rx
      << (float)ge_stat.global_fid.num_obj/b._n_tx << ", "
      // flector by speed
      << b._n_fast << ", "
      << b._n_slow << ", "
      << b._n_fixed << ", "
      << b.platforms() << ", "
      // Wafer Level Statistic
      << w_stats.num_wafer <<", "
      << w_stats.num_ppu_chiplet << ", "
      << w_stats.num_ge_chiplet << ", "
      // Clutter
      << b._frac_clutter << ", "
      // Size
      << b._avg_coef_per_object << ", "
      // Range
      << b._range << ", "
      // GE total area
      << ge_stat.total_area << ", "
      << ge_stat.compute_area << ", "
      << ge_stat.mem_area << ", "
      // GE fidelity
      << ge_stat.coordinate_trans.ta1_scene_upd_rate << ", "
      << ge_stat.nr_engine.interpolation_ord << ", "
      << ge_stat.nr_engine.conv_fed << ", "
      << ge_stat.antenna_gain.order << ", "
      << ge_stat.antenna_gain.num_antenna << ", "
      << ge_stat.antenna_gain.res_angle << ", "
      << ge_stat.antenna_gain.dict_dim << ", "
      << ge_stat.rcs.order << ", "
      << ge_stat.rcs.points << ", "
      << ge_stat.rcs.angle << ", "
      << ge_stat.rcs.freq << ", "
      << ge_stat.rcs.plzn << ", "
      // Coordinate
      << ge_stat.coordinate_trans.compute << ", "
      << ge_stat.coordinate_trans.memory << ", "
      << ge_stat.coordinate_trans.bandwidth << ", "
      << ge_stat.coordinate_trans.latency << ", "
      // NR Engine
      << ge_stat.nr_engine.compute << ", "
      << ge_stat.nr_engine.memory << ", "
      << ge_stat.nr_engine.bandwidth << ", "
      << ge_stat.nr_engine.latency << ", "
      // Relative Orientation
      << ge_stat.relative_orientation.compute << ", "
      << ge_stat.relative_orientation.memory << ", "
      << ge_stat.relative_orientation.bandwidth << ", "
      << ge_stat.relative_orientation.latency << ", "
      // Antenna
      << ge_stat.antenna_gain.compute << ", "
      << ge_stat.antenna_gain.memory << ", "
      << ge_stat.antenna_gain.bandwidth << ", "
      << ge_stat.antenna_gain.latency << ", "
      // Path gain and velocity
      << ge_stat.path_velocity.compute << ", "
      << ge_stat.path_velocity.memory << ", "
      << ge_stat.path_velocity.bandwidth << ", "
      << ge_stat.path_velocity.latency << ", "
      // RCS
      << ge_stat.rcs.compute << ", "
      << ge_stat.rcs.memory << ", "
      << ge_stat.rcs.bandwidth << ", "
      << ge_stat.rcs.latency << ", "  
      // Tu
      << ge_stat.tu.compute << ", "
      << ge_stat.tu.memory << ", "
      << ge_stat.tu.bandwidth << ", "
      << ge_stat.tu.latency << "\n";
    return ss.str();
  }

  ge_stat_per_band * ge_stat_vec;

};

class ScenarioGen {

  public:
  static void gen_scenarios(std::vector<Band>& scene_vec, bool direct_path, bool aidp,
      int num_scenarios, int fixed_platforms=0, int easy_scenario=0) {

    if(num_scenarios<=1) {
      Band b;
      int platforms=max_platforms();
      if(fixed_platforms) {
        platforms=fixed_platforms;
      }

      // Number of platforms of each type
      if(easy_scenario==0) { // Hard hard case
        b._n_fast  = platforms;
        b._n_slow  = platforms - b._n_fast;
        b._n_fixed = 0;

        b._high_update_period = 5000; //5us
        b._range=max_range();

        b._n_bands = 1;

        b._n_obj = 100;

        b._avg_frac_full_objects=0.04;
        b._avg_coef_per_object = 40;

        b._n_full_range_obj = b._n_fast;

      } else if (easy_scenario==1) { // easy case
        b._n_fast  = platforms;
        b._n_slow  = platforms - b._n_fast;
        b._n_fixed = 0;

        b._high_update_period = 5000; //5us
        b._range=max_range();

        b._n_bands = 1;

        b._n_obj = 20;

        b._avg_frac_full_objects=0.2; //needs to be higher to meet the link_complexity
        b._avg_coef_per_object = 60;

        b._n_full_range_obj = b._n_fast;

      } else {
        b._n_fast  = 0;
        b._n_slow  = platforms - b._n_fast;
        b._n_fixed = 0;

        b._high_update_period = 5000; //5us
        b._range=max_range();

        b._n_bands = 1;

        b._n_obj = 20;

        b._avg_frac_full_objects=0.2; //needs to be higher to meet the link_complexity
        b._avg_coef_per_object = 60;

        b._n_full_range_obj = 0;
      }

      b.recalculate_txrx();
        
      b._is_direct_path=direct_path;
      b._is_aidp=aidp;

      b._frac_clutter=0.1;

      scene_vec.emplace_back(b); 

      num_scenarios=0; //set to zero
    }

    // Scenario Generation -- this is for random scenarios only
    for(int i_scene = 0; i_scene < num_scenarios; i_scene++) {
      Band b;
      b._is_direct_path=direct_path;
      b._is_aidp=aidp;

      int platforms = rand_rng(min_platforms(),max_platforms());

      int x = rand_rng(0,platforms);
      int y = rand_rng(1,platforms);
      if(x>y) std::swap(x,y);
      
      b._n_fixed = x;
      b._n_slow = y;
      b._n_fast = platforms-y;

      b._n_bands = rand_rng(min_bands(),max_bands());

      b.recalculate_txrx();

      b._n_obj = rand_rng(min_objects(),max_objects());

      b._n_full_range_obj = rand_rng(0,b._n_obj); //just random is fine

      b._avg_coef_per_object = 40; //(rand_rng(min_coef_per_obj(),max_coef_per_obj()) / 10) * 10;
      b._avg_frac_full_objects=rand_rng(min_full_obj()*50,max_full_obj()*50)/50.0;

      b._range = (rand_rng(min_range(),max_range())/10)*10;
      
      b._high_update_period = 5000; //rand_rng(5000 /*5us*/,100000 /*100us*/);

      b._frac_clutter = rand_rng(min_clutter()*100,max_clutter()*100)/100.0f;

  /*
      printf("Stationary obj: %d, Slow Obj: %d, Fast Obj: %d, bands %d,  \
          coef_per_obj: %d, range: %d, update period: %f\n", 
          b._n_fixed, b._n_slow, b._n_fast, b._n_bands, 
          b._avg_coef_per_object, b._range, b._high_update_period);
  */
      if(scenario_is_between_threshold_and_objective(b)) {
        scene_vec.emplace_back(b);
      } else {
        i_scene--;
      }
    }
  }

  static bool scenario_is_between_threshold_and_objective(Band& b) {
    if(b.num_links() < min_links() || b.num_links() > max_links()) return false;

//    if(b.link_complexity() < min_coef_per_link()) return false;

    //The extra 60 is so that the pareto exploration will not stop short
    if(b.link_complexity() > (overprov_link_complexity() + 20)) {
      return false;
    }
    return true;
  }

  // Define the geometry engine fidelity
  static void set_fidelity(std::vector<Band>& scene_vec, GEStats & ge_stats, drbe_wafer & w, bool ge_hard){
    int i = 0;
    for(auto & b : scene_vec){
      // Global
      ge_stats.ge_stat_vec[i].global_fid.upd_rate = ge_hard ? 5e-6 : 50e-6; //CASE: nominal 50e-6 max 5e-6
      ge_stats.ge_stat_vec[i].global_fid.num_obj = b._n_obj;
      ge_stats.ge_stat_vec[i].global_fid.num_path = b.num_paths();
      // Coordinate Transformation
      ge_stats.ge_stat_vec[i].coordinate_trans.ta1_scene_upd_rate = 1e-3; 
      // NR Engine
      ge_stats.ge_stat_vec[i].nr_engine.interpolation_ord = ge_hard ? 6 : 4; //CASE: nominal 4 max 6
      ge_stats.ge_stat_vec[i].nr_engine.conv_fed = 2;
      // Antenna
      ge_stats.ge_stat_vec[i].antenna_gain.order = 5;
      ge_stats.ge_stat_vec[i].antenna_gain.general = 1;
      ge_stats.ge_stat_vec[i].antenna_gain.num_antenna = 16;
      ge_stats.ge_stat_vec[i].antenna_gain.res_angle = 2;
      ge_stats.ge_stat_vec[i].antenna_gain.dict_dim = ge_hard ? 800 : 400; //CASE: nominal 400 max 800
      // RCS
      ge_stats.ge_stat_vec[i].rcs.order = ge_hard ? 6 : 4; //CASE: nominal 4 max 6
      ge_stats.ge_stat_vec[i].rcs.points = ge_hard ? 2 : 1; //CASE: nominal 1 max 2
      ge_stats.ge_stat_vec[i].rcs.plates = 1;
      ge_stats.ge_stat_vec[i].rcs.pntAngle = 1;
      ge_stats.ge_stat_vec[i].rcs.angle = 4;
      ge_stats.ge_stat_vec[i].rcs.freq = 1;
      ge_stats.ge_stat_vec[i].rcs.plzn = ge_hard ? 2 : 1; //CASE: nominal 1 max 2
      i++;
    }
  }

  static int overprov_link_complexity() { 
    return max_coef_per_link() * LINK_COMPLEXITY_OVERPROV;
  }

  //const static int clutter_links=0;
  static int min_platforms() {return 80;}
  static int max_platforms() {return 200;}

  static int min_coef_per_link() {return 70;}
  static int max_coef_per_link() {return 200;}

  static int min_links() {return 6400;}
  static int max_links() {return 80000;}

  static int min_objects() {return 10;} //threshold (obj/link)
  static int max_objects() {return 100;} //objective

  static int min_bands() {return 1;}
  static int max_bands() {return 4;}

  static int min_coef_per_obj() {return  20;}
  static int max_coef_per_obj() {return 100;}

  static int min_range() {return 200;} // in km
  static int max_range() {return 500;}

  static float min_clutter() {return 0;}
  static float max_clutter() {return 0.1;}

  static float min_full_obj() {return 0;}
  static float max_full_obj() {return 0.3;}
};

