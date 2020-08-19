#include <assert.h>
#include <getopt.h>
#include <signal.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <dsp.h>
#include "hw.h"
#include "util.h"

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
  float samples = 0.0;
};

struct st_tu : cmbl {};

struct single_band_ge_stats {
  st_global_fid global_fid;
  st_coordinate_trans coordinate_trans;
  st_nr_engine nr_engine;
  st_relative_orientation relative_orientation;
  st_antenna_gain antenna_gain;
  st_path_velocity path_velocity;
  st_rcs rcs;
  st_tu tu;
};




//Note for later: Sumeet says that if either tx or rx is fast, then we need variable doppler
//yikes

class Band {
  public:

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
  float max_systems_per_wafer(int wafer_inputs) { 
    if(_is_direct_path) {
      return wafer_inputs / (float)(max(_n_obj-1,1)); 
    } else if (_is_aidp) {
      return wafer_inputs / _k_rcs_points;
    } else {
      return wafer_inputs;
    }
  }


  float calc_links_per_wafer(float wafer_unconstrained_ppus, path_proc_unit* ppu, drbe_wafer& w, bool verbose = false) {

    //We're goint to just cut the biggest square we can out, and assume that's what we'll
    //use all over.

    float num_wafers = wafer_unconstrained_ppus / (w.num_units() * ppu->_ppus_per_chiplet);
    
    float links_per_wafer = num_links() / num_wafers;

    if(w.limit_wafer_io()==false) {
      return links_per_wafer;
    }

    float systems_per_wafer = sqrt(links_per_wafer); 

    systems_per_wafer = min(systems_per_wafer,max_systems_per_wafer(w.max_input()));

    links_per_wafer = systems_per_wafer * systems_per_wafer;

    return links_per_wafer;
  }

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

    //iterate over combinations of fast and slow ppus
    for(int speed = 0; speed < 3; ++speed) {
      for(int clutter = 0; clutter < 2; ++clutter) {
        float frac_speed = frac_speed_vec[speed];
        float frac_clutter = clutter==0 ? (1-native_frac_clutter) : native_frac_clutter;
        float frac= frac_speed * frac_clutter;

        assert(frac_speed >=0 && frac_clutter >=0 && frac >=0);
        assert(frac_speed <=1 && frac_clutter <=1 && frac <=1);
        if(frac==0) continue; //save some time : )

        if(verbose) {
          print_clutter_impact_statement=true;
          printf("MAPPING WITH PARAMS clutter: %d, speed %d, frac %f\n",clutter,speed, frac);
        }

        float links_per_ppu = calc_links_per_ppu(ppu,w,_is_direct_path, _is_aidp,                                                                   speed,clutter,false /*just clutter*/,verbose);
        if(links_per_ppu==0) {
          failures+=1;
          continue;
        }

        if(verbose) printf("\n");


        float ppus = ceil(num_links() * frac / links_per_ppu);

        if(clutter && speed !=0) {
          ppus_affected_by_var_doppler_clutter+=ppus;
        }
        //printf("%f %f %f %f; ", frac_speed, frac_clutter,frac,ppus);
        total_ppus +=ppus;
      }
    }
    //printf(" PPUs\n");

    //Need to pay for FULL tapped delay lines in direct-path models
    if(_is_direct_path || _is_aidp) {
      for(int speed = 0; speed < 3; ++speed) {
        float frac_speed = frac_speed_vec[speed];
        float frac= frac_speed * _frac_clutter;
        assert(frac_speed >=0 && frac >=0);
        assert(frac_speed <=1 && frac <=1);
        if(frac==0) continue; //save some time : )

        if(verbose) {
          print_clutter_impact_statement=true;
          printf("Separate TDL MAPPING WITH PARAMS for clutter: speed %d, frac %f\n",speed, frac);
        }

        float links_per_ppu = calc_links_per_ppu(ppu,w,false,false/*tap-delay*/,
                                                 speed,true/*clutter*/,true/*just clutter*/,
                                                 verbose);
        if(links_per_ppu==0) {
          failures+=1;
          continue;
        }

        if(verbose) printf("\n");

        float ppus = ceil(num_links() * frac / links_per_ppu);

        if(speed !=0) {
          ppus_affected_by_var_doppler_clutter+=ppus;
        }
        //printf("%f %f %f %f; ", frac_speed, frac_clutter,frac,ppus);
        total_ppus +=ppus;
      }
    }

    if(verbose && print_clutter_impact_statement) {
      printf("PPUS: %f, Var Doppler PPUS: %f, var_doppler clutter: %f\n", 
          total_ppus, ppus_affected_by_var_doppler_clutter,
          ppus_affected_by_var_doppler_clutter / total_ppus);
    }
    if(verbose) printf("\n");
    return total_ppus;
  }
              

  //fast==0 (fixed), fast==1
  float calc_links_per_ppu(path_proc_unit& ppu, drbe_wafer& wafer, 
                     bool is_direct_path, bool is_aidp,
                     int speed_txrx, bool clutter, bool just_clutter,
                     bool verbose = false) {
    if(ppu._num_clusters == 0) {
      return 0;
    }

    // Can't have clutter and direct path or scatter seperable
    assert(!(clutter && (is_direct_path || is_aidp)));
    assert(!(just_clutter && !clutter));

    // First compute degree of multi-system required
    int total_buffer_needed = 3300000.0 * _range / 500.0;  

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

    int objects_contributed_per_ppu;
    if(is_direct_path) {
      // also need to expand ppus_in_full_range due to i/o limitations
      // we're going to get _n_obj inputs
      int io_ppus_per_side = ceil((float)_n_obj / (float)ppu._output_router._in_degree);
      ppus_in_full_range = max(ppus_in_full_range, io_ppus_per_side * io_ppus_per_side);

      objects_contributed_per_ppu = ceil(_n_obj/(float)ppus_in_full_range);
    } else if(is_aidp) {
      // Similar constraint here due to _k_rcs_points to aggregate
      int io_ppus_per_side = ceil(_k_rcs_points / (float)ppu._output_router._in_degree);
      ppus_in_full_range = max(ppus_in_full_range, io_ppus_per_side * io_ppus_per_side);

      // FOR ojects/ppu, THIS IS A BIT OF A WEIRD MODEL, 
      // SINCE COMPUTATION DOES NOT HAPPEN IN TERMS OF OBJECTS
      //instead, we will convert to effective number of object
      float coef_required = _k_rcs_points * _coef_per_rcs_point;
      //compute the equivalent objects:
      objects_contributed_per_ppu = ceil(coef_required/_avg_coef_per_object);

      // this will get multiplied out correctly later
      // Note also taht we don't both to split these up among ppus since we expect the
      // total clusters to be low anyways

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

    int clutter_objects=0;


    if(clutter) {
      clutter_objects=objects_contributed_per_ppu;
      
      if(!just_clutter) { //if its not just clutter, then add these objects
        objects_contributed_per_ppu+=clutter_objects;
      }
      if(verbose) {
        printf("Adding clutter; objects: %d, clutter objects: %d\n",
            objects_contributed_per_ppu,clutter_objects);
      }
    }

    

    // ciel: b/c coef_per_clust may be too large
    int clusters_required_per_ppu_per_link = objects_contributed_per_ppu * 
                        ceil((float)_avg_coef_per_object / (float)ppu._coef_per_cluster);  

    // How many replicas do we need to make
    int obj_input_replication_factor = ceil(clusters_required_per_ppu_per_link 
                                        / (float)ppu._num_clusters);



    //how many flexible cluters do we require
    //1. if we are fast or slow, then we require clutter clusters to be flexible
    //2. if we are fixed, no clutter clusters need to be flexible

    int clutter_input_replication_factor = 0;
  
    int clutter_clusters_required_per_ppu_per_link=0;
    if(clutter && speed_txrx!=0) { //if not fixed txrx
      // We're goint to need flexible clusters...
      clutter_clusters_required_per_ppu_per_link = clutter_objects * 
                          ceil((float)_avg_coef_per_object / (float)ppu._coef_per_cluster);  

      // How many replicas do we need to make
      clutter_input_replication_factor = ceil(clutter_clusters_required_per_ppu_per_link 
                                          / (float)ppu.num_flexible_clusters());
    }


    int input_replication_factor = max(obj_input_replication_factor,
                                     clutter_input_replication_factor);


    if(verbose) {
      printf("obj/ppu/link %d, clusters/ppu/link : %d, obj input repl: %d\n",
          objects_contributed_per_ppu, clusters_required_per_ppu_per_link, obj_input_replication_factor);
      printf("h/w param -- clusters %d, flex clusters %d\n",
          ppu._num_clusters, ppu._num_flexible_clusters);
      printf("CLUTTER: obj/ppu/link %d, clutter input repl: %d\n", clutter_objects, 
          clutter_input_replication_factor);
    }


    assert(input_replication_factor > 0);

    //if(clusters_required_per_ppu_per_link > ppu._num_clusters) return 0;

    //Redistributed the links in integer number of PPUS
    int obj_comp_constrainted_links_per_ppu = 
      input_replication_factor * ppu._num_clusters / clusters_required_per_ppu_per_link;

    int comp_constrainted_links_per_ppu;

    if(clutter && speed_txrx!=0) {
      int clutter_comp_constrained_links_per_ppu = 
        clutter_input_replication_factor * 
        ppu._num_flexible_clusters / clutter_clusters_required_per_ppu_per_link;

        comp_constrainted_links_per_ppu = min(obj_comp_constrainted_links_per_ppu,
                                                  clutter_comp_constrained_links_per_ppu);
    } else {
        comp_constrainted_links_per_ppu = obj_comp_constrainted_links_per_ppu;
    }
   
    assert(comp_constrainted_links_per_ppu >= 1);
      
      
        


    //TODO: for now i don't think its worth it to consider multiple links sharing the same
    //PPU for a direct-path model, b/c the I/O will almost always be the constraint.
    //similar reasoning for scatter sep, except the limitation is just memory
    int max_sharing_factor = (is_direct_path || is_aidp) ? 1 : ppu._output_router._in_degree;

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
    if(speed_txrx==2) { //FAST sees all as fast
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
      printf("num_links_sharing_ppu %d, ppus per link %f\n",num_links_sharing_ppu,ppus_per_link);
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
    if(_is_direct_path || _is_aidp) {
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

  int num_paths(){
    return _n_tx * _n_rx * reflectors();
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

  // ---------------- Geometry Engine -----------------
  // Coordinate Translation
  void coordinate_trans_cmbl(single_band_ge_stats & fed_cmbl);
  // NR Engine
  void nr_engine_cmbl(single_band_ge_stats & fed_cmbl);
  // Relative Orientation
  void relative_orientation_cmbl(single_band_ge_stats & fed_cmbl);
  // Antenna Gain
  void antenna_gain_cmbl(single_band_ge_stats & fed_cmbl);
  // Path Gain and Velocity
  void path_gain_cmbl(single_band_ge_stats & fed_cmbl);
  // RCS
  void rcs_cmbl(single_band_ge_stats & fed_cmbl);
  // Update Time
  void tu_cmbl(single_band_ge_stats & fed_cmbl);
  // All
  void get_ge_cmbl(single_band_ge_stats & fed_cmbl);

  bool _is_direct_path=false;
  bool _is_aidp=false;

  int _n_tx;
  int _n_rx;
  int _n_obj; //multipath
  int _n_bands;

  int _k_rcs_points=20; //number of points in RCS model
  int _coef_per_rcs_point=20; //number of points in RCS model
 

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
  float _frac_clutter = 0.0;

  // Geometry Engine Tradeoff
  single_band_ge_stats ge_stat;

  std::vector<float> _norm_features;
};


class ScenarioGen {

  public:

  static void gen_scenarios(std::vector<Band>& scene_vec, bool direct_path, bool aidp,
      bool challenge_scenario, int num_scenarios) {

    if(challenge_scenario) {
      Band b;
      b._n_fixed = max_platforms() -  max_reflectors();
      b._n_slow = 0;
      b._n_fast = max_reflectors();
      b._n_bands = 1;
      b.recalculate_txrx();
      b._avg_coef_per_object = 60;
      b._n_full_range_obj = b._n_fast;
      b._range=max_range();
      b._high_update_period = 10000; //10us
      b._is_direct_path=direct_path;
      b._is_aidp=aidp;
      scene_vec.emplace_back(b); 

      num_scenarios=0; //set to zero
    }

    // Scenario Generation
    for(int i_scene = 0; i_scene < num_scenarios; i_scene++) {
      Band b;
      int platforms = rand_rng(min_platforms(),max_platforms());
      int max_rand_reflect = min(max_reflectors(),platforms);
      int num_reflectors = rand_rng(1,max_rand_reflect);

      int slow = rand_rng(0,num_reflectors);

      b._is_direct_path=direct_path;
      b._is_aidp=aidp;

      b._n_fixed = platforms-num_reflectors;
      assert(b._n_fixed>=0);
      b._n_slow  = slow;
      b._n_fast = num_reflectors-slow;
      //int half_fast = b._n_fast/2;
      //b._n_fast -= half_fast;
      //b._n_slow += half_fast;
  
      b._n_bands = rand_rng(min_bands(),max_bands());

      b.recalculate_txrx();
  
      b._avg_coef_per_object = (rand_rng(min_coef_per_obj(),max_coef_per_obj()) / 4) * 4;
      b._n_full_range_obj = b._n_fast;
      b._range = (rand_rng(min_range(),max_range())/10)*10;
      
      b._high_update_period = rand_rng(10000 /*10us*/,100000 /*100us*/);

      b._frac_clutter=rand_rng(min_clutter(),max_clutter())/100.0f;

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

  //const static int clutter_links=0;
  static int min_platforms   () {return min(max_platforms(),80);}
  static int max_platforms   () {return 200;}
  static int max_reflectors  () {return min(max_platforms(),100);}
  static int min_bands       () {return 1;}
  static int max_bands       () {return 4;}
  static int min_coef_per_obj() {return 20;}
  static int max_coef_per_obj() {return 60;}
  static int min_range       () {return 50;}
  static int max_range       () {return 500;}
  static int min_clutter     () {return 0;}
  static int max_clutter     () {return 30;}
};

