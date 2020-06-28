#include "drbe.h"


using namespace std;

//using sec = chrono::seconds;
//using get_time = chrono::steady_clock;
//
//static struct option long_options[] = {
//    {"algorithm",      required_argument, nullptr, 'a',},
//    {"verbose",        no_argument,       nullptr, 'v',},
//    {0, 0, 0, 0,},
//};



float Band::normalized_distance_to(Band &other) {
  auto& my_vec = normalized_vec();
  auto& other_vec = other.normalized_vec();

  float dist = 0;
  for(unsigned i = 0; i < my_vec.size(); ++i) {
    float absdiff = abs(my_vec[i] - other_vec[i]);
    dist += absdiff*absdiff;
  }
  return dist;
}

bool Band::increase_difficulty() {
  int which = rand_bt(0,6);
  for(int i = 0; i < 6; ++i,which=(which+1)%6) {

    switch(which) {
      case 0: { //add a fixed object
        if(platforms() + 1 >ScenarioGen::max_platforms) break;
        _n_fixed+=1; //one object per band
        recalculate_txrx();
        return true;
      }
      case 1: { //upgrade a fixed to a slow object
        if(_n_fixed - 1 < 0) break;
        _n_fixed -= 1;
        _n_slow += 1;
        recalculate_txrx();
        return true;
      }
      case 2: { //upgrade a slow to a fast object
        if(_n_slow - 1 < 0) break;
        _n_slow -= 1;
        _n_fast += 1;
        recalculate_txrx();
        return true;
      }
      case 3: { //subtract a band
        if(_n_bands==1) break;
        _n_bands--;
        recalculate_txrx();
        return true;
      }
      case 4: { //avg_coef_per_obj
        if(_avg_coef_per_object +4 > ScenarioGen::max_coef_per_obj) break;
        _avg_coef_per_object+=4;
        return true;
      }
      case 5: { //range
        if(_range + 10 > ScenarioGen::max_range) break;
        _range+=10;
        return true;
      }
    }
  }
  return false;
}

  void Band::print_norm_csv() {
    printf("%0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f", 
        platforms() / (float) ScenarioGen::max_platforms, 
        reflectors() / (float)ScenarioGen::max_platforms,
        _n_fast / (float)ScenarioGen::max_platforms, 
        1.0/(float)_n_bands, 
        _avg_coef_per_object / (float)ScenarioGen::max_coef_per_obj, 
        _range / (float)ScenarioGen::max_range);
  }


std::vector<float>& Band::normalized_vec() {
  if(_norm_features.size()==0) {
    _norm_features.push_back(platforms() / (float)ScenarioGen::max_platforms);
    _norm_features.push_back(reflectors() / (float)ScenarioGen::max_platforms);
    _norm_features.push_back(1.0 / (float)_n_bands);
    _norm_features.push_back(_avg_coef_per_object / (float)ScenarioGen::max_coef_per_obj);
    _norm_features.push_back(_n_fast / (float)ScenarioGen::max_platforms);
    _norm_features.push_back(_range / (float)ScenarioGen::max_range);
  }
  return _norm_features;
}





int get_best_ppu() {

  return 0;
}

struct PPUStats {
  float avg_coef_per_ppu=0;
  float avg_coef_per_mm2=0;
  int total_failures=0;
  float avg_wafers=0;
  float percent_in_target_wafer=0;

  float avg_ppus_per_link=0;
  int wafer_histo[100000];
  int coef_per_mm2[100000]; 

  void print_wafer_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",wafer_histo[i]);
    }
  }
};

bool model_succesful(path_proc_unit* ppu, Band& band, drbe_wafer& w, int max_wafers) {
  float ppus_per_link=0;
  int coef_per_ppu = band.coef_per_ppu(*ppu,w,ppus_per_link,false/*verbose*/);
  if(coef_per_ppu==0) return false; // don't include this scenario if failed

  float ppus = band.num_links() * ppus_per_link;
  float num_wafers = ppus / w.num_units();

  return ceil(num_wafers) <= max_wafers;
}

std::vector<Band> top_k_pareto_scenarios(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, int k, int max_wafers) {
  std::vector<Band> top_k_scenarios;
  std::vector<Band> successful_pareto_bands;

  top_k_scenarios.resize(1);

  for(auto& band : scene_vec) {
    bool success = model_succesful(ppu,band,w,max_wafers);
    if(!success) continue;

    // Make sure we are harder than any currently successful band
    bool found_categorically_harder_band=false;
    for(auto& other : successful_pareto_bands) {
      if(!band.could_be_harder_than(other)) {
        //other.print_csv();
        //printf("   --dominates--   ");
        //band.print_csv();
        //printf("\n");
        found_categorically_harder_band=true;
        break;
      }
    }
  
    if(found_categorically_harder_band) continue;

    // Now lets increase the difficulty as much as possible:
    bool could_increase=true;
    Band new_band;

    int num_tries=0;

    while(could_increase && num_tries < 100) {
      new_band=band;
      
      if(new_band.increase_difficulty()) {
        if(model_succesful(ppu,new_band,w,max_wafers)) {
          std::swap(band,new_band); //swap in the new band
          num_tries=0;
        } else {
          num_tries++;
        }
      } else {
        could_increase=false;
      }
    }

    //while we do everything, check and see if we found a higher # platforms
    if(top_k_scenarios[0].platforms() < band.platforms()) {
      top_k_scenarios[0] = band;
    }

    //get rid of all the other bands we previously included but aren't necessary
    successful_pareto_bands.erase(
        std::remove_if(successful_pareto_bands.begin(), 
                       successful_pareto_bands.end(), 
                       [&band](Band& other) { return !other.could_be_harder_than(band); }), 
                       successful_pareto_bands.end());

    successful_pareto_bands.push_back(band);
  }

  //for(auto& band : successful_pareto_bands) {
  //  band.print_csv();
  //  printf("\n");
  //}
  //printf("num succ: %d \n", (int)successful_pareto_bands.size());


  //get the most different bands
  for(int i = 0; i < k; ++i) {
    float greatest_min_distance=0;
    Band chosen_band = successful_pareto_bands[0];
    for(auto& band : successful_pareto_bands) {
      float min_distance=100000000000;
      for(auto& top_k_band : top_k_scenarios) {
        float dist = band.normalized_distance_to(top_k_band);
        if(dist < min_distance) {
          min_distance = dist;
        }
      }
      if(min_distance > greatest_min_distance) {
        greatest_min_distance=min_distance;
        chosen_band=band;
      }
    }
    top_k_scenarios.push_back(chosen_band);
  }

  for(auto& band : top_k_scenarios) {
    band.print_norm_csv();
    printf("\n");
  }
  printf("num succ: %d \n", (int)successful_pareto_bands.size());


  return top_k_scenarios;
}

void evaluate_ppu(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, 
                   PPUStats& stats, int num_wafers_target, bool verbose = false) {
  stats.total_failures=0;

  float total_coef_per_ppu=0;
  float wafers=0, total_ppus_per_link=0;
  int total_in_target_wafer=0;

  for(auto& band : scene_vec) {
    float ppus_per_link=0;
    int coef_per_ppu = band.coef_per_ppu(*ppu,w,ppus_per_link,verbose);
    if(coef_per_ppu==0) stats.total_failures+=1;
    total_ppus_per_link += ppus_per_link;

    total_coef_per_ppu += coef_per_ppu;

    float ppus = band.num_links() * ppus_per_link;
    float num_wafers = ppus / w.num_units();
    int int_num_wafers = ceil(num_wafers);
    if(int_num_wafers<=num_wafers_target) total_in_target_wafer++;
    wafers += int_num_wafers;
    assert(int_num_wafers<100000);
    stats.wafer_histo[int_num_wafers]++;
    int int_avg_coef_per_mm2 = stats.avg_coef_per_ppu / w.chiplet_area();
    assert(int_avg_coef_per_mm2 < 100000);
    stats.coef_per_mm2[int_num_wafers]++;
    //printf("Avg wafer: %d, PPUs: %f, band.num_links():%d, ppus_per_link:%f\n",
    //    int_num_wafers, ppus, band.num_links(), ppus_per_link);
  }
  stats.percent_in_target_wafer = total_in_target_wafer / (float) scene_vec.size();
  stats.avg_ppus_per_link = total_ppus_per_link / scene_vec.size();
  stats.avg_wafers = wafers / scene_vec.size();
  stats.avg_coef_per_ppu = total_coef_per_ppu / scene_vec.size();
  stats.avg_coef_per_mm2 = stats.avg_coef_per_ppu / w.chiplet_area();
}

path_proc_unit* design_ppu_for_scenarios(std::vector<Band>& scene_vec, drbe_wafer& w) {
  PPUStats best_stats;
  path_proc_unit* best_ppu = NULL; //needs to be null so we don't delete a non-pointer
  float ppu_area = w.chiplet_area();
  tech_params& t = *w._t;

  // AggNet computaitons
  for(int agg_network = 1; agg_network < 20; agg_network+=1) {
    float side_length=sqrt((float)ppu_area);
    //side_length -= 1; //leave one mm out
    float bits_per_side=t._chiplet_io_bits_per_mm2*side_length;
    int total_ins_per_side = ceil(bits_per_side/32.0); //divide by 32 bit width
    int remain=total_ins_per_side - (2*2 + agg_network*2);
    if(remain < 2) break;

    // Loop over coefficients per cluster
    for(int cpc = 10; cpc < 40; cpc+=5) {

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
 
        PPUStats stats; 
        evaluate_ppu(ppu,scene_vec,w,stats,1);

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

  std::vector<Band> scene_vec;
  ScenarioGen::gen_scenarios(scene_vec);



  
  printf("\n  --- DSE Summary --- \n");

  //for(int ppu_area = 20; ppu_area < 21; ppu_area+=5) {
  int ppu_area=20;

  float fast_update_period=10000;
  //for(fast_update_period = 1000; fast_update_period < 1000000; 
  //    fast_update_period*=1.2589254117941672104239541063958) {

     //for(auto& band : scene_vec) {
     //  band._high_update_period=fast_update_period;
     //}

    drbe_wafer w(&t,300,(float)ppu_area);

  //for(int spmm = 1; spmm < 200; ++spmm) {
    //t._sram_Mb_per_mm2=spmm;

    //for(int agg_network = 1; agg_network < 200; agg_network+=1) {
    //agg_network=11

    path_proc_unit* best_ppu = design_ppu_for_scenarios(scene_vec,w);

    PPUStats stats;
    int num_wafers_target=1;
    evaluate_ppu(best_ppu,scene_vec,w,stats,num_wafers_target,false /*verbose*/);

    top_k_pareto_scenarios(best_ppu,scene_vec,w,5,num_wafers_target);

    printf("Fast Update Period: %0.0fus, ", fast_update_period/1000);
    printf("%dmm^2 PPU (%0.2f), in-MB: %0.2f, clust: %d, coef_per_clust %d, "\
           "Agg nets: %d, Mem Ratio: %d, Coef per mm2: %0.2f, "\
           "avg_wafers: %0.2f, avg_ppus_per_link: %0.2f, per_1_wafer:%0.1f, fails: %d\n",
        ppu_area, 
        best_ppu->area(), 
        best_ppu->input_buf_area()*t._sram_Mb_per_mm2/8,
        best_ppu->_num_clusters, best_ppu->_coef_per_cluster,
        best_ppu->_output_router._in_degree, 
        (int)best_ppu->_mem_ratio, 
        stats.avg_coef_per_mm2,
        stats.avg_wafers,
        stats.avg_ppus_per_link,
        stats.percent_in_target_wafer*100,
        stats.total_failures);
    //printf("Wafers Needed Histogram: ");
    //stats.print_wafer_histo();
    //printf("\n");

    //best_ppu->print_params();
    //best_ppu->print_area_breakdown();
  //}


  return 0;
}

