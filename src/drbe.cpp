#include "drbe.h"
using namespace std;

static struct option long_options[] = {
    {"direct-path",        no_argument,       nullptr, 'd',},
    {"kp-aidp",            no_argument,       nullptr, 'k',},
    {"dynamic-reconfig",   no_argument,       nullptr, 'r',},
    {"print-pareto",       no_argument,       nullptr, 'p',},
    {"challenge-scenario", no_argument,       nullptr, 'c',},
    {"num-scenarios",      required_argument, nullptr, 'n',},
    {"verbose",            no_argument,       nullptr, 'v',},
    {"verbose-ge-info",    no_argument,       nullptr, 'g',},
    {"wafer-io",           no_argument,       nullptr, 'w',},
    {0, 0, 0, 0,},
};

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

bool Band::increase_difficulty(int kind=-1) {

  int which = rand_bt(0,6);
  int max_types=6;

  if(kind!=-1) {
    which=kind;
    max_types=1;
  }

  for(int i = 0; i < max_types; ++i,which=(which+1)%6) {

    switch(which) {
      case 0: { //add a fixed object
        if(platforms() + 1 >ScenarioGen::max_platforms()) break;
        _n_fixed+=1; //one object per band
        recalculate_txrx();
        return true;
      }
      case 1: { //upgrade a fixed to a slow object
        if(_n_fixed - 1 < 0) break;
        if(reflectors() +1 >= ScenarioGen::max_reflectors()) break;
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
        if(_n_bands==1 || rand_bt(0,20)!=0) break; //rand check to make it less likely
        _n_bands--;
        recalculate_txrx();
        return true;
      }
      case 4: { //avg_coef_per_obj
        if(_avg_coef_per_object +1 > ScenarioGen::max_coef_per_obj()) break;
        _avg_coef_per_object+=1;
        return true;
      }
      case 5: { //range
        if(_range + 10 > ScenarioGen::max_range()) break;
        _range+=10;
        return true;
      }
    }
  }
  return false;
}

  void Band::print_norm_csv() {
    printf("%0.3f, %0.3f, %0.3f, %0.3f, %0.3f, %0.3f", 
        platforms() / (float) ScenarioGen::max_platforms(), 
        reflectors() / (float)ScenarioGen::max_platforms(),
        _n_fast / (float)ScenarioGen::max_platforms(), 
        1.0/(float)_n_bands, 
        _avg_coef_per_object / (float)ScenarioGen::max_coef_per_obj(), 
        _range / (float)ScenarioGen::max_range());
  }



std::vector<float>& Band::normalized_vec() {
  if(_norm_features.size()==0) {
    _norm_features.push_back(platforms() / (float)ScenarioGen::max_platforms());
    _norm_features.push_back(reflectors() / (float)ScenarioGen::max_platforms());
    _norm_features.push_back(1.0 / (float)_n_bands);
    _norm_features.push_back(_avg_coef_per_object / (float)ScenarioGen::max_coef_per_obj());
    _norm_features.push_back(_n_fast / (float)ScenarioGen::max_platforms());
    _norm_features.push_back(_range / (float)ScenarioGen::max_range());
  }
  return _norm_features;
}

  // ---------------- Geometry Engine -----------------
  // Coordinate Translation
  void Band::coordinate_trans_cmbl(single_band_ge_stats & fed_cmbl){
    float NO = fed_cmbl.global_fid.num_obj;
    float SU = fed_cmbl.coordinate_trans.ta1_scene_upd_rate;
    float C = NO*(13 + ComputeTable(Cos,Corl_C)+
                      ComputeTable(Sin,Corl_C)+
                     ComputeTable(Sqrt,Corl_C) +
                     ComputeTable(Div,Corl_C))/SU; // Compute - in terms of 64bit FP MACs
    float M = 2+3*NO; // Memory in terms of 64 bit FP Units(reused varable?)
                      //considering reused varable
                      //M= 2+6*NO
    float B = 3*NO/SU; // Memory Bandwidth in terms of 64bit FP I/O operations required
    float L = 0; // Latency in terms of number of clock cycles required for this step
    fed_cmbl.coordinate_trans.set_cmbl(C, M, B, L);
  }
  // NR Engine
  void Band::nr_engine_cmbl(single_band_ge_stats & fed_cmbl){
    // Given Fidelites (Change these according to your blocks inputs)
    float IO = fed_cmbl.nr_engine.interpolation_ord;
    float SU = fed_cmbl.coordinate_trans.ta1_scene_upd_rate;
    float TU = fed_cmbl.global_fid.upd_rate;
    float CF = fed_cmbl.nr_engine.conv_fed;
    float NP = fed_cmbl.global_fid.num_path;
    float NO = fed_cmbl.global_fid.num_obj;
    // compute table comtains the equivilent MAC operational values of the
    // Sine/exponent/division functions etc.
    // For this block we need Division which is located at ComputeTable('div_HF','C')
    
    float C = 0; // Compute - in terms of 64bit FP MACs
    float M = 0; // Memory in terms of 64 bit FP Units
    float B = 0; // Memory Bandwidth in terms of 64bit FP I/O operations required
    float L = 0; // Latency in terms of number of clock cycles required for this step

    // intial scenario extrapolation equation 6th order usually if 6 vaues
    // are given
    C = C + (3*IO*IO + IO + 2)*(1/SU)*NO;
    M = M + (2*IO + 3)*(1/SU)*NO;
    B = B + IO*(1/SU)*NO;
    L = L + ceil(log2(IO*(1/SU)));
    // Each cycle equation parameter update with t0 correction
    C = C + (3*3*IO)*NP*(1/TU);
    L = L + ceil(log2(IO));
    // NR_ iteration caluclations
    C = C + ((3*2*IO) + ComputeTable(Div_HF,Corl_C))*CF*NP*(1/TU);
    M = M + 2*IO*CF*NP;
    L = L + CF*(ceil(log2(2*IO)) + ComputeTable(Div_HF,Corl_L));
    B = B + 2*NP*(1/TU);
    // Post Tn computaion of devived position velocity and acceleration
    C = C + 3*6*IO*NP*(1/TU);
    M = M + 9*NP;
    B = B + 3*3*IO*NP*(1/TU);
    L = L + ceil(log2(IO));
    
    // THis calculation has to be done 2X per path
    C = 2*C;
    M = 2*M;
    B = 2*B;
    L = 2*L;
    fed_cmbl.nr_engine.set_cmbl(C, M, B, L);
  }
  // Relative Orientation
  void Band::relative_orientation_cmbl(single_band_ge_stats & fed_cmbl){
    //NO = NumberofObjects; // don"t think this is correct should be paths
    float NO = fed_cmbl.global_fid.num_obj;
    float TU = fed_cmbl.global_fid.upd_rate;
    
    //transformation matrix composition 4x4x4 =64
    //each coordinate transformation 4x4=16  
    //incident and radiating angle calculation ?(under development)    
    //velocity transformation matrix 3*3=9
    //velocity transformation 4*4
    //acceleration similar 9+16=25
    
    //theta phi calculaiton x2
    //relative pos 3
    //theta= arctan(y/x)  arctan+div
    //phi  2+sqrt+div+arctan
    float C1 = (64+16+25*2+2*(3+2+2*ComputeTable(Div,Corl_C)+
        2*ComputeTable(Sqrt,Corl_L)+
        2*ComputeTable(Arctan,Corl_C)));
    float C = NO*C1/TU; // Compute - in terms of 64bit FP MACs
    
    
    //input orientation 3*3=9
    //relative orientation 3*3 = 9
    //Position, Velocities and Accelerations in&out      9*2=18
    //output angles 4
    
    //M = NO*(40); // Memory in terms of 64 bit FP Units
    float M = 0; // dont think you need to store the rest of the data that star thinks needs to be stored

    
    //B = NO*(40)/TU; // Memory Bandwidth in terms of 64bit FP I/O operations required
    float B = NO*(9+18)/TU; // input orientation and output andgeles dont need to be accesed
    
    
    float L = 0; // Latency in terms of number of clock cycles required for this step
    L = ceil(log2(C1));
    fed_cmbl.relative_orientation.set_cmbl(C, M, B, L);
  }
  // Antenna Gain
  void Band::antenna_gain_cmbl(single_band_ge_stats & fed_cmbl){
    float no_antenna=fed_cmbl.antenna_gain.num_antenna;
    float o=fed_cmbl.antenna_gain.order;
    float Tu=fed_cmbl.global_fid.upd_rate;
    float K=fed_cmbl.antenna_gain.res_angle;
    float no_paths=fed_cmbl.global_fid.num_path;
    float no_Tx=_n_tx;
    float no_Rx=_n_rx;
    float dict_size=fed_cmbl.antenna_gain.dict_dim;
    
    float C=0; //Compute in 64-bit MAC
    float M=0; //Memory in terms of 64 bit FP Units
    float B=0; //Memory Bandwidth in terms of 64bit FP I/O operations required
    float L=0; //Latency in terms of number of clock cycles required for this step
    
    if(o==0){
        C=0;
        M=M+1;
        B=B+(no_paths*2/Tu); //getting gain value per update time
        L=20; //not sure about latency
    }else if(o==1){
        C=0;
        M=M+(no_Tx + no_Rx); //gains for each antenna
        B=B+((no_paths*2)/Tu); // getting each antenna's gain per update time
        L=20;
    }else if(o==2){
        C=C+2*(8)+1;
        C=no_paths*(C/Tu);
        M=M+(5*(no_Tx + no_Rx)+1);
        B=B+((5*(no_paths*2)+1)/Tu); //getting each antenna's gains and beamwidths per update time
        L=21;
    }else if(0==3){
        C=C+2*(8)+1;
        C=no_paths*(C/Tu);
        M=M+(6*(no_Tx + no_Rx));
        B=B+((6*(no_paths*2))/Tu); //getting each antenna's gains and beamwidths per update time
        L=21;
    }else if(o==4){
        K=pow(2, (ceil(log2(180/K))));
        C=C+2*(8)+1; //to get 2nd/3rd order gains
        float C1=K*(2*ComputeTable(Cos,Corl_C)+3)+no_antenna+(no_antenna/2)*log2(no_antenna)+no_antenna+K+(K/2)*log2(K)+2*K;
        C=C+C1;
        C=2*C+1;
        C=no_paths*(C/Tu);
        M=M+(6*(no_Tx + no_Rx))+(no_Tx + no_Rx)*(4*K+2*no_antenna); //(2K+N)2 (bytes to bits)
        B=B+((6*(no_paths*2))/Tu)+(4*K+2*no_antenna)*no_paths/Tu;
        L=21+ceil(log2(C1));
    }else if(o==5){
        float N=(180/K)*(90/K)*(90/K);
        float n=(180/K);
        M=(no_Tx + no_Rx)*(n*dict_size+dict_size*N);
        C=C+1620+dict_size; //assuming 540 MACs for absolute value calc
        C=no_paths*(C/Tu);
        B=B+((dict_size*2)*(no_paths*2)/Tu); //getting 3 numbers from memory per update time;
        L=20+ceil(log2(1620))+ceil(log2(dict_size));
    }

    fed_cmbl.antenna_gain.set_cmbl(C, M, B, L);
  }
  // Path Gain and Velocity
  void Band::path_gain_cmbl(single_band_ge_stats & fed_cmbl){
    // Given Fidelites (Change these according to your blocks inputs)
    float SU = fed_cmbl.coordinate_trans.ta1_scene_upd_rate;
    float TU = fed_cmbl.global_fid.upd_rate;
    float NF = num_paths() * 2;
    
    // compute table comtains the equivilent MAC operational values of the
    // Sine/exponent/division functions etc.
    // For this block we need Division which is located at ComputeTable('div_HF','C')
    
    float C = 0; // Compute - in terms of 64bit FP MACs
    float M = 0; // Memory in terms of 64 bit FP Units
    float B = 0; // Memory Bandwidth in terms of 64bit FP I/O operations required
    float L = 0; // Latency in terms of number of clock cycles required for this step

// ------ Path distance compuation ------
//     Path_distance is used to calculate the coefficients of fractional delay filter, path gain.
//     The path_distance at each scinero update is calculated with equation below
//     relative_location = {(x1-x2), (y1-y2), (z1-z2)};
//     Path_distance = ((x1-x2)^2 + (y1-y2)^2 +(z1-z2)^2)^0.5
//     The path_distance at each update is calculated by first update the relative_location with relative_location increment,
//     then do the distance calcualtion with new relative_location
//     relative_location = relative_location + relative_location_increment;
//     relative_location_increment = v * t + 1/2 * a * t^2;
//     v: path velocity, the detailed equation is givien in the next section.
//     a: path accelelation, the detailed euqation is given in the next section. 
//     t: the time interval between two updates.
    C = C + (3 / TU) + 1 / TU + (3 / TU) + ComputeTable(Sqrt,Corl_C); // 5 * TU to get the relative_location at each update(including scenarioUpdate),  3 to get sum of square, 
    // 1 at updating rate for the distance increment, we get new v at each update, yet the a get updated only at scenario update.
    M = M + 3; // 3 to store relative_location.
    B = B + 6; // Bandwidth to transmit 2 set of coordinates
    L = L + 2 + ComputeTable(Sqrt,Corl_L); // 1 for relative_location update, 1 for sum of squares, 1 for square root 
//% End of the path distance compuation


    //% Path veloctiy compuation
    //     Path_velocity is calculated by first calculate the path velocity between two object at each scienro update
    //     second we calculate the path acceleration at scienro update.
    //     Relative_velocity = {v1_x - v2_x, v1_y - v2_y, v1_z - v2_z};
    //     Relative_acceleration = t * {a1_x - a2_x, a1_y - a2_y, a1_z - a2_z};
    //     path_velocity = sqrt(relative_velocity.^2);
    //     the relative speed at each update is calculated with equation below
    //     relative_velocity = relative_velocity + relative acceleration
    //     t: the time interval between two updates
    C = C + 6 / SU + 3 / TU + (3 + ComputeTable(Sqrt,Corl_C)) / TU; // 6 floating add at SU for relative velocity and acceleration,
                                                                                                  // 3 floating addition for relative_velocity, 3 mac and 1 sqrt for path_velocity, 
    M = M + 6; // 3 for path_velocity, 3 for path_acceleration
    L = L + 3; // 1 for floating additiona, 1 for floating multiplication, 1 for sqrt
    B = B + 12 / SU;// 6 for two sets of path_velocity, 6 for two sets of path acceleration
    //% end of the Path veloctiy compuation


    //% Path gain computation
    //     Path_gain is calculated by first evluate the path distance at scenario update, the compelete equation for the gain is the  
    //     Path_gain = 1/(16pi^2) * 1/path_distance^2 * 1/f0^2
    //     f0: the carrier frequency of the channel 

    C = C + ComputeTable(Div_HF,Corl_C) / TU; // path_distance^2 is available from the path distance computation, only one floating division.
    M = M + 1; // 1 for constant;
    L = L + ComputeTable(Div_HF,Corl_C); // 1 for foating division;
    B = B + 1/TU; // 1 for path_distance^2; // added to wenhaos code the dependence of bandwidth on TU
    //% end of the path gain computation


    // Path delay computation
    //     Path_delay is computed by divide the path distance by speed of light
    //     path_delay = path_distance / c;
    C = C + ComputeTable(Div_HF,Corl_C) / TU; // path_distance2 is available from the path distance computation, only one floating division to get delay.
    M = M + 1; // 1 for speed of light;
    L = L + ComputeTable(Div_HF,Corl_L); // 1 for foating division;
    B = B + 1/TU; // 1 for distance; // added to wenhaos code the dependence of bandwidth on TU
    // end of the path delay computation

    C = C * NF;
    M = M * NF;
    L = L;
    B = B * NF;
  }

  // RCS
  void Band::rcs_cmbl(single_band_ge_stats & fed_cmbl){
    // TU - Update time
    // PT - Number of paths
    // OB - Number of objects
    
    float TU = fed_cmbl.global_fid.upd_rate;
    float PT = fed_cmbl.global_fid.num_path;
    float OB = fed_cmbl.global_fid.num_obj;

    float C = 0; // Compute - in terms of 64bit FP MACs - Total number of MAC units accessed per sec
    float M = 0; // Memory in terms of 64 bit FP Units - Total memory needed for all objects
    float B = 0; // Memory Bandwidth in terms of 64bit FP I/O operations required - Total memory access per second
    float L = 0; // Latency in terms of number of clock cycles required for this step

    //// RCS fidelity
    // Order 0
    if (fed_cmbl.rcs.order == 0){
        C = 0;
        M = 1;
        B = 1;
        L = 0;
    }
    // Order 1
    else if (fed_cmbl.rcs.order == 1){
        std::cout << 'RCS fidelity order not supported.';
        C = 0;
        M = 0;
        B = 0;
        L = 0;
    }

    // Order 2
    else if (fed_cmbl.rcs.order == 2){
        C = PT*(5 + 3*fed_cmbl.rcs.points + 2*ComputeTable(Cos,Corl_L) + 3*ComputeTable(Sin,Corl_C))/TU;
        M = OB*4*fed_cmbl.rcs.points;
        B = PT*4*fed_cmbl.rcs.points/TU;
        L = 0;
    }

    // Order 3
    else if (fed_cmbl.rcs.order == 3){
        C = PT*(5 + 7*fed_cmbl.rcs.points + 2*ComputeTable(Cos,Corl_L) + 3*ComputeTable(Sin,Corl_C))/TU;
        M = OB*7*fed_cmbl.rcs.points;
        B = PT*7*fed_cmbl.rcs.points/TU;
        L = 0;
    }

    // Order 4
    else if (fed_cmbl.rcs.order == 4){
        C = PT*(5 + 12*fed_cmbl.rcs.points + 2*ComputeTable(Cos,Corl_L) + 3*ComputeTable(Sin,Corl_C))/TU;
        M = OB*8*fed_cmbl.rcs.points;
        B = PT*8*fed_cmbl.rcs.points/TU;
        L = 0;
    }

    // Order 5
    else if (fed_cmbl.rcs.order == 5){
        std::cout << 'RCS fidelity order not supported.';
        C = 0;
        M = 0;
        B = 0;
        L = 0;
    }
    // Order 6
    else if (fed_cmbl.rcs.order == 6){
        C = 0;
        M = OB*(360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*(360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*fed_cmbl.rcs.freq*fed_cmbl.rcs.plzn*fed_cmbl.rcs.samples;
        B = PT*fed_cmbl.rcs.samples/TU;
        L = 0;
    }
    else{
        std::cout << 'RCS fidelity order not supported.';
        C = 0;
        M = 0;
        B = 0;
        L = 0;
    }

    fed_cmbl.rcs.set_cmbl(C, M, B, L);
  }

  // Time Update
  void Band::tu_cmbl(single_band_ge_stats & fed_cmbl){
    // Given Fidelites (Change these according to your blocks inputs)
    float TU = fed_cmbl.global_fid.upd_rate;
    float NF = fed_cmbl.global_fid.num_path; // no need to 2X factor
  
    float C = 0; // Compute - in terms of 64bit FP MACs
    float M = 0; // Memory in terms of 64 bit FP Units
    float B = 0; // Memory Bandwidth in terms of 64bit FP I/O operations required
    float L = 0; // Latency in terms of number of clock cycles required for this step



    //% UpdateRate computation
    //     UpdateRate is inversly propotional to path speed for a fixed noise threshold
    //     TU = C / v;
    //     C: some constant including noise threshold
    //     v: path velocity from other block

    C = C + ComputeTable(Div_HF, Corl_C) / TU; // path_distance^2 is available from the path distance computation, only one floating division.
    M = M + 1; // 1 for constant;
    L = L + ComputeTable(Div_HF, Corl_C); // 1 for foating division;
    B = B + 1/TU; // 1 for path velocity
    //% end of the path gain computation

    C = C * NF;
    M = M * NF;
    L = L;
    B = B * NF;

    fed_cmbl.tu.set_cmbl(C, M, B, L);
  }

  void Band::get_ge_cmbl(single_band_ge_stats & fed_cmbl){
    coordinate_trans_cmbl(fed_cmbl);
    nr_engine_cmbl(fed_cmbl);
    relative_orientation_cmbl(fed_cmbl);
    antenna_gain_cmbl(fed_cmbl);
    path_gain_cmbl(fed_cmbl);
    rcs_cmbl(fed_cmbl);
    tu_cmbl(fed_cmbl);
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

  float avg_links_per_mm2=0;

  // Histogram
  int wafer_histo[100000];
  int wafer_ppu_ge_histo[100000];
  void print_wafer_histo() {
    for(int i = 1; i < 400; ++i) {
      printf("%d ",wafer_histo[i]);
    }
  }

  void print_wafer_ge_ppu_histo(){
    for(int i = 1; i < 400; ++i) {
      printf("%d ",wafer_ppu_ge_histo[i]);
    }
  }
};

void evaluate_ge(std::vector<Band>& scene_vec){  
  for(auto & b : scene_vec) {
    b.get_ge_cmbl(b.ge_stat);
  }
}

ge_core * design_ge_core_for_scenario(std::vector<Band>& scene_vec, drbe_wafer& w) {
  ge_core* ge = new ge_core(&*w._t); //needs to be null so we don't delete a non-pointer

  evaluate_ge(scene_vec);

  // TODO: how to design the ge based on the statistic of all scenarios
  // ......

  return ge;
}

bool model_succesful(path_proc_unit* ppu, Band& band, drbe_wafer& w, int max_wafers) {
  int failures=0;
  int ppus = band.ppus_per_band(*ppu,w,failures,true);

  if(failures!=0) return false; // don't include this scenario if failed

  float num_wafers = ppus / w.num_units();

  return ceil(num_wafers) <= max_wafers;
}

std::vector<Band> top_k_pareto_scenarios(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, int k, int max_wafers) {
  std::vector<Band> top_k_scenarios;
  std::vector<Band> successful_pareto_bands;

  top_k_scenarios.resize(2);

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
    Band try_band;
    Band hardest_possible_band=band;

    int num_tries=0;

    while(could_increase && num_tries < 20) {
      try_band=hardest_possible_band;
      
      if(try_band.increase_difficulty()) {
        if(model_succesful(ppu,try_band,w,max_wafers)) {
          hardest_possible_band=try_band;
        } else {
          num_tries++;
        }
      } else {
        could_increase=false;
      }
    }

    //while we do everything, check and see if we found a higher # platforms
    if(top_k_scenarios[0].platforms() < hardest_possible_band.platforms()) {
      top_k_scenarios[0] = hardest_possible_band;
    }

    //while we do everything, check and see if we found a higher # platforms
    if(top_k_scenarios[1].reflectors() < hardest_possible_band.reflectors()) {
      top_k_scenarios[1] = hardest_possible_band;
    }

    //get rid of all the other bands we previously included but aren't necessary
    successful_pareto_bands.erase(
        std::remove_if(successful_pareto_bands.begin(), 
                       successful_pareto_bands.end(), 
                       [&hardest_possible_band](Band& other) { 
                       return !other.could_be_harder_than(hardest_possible_band); }), 
                       successful_pareto_bands.end());

    successful_pareto_bands.push_back(hardest_possible_band);
  }

  //for(auto& band : successful_pareto_bands) {
  //  band.print_csv();
  //  printf("\n");
  //}
  //printf("num succ: %d \n", (int)successful_pareto_bands.size());

  if(top_k_scenarios[0].reflectors() > ScenarioGen::max_platforms() * .9) {
    top_k_scenarios.resize(1);
  }

  //get the most different bands
  for(int i = top_k_scenarios.size(); i < k; ++i) {
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
  float wafers=0;
  int total_in_target_wafer=0;

  float total_links = 0;
  float total_ppus = 0;
  float total_coef = 0;

  for(auto & b : scene_vec) {

    float wafer_unconstrained_ppus = b.ppus_per_band(*ppu,w,stats.total_failures,verbose);
    total_ppus += wafer_unconstrained_ppus;

    float links_per_wafer = b.calc_links_per_wafer(wafer_unconstrained_ppus,ppu,w,verbose);
    float num_wafers = b.num_links() / links_per_wafer;

    assert(num_wafers>0);

    float wafer_constrained_ppus = num_wafers * w.num_units() * ppu->_ppus_per_chiplet;


    // For status purposes, I need to use algorithmic links to make a fair comparison
    // between channel models

    total_coef += b._n_obj * b._avg_coef_per_object * b.algorithmic_num_links();
    total_links += b.algorithmic_num_links();


    //printf("%d chiplets (PPU) are needed for this scenario\n", int_num_ppus);
    b.num_chiplet    += wafer_constrained_ppus / ppu->_ppus_per_chiplet;
    b.num_ppu_chiplet = wafer_constrained_ppus / ppu->_ppus_per_chiplet; //FIXME: Why two of these?
    //printf("after count the chiplet from PPU, there are %d chiplets\n", b.num_chiplet);


    int int_num_wafers = ceil(num_wafers);
    if(int_num_wafers<=num_wafers_target) total_in_target_wafer++;
    wafers += int_num_wafers;
    int max_histo_elem=sizeof(stats.wafer_histo)/sizeof(stats.wafer_histo[0]);
    int_num_wafers=max(0,min(max_histo_elem-1,int_num_wafers));
    stats.wafer_histo[int_num_wafers]++;
    
  }

  stats.percent_in_target_wafer = total_in_target_wafer / (float) scene_vec.size();
  stats.avg_links_per_mm2 = total_links / (total_ppus * ppu->area());
  stats.avg_wafers = wafers / scene_vec.size();
  stats.avg_coef_per_ppu = total_coef / total_ppus;
  stats.avg_coef_per_mm2 = total_coef / (total_ppus * ppu->area());
}

path_proc_unit* design_ppu_for_scenarios(std::vector<Band>& scene_vec, drbe_wafer& w,
                                         bool dynamic_reconfig) {
  PPUStats best_stats;
  path_proc_unit* best_ppu = NULL; //needs to be null so we don't delete a non-pointer
  float chiplet_area = w.chiplet_area();
  tech_params& t = *w._t;

  bool seen_aidp=false;
  // Preprocess scene vec to narrow params
  for(Band& b : scene_vec) {
    seen_aidp |= b._is_aidp;
  }

  int max_dielets_per_side=1;
  if(seen_aidp) {
    max_dielets_per_side=2;
  } 

  for(int dielets_per_side = 1; dielets_per_side <= max_dielets_per_side; ++dielets_per_side) {
    float ppus_per_chiplet = dielets_per_side * dielets_per_side;
    float ppu_area = chiplet_area / ppus_per_chiplet;

    // AggNet computaitons
    for(int agg_network = 1; agg_network < 80; agg_network+=8) {
      float side_length=sqrt((float)ppu_area);
      //for now, round down to 4.4 (round to nearst tenth mm), so that numbers match
      side_length = ((float)((int)(side_length*10)))/10.0;

      //side_length -= 1; //leave one mm out
      int bits_per_side=t._chiplet_io_bits_per_mm2*side_length;    
      int total_ins_per_side = (int)(bits_per_side/32.0); //divide by 32 bit width

      int remain=total_ins_per_side - (2*2 + agg_network*2);
      if(remain < 2) break;


      // Loop over coefficients per cluster
      for(int cpc = 20; cpc < 21; cpc+=5) { //FIXME: turn off for speed


        for(int clutter_index = 0; clutter_index <= 100; clutter_index+=10) {
          int clutter_ratio=clutter_index;
          if(clutter_index==0) {
            clutter_ratio=1;
          }

          for(int mem_ratio =4; mem_ratio < 60; mem_ratio+=2) {
            path_proc_unit* ppu =new path_proc_unit(&t);

            ppu->_is_dyn_reconfig=dynamic_reconfig;
            ppu->_coef_per_cluster=cpc;
            ppu->_input_router._in_degree=2;
            ppu->_input_router._out_degree=2;
            ppu->_output_router._in_degree=agg_network;
            ppu->_output_router._out_degree=agg_network;
            ppu->_coef_router._in_degree=remain/2;
            ppu->_coef_router._out_degree=remain/2;
            ppu->_ppus_per_chiplet=ppus_per_chiplet;

            ppu->set_params_by_mem_ratio(mem_ratio,clutter_ratio,ppu_area);
            if(ppu->_num_clusters==0) break;
            //ppu_vec[i].print_area_breakdown();
 
            PPUStats stats; 
            evaluate_ppu(ppu,scene_vec,w,stats,1 /*target wafers*/);

            //Optimize the average coefficient per mm2
            //if((stats.percent_in_target_wafer > best_stats.percent_in_target_wafer) 
            //    && stats.total_failures==0) {
            if((stats.avg_coef_per_mm2 > best_stats.avg_coef_per_mm2) && stats.total_failures==0) {
              path_proc_unit* old_best = best_ppu;
              best_stats = stats;
              best_ppu=ppu;
              if(old_best) delete old_best;
            } else {
              delete ppu;
            }

            //printf("Mem Ratio: %d, Clusters/PPU: %d, Coef/PPU/mm2: %f, coverage: %f\n",
            //    mem_ratio,ppu->_num_clusters,
            //    best_stats.avg_coef_per_mm2,
            //    best_stats.percent_in_target_wafer);
          }
        }
      }
    }
  }

  assert(best_ppu);
  return best_ppu;
}



int main(int argc, char* argv[]) {
  bool verbose = false;
  bool direct_path = false;
  bool print_pareto = false;
  bool dynamic_reconfig = false;
  bool verbose_ge_info = false;
  bool aidp = false;
  bool challenge_scenario=false;
  int num_scenarios=1000;
  bool limit_wafer_io=false;
  tech_params t;

  int opt;
  while ((opt = getopt_long(argc, argv, "vgkdprcn:w", long_options, nullptr)) != -1) {
    switch (opt) {
      case 'd': direct_path=true; break;
      case 'c': challenge_scenario=true; break;
      case 'k': aidp=true; break;
      case 'r': dynamic_reconfig=true; break;
      case 'v': verbose = true; break;
      case 'g': verbose_ge_info = true; break;
      case 'p': print_pareto = true; break;
      case 'n': num_scenarios = atoi(optarg); break;
      case 'w': limit_wafer_io = true; break;
      default: exit(1);
    }
  }

  argc -= optind;
  argv += optind;

  if(verbose) {
    printf("Verbose Mode On\n");
  }

  if (argc != 0) {
    cerr << "Usage: drbe_model_hpc [FLAGS]\n";
    exit(1);
  }

  cout << "Assuming ";
  if(direct_path) cout << "direct path";
  else if(aidp) cout << "aidp";
  else cout << "tapped delay";
  cout << " channel model\n";


  std::vector<Band> scene_vec;
  ScenarioGen::gen_scenarios(scene_vec,direct_path,aidp,
                             challenge_scenario, num_scenarios);

  // You can run different experiments by uncommenting different loops here ... we should make this
  // more configurable in the future. : )
  

  //for(int ppu_area = 20; ppu_area < 21; ppu_area+=5) {
  int ppu_area=20;


  drbe_wafer w(&t,300,(float)ppu_area);
  w.set_limit_wafer_io(limit_wafer_io);

  //float fast_update_period=10000;
  
  //for(fast_update_period = 1000; fast_update_period < 1000000; 
  //    fast_update_period*=1.2589254117941672104239541063958) {
  
  //for(float frac_clutter = 0; frac_clutter <= 100; frac_clutter += 20) {

  //  for(auto& b : scene_vec) {
  //    //b._high_update_period=fast_update_period;
  //    b._frac_clutter = frac_clutter/100.0f;
  //  }

  //printf("\nTechnology Density Experiment: Increase the Density (factor \"v\" below)\n");
  float old = w.wafer_io();//t.area_multiplier();
  float v_range = 100;
  float factor = 1.0905077326652576;
  //factor = factor * factor;
  factor = factor * factor * factor * factor;
  for(float v = old; v < old*v_range+0.01; v*=factor) {
    //t.set_area_multiplier(v);
    w.set_wafer_io(v);
  


  //for(int spmm = 1; spmm < 200; ++spmm) {
    //t._sram_Mb_per_mm2=spmm;

    //for(int agg_network = 1; agg_network < 200; agg_network+=1) {
    //agg_network=11

    path_proc_unit* best_ppu = design_ppu_for_scenarios(scene_vec,w,dynamic_reconfig);

    PPUStats stats;
    int num_wafers_target=1;


    // clear the number of chiplet of each band from previous evaluation
    for(auto & b : scene_vec){
      b.num_chiplet = 0;
    }

    evaluate_ppu(best_ppu,scene_vec,w,stats,num_wafers_target,verbose /*verbose*/);

    if(print_pareto) {
      top_k_pareto_scenarios(best_ppu,scene_vec,w,6,num_wafers_target);
    }

    float average_clutter=0;
    for(auto & b : scene_vec){
      average_clutter+=b._frac_clutter;
    }
    average_clutter/=scene_vec.size();



    //printf("Fast Update Period: %0.0fus, ", fast_update_period/1000);
    
    printf("v: %0.3f, ", v);
    
    printf("avg_clut: %f, "\
           "%dmm^2 PPU (%0.2f), in-MB: %0.2f, clust: %d, flex_clust: %d, coef/clust %d, "\
           "In: %d/%d Agg: %d/%d, Coef: %d/%d, Mem Ratio: %d, "\
           "ppus/die: %d, " \
           "Coef/mm2: %0.2f, links/cm2: %0.2f, "\
           "avg_waf: %0.2f, links/ppu: %0.2f, per_targ_waf: %0.1f, fail: %d\n",
        average_clutter,
        ppu_area, 
        best_ppu->area(), 
        best_ppu->input_buf_area()*t.sram_Mb_per_mm2()/8,
        best_ppu->_num_clusters, 
        best_ppu->_num_flexible_clusters,
        best_ppu->_coef_per_cluster,
        best_ppu->_input_router._in_degree, 
        best_ppu->_input_router._out_degree, 
        best_ppu->_output_router._in_degree, 
        best_ppu->_output_router._out_degree, 
        best_ppu->_coef_router._in_degree, 
        best_ppu->_coef_router._out_degree, 
        (int)best_ppu->_mem_ratio, 
        (int)best_ppu->_ppus_per_chiplet,
        stats.avg_coef_per_mm2,
        stats.avg_links_per_mm2 * 100,
        stats.avg_wafers,
        stats.avg_links_per_mm2 * best_ppu->area(),
        stats.percent_in_target_wafer*100,
        stats.total_failures);

    //printf("wafer_sram_MB: %0.2f, wafer_macc/s: %0.2f\n", 
    //    2025 * (best_ppu->input_buf_area() + best_ppu->_input_tap_fifo.area() 
    //           + best_ppu->_coef_storage.area())*t.sram_Mb_per_mm2()/8,
    //    2025 * best_ppu->_num_clusters*best_ppu->_coef_per_cluster*t._macc_per_complex_op/1024.0);


    //printf("Wafers Needed Histogram: ");
    //stats.print_wafer_histo();
    //printf("\n");


    //printf("PPU Wafer Needed (No GE) Histogram: ");
    //stats.print_wafer_histo();
    //printf("\n");
    
    //best_ppu->print_params();
    //best_ppu->print_area_breakdown();

  } PPUStats stats;


  // -------------------------- Design the GE core ---------------------------

  // Define the Fidelity
  for(auto & b : scene_vec){
    // Global
    b.ge_stat.global_fid.upd_rate = w._t->ge_freq();
    b.ge_stat.global_fid.num_obj = b._n_obj;
    b.ge_stat.global_fid.num_path = b.num_paths();
    // Coordinate Transformation
    b.ge_stat.coordinate_trans.ta1_scene_upd_rate = 1e-3;
    // NR Engine
    b.ge_stat.nr_engine.interpolation_ord = 6;
    b.ge_stat.nr_engine.conv_fed = 2;
    // Antenna
    b.ge_stat.antenna_gain.order = 5;
    b.ge_stat.antenna_gain.num_antenna = 16;
    b.ge_stat.antenna_gain.res_angle = 1;
    b.ge_stat.antenna_gain.dict_dim = 800;
    // RCS
    b.ge_stat.rcs.order = 6;
    b.ge_stat.rcs.points = 20;
    b.ge_stat.rcs.angle = 2;
    b.ge_stat.rcs.freq = 1;
    b.ge_stat.rcs.plzn = 4;
    b.ge_stat.rcs.samples = 10;
  }

  // ------- Design the Geometry Engine ------
  ge_core * ge = design_ge_core_for_scenario(scene_vec, w);

  // Print out all tradeoff of all scenario
  std::ofstream ge_tradeoff;
  ge_tradeoff.open("ge_tradeoff.csv");
  printf("Start to write GE tradeoff to `ge_tradeoff.csv`\n");
  // Print the header
  ge_tradeoff << "update-rate, number-object, number-path, "
              << "coor-trans-compute, coor-trans-memory, coor-trans-bandwidth, coor-trans-latency, "
              << "nr-engine-compute, nr-engine-memory, nr-engine-bandwidth, nr-engine-latency, "
              << "relative-orientation-compute, relative-orientation-memory, relative-orientation-bandwidth, relative-orientation-latency, "
              << "antenna-compute, antenna-memory, antenna-bandwidth, antenna-latency, "
              << "path-gain-compute, path-gain-memory, path-gain-bandwidth, path-gain-latency, "
              << "rcs-compute, rcs-memory, rcs-bandwidth, rcs-latency, "
              << "tu-compute, tu-memory, tu-bandwidth, tu-latency\n";
  for(auto & b : scene_vec){
    ge_tradeoff << b.ge_stat.global_fid.upd_rate << ", "
                << b.ge_stat.global_fid.num_obj << ", "
                << b.ge_stat.global_fid.num_path << ", "
                // Coordinate
                << b.ge_stat.coordinate_trans.compute << ", "
                << b.ge_stat.coordinate_trans.memory << ", "
                << b.ge_stat.coordinate_trans.bandwidth << ", "
                << b.ge_stat.coordinate_trans.latency << ", "
                // NR Engine
                << b.ge_stat.nr_engine.compute << ", "
                << b.ge_stat.nr_engine.memory << ", "
                << b.ge_stat.nr_engine.bandwidth << ", "
                << b.ge_stat.nr_engine.latency << ", "
                // Relative Orientation
                << b.ge_stat.relative_orientation.compute << ", "
                << b.ge_stat.relative_orientation.memory << ", "
                << b.ge_stat.relative_orientation.bandwidth << ", "
                << b.ge_stat.relative_orientation.latency << ", "
                // Antenna
                << b.ge_stat.antenna_gain.compute << ", "
                << b.ge_stat.antenna_gain.memory << ", "
                << b.ge_stat.antenna_gain.bandwidth << ", "
                << b.ge_stat.antenna_gain.latency << ", "
                // Path gain and velocity
                << b.ge_stat.path_velocity.compute << ", "
                << b.ge_stat.path_velocity.memory << ", "
                << b.ge_stat.path_velocity.bandwidth << ", "
                << b.ge_stat.path_velocity.latency << ", "
                // RCS
                << b.ge_stat.rcs.compute << ", "
                << b.ge_stat.rcs.memory << ", "
                << b.ge_stat.rcs.bandwidth << ", "
                << b.ge_stat.rcs.latency << ", "  
                // Tu
                << b.ge_stat.tu.compute << ", "
                << b.ge_stat.tu.memory << ", "
                << b.ge_stat.tu.bandwidth << ", "
                << b.ge_stat.tu.latency << "\n";
  }
  ge_tradeoff.close();
  // ------- PPU vs. GE DSE ------

  // ------- PPU vs. GE + GM DSE ------

  return 0;
}

