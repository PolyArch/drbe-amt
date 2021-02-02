#include "drbe.h"
using namespace std;

enum flags {
  ANALYSIS_PARETO=128,
  ANALYSIS_WAFER_SCALING,
  ANALYSIS_HW_CONFIG,
  ANALYSIS_PARETO_AFTER_GE,
  SENSITIVITY_TECH_SCALING,
  SENSITIVITY_WAFER_IO,
  SENSITIVITY_TX_SPARSITY,
  PARAM_PULSED_TX_DUTY,
  PARAM_DIELET_AREA,
  PARAM_TECH_SCALING,
  PARAM_DIELETS_PER_WAFER,
  GE_CPU,
  GE_ASIC,
  GE_CGRA,
  GE_HARD
};

static struct option long_options[] = {
    {"direct-path",          no_argument,       nullptr, 'd',},
    {"kp-aidp",              no_argument,       nullptr, 'k',},
    {"dynamic-reconfig",     no_argument,       nullptr, 'r',},
    {"challenge-scenario",   no_argument,       nullptr, 'c',},
    {"num-scenarios",        required_argument, nullptr, 'n',},
    {"target-wafers",        required_argument, nullptr, 't',},
    {"final-target-wafers",  required_argument, nullptr, 'u',},
    {"easy-scenario",        no_argument,       nullptr, 'e',},
    {"very-easy-scenario",   no_argument,       nullptr, 'z',},
    {"verbose",              no_argument,       nullptr, 'v',},
    {"verbose-ge-info",      no_argument,       nullptr, 'g',},
    {"wafer-io",             no_argument,       nullptr, 'w',},
    {"layer",                required_argument, nullptr, 'l',},
    {"dump-file-name",       required_argument, nullptr, 'f',},


    {"pre-ge-summary",       no_argument,       nullptr, 'p',},

    {"pulsed-duty-cycle",    required_argument, nullptr, PARAM_PULSED_TX_DUTY},

    {"dielet-area",          required_argument, nullptr, PARAM_DIELET_AREA},
    {"tech-scale",           required_argument, nullptr, PARAM_TECH_SCALING},
    {"dielets-per-wafer",    required_argument, nullptr, PARAM_DIELETS_PER_WAFER},

    {"ge-cpu",               no_argument,       nullptr, GE_CPU},
    {"ge-asic",              no_argument,       nullptr, GE_ASIC},
    {"ge-cgra",              no_argument,       nullptr, GE_CGRA},
    {"ge-hard",              no_argument,       nullptr, GE_HARD},

    {"sense-tech-scaling",   no_argument,       nullptr, SENSITIVITY_TECH_SCALING},
    {"sense-wafer-io",       no_argument,       nullptr, SENSITIVITY_WAFER_IO},
    {"sense-tx-sparsity",    no_argument,       nullptr, SENSITIVITY_TX_SPARSITY},

    {"print-pareto",         no_argument,       nullptr, ANALYSIS_PARETO},
    {"print-pareto-after-ge",no_argument,       nullptr, ANALYSIS_PARETO_AFTER_GE},
    {"print-wafer-scaling",  no_argument,       nullptr, ANALYSIS_WAFER_SCALING},
    {"print-hw-config",      no_argument,       nullptr, ANALYSIS_HW_CONFIG},
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

int Band::increase_difficulty(int kind=-1) {

  int max_types=10;
  int which = rand_bt(0,max_types);

  if(kind!=-1) {
    which=kind;
    max_types=1;
  }

  assert(num_links() <= ScenarioGen::max_links());

  //cout << "start\n";
  for(int i = 0; i < max_types; ++i,which=(which+1)%max_types) {
    //cout << "which:" << which <<"\n";
    switch(which) {
      case 0: { //add a fixed platform
        //if(platforms() + 1 >ScenarioGen::max_platforms()) break;
        int change_amount = rand_bt(1,5);
        _n_fixed+=change_amount; //one object per band
        recalculate_txrx();

        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _n_fixed -= change_amount;
          recalculate_txrx();
        } 

      }
      case 1: { //upgrade a fixed platform to slow platform
        if(_n_fixed == 0) break;
        int change_amount = std::min(_n_fixed,20);
        _n_fixed -= change_amount;
        _n_slow  += change_amount;
        recalculate_txrx();
        return which+1;
      }
      case 2: { //upgrade a slow platform to fast platform
        if(_n_slow == 0) break;
        int change_amount = std::min(_n_slow,20);
        _n_slow -= change_amount;
        _n_fast += change_amount;
        recalculate_txrx();
        return which+1;
      }
      case 3: { //subtract a band
        if(_n_bands==1 || rand_bt(0,2)!=0) break; //rand check to make it less likely

        _n_bands--;
        recalculate_txrx();

        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _n_bands++;
          recalculate_txrx();
        } 
      }
      case 4: { //avg_coef_per_obj
        break;
        if(_avg_coef_per_object + OBJ_SIZE_GRAN > ScenarioGen::max_coef_per_obj()) break;
        _avg_coef_per_object+=OBJ_SIZE_GRAN;
        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _avg_coef_per_object-=1;
        } 
      }
      case 5: { //range
        if(_range + 1 > ScenarioGen::max_range()) break;
        _range+=1;
        return which+1;
      }
      case 6: { // more objects
        if(_n_obj + 1 > ScenarioGen::max_objects()) break;
        _n_obj += 1;
        _n_full_range_obj += 1;

        //Remove full objects until we're within the thresh/objective
        while(!ScenarioGen::scenario_is_between_threshold_and_objective(*this) &&
               _avg_frac_full_objects>0.001) {
          _avg_frac_full_objects -= 0.001;
        }

        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          //Failed to add a new object b/c link complexity is too high
          _n_obj -=1;
          _n_full_range_obj -= 1;
        } 
      }
      case 7: { // more clutter
        if(_frac_clutter + 0.01 > ScenarioGen::max_clutter()) break;
        _frac_clutter += 0.01;
        return which+1;
      }
      case 8: { // more link complexity
        _avg_frac_full_objects += 0.02;
        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _avg_frac_full_objects -= 0.02;
        } 
        assert(num_links() <= ScenarioGen::max_links());
      }
      case 9: { //increase the number of full range objects
        if(_n_full_range_obj < _n_obj) {
          _n_full_range_obj += 1;
          return which+1;
        }
      }
    }
  }
  return false;
}

float Band::frac_max_links() {
  float val = algorithmic_num_links() / (float)ScenarioGen::max_links();
  float rounded_val = ceil(val*100)/100.0f;
  return rounded_val;
}
float Band::frac_max_link_complexity() {
  float val = link_complexity() / (float)ScenarioGen::overprov_link_complexity();
  float rounded_val = ceil(val*100)/100.0f;
  return std::min(1.0f,rounded_val);
}

void Band::print_norm_csv() {
  static bool printed = false;
  if(!printed) {
    printf("%-8s %-8s %-8s %-8s %-8s %-8s %-8s \n", 
           "Links",   "Objects",   "Tap/Link",   
           "Fast", "RngOvrProv", "%%Clutter", "Range");
    printed=true;
  }

  printf("%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f", 
    frac_max_links(),
    _n_obj/(float)ScenarioGen::max_objects(),
    frac_max_link_complexity(),
    _n_fast / (float)platforms(),
    frac_full_range(),
    _frac_clutter / (float)ScenarioGen::max_clutter(),
    _range / (float)ScenarioGen::max_range()
  );
//  printf("%d %d ", _n_obj, _n_full_range_obj);
}

void Band::print_csv() {
  printf("%3d, %d, %0.3f, %0.3f, %0.3f, %0.3f, %3d", 
    num_links(),_n_obj,link_complexity(),_n_fast/(float)ScenarioGen::max_platforms(),
    frac_full_range(),_frac_clutter, _range);
}

std::vector<float>& Band::normalized_vec() {
  if(_norm_features.size()==0) {
    _norm_features.push_back(frac_max_links());
    _norm_features.push_back(_n_obj/(float)ScenarioGen::max_objects());
    _norm_features.push_back(frac_max_link_complexity());
    _norm_features.push_back(_n_fast / (float)platforms());
    _norm_features.push_back(frac_full_range());
    _norm_features.push_back(_frac_clutter / (float)ScenarioGen::max_clutter());
    _norm_features.push_back(_range / (float)ScenarioGen::max_range());
  }
  return _norm_features;
}

  // ---------------- Geometry Engine -----------------
  // Coordinate Translation
  void coordinate_trans_cmbl(ge_stat_per_band & fed_cmbl){
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
  void nr_engine_cmbl(ge_stat_per_band & fed_cmbl){
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
  void relative_orientation_cmbl(Band & b, ge_stat_per_band & fed_cmbl){
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
void antenna_gain_cmbl(Band & b, ge_stat_per_band & fed_cmbl){
  float no_antenna=fed_cmbl.antenna_gain.num_antenna;
  float o=fed_cmbl.antenna_gain.order;
  float Tu=fed_cmbl.global_fid.upd_rate;
  float K=fed_cmbl.antenna_gain.res_angle;
  float no_paths=fed_cmbl.global_fid.num_path;
  float no_Tx=b._n_tx * b._n_bands; //_n_tx and _n_tx
  float no_Rx=b._n_rx * b._n_bands;
  float dict_size=fed_cmbl.antenna_gain.dict_dim;
  int general = fed_cmbl.antenna_gain.general;
  
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
      C=no_paths*(2*C/Tu);
      M=M+(3*(no_Tx + no_Rx)+1);
      B=B+((3*(no_paths*2)+1)/Tu); 
      L=21;
  }else if(0==3){
        C=C+2*(8)+1;
        C=no_paths*(2*C/Tu);
        M=M+(4*(no_Tx + no_Rx));
        B=B+((4*(no_paths*2))/Tu); //getting each antenna's gains and beamwidths per update time
        L=21;
  }else if(o==4){

        float no_bins=180/K;
        float bins_for_ifft=pow(2, ceil(log2(no_bins)));

        C=C+1+ComputeTable(mod, Corl_C/*Not used*/)+2+no_antenna+2+no_antenna+2+no_antenna+no_antenna+2*no_antenna; 
          // to convert to AWV domain

        C=C+6*no_antenna; // complex point-wise multiplication
        C=C+6*(no_antenna-1)+(bins_for_ifft/2)*log2(bins_for_ifft)+bins_for_ifft; // to convert to angular
        C=C+ComputeTable(division, Corl_C/*Not used*/)+(bins_for_ifft+1); // to normalize gain
        C=C+ComputeTable(division, Corl_C/*Not used*/);
        
        C=2*C; // For elevation dimension
        C=C+6; // Multiplication of elevation and azimuth gains;
        C=no_paths*(2*C/Tu);
        int no_inputs=4;
        M=M+no_inputs*(no_Tx + no_Rx)+no_antenna; // window precompute and store.
        B=B+((no_inputs*no_paths*2)+(no_antenna*no_paths*2))/Tu;

  }else if(o==5 && general == 0){
        int To = 3;
        float N=floor( pow((180/K)*(90/K), 2) ); // needs to be integer
        float n=floor((180/K)); // needs to be integer
        float x=floor(90/K);
        float y=floor(180/K);
        M=(no_Tx + no_Rx)*(n*dict_size+2*To*N); //n*dict_size+2*T0*N
        C=C+2*(((2*x)+(2*x))+2*((2*y)+(2*y))+2*((x)-1)+2*((y)-1)+2*To-1+6); //11 (dependent on To)
        C=no_paths*(C/Tu);
        B=B+(To*(no_paths*2)/Tu); //getting To numbers from memory for one path per update time;
  }else if(o==5 && general == 1){
        int N=floor((180/K)* pow((90/K), 2) ); // needs to be integer
        int n=floor((180/K)); // needs to be integer
        M=(no_Tx + no_Rx)*(n*dict_size+2*3*N); //n*dict_size+2*T0*N
        C=C+2*((4*ComputeTable(Div, Corl_C)+4)+(2*3-1+6));
        C=no_paths*(C/Tu);
        B=B+(3*(no_paths*2)/Tu); //getting To numbers from memory for one path per update time;
    
  }

  fed_cmbl.antenna_gain.set_cmbl(C, M, B, L);
}
// Path Gain and Velocity
void path_gain_cmbl(Band & b, ge_stat_per_band & fed_cmbl){
  // Given Fidelites (Change these according to your blocks inputs)
  float SU = fed_cmbl.coordinate_trans.ta1_scene_upd_rate;
  float TU = fed_cmbl.global_fid.upd_rate;
  float NF = b.num_paths() * 2;
  
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
//  C = C + (3 / TU) + 1 / TU + (3 / TU) + ComputeTable(Sqrt,Corl_C); // 5 * TU to get the relative_location at each update(including scenarioUpdate),  3 to get sum of square, 
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
  fed_cmbl.path_velocity.set_cmbl(C, M, B, L);
}

// RCS
void rcs_cmbl(ge_stat_per_band & fed_cmbl){
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
      M = OB*1;
      B = (PT/TU)*1;
      L = 20;
  }
  // Order 1
  else if (fed_cmbl.rcs.order == 1){
      C = (PT/TU)*(7 + 2*ComputeTable(Cos,Corl_C) + 3*ComputeTable(Sin,Corl_C) + 1*ComputeTable(Cosi,Corl_C));
      M = OB*250e3;
      B = (PT/TU)*3;
      L = 20 + ComputeTable(Cos,Corl_L) + ComputeTable(Cosi,Corl_L) + ceil(log2(7));
  }

  // Order 2
  else if (fed_cmbl.rcs.order == 2){
      C = (PT/TU)*(4*fed_cmbl.rcs.points + 7 + 4*ComputeTable(Cos,Corl_C) + 6*ComputeTable(Sin,Corl_C));
      M = OB*4*fed_cmbl.rcs.points;
      B = (PT/TU)*4*fed_cmbl.rcs.points;
      L = 20 + ComputeTable(Cos,Corl_L) + ceil(log2(4*fed_cmbl.rcs.points + 7));
  }

  // Order 3
  else if (fed_cmbl.rcs.order == 3){
      C = (PT/TU)*(8*fed_cmbl.rcs.points + 7 + 4*ComputeTable(Cos,Corl_C) + 6*ComputeTable(Sin, Corl_C));
        M = OB*7*fed_cmbl.rcs.points;
        B = (PT/TU)*7*fed_cmbl.rcs.points;
        L = 20 + ComputeTable(Cos,Corl_L) + ceil(log2(7 + 7*fed_cmbl.rcs.points));
  }

  // Order 4
  else if (fed_cmbl.rcs.order == 4){
        C = (PT/TU)*((4 + 4*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points + 7 + 4*ComputeTable(Cos,Corl_C) + 6*ComputeTable(Sin,Corl_C));
        M = OB*(5 + 2*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points;
        B = (PT/TU)*(5 + 2*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points;
        L = 20 + ComputeTable(Cos,Corl_L) + ceil(log2(7 + (4 + 4*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points));
  }

  // Order 5
  else if (fed_cmbl.rcs.order == 5){
        C = (PT/TU)*(((4 + 4*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points + 7 + 4*ComputeTable(Cos,Corl_C) + 6*ComputeTable(Sin,Corl_C))
            + (15*fed_cmbl.rcs.plates + 4 + 4*ComputeTable(Cos,Corl_C) + 6*ComputeTable(Sin,Corl_C)));
        M = OB*((5 + 2*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points + 8*fed_cmbl.rcs.plates);
        B = (PT/TU)*((5 + 2*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points + 8*fed_cmbl.rcs.plates);
        L = 20 + ComputeTable(Cos,Corl_L) + ceil(log2(((4 + 4*fed_cmbl.rcs.pntAngle)*fed_cmbl.rcs.points)+7) + (15*fed_cmbl.rcs.plates + 4));
  }
  // Order 6
  else if (fed_cmbl.rcs.order == 6){
        C = 0.1; //Dummy
        M = OB * 
              (360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*
              (360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*
              fed_cmbl.rcs.freq*fed_cmbl.rcs.plzn*fed_cmbl.rcs.points*2;
        B = (PT/TU)*fed_cmbl.rcs.points*2;
        L = 20;
  }
  else{
      std::cout << "RCS fidelity order not supported.";
      C = 0;
      M = 0;
      B = 0;
      L = 0;
  }

  fed_cmbl.rcs.set_cmbl(C, M, B, L);
}

// Time Update
void tu_cmbl(ge_stat_per_band & fed_cmbl){
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

float Band::average_clutter(std::vector<Band> scene_vec) {
  float average_clutter=0;
  for(auto & b : scene_vec){
    average_clutter+=b._frac_clutter;
  }
  average_clutter/=scene_vec.size();
  return average_clutter;
}

void get_ge_cmbl(Band & b, ge_stat_per_band & fed_cmbl, drbe_wafer& w, bool ge_cpu){
  // All of the calculation per band
  float utilization_fp_macc = 1;
  if(ge_cpu) utilization_fp_macc = 0.1;

  coordinate_trans_cmbl(fed_cmbl);
  nr_engine_cmbl(fed_cmbl);
  relative_orientation_cmbl(b, fed_cmbl);
  antenna_gain_cmbl(b, fed_cmbl);
  path_gain_cmbl(b, fed_cmbl);
  rcs_cmbl(fed_cmbl);
  tu_cmbl(fed_cmbl);

  // Calculate the area for compute
  // Total 64-bit MACC operation need per second
  // Unit is 64-bit MACC / sec
  float total_compute_per_sec = (fed_cmbl.coordinate_trans.compute + 
                        fed_cmbl.nr_engine.compute +
                        fed_cmbl.relative_orientation.compute +
                        fed_cmbl.antenna_gain.compute +
                        fed_cmbl.path_velocity.compute +
                        fed_cmbl.rcs.compute +
                        fed_cmbl.tu.compute)  / utilization_fp_macc;
  // Total compute density can provide per second
  // Unit is 64-bit MACC / mm2 / second
  float compute_density_per_sec = w._t->ge_comp_density() * w._t->ge_freq();

  // Unit is 64-bit
  float total_mem = fed_cmbl.coordinate_trans.memory + 
                    fed_cmbl.nr_engine.memory +
                    fed_cmbl.relative_orientation.memory +
                    fed_cmbl.antenna_gain.memory +
                    fed_cmbl.path_velocity.memory +
                    fed_cmbl.rcs.memory +
                    fed_cmbl.tu.memory;
  float mem_density = w._t->ge_dram_MB() * 1024 * 1024 * 8 / 64;

  // Record the compute and memory requirement
  fed_cmbl.total_compute = total_compute_per_sec;
  fed_cmbl.total_memory = total_mem;

  // Record the Compute and Memory Area
  fed_cmbl.compute_area = total_compute_per_sec / compute_density_per_sec;
  fed_cmbl.mem_area = total_mem / mem_density;

  // Record the total area
  fed_cmbl.total_area = fed_cmbl.compute_area + fed_cmbl.mem_area;
}

int get_best_ppu() {
  return 0;
}

void dump_ge_tradeoff(ofstream & ge_tradeoff, GEStats & ge_stats, std::vector<Band>& scene_vec, WaferStats & w_stats){
  for(unsigned i = 0; i < scene_vec.size(); i++){
    ge_tradeoff << ge_stats.print_ge_tradeoff(scene_vec[i],ge_stats.ge_stat_vec[i], w_stats);
  }
}

void evaluate_ge(std::vector<Band>& scene_vec, GEStats & ge_stats, drbe_wafer& w, bool ge_cpu){
  int i = 0;
  for(auto & b : scene_vec) {
    get_ge_cmbl(b, ge_stats.ge_stat_vec[i++], w, ge_cpu);
  }
}

void print_ge_breakdown(ge_core & ge){
  printf(" ------ ASIC Compute Area ------\n");
  printf("\
          coord_trans = %.6f um2\n \
          NR engine = %.2f um2\n \
          relative orientation = %.6f um2\n \
          antenna = %.2f um2\n \
          path gain = %.2f um2\n \
          RCS = %.2f um2\n \
          Tu = %.2f um2\n",
          ge.coord_trans_compute_area,
          ge.nr_engine_compute_area,
          ge.relative_orient_compute_area,
          ge.antenna_compute_area,
          ge.path_gain_compute_area,
          ge.rcs_compute_area,
          ge.tu_compute_area);
  printf(" ------ ASIC Memory Area ------\n");
    printf("\
          coord_trans = %.6f um2\n \
          NR engine = %.2f um2\n \
          relative orientation = %.6f um2\n \
          antenna = %.2f um2\n \
          path gain = %.2f um2\n \
          RCS = %.2f um2\n \
          Tu = %.2f um2\n",
          ge.coord_trans_memory_area,
          ge.nr_engine_memory_area,
          ge.relative_orient_memory_area,
          ge.antenna_memory_area,
          ge.path_gain_memory_area,
          ge.rcs_memory_area,
          ge.tu_memory_area);
}

int get_num_supported_scenario(
    ge_core & ge, 
    drbe_wafer& w, 
    int num_scenarios,
    int num_ppu_chiplet,
    GEStats & ge_stats,
    PPUStats & ppu_stats){
      float compute_density_per_sec = w._t->ge_comp_density() * w._t->ge_freq();
      float mem_density = w._t->ge_dram_MB() * 1024 * 1024 * 8 / 64;
      int num_supported_scenario = 0;
      for(int band_idx = 0; band_idx < num_scenarios;band_idx ++){
        // memory and compute resource need satisfied
        bool can_support_this =
        /* PPU */ 
            ppu_stats.ppu_stat_vec[band_idx].num_ppu_chiplet <= num_ppu_chiplet
        /* GE Compute */ 
        &&  ge_stats.ge_stat_vec[band_idx].coordinate_trans.compute / compute_density_per_sec <= ge.coord_trans_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].nr_engine.compute / compute_density_per_sec <= ge.nr_engine_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].relative_orientation.compute / compute_density_per_sec <= ge.relative_orient_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].antenna_gain.compute / compute_density_per_sec <= ge.antenna_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].path_velocity.compute / compute_density_per_sec <= ge.path_gain_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].rcs.compute / compute_density_per_sec <= ge.rcs_compute_area
        &&  ge_stats.ge_stat_vec[band_idx].tu.compute / compute_density_per_sec <= ge.tu_compute_area
        /* GE Memory */
        &&  ge_stats.ge_stat_vec[band_idx].coordinate_trans.memory / mem_density <= ge.coord_trans_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].nr_engine.memory / mem_density <= ge.nr_engine_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].relative_orientation.memory / mem_density <= ge.relative_orient_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].antenna_gain.memory / mem_density <= ge.antenna_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].path_velocity.memory / mem_density <= ge.path_gain_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].rcs.memory / mem_density <= ge.rcs_memory_area
        &&  ge_stats.ge_stat_vec[band_idx].tu.memory / mem_density <= ge.tu_memory_area;

        if(can_support_this){
          num_supported_scenario ++;
        }
      }
      return num_supported_scenario;
  }

ge_core * design_ge_core_for_scenario(path_proc_unit * ppu, std::vector<Band>& scene_vec, drbe_wafer& w, WaferStats & w_stats, 
                                      GEStats & ge_stats, PPUStats & ppu_stats, bool ge_cpu, bool ge_asic, bool ge_cgra) {
  ge_core* ge = new ge_core(&*w._t); //needs to be null so we don't delete a non-pointer

  // This function calculate the compute/memory/latency/bandwidth required by each block
  evaluate_ge(scene_vec, ge_stats, w, ge_cpu);

  /*
  We use all scenarios to find the optimal point of CPU/CGRA
  */    
  // Get the total amount of chiplet per wafer
  float total_chiplets = w.num_units() * w_stats.num_wafer;
  int most_num_scenario_supported = 0;
  int most_scenario_ge_chiplet = 0;
  int most_scenario_gc_chiplet = 0;
  int most_scenario_gm_chiplet = 0;
  // Find: GE chiplet vs. PPU Chiplet
  for(float ge_chiplet_ratio = 0.5; ge_chiplet_ratio < 100; ge_chiplet_ratio += 0.5 ){
    int ge_chiplet = (float)ge_chiplet_ratio / 100 * total_chiplets;
    // This is just the area for Geometry Engine, mixing Compute and Memory
    float ppu_chiplet = total_chiplets - ge_chiplet;

    // Find: GC chiplet vs. GM chiplet
    for(float gm_chiplet_ratio = 0.5; gm_chiplet_ratio < 100; gm_chiplet_ratio += 0.5){
      int gm_chiplet = (float)gm_chiplet_ratio / 100 * ge_chiplet;
      int gc_chiplet = ge_chiplet - gm_chiplet; // calculate the number of geometry compute chiplet
      float gc_area = gc_chiplet * 20;// calculate the geometry compute area
      float gm_area = gm_chiplet * 20;// calculate the geometry memory area
      // Loop over all scenarios
      int num_support_scenario = 0;
      for(unsigned band_idx = 0; band_idx < scene_vec.size(); ++band_idx) {
        // Get the PPU area required by this scenario
        float current_band_ppu_chiplet = ppu_stats.ppu_stat_vec[band_idx].num_ppu_chiplet;
        // Get the geometry compute area required by this scenario
        float current_band_gc_area = ge_stats.ge_stat_vec[band_idx].compute_area;
        // Get the geometry memory area required by this scenario
        float current_band_gm_area = ge_stats.ge_stat_vec[band_idx].mem_area;
        // debug
        //printf("Band[%d] needs %f compute area and %f memory area\n", band_idx, current_band_gc_area, current_band_gm_area);
        //return ge;
        //printf("The memory of RCS is %f\n (Unit is 64-bit) ",  ge_stats.ge_stat_vec[band_idx].rcs.memory);
        // We first check whether the memory and ppu requirements are met
        if(current_band_gm_area <= gm_area && /* Geometry Memory Area Satisfied*/
          current_band_ppu_chiplet <= ppu_chiplet/* PPU Chiplets count Satisfied*/){
            if(current_band_gc_area <= gc_area){
              num_support_scenario++;
            }
        }
      }// End of All scenarios
      if(num_support_scenario > most_num_scenario_supported){
        most_num_scenario_supported = num_support_scenario;
        most_scenario_ge_chiplet = ge_chiplet;
        most_scenario_gc_chiplet = gc_chiplet;
        most_scenario_gm_chiplet = gm_chiplet;
      }
    }// End of loop over from 0% of GM to 100% of GM (100% GC to 0% GC)
  }// End of loop over from 0% of GE to 100 % of GE
  // debug
  printf("%d out of %d scenarios (%.2f) are support by CGRA/CPU design\n", 
    most_num_scenario_supported, (int)scene_vec.size(), 100 * (float)most_num_scenario_supported / (float)scene_vec.size() );

  w_stats.num_ge_chiplet = most_scenario_ge_chiplet;
  w_stats.num_ppu_chiplet = total_chiplets - most_scenario_ge_chiplet;
  printf("Total #chiplet is %f, most scenario GE #chiplet is %d\n", total_chiplets, most_scenario_ge_chiplet);

  w_stats.num_gc_chiplet = most_scenario_gc_chiplet;
  w_stats.num_gm_chiplet = most_scenario_gm_chiplet;
  /*
  If we are using ASIC, then we use the average of all scenarios as start point,
  the point of doing ASIC design, is that we want to see how it perform compared against CGRA
  */
  float compute_density_per_sec = w._t->ge_comp_density() * w._t->ge_freq();
  float mem_density = w._t->ge_dram_MB() * 1024 * 1024 * 8 / 64;
  if(ge_asic){
    // In order to make the comparison fair, we use the same amout of GE chiplet
    int num_asic_ge_chiplet = w_stats.num_ge_chiplet;
    // Let us use the average Compute and Memory to distribute the ASIC
    // Calculate the average compute/memory area across all scenarios
    float all_scenario_compute_area = 0.0; 
    float all_scenario_memory_area = 0.0;
    // Compute Area breakdown of all scenarios
    float all_scenario_coord_trans_compute_area = 0.0;
    float all_scenario_nr_engine_compute_area = 0.0;
    float all_scenario_relative_orient_compute_area = 0.0;
    float all_scenario_antenna_compute_area = 0.0;
    float all_scenario_path_gain_compute_area = 0.0;
    float all_scenario_rcs_compute_area = 0.0;
    float all_scenario_tu_compute_area = 0.0;
    // memory Area breakdown of all scenarios
    float all_scenario_coord_trans_memory_area = 0.0;
    float all_scenario_nr_engine_memory_area = 0.0;
    float all_scenario_relative_orient_memory_area = 0.0;
    float all_scenario_antenna_memory_area = 0.0;
    float all_scenario_path_gain_memory_area = 0.0;
    float all_scenario_rcs_memory_area = 0.0;
    float all_scenario_tu_memory_area = 0.0;
    // Start gathering the statistics
    for(unsigned band_idx=0; band_idx < scene_vec.size(); ++band_idx) {
      // Overall statistic
      all_scenario_memory_area += ge_stats.ge_stat_vec[band_idx].mem_area;
      all_scenario_compute_area += ge_stats.ge_stat_vec[band_idx].compute_area;
      // Compute
      all_scenario_coord_trans_compute_area += ge_stats.ge_stat_vec[band_idx].coordinate_trans.compute / compute_density_per_sec;
      all_scenario_nr_engine_compute_area += ge_stats.ge_stat_vec[band_idx].nr_engine.compute / compute_density_per_sec;
      all_scenario_relative_orient_compute_area += ge_stats.ge_stat_vec[band_idx].relative_orientation.compute / compute_density_per_sec;
      all_scenario_antenna_compute_area += ge_stats.ge_stat_vec[band_idx].antenna_gain.compute / compute_density_per_sec;
      all_scenario_path_gain_compute_area += ge_stats.ge_stat_vec[band_idx].path_velocity.compute / compute_density_per_sec;
      all_scenario_rcs_compute_area += ge_stats.ge_stat_vec[band_idx].rcs.compute / compute_density_per_sec;
      all_scenario_tu_compute_area += ge_stats.ge_stat_vec[band_idx].tu.compute / compute_density_per_sec;
      // Memory
      all_scenario_coord_trans_memory_area += ge_stats.ge_stat_vec[band_idx].coordinate_trans.memory / mem_density;
      all_scenario_nr_engine_memory_area += ge_stats.ge_stat_vec[band_idx].nr_engine.memory / mem_density;
      all_scenario_relative_orient_memory_area += ge_stats.ge_stat_vec[band_idx].relative_orientation.memory / mem_density;
      all_scenario_antenna_memory_area += ge_stats.ge_stat_vec[band_idx].antenna_gain.memory / mem_density;
      all_scenario_path_gain_memory_area += ge_stats.ge_stat_vec[band_idx].path_velocity.memory / mem_density;
      all_scenario_rcs_memory_area += ge_stats.ge_stat_vec[band_idx].rcs.memory / mem_density;
      all_scenario_tu_memory_area += ge_stats.ge_stat_vec[band_idx].tu.memory / mem_density;
    }
    float average_memory_area = all_scenario_memory_area / (float) scene_vec.size();
    float average_compute_area = all_scenario_compute_area / (float) scene_vec.size();
    float average_total_area = average_memory_area + average_compute_area;
    // Calculate the initial number of GC/GM chiplet
    int num_asic_gc_chiplet = ceil((float)num_asic_ge_chiplet * 
      (average_compute_area / average_total_area));
    float gc_chiplet_area = num_asic_gc_chiplet * 20;// each chiplet is 20 um2
    int num_asic_gm_chiplet = floor((float)num_asic_ge_chiplet * 
      (average_memory_area / average_total_area));
    float gm_chiplet_area = num_asic_gm_chiplet * 20;// each chiplet is 20 um2
    // Calculate the initial breakdown of each block (using the average breakdown)
    ge->coord_trans_compute_area = gc_chiplet_area * all_scenario_coord_trans_compute_area / all_scenario_compute_area;
    ge->nr_engine_compute_area = gc_chiplet_area * all_scenario_nr_engine_compute_area / all_scenario_compute_area;
    ge->relative_orient_compute_area = gc_chiplet_area * all_scenario_relative_orient_compute_area / all_scenario_compute_area;
    ge->antenna_compute_area = gc_chiplet_area * all_scenario_antenna_compute_area / all_scenario_compute_area;
    ge->path_gain_compute_area = gc_chiplet_area * all_scenario_path_gain_compute_area / all_scenario_compute_area;
    ge->rcs_compute_area = gc_chiplet_area * all_scenario_rcs_compute_area / all_scenario_compute_area;
    ge->tu_compute_area = gc_chiplet_area * all_scenario_tu_compute_area / all_scenario_compute_area;
    // Calculate the initial breakdown of each block (using the average breakdown)
    ge->coord_trans_memory_area = gm_chiplet_area * all_scenario_coord_trans_memory_area / all_scenario_memory_area;
    ge->nr_engine_memory_area = gm_chiplet_area * all_scenario_nr_engine_memory_area / all_scenario_memory_area;
    ge->relative_orient_memory_area = gm_chiplet_area * all_scenario_relative_orient_memory_area / all_scenario_memory_area;
    ge->antenna_memory_area = gm_chiplet_area * all_scenario_antenna_memory_area / all_scenario_memory_area;
    ge->path_gain_memory_area = gm_chiplet_area * all_scenario_path_gain_memory_area / all_scenario_memory_area;
    ge->rcs_memory_area = gm_chiplet_area * all_scenario_rcs_memory_area / all_scenario_memory_area;
    ge->tu_memory_area = gm_chiplet_area * all_scenario_tu_memory_area / all_scenario_memory_area;
    /*
    Exploration Start
    */
    int initial_num_supported = get_num_supported_scenario(*ge, w, scene_vec.size(),w_stats.num_ppu_chiplet , ge_stats, ppu_stats);
    print_ge_breakdown(*ge);
    printf("%d out of %d scenarios (%.2f) are supportted by ASIC design\n", 
    initial_num_supported, (int)scene_vec.size(), 100 * (float)initial_num_supported / (float)scene_vec.size() );

  }

  
  
  // Assertion
  assert(w_stats.num_ge_chiplet == (w_stats.num_gc_chiplet + w_stats.num_gm_chiplet) && 
    "Requirement: num_ge_chiplet == num_gc_chiplet + num_gm_chiplet");
  return ge;
}

bool model_successful(path_proc_unit* ppu, Band& b, drbe_wafer& w, WaferStats & w_stats, int max_wafers) {

  int failures=0;
  bool verbose=false;
  PPUStats ppu_stats; //ignored for now
  float wafer_unconstrained_ppus = b.ppus_per_band(*ppu,w,ppu_stats,verbose);
  float links_per_wafer = b.calc_links_per_wafer(wafer_unconstrained_ppus,ppu,w,w_stats,verbose);
  float num_wafers = b.num_links() / links_per_wafer;

  if(failures) {
    printf("model mapping failures: %d\n", failures);
    assert(0);
  }

  //if(failures!=0) return false; // don't include this scenario if failed

  //if(num_wafers > max_wafers ){
  //  printf("model unsuccessful, no_io_ppus: %f, links_per_wafer: %f, num_wafers %f\n",
  //      wafer_unconstrained_ppus, links_per_wafer, num_wafers);
  //} 
  return ceil(num_wafers) <= max_wafers;
}

void make_harder(Band& band, path_proc_unit* ppu, drbe_wafer& w, 
                 WaferStats& w_stats, int max_wafers) {
  Band try_band;
  int num_tries=0;
  bool could_increase=true;

  while(could_increase && num_tries < 50) {
    try_band=band;
   
    assert(model_successful(ppu,try_band,w, w_stats, max_wafers));

    if( /*int which = */try_band.increase_difficulty()) {
      //cout << "increased" << which -1 << "\n";
      if(model_successful(ppu,try_band,w, w_stats, max_wafers)) {
        band=try_band;
        num_tries=0;
      } else {
        //cout << "could not increase " << which-1 <<"\n";
        num_tries++;
      }
    } else {
      could_increase=false;
    }
  }
}

std::vector<Band> top_k_pareto_scenarios(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, WaferStats & w_stats, int k, int max_wafers, int final_max_wafers) {
  std::vector<Band> top_k_scenarios;
  std::vector<Band> successful_pareto_bands;

  std::vector<vector<Band>> top_k_scenarios_per_max_wafers;

  top_k_scenarios.resize(2);

  for(auto& band : scene_vec) {
    bool success = model_successful(ppu,band,w,w_stats,max_wafers);
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
    Band hardest_possible_band=band;

    //cout << "BEFORE\n";
    //hardest_possible_band.print_csv();
    //cout << "\n";

    make_harder(hardest_possible_band,ppu,w,w_stats,max_wafers);

    //cout << "AFTER\n";
    //hardest_possible_band.print_csv();
    //cout << "\n";

    //while we do everything, check and see if we found a higher # platforms
    if(top_k_scenarios[0].num_links() < hardest_possible_band.num_links()) {
      top_k_scenarios[0] = hardest_possible_band;
    }

    if(top_k_scenarios[1].link_complexity() < hardest_possible_band.link_complexity()) {
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

  //if(top_k_scenarios[0].reflectors() > ScenarioGen::max_platforms() * .9) {
  //  top_k_scenarios.resize(1);
  //}

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
  for(auto& band : top_k_scenarios) {
    band.print_csv();
    printf("\n");
  }

  for(auto& band : top_k_scenarios) {
    band.print_detailed();
    PPUStats ppu_stats;

    int ppus = band.ppus_per_band(*ppu,w,ppu_stats,false);

    printf(" PPUs: %d \n",ppus);
  }

  printf("num succ: %d \n", (int)successful_pareto_bands.size());


  // ----------------------------------------------------------------
  // Now we go nuts and try to improve each scenario
  top_k_scenarios_per_max_wafers.push_back(top_k_scenarios);

  for(int w_iter=1; w_iter < final_max_wafers-max_wafers; ++w_iter) {
    //Copy the scenarios from the previous iteration
    top_k_scenarios_per_max_wafers.push_back(top_k_scenarios_per_max_wafers[w_iter-1]);

    //Now make each one as hard as possible with the new number of wafers
    int cur_max_wafers = w_iter + max_wafers;
    for(int scene_ind = 0; scene_ind < k; ++scene_ind) {
      make_harder(top_k_scenarios_per_max_wafers[w_iter][scene_ind],
                  ppu,w,w_stats,cur_max_wafers);
    }
  }
  
  for(int scene_ind = 0; scene_ind < k; ++scene_ind) {
    printf("Scene %d\n",scene_ind+1);
    for(int w_iter=0; w_iter < top_k_scenarios_per_max_wafers.size(); ++w_iter) {

      top_k_scenarios_per_max_wafers[w_iter][scene_ind].print_norm_csv();
      printf("\n");
    }
  }

  printf("num_wafers,");
  Band::print_detailed_header();
  printf("\n");

  for(int scene_ind = 0; scene_ind < k; ++scene_ind) {
    printf("Scene %d\n",scene_ind+1);
    for(int w_iter=0; w_iter < top_k_scenarios_per_max_wafers.size(); ++w_iter) {
      printf("%d,",w_iter+max_wafers);
      top_k_scenarios_per_max_wafers[w_iter][scene_ind].print_detailed_body();
      printf("\n");
    }
  }



  return top_k_scenarios;
}

void evaluate_ppu(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, 
                   PPUStats& ppu_stats, WaferStats & w_stats, int num_wafers_target, bool verbose = false) {
  ppu_stats.total_failures=0;
  float wafers=0;
  int total_in_target_wafer=0;

  float total_links = 0;
  float total_ppus = 0;
  float total_coef = 0;

  int i = 0; // band index
  for(auto & b : scene_vec) {

    float wafer_unconstrained_ppus = b.ppus_per_band(*ppu,w, ppu_stats,verbose);
    total_ppus += wafer_unconstrained_ppus;
    ppu_stats.ppu_stat_vec[i].num_ppu_chiplet = wafer_unconstrained_ppus / (ppu->_ppus_per_chiplet);

    float links_per_wafer = b.calc_links_per_wafer(wafer_unconstrained_ppus,ppu,w,w_stats,verbose);

    float num_wafers = b.num_links() / links_per_wafer;

    assert(num_wafers>0);

    float wafer_constrained_ppus = num_wafers * w.num_units() * ppu->_ppus_per_chiplet;

    // For status purposes, I need to use algorithmic links to make a fair comparison
    // between channel models

    total_coef += b._n_obj * b._avg_coef_per_object * b.algorithmic_num_links();
    total_links += b.algorithmic_num_links();
    ppu_stats.ppu_stat_vec[i].num_ppu_chiplet = wafer_constrained_ppus / ppu->_ppus_per_chiplet; 
    //printf("after count the chiplet from PPU, there are %d chiplets\n", b.num_ppu_chiplet);


    int int_num_wafers = ceil(num_wafers);
    //w_stats.num_wafer = int_num_wafers;
    if(int_num_wafers<=num_wafers_target) total_in_target_wafer++;
    //wafers += int_num_wafers;
    wafers += num_wafers;
    int max_histo_elem=sizeof(ppu_stats.wafer_histo)/sizeof(ppu_stats.wafer_histo[0]);
    int_num_wafers=max(0,min(max_histo_elem-1,int_num_wafers));
    ppu_stats.wafer_histo[int_num_wafers]++;
    i++;
  }

  ppu_stats.percent_in_target_wafer = total_in_target_wafer / (float) scene_vec.size();
  ppu_stats.avg_links_per_wafer = total_links / wafers;
  ppu_stats.avg_links_per_mm2 = total_links / (total_ppus * ppu->area());
  ppu_stats.avg_wafers = wafers / scene_vec.size();
  ppu_stats.avg_coef_per_ppu = total_coef / total_ppus;
  ppu_stats.avg_coef_per_mm2 = total_coef / (total_ppus * ppu->area());
}

//path_proc_unit* design_ppu_for_max_capability(drbe_wafer& w, bool dynamic_reconfig,
//                                              bool direct_path, bool aidp, bool easy_scenario) {
//
//  for(int fixed_platforms = ScenarioGen::min_platforms() - 20; 
//          fixed_platforms <= 300; fixed_platforms +=10){
//    std::vector<Band> scene_vec;
//    // Generate Scenarios
//    ScenarioGen::gen_scenarios(scene_vec, direct_path, aidp,
//                               1 /*num scenarios=1*/, fixed_platforms, !easy_scenario);
//  }
//
//
//
//}



path_proc_unit* design_ppu_for_scenarios(std::vector<Band>& scene_vec, drbe_wafer& w, WaferStats & w_stats,
                                         bool dynamic_reconfig) {

  int num_scenario = scene_vec.size();
  PPUStats best_stats;
  best_stats.ppu_stat_vec = new ppu_stat_per_band[num_scenario];
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

      float ratio_of_input_to_output = (float) OUTPUT_BITWIDTH / (float) INPUT_BITWIDTH;

      // we multiply both of these expressions by 2 because we need both inputs and outputs
      int remain=total_ins_per_side - (2*2 + agg_network*2*ratio_of_input_to_output);
      if(remain < 2) break;


      // Loop over coefficients per cluster
      for(int cpc = 20; cpc < 21; cpc+=5) { //FIXME: turn off for speed


        for(int clutter_index = 0; clutter_index <= 100; clutter_index+=5) {
          int clutter_ratio=clutter_index;

          for(int mem_ratio =1; mem_ratio < 99; mem_ratio+=2) {
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
 
            PPUStats ppu_stats;
            ppu_stats.ppu_stat_vec = new ppu_stat_per_band[scene_vec.size()];
            evaluate_ppu(ppu,scene_vec,w, ppu_stats,w_stats,1 /*target wafers*/);

            //Optimize the average coefficient per mm2
            //if((stats.percent_in_target_wafer > best_stats.percent_in_target_wafer) 
            //    && ppu_stats.total_failures==0) {
            
            if(ppu_stats.total_failures!=0){
              std::cout << "scenario failed\n";
            }
            /*
            std::cout << "avg_corf: " << ppu_stats.avg_coef_per_mm2
                      << "best_coef: " << best_stats.avg_coef_per_mm2
                      << "fail: " << ppu_stats.total_failures<< endl;
            */
            if((ppu_stats.avg_coef_per_mm2 > best_stats.avg_coef_per_mm2) && ppu_stats.total_failures==0) {
              path_proc_unit* old_best = best_ppu;
              best_stats = ppu_stats;
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

void print_ppu_overheads(PPUStats& ppu_stats) {
  for(int i = DatapathOverheads::NUM_ENTRIES-1; i >=0; --i) {
    printf("%s: %0.1f,", DatapathOverheads::nameof(i+1), 
        ppu_stats.overhead.reason[i]/ppu_stats.num_scenarios);
  }
}

void print_performance_summary(float v, int dielet_area, path_proc_unit* best_ppu, 
                               tech_params& t,
                               std::vector<Band>& scene_vec, PPUStats& ppu_stats,
                               int num_wafers_target) {
    printf("v: %0.3f, ", v);
    
    printf("avg_clut: %f, "\
          "%dmm^2 die (%0.2f), in-MB: %0.2f, clust: %d, flex_clust: %d, coef/clust %d, "\
          "In: %d/%d Agg: %d/%d, Coef: %d/%d, Mem Ratio: %d, "\
          "ppus/die: %d, " \
          "Coef/mm2: %0.2f, links/cm2: %0.2f, "\
          "avg_waf: %0.2f, links/ppu: %0.2f, links/wafer: %0.0f, per_targ_waf: %0.1f, targ_waf: %d, fail: %d ",
        Band::average_clutter(scene_vec),
        dielet_area, 
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
        ppu_stats.avg_coef_per_mm2,
        ppu_stats.avg_links_per_mm2 * 100,
        ppu_stats.avg_wafers,
        ppu_stats.avg_links_per_mm2 * best_ppu->area(),
        ppu_stats.avg_links_per_wafer,
        ppu_stats.percent_in_target_wafer*100,
        num_wafers_target,
        ppu_stats.total_failures);
    print_ppu_overheads(ppu_stats);
}


void print_wafer_tradeoff(path_proc_unit& ppu, drbe_wafer& w, WaferStats & w_stats,
                          bool direct_path, bool aidp, int easy, 
                          bool ge_cpu, bool ge_asic, bool ge_cgra, bool ge_hard, bool pre_ge) {

    int fixed_platforms = 8;
 
    static const int MAX_WAFERS=21;
    // #platform, #links, mem_ratio
    std::vector<std::tuple<int,int,float>> metric;
    metric.resize(MAX_WAFERS);

    w.set_limit_wafer_io(false);

    for(int num_wafers=1; num_wafers < MAX_WAFERS; ++num_wafers) {
      if(num_wafers ==2) {
        w.set_limit_wafer_io(true);
      }
      w_stats.num_wafer = num_wafers;

      for(; fixed_platforms <= 10000; fixed_platforms +=1){
      
        std::vector<Band> scene_vec;
        ScenarioGen::gen_scenarios(scene_vec, direct_path, aidp,1/*num_scene=1*/, 
                                   fixed_platforms, easy);
        PPUStats ppu_stats;
        ppu_stats.ppu_stat_vec = new ppu_stat_per_band[scene_vec.size()];
        GEStats ge_stats;
        ge_stats.ge_stat_vec = new ge_stat_per_band[scene_vec.size()];
        ScenarioGen::set_fidelity(scene_vec, ge_stats, w, ge_hard);
        
        // Design PPU and evaluate it
        path_proc_unit* new_ppu = design_ppu_for_scenarios(scene_vec,w,w_stats, false);
        evaluate_ppu(new_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers,false /*verbose*/);  
       
        if(!pre_ge) { 
          // Design GE and evaluate it
          ge_core* ge = design_ge_core_for_scenario(new_ppu,scene_vec,w,w_stats,ge_stats, ppu_stats, ge_cpu, ge_asic, ge_cgra);

          // We need to redesign the PPU since some resources are used by GE
          new_ppu = design_ppu_for_scenarios(scene_vec,w,w_stats, false);
          evaluate_ppu(new_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers,false /*verbose*/);  
        }
        
        if(ppu_stats.avg_wafers <= num_wafers && (w_stats.num_ge_chiplet > 0 || pre_ge)) {
          //we succeeded, record fixed platforms
          //printf("#wafer = %d, #platform = %d, #ge_chiplet = %d (#gc_chiplet = %d, #gm_chiplet = %d)\n",num_wafers, fixed_platforms, w_stats.num_ge_chiplet, w_stats.num_gc_chiplet, w_stats.num_gm_chiplet);
        } else {
          //printf("#wafer = %d, #platform = %d, #ge_chiplet = %d (#gc_chiplet = %d, #gm_chiplet = %d)\n",num_wafers, fixed_platforms, w_stats.num_ge_chiplet, w_stats.num_gc_chiplet, w_stats.num_gm_chiplet);
          metric[num_wafers]=std::make_tuple(scene_vec[0].platforms(),
                                             scene_vec[0].num_links(),
                                             new_ppu->_mem_ratio);
          // we failed, break out and try increasing num_wafers
          break;
        }
     }
   }

   printf("num_wafers num_platforms num_links, mem_ratio\n");
   for(int num_wafers=1; num_wafers < MAX_WAFERS; ++num_wafers) {
     printf("%8d %8d %8d %8.2f\n", num_wafers, std::get<0>(metric[num_wafers]),
                                            std::get<1>(metric[num_wafers]),
                                            std::get<2>(metric[num_wafers]));

   }
   printf("\n");
}

void print_graph(std::vector<float> x, std::vector<float> y) {
  for(float v : x) {
    cout << v << " ";
  }
  cout << "\n";
  for(float v : y) {
    cout << v << " ";
  }
  cout << "\n";
}

int main(int argc, char* argv[]) {
  bool verbose = false;
  bool direct_path = false;
  bool dynamic_reconfig = false;
  bool verbose_ge_info = false;
  bool aidp = false;
  bool challenge_scenario=false;
  int num_scenarios=1000;
  bool limit_wafer_io=true; //DEFAULT IS TRUE
  int chiplet_io_layer=4;
  int easy_scenario=0; //false

  bool print_pareto = false;
  bool print_wafer_scaling = false;
  bool print_hw_config = false;
  bool print_pareto_after_ge = false;

  int num_wafers_target=1;
  int final_wafers_target=1;

  string dump_file_name="";
  tech_params t;

  bool ge_cpu = false;
  bool ge_cgra= true;
  bool ge_asic= false;
  bool ge_hard = false;

  bool sense=false;
  bool sense_tech_scaling=false;
  bool sense_wafer_io=false;
  bool sense_tx_sparsity=false;

  bool print_pre_ge_summary=false;

  float pulsed_duty_cycle = 0.5f;

  float dielet_area = 20; //mm^2
  float tech_scale=1;

  
  int dielets_per_wafer = -1; //units

  int opt;
  while ((opt = getopt_long(argc, argv, "vgkdprcn:wl:f:ezt:", long_options, nullptr)) != -1) {
    switch (opt) {
      case 'd': direct_path=true; break;
      case 'c': challenge_scenario=true; break;
      case 'k': aidp=true; break;
      case 'r': dynamic_reconfig=true; break;
      case 'v': verbose = true; break;
      case 'g': verbose_ge_info = true; break;
      case 'n': num_scenarios = atoi(optarg); break;
      case 'w': limit_wafer_io = false; break;
      case 'e': easy_scenario = 1; break;
      case 'z': easy_scenario = 2; break;
      case 'l': chiplet_io_layer = atoi(optarg); break;
      case 't': num_wafers_target = atoi(optarg); break;
      case 'u': final_wafers_target = atoi(optarg); break;
      case 'f': dump_file_name = optarg; break;
      case 'p': print_pre_ge_summary = true; break;

      case PARAM_PULSED_TX_DUTY: pulsed_duty_cycle = atof(optarg); break;

      case PARAM_DIELET_AREA: dielet_area = atof(optarg); break;
      case PARAM_TECH_SCALING: tech_scale = atof(optarg); break;
      case PARAM_DIELETS_PER_WAFER: dielets_per_wafer = atoi(optarg); break;

      case GE_CPU: {ge_cpu = true;ge_asic = false; ge_cgra = false; break;}
      case GE_ASIC: {ge_cpu = false;ge_asic = true; ge_cgra = false; break;}
      case GE_CGRA: {ge_cpu = false;ge_asic = false; ge_cgra = true; break;}
      case GE_HARD: {ge_hard = true; break;}

      case SENSITIVITY_TECH_SCALING: {sense=true; sense_tech_scaling = true; break;}
      case SENSITIVITY_WAFER_IO:     {sense=true; sense_wafer_io = true; break;}
      case SENSITIVITY_TX_SPARSITY:  {sense=true; sense_tx_sparsity = true; break;}

      case ANALYSIS_PARETO:        print_pareto = true; break;
      case ANALYSIS_WAFER_SCALING: print_wafer_scaling = true; break;
      case ANALYSIS_HW_CONFIG:     print_hw_config = true; break;
      case ANALYSIS_PARETO_AFTER_GE: print_pareto_after_ge = true; break;

      default: exit(1);
    }
  }
  // We use same easy/hard for ppu and GE
 
  printf("Num Scenarios: %d\n", num_scenarios); 

  if(limit_wafer_io==false) {
    cout << "WARNING: Unlimited Wafer IO.  (this only pertains to sensitivity analysis)\n";
  }

  argc -= optind;
  argv += optind;

  if(pulsed_duty_cycle != 1.0) {
    cout << "Pulsed TX Duty Cycle: " << pulsed_duty_cycle << "\n";
  }

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

  if(ge_cpu) printf("Assuming GE CPU model\n");
  if(ge_asic) printf("Assuming GE ASIC model\n");
  if(ge_cgra) printf("Assuming GE CGRA model\n");
  if(ge_hard) printf("GE Hard Case\n"); else printf("GE Easy Case\n");


  t.set_area_multiplier(tech_scale);

  std::vector<Band> scene_vec;
  // Generate Scenarios
  // Note that "easy_scenario" is ignored unless num_scenarios == 1
  ScenarioGen::gen_scenarios(scene_vec, direct_path, aidp,
                             num_scenarios, 0, easy_scenario);

  for(auto& b : scene_vec) {
    b._pulsed_duty_cycle = pulsed_duty_cycle;
  }



  // Print the header for ge_trade
  std::ofstream ge_tradeoff;
  if(dump_file_name != "") ge_tradeoff = print_ge_tradeoff(dump_file_name);

  // You can run different experiments by uncommenting different loops here ... we should make this
  // more configurable in the future. : )

  drbe_wafer w(&t,300,(float)dielet_area);

  if(dielets_per_wafer!=-1) {
    w.set_num_units(dielets_per_wafer);
  }

  //float fast_update_period=10000;

  //for(fast_update_period = 1000; fast_update_period < 1000000; 
  //    fast_update_period*=1.2589254117941672104239541063958) {
  
  //for(float frac_clutter = 0; frac_clutter <= 100; frac_clutter += 20) {

  //  for(auto& b : scene_vec) {
  //    //b._high_update_period=fast_update_period;
  //    b._frac_clutter = frac_clutter/100.0f;
  //  }

  float start_v = 2;
  float end_v = 2.01; 
  float v_factor = 1.4142857;
  bool use_plus=false;
  float v_increment=1;


  if(sense_tech_scaling) {
    printf("\nSensitivity to Technology Density (factor \"v\" below)\n");
    start_v = t.area_multiplier();
    end_v = start_v * 4+0.001;
  } else if(sense_wafer_io) {
    printf("\nSensitivity to Wafer IO (factor \"v\" below)\n");
    start_v = 10; //w.io_tb_per_sec();
    end_v = 200+0.001; 
  } else if(sense_tx_sparsity) {
    printf("\nSensitivity to TX Sparsity (factor \"v\" below)\n");
    start_v = 0; //w.io_tb_per_sec();
    end_v = 1.0001; 
    v_increment=0.1;
    use_plus=true;
  } else  {
     
  }

  //factor = factor * factor;
  for(float v = start_v; v <= end_v; v= use_plus? v+v_increment : v*v_factor) {
    if(sense_tech_scaling) {
      t.set_area_multiplier(v);
    } else if(sense_wafer_io) {
      w.set_io_tb_per_sec(v);
    } else if(sense_tx_sparsity) {
      for(auto& b : scene_vec) {
        b._frac_pulsed=std::min(v,1.0f);
      }
    }
    //int num_wafers_target=v;
    
    //if(num_wafers_target == 1){
    //  limit_wafer_io = false;
    //}else{
    //  limit_wafer_io = true;
    //}
    
    w.set_limit_wafer_io(limit_wafer_io);
    w.set_chiplet_io_layer(chiplet_io_layer);

    // Initialize the recording statistic structure
    // Sihao -- you are leaking memory everywhere -- FIXME : )
    PPUStats ppu_stats;
    ppu_stats.ppu_stat_vec = new ppu_stat_per_band[num_scenarios];
    GEStats ge_stats;
    ge_stats.ge_stat_vec = new ge_stat_per_band[num_scenarios];
    WaferStats w_stats;
    w_stats.num_wafer = num_wafers_target;
    // Set GE fidelity
    ScenarioGen::set_fidelity(scene_vec, ge_stats, w, ge_hard);
    path_proc_unit* best_ppu = design_ppu_for_scenarios(scene_vec,w,w_stats,dynamic_reconfig);


    if(best_ppu == NULL){// if cannot design ppu continue
      continue;
    }

    evaluate_ppu(best_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);

    if(print_pre_ge_summary) {
      print_performance_summary(v, dielet_area, best_ppu, t, scene_vec, ppu_stats, num_wafers_target);
      printf("\n");
    }

    if(print_pareto) {
      if(final_wafers_target < num_wafers_target) {
        final_wafers_target = num_wafers_target;  //just the one
      }

      std::vector<Band> pband = 
        top_k_pareto_scenarios(best_ppu,scene_vec,w,w_stats, 6,num_wafers_target,final_wafers_target);

      // cout << "\nNum Obj\n";
      //for(Band band : pband) {
      //  std::vector<float> x;
      //  std::vector<float> y;
      //  for(int i = 50; i <= 150; i+=5) {
      //    std::vector<Band> tband;
      //    tband.push_back(band);
      //    tband[0]._n_obj=i * tband[0]._n_obj/100;
      //    evaluate_ppu(best_ppu,tband,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);
      //    x.push_back(i);
      //    y.push_back(ppu_stats.avg_wafers);
      //  }
      //  print_graph(x,y);
      //}
      
      cout << "\nLink Complexity\n";
      for(Band band : pband) {
        std::vector<float> x;
        std::vector<float> y;
        for(int i = 0; i <= 20; i+=1) {
          std::vector<Band> tband;
          tband.push_back(band);
          tband[0]._avg_frac_full_objects=i/100.0;
          if(tband[0]._avg_frac_full_objects<0){
            continue;
          }

          evaluate_ppu(best_ppu,tband,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);
          x.push_back(tband[0].link_complexity());
          y.push_back(ppu_stats.avg_wafers);
        }
        print_graph(x,y);
      }

      cout << "\nRange Overprov\n";
      for(Band band : pband) {
        std::vector<float> x;
        std::vector<float> y;
        for(int i = 0; i <= 100; i+=5) {
          std::vector<Band> tband;
          tband.push_back(band);
          tband[0]._n_full_range_obj= i/100.0 * tband[0]._n_obj;
          if(tband[0]._n_full_range_obj<0 || tband[0]._n_full_range_obj>tband[0]._n_obj){
            continue;
          }

          evaluate_ppu(best_ppu,tband,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);
          x.push_back(i);
          y.push_back(ppu_stats.avg_wafers);
        }
        print_graph(x,y);
      }

      cout << "\nRange\n";
      for(Band band : pband) {
        std::vector<float> x;
        std::vector<float> y;
        for(int i = 200; i <= 500; i+=20) {
          std::vector<Band> tband;
          tband.push_back(band);
          tband[0]._range=i;
          if(tband[0]._range<200 || tband[0]._range > 500){
            continue;
          }

          evaluate_ppu(best_ppu,tband,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);
          x.push_back(i);
          y.push_back(ppu_stats.avg_wafers);
        }
        print_graph(x,y);
      }

    }

    if(print_wafer_scaling) {
      print_wafer_tradeoff(*best_ppu, w, w_stats, direct_path, aidp, easy_scenario, 
      ge_cpu, ge_asic, ge_cgra, ge_hard, print_pre_ge_summary);
    }


    if(print_hw_config) {
      best_ppu->print_params();
      best_ppu->print_area_breakdown();


      printf("%d,%d,%d,%d\n",w.num_units(), best_ppu->_ppus_per_chiplet,
          (best_ppu->num_full_clusters()*best_ppu->_coef_per_cluster +
          best_ppu->num_point_clusters()), SCALAR_MACC_PER_COMPLEX_MACC);
      printf("wafer_sram_MB: %0.2f, wafer_macc: %0.2f P-MACC/s, wafer_SRAM_bw:%0.2f PB/s\n", 
          w.num_units() * best_ppu->_ppus_per_chiplet *
           (best_ppu->input_buf_area() + best_ppu->_input_tap_fifo.area() 
                 + best_ppu->_coef_storage.area())*t.sram_Mb_per_mm2()/8,
          w.num_units() * best_ppu->_ppus_per_chiplet *
          (best_ppu->num_full_clusters()*best_ppu->_coef_per_cluster +
          best_ppu->num_point_clusters()) *
          SCALAR_MACC_PER_COMPLEX_MACC * 1000000000.0 / 1024.0 / 1024.0 / 1024.0 / 1024.0 / 1024.0,
          w.num_units() * best_ppu->_ppus_per_chiplet *          
          best_ppu->_num_clusters * best_ppu->_input_bitwidth 
          * 1000000000.0 / 1024.0 / 1024.0 / 1024.0 / 1024.0 / 1024.0 / 8.0
          );
    }

    //  --------------------- Find the optimal ratio of GE chiplet ---------------- 
    ge_core * ge = design_ge_core_for_scenario(best_ppu, scene_vec, w, w_stats, ge_stats, ppu_stats, ge_cpu, ge_asic, ge_cgra);


    if(dump_file_name != ""){
      dump_ge_tradeoff(ge_tradeoff, ge_stats, scene_vec, w_stats);
      if(ppu_stats.percent_in_target_wafer*100 < 1){
        continue;
      }
    }

    // Evaluate PPU again so that the number of ppu chiplet will be substracted by GE chiplet number
    evaluate_ppu(best_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers_target,verbose /*verbose*/);

    if(!print_pre_ge_summary) {
      print_performance_summary(v, dielet_area, best_ppu, t, scene_vec, ppu_stats, num_wafers_target);
      printf(", #wafer = %d, #ge_chiplet = %d (#gc_chiplet = %d, #gm_chiplet = %d), #ppu_chiplet = %d\n", 
              w_stats.num_wafer, w_stats.num_ge_chiplet, w_stats.num_gc_chiplet, w_stats.num_gm_chiplet,
              w_stats.num_ppu_chiplet);
    }

    if(print_pareto_after_ge) {
      top_k_pareto_scenarios(best_ppu,scene_vec,w,w_stats, 6,num_wafers_target,final_wafers_target);
    }

  }
  if(dump_file_name != "") ge_tradeoff.close();

  return 0;
}

