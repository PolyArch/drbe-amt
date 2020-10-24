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
  GE_CPU,
  GE_ASIC,
  GE_CGRA,
};

static struct option long_options[] = {
    {"direct-path",          no_argument,       nullptr, 'd',},
    {"kp-aidp",              no_argument,       nullptr, 'k',},
    {"dynamic-reconfig",     no_argument,       nullptr, 'r',},
    {"challenge-scenario",   no_argument,       nullptr, 'c',},
    {"num-scenarios",        required_argument, nullptr, 'n',},
    {"target-wafers",        required_argument, nullptr, 't',},
    {"easy-scenario",        no_argument,       nullptr, 'e',},
    {"verbose",              no_argument,       nullptr, 'v',},
    {"verbose-ge-info",      no_argument,       nullptr, 'g',},
    {"wafer-io",             no_argument,       nullptr, 'w',},
    {"layer",                required_argument, nullptr, 'l',},
    {"dump-file-name",       required_argument, nullptr, 'f',},

    {"pre-ge-summary",       no_argument,       nullptr, 'p',},

    {"pulsed-duty-cycle",    required_argument, nullptr, PARAM_PULSED_TX_DUTY},

    {"ge-cpu",               no_argument,       nullptr, GE_CPU},
    {"ge-asic",              no_argument,       nullptr, GE_ASIC},
    {"ge-cgra",              no_argument,       nullptr, GE_CGRA},

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

  int max_types=9;
  int which = rand_bt(0,max_types);

  if(kind!=-1) {
    which=kind;
    max_types=1;
  }

  assert(num_links() <= ScenarioGen::max_links());

  for(int i = 0; i < max_types; ++i,which=(which+1)%max_types) {
    switch(which) {
      case 0: { //add a fixed platform
        if(platforms() + 1 >ScenarioGen::max_platforms()) break;
        _n_fixed+=1; //one object per band
        recalculate_txrx();

        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _n_fixed -=1;
          recalculate_txrx();
        } 

      }
      case 1: { //upgrade a fixed platform to slow platform
        if(_n_fixed - 1 < 0) break;
        _n_fixed -= 1;
        _n_slow += 1;
        recalculate_txrx();
        return which+1;
      }
      case 2: { //upgrade a slow platform to fast platform
        if(_n_slow - 1 < 0) break;
        _n_slow -= 1;
        _n_fast += 1;
        recalculate_txrx();
        return which+1;
      }
      case 3: { //subtract a band
        if(_n_bands==1 || rand_bt(0,20)!=0) break; //rand check to make it less likely

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
        if(_avg_coef_per_object +1 > ScenarioGen::max_coef_per_obj()) break;
        _avg_coef_per_object+=1;
        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _avg_coef_per_object-=1;
        } 
      }
      case 5: { //range
        if(_range + 10 > ScenarioGen::max_range()) break;
        _range+=10;
        return which+1;
      }
      case 6: { // more objects
        if(_n_obj + 1 > ScenarioGen::max_objects()) break;
        _n_obj += 1;
        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _n_obj -=1;
        } 
      }
      case 7: { // more clutter
        if(_frac_clutter + 0.01 > ScenarioGen::max_clutter()) break;
        _frac_clutter += 0.01;
        return which+1;
      }
      case 8: { // more link complexity
        _avg_frac_full_objects += 0.001;
        if(ScenarioGen::scenario_is_between_threshold_and_objective(*this)) {
          return which+1;
        } else {
          _avg_frac_full_objects -= 0.001;
        } 
      }
      assert(num_links() <= ScenarioGen::max_links());
    }
  }
  return false;
}

void Band::print_norm_csv() {
  static bool printed = false;
  if(!printed) {
    printf("Links, Tap/Link, Fast, %%Clutter, Range\n");
    printed=true;
  }

  printf("%0.3f, %0.3f, %0.3f, %0.3f, %0.3f", 
    num_links() / (float)ScenarioGen::max_links(),
    link_complexity() / (float)ScenarioGen::overprov_link_complexity(),
    _n_fast / (float)platforms(),
    _frac_clutter / (float)ScenarioGen::max_clutter(),
    _range / (float)ScenarioGen::max_range()
  );
}

void Band::print_csv() {
  printf("%3d, %0.3f, %0.3f, %0.3f, %3d", 
    num_links(),link_complexity(),_n_fast/(float)ScenarioGen::max_platforms(),_frac_clutter, _range);
}

std::vector<float>& Band::normalized_vec() {
  if(_norm_features.size()==0) {
    _norm_features.push_back(num_links() / (float)ScenarioGen::max_links());
    _norm_features.push_back(link_complexity() / (float)ScenarioGen::overprov_link_complexity());
    _norm_features.push_back(_n_fast / (float)platforms());
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
        M = OB*(360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*(360/fed_cmbl.rcs.angle)*(180/fed_cmbl.rcs.angle)*fed_cmbl.rcs.freq*fed_cmbl.rcs.plzn*fed_cmbl.rcs.points*2;
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

struct ppu_stat_per_band {
  float num_ppu = 0.0;
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

  // Per Band Statistic
  ppu_stat_per_band * ppu_stat_vec;
};

void dump_ge_tradeoff(ofstream & ge_tradeoff, GEStats & ge_stats, std::vector<Band>& scene_vec, WaferStats & w_stats){
  for(int i = 0; i < scene_vec.size(); i++){
    ge_tradeoff << ge_stats.print_ge_tradeoff(scene_vec[i],ge_stats.ge_stat_vec[i], w_stats);
  }
}

void evaluate_ge(std::vector<Band>& scene_vec, GEStats & ge_stats, drbe_wafer& w, bool ge_cpu){
  int i = 0;
  for(auto & b : scene_vec) {
    get_ge_cmbl(b, ge_stats.ge_stat_vec[i++], w, ge_cpu);
  }
}

ge_core * design_ge_core_for_scenario(path_proc_unit * ppu, std::vector<Band>& scene_vec, drbe_wafer& w, WaferStats & w_stats, 
                                      GEStats & ge_stats, PPUStats & ppu_stats, bool ge_cpu, bool ge_asic, bool ge_cgra) {
  ge_core* ge = new ge_core(&*w._t); //needs to be null so we don't delete a non-pointer

  evaluate_ge(scene_vec, ge_stats, w, ge_cpu);

  // Get the total amount of chiplet per wafer
  float total_chiplets = w.num_units() * w_stats.num_wafer;
  int most_num_scenario_supported = 0;
  int most_scenario_ge_chiplet = 0;
  int most_scenario_gc_chiplet = 0;
  int most_scenario_gm_chiplet = 0;
  // Find: GE chiplet vs. PPU Chiplet
  for(int ge_chiplet_ratio = 1; ge_chiplet_ratio < 100; ge_chiplet_ratio ++){
    int ge_chiplet = (float)ge_chiplet_ratio / 100 * total_chiplets;
    // This is just the area for Geometry Engine, mixing Compute and Memory
    float ge_area = ge_chiplet * 20; // each chiplet is 20 mm2
    float ppu_chiplet = total_chiplets - ge_chiplet;

    // Find: GC chiplet vs. GM chiplet
    for(int gm_chiplet_ratio = 1; gm_chiplet_ratio < 100; gm_chiplet_ratio ++){
      int gm_chiplet = (float)gm_chiplet_ratio / 100 * ge_chiplet;
      int gc_chiplet = ge_chiplet - gm_chiplet; // calculate the number of geometry compute chiplet
      float gc_area = gc_chiplet * 20;// calculate the geometry compute area
      float gm_area = gm_chiplet * 20;// calculate the geometry memory area
      // Loop over all scenarios
      int band_idx = 0;
      int num_support_scenario = 0;
      for(auto & b : scene_vec){
        // Get the PPU area required by this scenario
        float current_band_ppu_chiplet = ppu_stats.ppu_stat_vec[band_idx].num_ppu;
        // Get the geometry compute area required by this scenario
        float current_band_gc_area = ge_stats.ge_stat_vec[band_idx].compute_area;
        // Get the geometry memory area required by this scenario
        float current_band_gm_area = ge_stats.ge_stat_vec[band_idx].mem_area;
        // We first check whether the memory and ppu requirements are met
        if(current_band_gm_area <= gm_area && /* Geometry Memory Area Satisfied*/
          current_band_ppu_chiplet <= ppu_chiplet/* PPU Chiplets count Satisfied*/){
          if(!ge_asic){
            if(current_band_gc_area <= gc_area){
              num_support_scenario++;
            }
          }else{
              // ASIC design place holder
          }
        }

        band_idx ++;
      }// End of All scenarios
      if(num_support_scenario > most_num_scenario_supported){
        most_num_scenario_supported = num_support_scenario;
        most_scenario_ge_chiplet = ge_chiplet;
        most_scenario_gc_chiplet = gc_chiplet;
        most_scenario_gm_chiplet = gm_chiplet;
      }
    }// End of loop over from 0% of GM to 100% of GM (100% GC to 0% GC)
  }// End of loop over from 0% of GE to 100 % of GE

  w_stats.num_ge_chiplet = most_scenario_ge_chiplet;
  w_stats.num_ppu_chiplet = total_chiplets - most_scenario_ge_chiplet;
  w_stats.num_gc_chiplet = most_scenario_gc_chiplet;
  w_stats.num_gm_chiplet = most_scenario_gm_chiplet;
  
  
  // Assertion
  assert(w_stats.num_ge_chiplet == (w_stats.num_gc_chiplet + w_stats.num_gm_chiplet) && 
    "Requirement: num_ge_chiplet == num_gc_chiplet + num_gm_chiplet");
  return ge;
}

bool model_succesful(path_proc_unit* ppu, Band& b, drbe_wafer& w, WaferStats & w_stats, int max_wafers) {

  int failures=0;
  bool verbose=false;

  float wafer_unconstrained_ppus = b.ppus_per_band(*ppu,w,failures,verbose);
  float links_per_wafer = b.calc_links_per_wafer(wafer_unconstrained_ppus,ppu,w,w_stats,verbose);
  float num_wafers = b.num_links() / links_per_wafer;

  if(failures!=0) return false; // don't include this scenario if failed

  return ceil(num_wafers) <= max_wafers;
}

std::vector<Band> top_k_pareto_scenarios(path_proc_unit* ppu, std::vector<Band>& scene_vec, drbe_wafer& w, WaferStats & w_stats, int k, int max_wafers) {
  std::vector<Band> top_k_scenarios;
  std::vector<Band> successful_pareto_bands;

  top_k_scenarios.resize(2);

  for(auto& band : scene_vec) {
    bool success = model_succesful(ppu,band,w,w_stats,max_wafers);
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

    cout << "BEFORE\n";
    hardest_possible_band.print_csv();
    cout << "\n";

    while(could_increase && num_tries < 20) {
      try_band=hardest_possible_band;
      
      if(int which = try_band.increase_difficulty()) {
        cout << "increased" << which -1 << "\n";
        if(model_succesful(ppu,try_band,w, w_stats, max_wafers)) {
          hardest_possible_band=try_band;
        } else {
          num_tries++;
        }
      } else {
        could_increase=false;
      }
    }

    cout << "AFTER\n";
    hardest_possible_band.print_csv();
    cout << "\n";


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
  printf("num succ: %d \n", (int)successful_pareto_bands.size());


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

    float wafer_unconstrained_ppus = b.ppus_per_band(*ppu,w, ppu_stats.total_failures,verbose);
    total_ppus += wafer_unconstrained_ppus;
    ppu_stats.ppu_stat_vec[i].num_ppu = wafer_unconstrained_ppus;
    ppu_stats.ppu_stat_vec[i].num_ppu_chiplet = wafer_unconstrained_ppus / (ppu->_ppus_per_chiplet);

    float links_per_wafer = b.calc_links_per_wafer(wafer_unconstrained_ppus,ppu,w,w_stats,verbose);

    float num_wafers = b.num_links() / links_per_wafer;

    assert(num_wafers>0);

    float wafer_constrained_ppus = num_wafers * w.num_units() * ppu->_ppus_per_chiplet;

    // For status purposes, I need to use algorithmic links to make a fair comparison
    // between channel models

    total_coef += b._n_obj * b._avg_coef_per_object * b.algorithmic_num_links();
    total_links += b.algorithmic_num_links();
    w_stats.num_ppu_chiplet = wafer_constrained_ppus / ppu->_ppus_per_chiplet; 
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


void print_performance_summary(float v, int ppu_area, path_proc_unit* best_ppu, 
                               tech_params& t,
                               std::vector<Band>& scene_vec, PPUStats& ppu_stats,
                               int num_wafers_target) {
    printf("v: %0.3f, ", v);
    
    printf("avg_clut: %f, "\
          "%dmm^2 PPU (%0.2f), in-MB: %0.2f, clust: %d, flex_clust: %d, coef/clust %d, "\
          "In: %d/%d Agg: %d/%d, Coef: %d/%d, Mem Ratio: %d, "\
          "ppus/die: %d, " \
          "Coef/mm2: %0.2f, links/cm2: %0.2f, "\
          "avg_waf: %0.2f, links/ppu: %0.2f, links/wafer: %0.0f, per_targ_waf: %0.1f, targ_waf: %d, fail: %d",
        Band::average_clutter(scene_vec),
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
        ppu_stats.avg_coef_per_mm2,
        ppu_stats.avg_links_per_mm2 * 100,
        ppu_stats.avg_wafers,
        ppu_stats.avg_links_per_mm2 * best_ppu->area(),
        ppu_stats.avg_links_per_wafer,
        ppu_stats.percent_in_target_wafer*100,
        num_wafers_target,
        ppu_stats.total_failures);
}


void print_wafer_tradeoff(path_proc_unit& ppu, drbe_wafer& w, WaferStats & w_stats,
                          bool direct_path, bool aidp, bool easy, 
                          bool ge_cpu, bool ge_asic, bool ge_cgra) {

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
                                   fixed_platforms, !easy);
        PPUStats ppu_stats;
        ppu_stats.ppu_stat_vec = new ppu_stat_per_band[scene_vec.size()];
        GEStats ge_stats;
        ge_stats.ge_stat_vec = new ge_stat_per_band[scene_vec.size()];
        ScenarioGen::set_fidelity(scene_vec, ge_stats, w, !easy);
        
        // Design PPU and evaluate it
        path_proc_unit* new_ppu = design_ppu_for_scenarios(scene_vec,w,w_stats, false);
        evaluate_ppu(new_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers,false /*verbose*/);  
        
        // Design GE and evaluate it
        ge_core* ge = design_ge_core_for_scenario(new_ppu,scene_vec,w,w_stats,ge_stats, ppu_stats, ge_cpu, ge_asic, ge_cgra);

        // We need to redesign the PPU since some resources are used by GE
        new_ppu = design_ppu_for_scenarios(scene_vec,w,w_stats, false);
        evaluate_ppu(new_ppu,scene_vec,w, ppu_stats,w_stats,num_wafers,false /*verbose*/);  
        
        if(ppu_stats.avg_wafers <= num_wafers && w_stats.num_ge_chiplet > 0) {
          //we succeeded, record fixed platforms
        } else {
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
  bool easy_scenario=false;

  bool print_pareto = false;
  bool print_wafer_scaling = false;
  bool print_hw_config = false;
  bool print_pareto_after_ge = false;

  int num_wafers_target=1;

  string dump_file_name="";
  tech_params t;

  bool ge_cpu = false;
  bool ge_cgra= true;
  bool ge_asic= false;

  bool sense=false;
  bool sense_tech_scaling=false;
  bool sense_wafer_io=false;
  bool sense_tx_sparsity=false;

  bool print_pre_ge_summary=false;

  float pulsed_duty_cycle = 0.5f;

  int opt;
  while ((opt = getopt_long(argc, argv, "vgkdprcn:wl:f:et:", long_options, nullptr)) != -1) {
    switch (opt) {
      case 'd': direct_path=true; break;
      case 'c': challenge_scenario=true; break;
      case 'k': aidp=true; break;
      case 'r': dynamic_reconfig=true; break;
      case 'v': verbose = true; break;
      case 'g': verbose_ge_info = true; break;
      case 'n': num_scenarios = atoi(optarg); break;
      case 'w': limit_wafer_io = false; break;
      case 'e': easy_scenario = true; break;
      case 'l': chiplet_io_layer = atoi(optarg); break;
      case 't': num_wafers_target = atoi(optarg); break;
      case 'f': dump_file_name = optarg; break;
      case 'p': print_pre_ge_summary = true; break;

      case PARAM_PULSED_TX_DUTY: pulsed_duty_cycle = atof(optarg); break;

      case GE_CPU: {ge_cpu = true;ge_asic = false; ge_cgra = false; break;}
      case GE_ASIC: {ge_cpu = false;ge_asic = true; ge_cgra = false; break;}
      case GE_CGRA: {ge_cpu = false;ge_asic = false; ge_cgra = true; break;}

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
  bool ge_hard = !easy_scenario;

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

  std::vector<Band> scene_vec;
  // Generate Scenarios
  // Note that "easy_scenario" is ignored unless num_scenarios == 1
  ScenarioGen::gen_scenarios(scene_vec, direct_path, aidp,
                             num_scenarios, 0, !easy_scenario);

  for(auto& b : scene_vec) {
    b._pulsed_duty_cycle = pulsed_duty_cycle;
  }



  // Print the header for ge_trade
  std::ofstream ge_tradeoff;
  if(dump_file_name != "") ge_tradeoff = print_ge_tradeoff(dump_file_name);

  // You can run different experiments by uncommenting different loops here ... we should make this
  // more configurable in the future. : )

  //for(int ppu_area = 20; ppu_area < 21; ppu_area+=5) {
  int ppu_area=20;
  drbe_wafer w(&t,300,(float)ppu_area);
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

    if(print_pareto) {
      top_k_pareto_scenarios(best_ppu,scene_vec,w,w_stats, 6,num_wafers_target);
    }

    if(print_pre_ge_summary) {
      print_performance_summary(v, ppu_area, best_ppu, t, scene_vec, ppu_stats, num_wafers_target);
      printf("\n");
    }

    if(print_wafer_scaling) {
      print_wafer_tradeoff(*best_ppu, w, w_stats, direct_path, aidp, easy_scenario, 
      ge_cpu, ge_asic, ge_cgra);
    }


    if(print_hw_config) {
      best_ppu->print_params();
      best_ppu->print_area_breakdown();

      printf("wafer_sram_MB: %0.2f, wafer_macc: %0.2f E-MACC/s, wafer_SRAM_bw:%0.2f EB/s\n", 
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
      print_performance_summary(v, ppu_area, best_ppu, t, scene_vec, ppu_stats, num_wafers_target);
      printf(", #wafer = %d, #ge_chiplet = %d (#gc_chiplet = %d, #gm_chiplet = %d), #ppu_chiplet = %d\n", 
              w_stats.num_wafer, w_stats.num_ge_chiplet, w_stats.num_gc_chiplet, w_stats.num_gm_chiplet,
              w_stats.num_ppu_chiplet);
    }

    if(print_pareto_after_ge) {
      top_k_pareto_scenarios(best_ppu,scene_vec,w,w_stats, 6,num_wafers_target);
    }

  }
  if(dump_file_name != "") ge_tradeoff.close();

  return 0;
}

