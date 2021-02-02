#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

// Math Constants
#define SCALAR_MACC_PER_COMPLEX_MACC (4)
#define SCALAR_MACC_PER_COMPLEX_ADD (1)

#define INPUT_BITWIDTH  (32)
#define OUTPUT_BITWIDTH (64)
#define COEF_BITWIDTH   (32)

#define REG_OVER_SRAM (4)

#define EXTRA_MEM_BW (1.05)

using namespace std;

class tech_params {
public:
  void set_area_multiplier(float m) {_area_multiplier = m;}
  float area_multiplier() {return _area_multiplier;}
  float fp_macc_per_mm2() {return _area_multiplier * _fp_macc_per_mm2;}
  float int_macc_per_mm2(){return _area_multiplier * _int_macc_per_mm2;}
  float sram_Mb_per_mm2() {return _area_multiplier * _sram_Mb_per_mm2;}
  float router_constant() {return _area_multiplier * _router_constant;}
  float mtr_per_mm2() {return _area_multiplier * 24000000;}
  float ge_comp_density(){return _area_multiplier * _ge_comp_density;}
  float ge_dram_MB(){return _area_multiplier * _ge_dram_MB_per_mm2;}
  float ge_freq(){return _ge_clock_rate;}
  float ppu_freq() {return _ppu_clock_rate;}

  //These are from uneeb's calculations
  //2.56 Tb/s per side, if each side is 4.4mm
  //Multiply by 2 for 4 layer I/O
  float _chiplet_io_bits_per_mm2=2.56/4.4*1024*1;
  float _router_constant=1.2244e-7; // area of a 1-bit, 1input, 1output router in mm^2

private:
  float _area_multiplier=1.0;
  float _fp_macc_per_mm2=175*4; //?
  float _int_macc_per_mm2=700*4;
  float _sram_Mb_per_mm2=8.1f;
  float _ppu_clock_rate = 1e9;// Clock of the PPU is 1 GHz

  // GE Params:
  float _ge_dram_MB_per_mm2 = 276; //MegaBytes per mm^2
  float _ge_clock_rate = 2e9;// Clock of the Geometry Engine is 2GHz
  float _ge_comp_density = 200; // 64 bit FP MACs per mm^2

};

class hw_unit {
public:
  virtual float area() = 0; //area in mm2
  //virtual float power() = 0; //area in Watts
  tech_params* _t;

  hw_unit(tech_params* t) : _t(t) {
  }
};

class router : hw_unit {
  static const int MaxDegree=16;
  static const int DelayBufferEntries=0;
public:
  router(tech_params* t) : hw_unit(t) {
  }
  virtual float area() {
    int num_networks = max(ceil(_in_degree/(float)MaxDegree), 
                           ceil(_out_degree/(float)MaxDegree));

    float area=0;
    area+=(_in_degree/num_networks*4)*(_out_degree/num_networks*4)*_bitwidth*_t->router_constant();
    area*=num_networks;

    // we also need some fifos -- also a weak model
    float delay_fifo_bits = DelayBufferEntries * _bitwidth * _in_degree * 4; //arbitrary constant
    float delay_fifo_Mbits = delay_fifo_bits/1024.0/1024.0;
    float delay_fifo_area = delay_fifo_Mbits * (1.0/_t->sram_Mb_per_mm2()) * REG_OVER_SRAM;
    area+=delay_fifo_area;

    return area;
  }
  int _in_degree=8;
  int _out_degree=8;
  int _bitwidth=INPUT_BITWIDTH;
};

class path_proc_unit;

class distribution_network : public hw_unit {
  public:
  distribution_network(tech_params* t, path_proc_unit* p) : hw_unit(t) {
    ppu = p;
  }
  virtual float area();
  path_proc_unit* ppu;
};


class input_tap_fifo : public hw_unit {
  public:
  input_tap_fifo(tech_params* t, path_proc_unit* p) : hw_unit(t) {
    ppu = p;
  }
  virtual float area();
  path_proc_unit* ppu;
};

class coef_storage : public hw_unit {
  public:
  coef_storage(tech_params* t, path_proc_unit* p) : hw_unit(t) {
    ppu = p;
  }
  virtual float area();
  path_proc_unit* ppu;
};

class path_proc_unit : public hw_unit {
public:
  virtual ~path_proc_unit() {
  }

  path_proc_unit(tech_params* t) : hw_unit(t), 
     _coef_router(t), _input_router(t), _output_router(t),
      _coef_storage(t,this), _input_tap_fifo(t,this), _dist_network(t,this)   {
    _input_router._in_degree=8;
    _input_router._out_degree=8; // don't really think we need much on input router
    _coef_router._in_degree=8;
    _coef_router._out_degree=8; // don't really think we need much on coef router
    _output_router._in_degree=25;
    _output_router._out_degree=25; //need much mroe on output router
    _output_router._bitwidth=OUTPUT_BITWIDTH;
  }  
 
  path_proc_unit(const path_proc_unit &from) : hw_unit(from._t), 
    _coef_router(from._coef_router), _input_router(from._input_router), 
    _output_router(from._output_router),
    _coef_storage(from._coef_storage),
    _input_tap_fifo(from._input_tap_fifo),
    _dist_network(from._dist_network) {
    _input_tap_fifo.ppu=this;
    _coef_storage.ppu=this;
  } 

  void set_params_by_mem_ratio(float mem_ratio, float clutter_ratio, float total_area) {
    _mem_ratio = mem_ratio;
    _num_clusters=100;
    _num_flexible_clusters=max(1.0f,_num_clusters*clutter_ratio/100);

    for(int i = 0; i < 4; ++i) {
      set_params_by_mem_ratio_once(mem_ratio,clutter_ratio,total_area);
      if(_num_clusters<=0) {
        _num_clusters=0;
        return;
      }
    }

    assert(_input_tap_fifo.area()>0);
    assert(_input_router.area()>0);
    assert(_input_router.area()>0);
    assert(_output_router.area()>0);
    assert(_coef_router.area()>0);
    assert(path_conv_area()>0);
    assert(_num_flexible_clusters>0); //need one flex cluster for clutter

  }

  void set_params_by_mem_ratio_once(float mem_ratio, float clutter_ratio, float total_area) {
    //First set the input buffer length
    
    //printf("mem ratio %f, total_area %f\n", mem_ratio, total_area);

    float area_for_ibuf = total_area * mem_ratio/100;
    float MBits = area_for_ibuf * _t->sram_Mb_per_mm2();
    int bits = MBits*1024*1024;
    _input_buffer_length=bits/_input_bitwidth;

    //_num_clusters=100; //nominal value for calculation (TODO: fix this)

    //Subtract off the other unit area
    float area_for_path_conv = total_area - area_for_ibuf - router_area();

    //Now set the amount of compute that scales with clusters
    float cluster_area = path_conv_area() + _coef_storage.area() +
                        _input_tap_fifo.area() + _dist_network.area();


    float area_per_cluster = cluster_area / _num_clusters;

    int num_clusters = area_for_path_conv / area_per_cluster;
    //printf("area per 1 cluster: %f\n",area_per_cluster);

    _num_clusters = num_clusters;
    _num_flexible_clusters = max(1.0f,num_clusters*clutter_ratio/100);
  }

  void init() {
  }

  int num_point_clusters() {
    int val = _num_clusters - _num_flexible_clusters;
    assert(val >= 0);
    return val;
  }

  int num_full_clusters() {
    return _num_flexible_clusters;
  }

  float path_conv_area() {
    float area=0;

    int fixed_macc = 0;
    // First, lets just add all the regular multipliers
    fixed_macc += num_point_clusters() * 1;
    fixed_macc += num_full_clusters()  * _coef_per_cluster; 

    // Now lets add in doppler
    fixed_macc += num_full_clusters();
    fixed_macc += num_point_clusters(); //don't multiply, cause only one coefficient

    // Multiply out to get the area
    area += fixed_macc * (1.0/_t->int_macc_per_mm2()) * SCALAR_MACC_PER_COMPLEX_MACC;

    // Now we add in some FP ops for aggregating clusters
    int fp_macc = _num_clusters;
    area += fp_macc * (1.0/_t->fp_macc_per_mm2()) * SCALAR_MACC_PER_COMPLEX_ADD;

    return area;
  }

  float input_buf_area() {
    float bits = _input_buffer_length * _input_bitwidth;
    float Mbits = bits/1024/1024;
    float area = Mbits * (1.0/_t->sram_Mb_per_mm2());
    return area;
  }

  float router_area() {
    return _coef_router.area() + _input_router.area() + _output_router.area();
  }

  float dist_network_area() {
    return _dist_network.area();
  }


  virtual float area() {
    return path_conv_area() + input_buf_area() + router_area() 
      + _input_tap_fifo.area() + _coef_storage.area() + _dist_network.area();
  }

  void print_params() {
    printf("Number of Clusters: %d\n", _num_clusters);
    printf("Coefficients per Cluster: %d\n", _coef_per_cluster);
    printf("Input Length: %d\n", _input_buffer_length);
  }

  void print_area_breakdown() {
    //printf("PPU Area: %f\n", area());
    printf("Conv Area: %f\n", path_conv_area());
    printf("Input Buf Area: %f\n", input_buf_area());
    printf("Router Area: %f\n", router_area());
    printf("Input Tap Fifo Area: %f\n", _input_tap_fifo.area());
    printf("Coef. Storage Area: %f\n", _coef_storage.area());
    printf("Distribution Network Area: %f\n", _dist_network.area());
  }

  // we want to avoid div by 0, so lets not return 0
  float num_flexible_clusters() {
    //if(_num_flexible_clusters==0) return 0.001;
    assert(_num_flexible_clusters>0); //ah screw it
    return _num_flexible_clusters;
  }

  int fifo_bitwidth () {
    return std::max(8,(int)(((float)_num_clusters / (float)_mem_banks) * EXTRA_MEM_BW));
  }

  router _coef_router;
  router _input_router;
  router _output_router;
  coef_storage _coef_storage;
  input_tap_fifo _input_tap_fifo;
  distribution_network _dist_network;
  int _ppus_per_chiplet=1;
  int _num_clusters=140;
  int _num_flexible_clusters=70; // subset of clusters, can be used for clutter
  int _coef_per_cluster=20;
  int _input_buffer_length=3300100;
  int _input_bitwidth=32;
  int _coef_bitwidth=32;
  float _mem_ratio=-1;

  int _compute_banks=4;
  int _mem_banks=4;

  bool _is_dyn_reconfig=false;

  bool _stat_max_platforms_found=0; // max platforms 
};

class drbe_wafer : public hw_unit {
public:
  drbe_wafer(tech_params* t, float diameter, float chiplet_area) : hw_unit(t) {
    _wafer_diameter=diameter;
    _chiplet_area=chiplet_area;
    set_io_tb_per_sec(_io_tb_per_sec);
  }

  virtual float area() {
    float side_len = _wafer_diameter/sqrt(2);
    float area = side_len * side_len;
    return area;
  }

  int depth() {
    float wafer_side_len = sqrt(area());
    float ppu_side_len = sqrt(_chiplet_area);
    return ceil(wafer_side_len / ppu_side_len);
  }
  int half_depth() {
    return ceil(depth()/2);
  }

  void set_num_units(int n) {
    _num_units=n;
  }
  int num_units() {
    if(_num_units==0) {
      //_num_units = area()/_chiplet_area;
     _num_units=ceil(2025 * 20.0 / _chiplet_area);
     // printf("num units: %d\n",_num_units);
    }
    return _num_units;
  }

  

  void reset_computed() {
    _num_units=0;
  }
  void set_chiplet_area(float area) {
    _chiplet_area=area;
    reset_computed();
  }

  float chiplet_area() {return _chiplet_area;}

  bool limit_wafer_io() {
    return _limit_wafer_io;
  }
  void set_limit_wafer_io(bool v) {
    _limit_wafer_io=v;
  }

  void set_chiplet_io_layer(int io_layer){
    switch (io_layer)
    {
    case 2:
      _t->_chiplet_io_bits_per_mm2 = 2.56/4.4*1024*1;
      break;
    case 4:
      _t->_chiplet_io_bits_per_mm2 = 2.56/4.4*1024*2;
      break;
    default:

      break;
    }
  }

  float io_tb_per_sec() {
    return _io_tb_per_sec;
  }

  void set_io_tb_per_sec(float v) {
    _io_tb_per_sec=v;

    // we need to decide i/o distribution
    int remaining_input_ios = (_io_tb_per_sec) / INPUT_BITWIDTH / _t->ppu_freq() 
                             * 1024 * 1024 * 1024 * 1024;

    _input_sig_io = min(200,remaining_input_ios/2);

    remaining_input_ios -= _input_sig_io;
    int remaining_agg_ios = remaining_input_ios * (float) INPUT_BITWIDTH / (float) OUTPUT_BITWIDTH;

    _agg_input_io  = remaining_agg_ios/2;
    _agg_output_io = remaining_agg_ios/2; 
    printf("input_sig %d, agg_input_io %d\n", _input_sig_io, _agg_input_io);
  }

  int input_sig_io()  {return _input_sig_io;}
  int agg_input_io()  {return _agg_input_io;}
  int agg_output_io() {return _agg_output_io;}

private:
  int _wafer_diameter=0; //mm^2
  float _chiplet_area=0;
  int _num_units=0;
  bool _limit_wafer_io=false;
  path_proc_unit* _ppu=0;

  float _io_tb_per_sec=100;
  int _input_sig_io=0;
  int _agg_input_io=0;
  int _agg_output_io=0;

};

class ge_core : public hw_unit{
  public:
  ge_core(tech_params * t) : hw_unit(t){
    tech = t;
  }
  
  // ----------- Compute AREA ----------- 
  // This is dedicated to ASIC design
  float coord_trans_compute_area = 0.0;
  float nr_engine_compute_area = 0.0;
  float relative_orient_compute_area = 0.0;
  float antenna_compute_area = 0.0;
  float path_gain_compute_area = 0.0;
  float rcs_compute_area = 0.0;
  float tu_compute_area = 0.0;

  // ----------- Memory AREA ----------- 
  // This is dedicated to ASIC design
  float coord_trans_memory_area = 0.0;
  float nr_engine_memory_area = 0.0;
  float relative_orient_memory_area = 0.0;
  float antenna_memory_area = 0.0;
  float path_gain_memory_area = 0.0;
  float rcs_memory_area = 0.0;
  float tu_memory_area = 0.0;

  virtual float area(){return 0;} //TODO:FIXME -- really?

  tech_params * tech;
};
