#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

// Math Constants
#define SCALAR_MACC_PER_COMPLEX_MACC (4)

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

  //These are from uneeb's calculations
  //2.56 Tb/s per side, if each side is 4.4mm
  //Multiply by 2 for 4 layer I/O
  float _chiplet_io_bits_per_mm2=2.56/4.4*1024*2;
  float _router_constant=1.2244e-7; // area of a 1-bit, 1input, 1output router in mm^2

private:
  float _area_multiplier=1.0;
  float _fp_macc_per_mm2=400; //?
  float _int_macc_per_mm2=700;
  float _sram_Mb_per_mm2=8.1f;
  
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
public:
  router(tech_params* t) : hw_unit(t) {
  }
  virtual float area() {
    int num_networks = max(ceil(_in_degree/(float)MaxDegree), 
                           ceil(_out_degree/(float)MaxDegree));

    float area=0;
    area+=(_in_degree/num_networks*4)*(_out_degree/num_networks*4)*_bitwidth*_t->_router_constant;
    area*=num_networks;
    return area;
  }
  int _in_degree=8;
  int _out_degree=8;
  int _bitwidth=32;
};

class path_proc_unit;

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

class ge_core : public hw_unit{
  public:
  ge_core(tech_params * t) : hw_unit(t){
    tech = t;
  }
  virtual float area(){};
  // ----------- Area ----------- 

  tech_params * tech;
};

class path_proc_unit : public hw_unit {
public:
  path_proc_unit(tech_params* t) : hw_unit(t), 
     _coef_router(t), _input_router(t), _output_router(t),
      _coef_storage(t,this), _input_tap_fifo(t,this)   {
    _input_router._in_degree=8;
    _input_router._out_degree=8; // don't really think we need much on input router
    _coef_router._in_degree=8;
    _coef_router._out_degree=8; // don't really think we need much on coef router
    _output_router._in_degree=25;
    _output_router._out_degree=25; //need much mroe on output router
  }  
 
  path_proc_unit(const path_proc_unit &from) : hw_unit(from._t), 
    _coef_router(from._coef_router), _input_router(from._input_router), 
    _output_router(from._output_router),
    _coef_storage(from._coef_storage),
     _input_tap_fifo(from._input_tap_fifo) {
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
    float cluster_area = path_conv_area() + _coef_storage.area();

    cluster_area+=_input_tap_fifo.area();

    float area_per_cluster = cluster_area / _num_clusters;

    int num_clusters = area_for_path_conv / area_per_cluster;
    //printf("area per 1 cluster: %f\n",area_per_cluster);

    _num_clusters = num_clusters;
    _num_flexible_clusters = max(1.0f,num_clusters*clutter_ratio/100);
  }

  void init() {
  }

  float path_conv_area() {
    float area=0;

    int fixed_macc = _num_clusters * _coef_per_cluster; 
    fixed_macc += _num_flexible_clusters * _coef_per_cluster;
    area += fixed_macc * (1.0/_t->int_macc_per_mm2()) * SCALAR_MACC_PER_COMPLEX_MACC;

    // These are FP ops for doppler + aggregating clusters
    // This is a bit conservative, since uneeb said that we don't need FP for this
    int fp_macc = _num_clusters;
    area += fp_macc * (1.0/_t->fp_macc_per_mm2()) * SCALAR_MACC_PER_COMPLEX_MACC;

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

  virtual float area() {
    return path_conv_area() + input_buf_area() + router_area() 
      + _input_tap_fifo.area() + _coef_storage.area();
  }

  void print_params() {
    printf("Number of Clusters: %d\n", _num_clusters);
    printf("Coefficients per Cluster: %d\n", _coef_per_cluster);
    printf("Input Length: %d\n", _input_buffer_length);
  }


  void print_area_breakdown() {
    printf("PPU Area: %f\n", area());
    printf("Conv Area: %f\n", path_conv_area());
    printf("Input Buf Area: %f\n", input_buf_area());
    printf("Router Area: %f\n", router_area());
    printf("Num clusters: %d\n", _num_clusters);
    printf("Input bitwidth: %d\n", _input_bitwidth);

    printf("Input Tap Fifo Area: %f\n", _input_tap_fifo.area());
    printf("Coef. Storage Area: %f\n", _coef_storage.area());
  }

  // we want to avoid div by 0, so lets not return 0
  float num_flexible_clusters() {
    //if(_num_flexible_clusters==0) return 0.001;
    assert(_num_flexible_clusters>0); //ah screw it
    return _num_flexible_clusters;
  }

  router _coef_router;
  router _input_router;
  router _output_router;
  coef_storage _coef_storage;
  input_tap_fifo _input_tap_fifo;
  int _ppus_per_chiplet=1;
  int _num_clusters=140;
  int _num_flexible_clusters=70; // subset of clusters, can be used for clutter
  int _coef_per_cluster=20;
  int _input_buffer_length=3300100;
  int _input_bitwidth=32;
  int _coef_bitwidth=32;
  float _mem_ratio=-1;

  bool _is_dyn_reconfig=false;
};

class drbe_wafer : public hw_unit {
public:
  drbe_wafer(tech_params* t, float diameter, float chiplet_area) : hw_unit(t) {
    _wafer_diameter=diameter;
    _chiplet_area=chiplet_area;
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
  int num_units() {
    if(_num_units==0) {
      //_num_units = area()/_chiplet_area;
     _num_units=2025;
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


  void set_wafer_io(int v) {
    _max_input=v;
    _max_output=v;
  }
  int wafer_io() {
    return _max_input;
  }


  int max_input()  {return _max_input;}
  int max_output() {return _max_output;}

private:
  int _wafer_diameter=0; //mm^2
  float _chiplet_area=0;
  int _num_units=0;
  bool _limit_wafer_io=false;
  path_proc_unit* _ppu=0;
  int _max_input=200;
  int _max_output=200;
};

