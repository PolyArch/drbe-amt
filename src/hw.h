#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

class tech_params {
public:
  float _fp_macc_per_mm2=400; //?
  float _int_macc_per_mm2=700;
  float _sram_Mb_per_mm2=8.1f;
  int _chiplet_io_bits_per_mm2=200;
  int _macc_per_complex_op=4;
  float _router_constant=1.2244e-7; // area of a 1-bit, 1input, 1output router in mm^2
};

class hw_unit {
public:
  virtual float area() = 0; //area in mm2
  //virtual float power() = 0; //area in Watts
  tech_params _t;

  hw_unit(tech_params& t) : _t(t) {
  }
};

class router : hw_unit {
public:
  router(tech_params& t) : hw_unit(t) {
  }
  virtual float area() {
    float area=0;
    area+=_in_degree*_out_degree*_bitwidth*_t._router_constant;
    return area;
  }
  int _in_degree=8;
  int _out_degree=8;
  int _bitwidth=32;
};

class path_proc_unit;

class input_tap_fifo : public hw_unit {
  public:
  input_tap_fifo(tech_params& t, path_proc_unit* p) : hw_unit(t) {
    ppu = p;
  }
  virtual float area();
  path_proc_unit* ppu;
};

class coef_storage : public hw_unit {
  public:
  coef_storage(tech_params& t, path_proc_unit* p) : hw_unit(t) {
    ppu = p;
  }
  virtual float area();
  path_proc_unit* ppu;
};



class path_proc_unit : public hw_unit {
public:
  path_proc_unit(tech_params& t) : hw_unit(t), 
     _coef_router(t), _input_router(t), _output_router(t),
     _input_tap_fifo(t,this), _coef_storage(t,this)   {
    _input_router._in_degree=8;
    _input_router._out_degree=8; // don't really think we need much on input router
    _coef_router._in_degree=8;
    _coef_router._out_degree=8; // don't really think we need much on coef router
    _output_router._in_degree=25;
    _output_router._out_degree=25; //need much mroe on output router
  }  
  
  void set_params_by_mem_ratio(float mem_ratio, float total_area) {
    //First set the input buffer length
    
    //printf("mem ratio %f, total_area %f\n", mem_ratio, total_area);

    float area_for_ibuf = total_area * mem_ratio/100;
    float MBits = area_for_ibuf * _t._sram_Mb_per_mm2;
    int bits = MBits*1024*1024;
    _input_buffer_length=bits/_input_bitwidth;

    //Subtract off the other unit area
    float area_for_path_conv = total_area - area_for_ibuf 
      - router_area() - _coef_router.area() - _input_tap_fifo.area();

    //Now set the amount of compute
    _num_clusters=100; //nominal value for calculation
    float area_per_cluster = path_conv_area() / _num_clusters;
    int num_clusters = area_for_path_conv / area_per_cluster;

    //printf("area per 1 cluster: %f\n",area_per_cluster);

    _num_clusters = num_clusters;

    assert(_num_clusters >= 0);
  }

  void init() {
  }

  float path_conv_area() {
    float area=0;

    int fixed_macc = _num_clusters * _coef_per_cluster; 
    area += fixed_macc * (1.0/_t._int_macc_per_mm2) * _t._macc_per_complex_op;

    int fp_macc = _num_clusters;
    area += fp_macc * (1.0/_t._fp_macc_per_mm2) * _t._macc_per_complex_op;

    return area;
  }

  float input_buf_area() {
    float bits = _input_buffer_length * _input_bitwidth;
    float Mbits = bits/1024/1024;
    float area = Mbits * (1.0/_t._sram_Mb_per_mm2);
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
    printf("Input Tap Fifo Area: %f\n", _input_tap_fifo.area());
    printf("Coef. Storage Area: %f\n", _coef_storage.area());
  }

  router _coef_router;
  router _input_router;
  router _output_router;
  coef_storage _coef_storage;
  input_tap_fifo _input_tap_fifo;
  int _num_clusters=140;
  int _coef_per_cluster=20;
  int _input_buffer_length=3300100;
  int _input_bitwidth=32;
  int _coef_bitwidth=32;
};

class drbe_wafer : public hw_unit {
public:
  drbe_wafer(tech_params& t) : hw_unit(t) {
  }

  virtual float area() {
    float side_len = 300/sqrt(2);
    float area = side_len * side_len;
    return area;
  }

  int _wafer_diameter; //mm^2
  int _num_units;
  path_proc_unit* _ppu;
};

