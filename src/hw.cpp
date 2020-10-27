#include "hw.h"

float input_tap_fifo::area() {
 
  // Shift register estimate, a pretty weak one

  int total_entries_per_cluster = ppu->fifo_bitwidth() * 2;
  float bits = total_entries_per_cluster * ppu->_num_clusters * ppu->_input_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->sram_Mb_per_mm2()) * REG_OVER_SRAM;  

  //printf("num_clusters %d, fifo_area %f, ",ppu->_num_clusters,area);

  return area;
}

float distribution_network::area() {
  //Crossbar -- this is a VERY agressive estimate, needs validation
  int crossbar_radix = ppu->_mem_banks * ppu->_compute_banks;
  float area = ppu->fifo_bitwidth() * crossbar_radix * _t->router_constant();

  // Distribution network estimate, needs validation 
  int dist_fanout=16;
  int bus_width=ppu->_num_clusters * ppu->_input_bitwidth;
  int stages = ceil(ppu->_num_clusters / (dist_fanout-1));
  int switches = bus_width * stages;
  int transistors = switches * 4;
  float dist_area = transistors / _t->mtr_per_mm2();

  //printf("dist_area %f\n",dist_area);
  area += dist_area;

  return area;
}

float coef_storage::area() {
  int total_entries = ppu->_coef_per_cluster * ppu->_num_flexible_clusters +
                                  ppu->num_point_clusters();
  float bits = total_entries * ppu->_coef_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->sram_Mb_per_mm2()) * REG_OVER_SRAM;
  return area;
}