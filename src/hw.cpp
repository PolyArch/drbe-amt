#include "hw.h"

float input_tap_fifo::area() {

  // Shift register estimate, a pretty weak one
  uint64_t total_entries_per_cluster = ppu->_num_clusters * 2;
  float bits = total_entries_per_cluster * ppu->_num_clusters * ppu->_input_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->sram_Mb_per_mm2());

  //printf("num_clusters %d, fifo_area %f, ",ppu->_num_clusters,area);


  // Distribution network estimate, also a pretty weak one
  int dist_fanout=16;
  int bus_width=ppu->_num_clusters * ppu->_input_bitwidth;
  int stages = ceil(ppu->_num_clusters / (dist_fanout-1));
  int switches = bus_width * stages;
  int transistors = switches * 4;
  float dist_area = transistors /_t->mtr_per_mm2();

  //printf("dist_area %f\n",dist_area);

  area += dist_area;
  return area;
}

float coef_storage::area() {
  int total_entries_per_cluster = ppu->_coef_per_cluster * 4; // some arb. contstant
  float bits = total_entries_per_cluster * ppu->_num_clusters * ppu->_coef_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->sram_Mb_per_mm2());
  return area;
}

