#include "hw.h"



float input_tap_fifo::area() {
  int total_entries_per_cluster = ppu->_num_clusters * 2;
  float bits = total_entries_per_cluster * ppu->_num_clusters * ppu->_input_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->_sram_Mb_per_mm2);
  return area;
}

float coef_storage::area() {
  int total_entries_per_cluster = ppu->_coef_per_cluster * 4; // some arb. contstant
  float bits = total_entries_per_cluster * ppu->_num_clusters * ppu->_coef_bitwidth;
  float Mbits = bits/1024.0/1024.0;
  float area = Mbits * (1.0/_t->_sram_Mb_per_mm2);
  return area;
}

