#include "util.h"
#include "assert.h"
#include <cstdlib>
#include <iostream>

int rand_bt(int s, int e) {
  assert(e > s && "bad range for rand_bt");
  return rand() % (e - s) + s;
}

int rand_rng(int s, int e) {
  e++;
  return rand_bt(s,e);
}

void print_ge_compare(float records[200][7]){
  std::ofstream compare;
  compare.open("ge-compare.csv");
  // print header
  compare << "PPU-limit, CPU-GE-limit, CGRA-GE-limit, ASIC-GE-limit, "
          << "PPU-CPU-GE-limit, PPU-CGRA-GE-limit, PPU-ASIC-GE-limit\n";
  for(int row = 0; row < 200; row ++){
    for(int col = 0; col < 7; col ++){
      compare << records[row][col];
      if(col == 6){
        compare << "\n";
      }else{
        compare << ", ";
      }
    }
  }
}

std::ofstream print_ge_tradeoff(std::string filename){
    std::ofstream ge_tradeoff;
    ge_tradeoff.open(filename);
    std::cout << "Write GE tradeoff to " << filename <<std::endl;
    // Print the header
    ge_tradeoff //Channel Model
              << "channel-model, "
              // Global Parameter
              << "update-time, "
              // Object(on sky, reflector), Transmitter, Receiver, Path and Bands
              << "number-object, number-tx-per-band, number-rx-per-band, number-path, number-band, "
              // Ratio: #Object/#Tx
              << "ratio-obj-tx/rx, "
              // Object by speed
              << "number-fast-obj, number-slow-obj, number-fixed-obj, number-platforms, "
              // Wafer Level 
              << "num_wafer, num_ppu_chiplet, num_ge_chiplet, "
              // Clutter
              << "frac-clutter, "
              // Size
              << "avg-coef-per-obj, "
              // Range
              << "range, "
              // GE Total Area
              << "ge-total-area, "
              << "ge-compute-area, "
              << "ge-mem-area, "
              // GE fidelity
              << "ta1-upd-rate, "
              << "nr-engine-interpolation_ord, "
              << "nr-engine-conv, "
              << "antenna-order, "
              << "antenna-number, "
              << "antenna-resolution-angle, "
              << "antenna-dict-dim, "
              << "rcs-order, "
              << "rcs-points, "
              << "rcs-angle, "
              << "rcs-freq, "
              << "rcs-plzn, "
              // Coordinate
              << "coor-trans-compute, coor-trans-memory, coor-trans-bandwidth, coor-trans-latency, "
              // NR Engine
              << "nr-engine-compute, nr-engine-memory, nr-engine-bandwidth, nr-engine-latency, "
              // Relative Orientation
              << "relative-orientation-compute, relative-orientation-memory, relative-orientation-bandwidth, relative-orientation-latency, "
              // Antenna
              << "antenna-compute, antenna-memory, antenna-bandwidth, antenna-latency, "
              // Path gain and velocity
              << "path-gain-compute, path-gain-memory, path-gain-bandwidth, path-gain-latency, "
              // RCS
              << "rcs-compute, rcs-memory, rcs-bandwidth, rcs-latency, "
              // Update
              << "tu-compute, tu-memory, tu-bandwidth, tu-latency\n";
    return ge_tradeoff;
}