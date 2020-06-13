#include "util.h"
#include "assert.h"
#include <cstdlib>

int rand_bt(int s, int e) {
  assert(e > s && "bad range for rand_bt");
  return rand() % (e - s) + s;
}

int rand_rng(int s, int e) {
  e++;
  return rand_bt(s,e);
}

