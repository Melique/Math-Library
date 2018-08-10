#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

double heron(const double alpha, double a, const int n){
  assert(alpha > 0);

  if(n ==0 ) return a;
  else{
    double a = 1/2(a + alpha/a);
    heron(alpha, a, --n);
  }

}
