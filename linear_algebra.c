#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "matrix_ADT.h"

//typedef struct Matrix *Matrix;

double dot_product(const struct Matrix *x, const struct Matrix *y, const int length){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));
  assert(sizeof_data(x) == sizeof_data(y));
  assert(length > 0);

  const int x_col = get_columns(x);
  const int x_row = get_rows(x);
  double re = 0;

  for(int i=0; i < length; ++i){
    int m = i%x_row;
    int n = i - m*x_col;

    re += (get_elem(x,m,n))*(get_elem(y,m,n));
  }
  return re;
}

// double magnitude(const double *x, const int length){
//   assert(x);
//   assert(length > 0);
//
//   double dot = dot_product(x, x, length);
//   return sqrt(dot);
// }
//
// bool is_unit_vector(const double *x, const int length){
//   assert(x);
//   assert(length > 0);
//
//   if(magnitude(x, length) == 1.0)
//     return 1;
//
//   return 0;
// }
//
// bool is_orthogonal(const double *x, const double *y, const int length){
//   assert(x);
//   assert(y);
//   assert(length > 0);
//
//   if(!(dot_product(x, y, length)))
//     return 1;
//
//   return 0;
// }
//
// double *cross_product(const double *x, const double *y){
//   assert(x);
//   assert(y);
//
//   double *cross = malloc(sizeof(double)*3);
//   cross[0] = (x[1]*y[2]) - (x[2]*y[1]);
//   cross[1] = (x[2]*y[0]) - (x[0]*y[2]);
//   cross[2] = (x[0]*y[1]) - (x[1]*y[0]);
//
//   return cross;
// }
//
// double *projection(const double *x, const double *y, const int length){
//   assert(x);
//   assert(y);
//   assert(length > 0);
//
//   double *re = malloc(sizeof(double)*length);
//   double dot = dot_product(x, y, length);
//   double mag = pow(magnitude(y, length), 2);
//
//   for(int i =0; i < length; ++i){
//     re[i] = (dot/mag)*y[i];
//   }
//
//   return re;
// }
//
// double *perpendicular(const double *x, const double *y, const int length){
//   assert(x);
//   assert(y);
//   assert(length > 0);
//
//   double *re = malloc(sizeof(double)*length);
//
//   double *proj = projection(x, y, length);
//
//   for(int i=0; i < length; ++i){
//     re[i] = x[i] - proj[i];
//   }
//
//   free(proj);
//   return re;
// }
