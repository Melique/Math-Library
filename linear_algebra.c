#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "matrix_ADT.h"

//typedef struct Matrix *Matrix;

double dot_product(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));
  assert(sizeof_data(x) == sizeof_data(y));
  assert((get_rows(x) == 1 && get_rows(y) == 1) ||
         (get_columns(x) == 1) && (get_columns(y) == 1));

  const int col = get_columns(x);
  const int row = get_rows(x);
  const double *data_x = get_data(x);
  const double *data_y = get_data(y);
  double re = 0;

  for(int i=0; i < col*row; ++i){
    //printf("x: %f, y: %f\n", data_x[i], data_y[i]);


    // int m = i%row;
    // int n = i/row + 1;

    re += data_x[i]*data_y[i]; //(get_elem(x,m,n))*(get_elem(y,m,n));
  }
  return re;
}

double magnitude(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  double dot = dot_product(x, x);
  return sqrt(dot);
}

bool is_unit_vector(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));
  assert(get_rows(x) == 1 || get_columns(x) == 1); //could be a column or row vector

  if (magnitude(x) == 1) {
    return 1;
  }

  return 0;
}


bool is_orthogonal(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));

  if(!dot_product(x,y)){
    return 1;
  }

  return 0;
}


struct Matrix *cross_product(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));
  assert(sizeof_data(x) == 3 && sizeof_data(y) == 3);
  //Ensuring the vector is in R^3
  assert((get_rows(x) == 3 && get_rows(y) == 3) ||
         (get_columns(x) == 3) && (get_columns(y) == 3));


  const int length = get_rows(x)*get_columns(x);

  const double *data_x =  get_data(x);
  const double *data_y =  get_data(y);

  double *data = malloc(sizeof(double)*length);

  double temp1 = (data_x[1]*data_y[2]) - (data_x[2]*data_y[1]);
  double temp2 = (data_x[2]*data_y[0]) - (data_x[0]*data_y[2]);
  double temp3 = (data_x[0]*data_y[1]) - (data_x[1]*data_y[0]);

  data[0] = temp1;
  data[1] = temp2;
  data[2] = temp3;

  struct Matrix *re = matrix_create(3,1,data);

  return re;
}

struct Matrix *projection(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));

  const int length = get_rows(x)*get_columns(x);

  double *data = malloc(sizeof(double)*length);
  double *data_y = get_data(y);
  double dot = dot_product(x, y);
  double mag = pow(magnitude(y), 2);

  for(int i =0; i < length; ++i){
    data[i] = (dot/mag)*data_y[i];
  }

  struct Matrix *re = matrix_create(get_rows(x), get_columns(y), data);

  return re;

}

struct Matrix *perpendicular(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));

  const int length = get_rows(x)*get_columns(x);
  const double *data_x = get_data(x);
  struct Matrix *proj = projection(x, y);

  double *data = malloc(sizeof(double)*length);

  for(int i=0; i < length; ++i){
    data[i] = data_x[i] - get_data(proj)[i];
  }

  struct Matrix *re = matrix_create(get_rows(x), get_columns(y), data);

  free(proj);
  return re;
}
