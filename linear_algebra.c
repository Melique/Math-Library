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
    re += data_x[i]*data_y[i];
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

static int arj_max(double *x, int a, int b, int ROWS, int COLS){
  assert(x);

  int track = a*COLS + b;
  double large = x[a*COLS + b];

  for(int i = a; i < ROWS; ++i){
    //printf("%f \n", x[i*COLS + b]);
    if(fabs(x[i*COLS + b]) > fabs(large)){
      large = x[i*COLS + b];
      track = i*COLS + b;
    }
  }
  // printf("MAX: %f\n", large);
  // printf("INDEX: %d\n", track);
  return track;
}


void swap(double *x, const int ROWS, const int COLS, int l, int m){
  assert(x);
  assert(l >= 0);
  assert(m >= 0);
  assert(l <= ROWS);
  assert(m <= ROWS);

  double *row_l = malloc(sizeof(double)*COLS);
  double *row_m = malloc(sizeof(double)*COLS);

  for(int i = 0; i < COLS; ++i){
    row_l[i] = x[l * COLS + i];
    row_m[i] = x[m * COLS + i];
  }

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      if(i == l){
        x[i*COLS + j] = row_m[j];
      }
      else if(i == m){
        x[i*COLS + j] = row_l[j];
      }
      else {
        x[i*COLS + j] = x[i*COLS + j];
      }
    }
  }

  free(row_l);
  free(row_m);

}

void GJE(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  const int COLS = get_columns(x);

  double *data = get_data(x);
  // printf("NO CHANGE: \n");
  // for(int p = 0; p < ROWS; ++p){
  //   for(int q = 0; q < COLS; ++q){
  //     printf("%f ", data[p*COLS + q]);
  //   }
  //   printf("\n");
  // }

  int row_track = 0;
  int col_track = 0;

  while(row_track < ROWS && col_track < COLS){ //dont know max{ROWS,COLS}
    int max_pivot = arj_max(data, row_track, col_track, ROWS, COLS); //try to find non-zerp
    if(max_pivot == 0){ // free var
      ++col_track;
    }

    swap(data, ROWS, COLS, (max_pivot/COLS), col_track);

    float temp_scalar = data[row_track*COLS];
    for(int i = col_track; i < COLS; ++i){
      data[row_track*COLS + i] = data[row_track*COLS + i] / temp_scalar;
    }

    // for(int j = row_track+1; j < ROWS; ++j){
    //   data[j*COLS + col_track] = 0;
    //
    // }

    for(int k = row_track + 1; k < ROWS; ++k){
      float shefali_scale = data[k*COLS + col_track];
      printf("SHEFALI SCALE: %f \n", shefali_scale);
      for(int l = col_track; l < COLS; ++l){
        data[k*COLS+l] = data[k*COLS+l] - shefali_scale*data[row_track*COLS+l];
      }
    }

    printf("\n");
    for(int p = 0; p < ROWS; ++p){
      for(int q = 0; q < COLS; ++q){
        printf("%f ", data[p*COLS + q]);
      }
      printf("\n");
    }
    ++row_track;
    ++col_track;
  }

}
