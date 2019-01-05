#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "matrix_ADT.h"

extern double det(const struct Matrix *x);
extern double cofactor(const struct Matrix *x, const int row, const int col);
// TODO: design for linear mappings and vector spaces
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

// static int arj_max(double *x, int a, int b, int ROWS, int COLS){
//   assert(x);
//
//   int track = a*COLS + b;
//   double large = x[a*COLS + b];
//
//   for(int i = a; i < ROWS; ++i){
//     //printf("%f \n", x[i*COLS + b]);
//     if(fabs(x[i*COLS + b]) > fabs(large)){
//       large = x[i*COLS + b];
//       track = i*COLS + b;
//     }
//   }
//   printf("INDEX: %d \nVALUE: %f", track, large);
//   return track;
// }
//
//
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

int rank(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const double *data = get_data(x);
  const int ROWS = get_rows(x);
  const int COLS = get_columns(x);
  int track = 0;

  for(int i = 0; i < ROWS; ++i){
    for(int j = i; j < COLS; ++j){
      if(data[i*COLS + j] != 0){
        ++track;
        break;
      }
    }
  }

  return track;
}

void GE(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  const int COLS = get_columns(x);

  double *data = get_data(x);
  int result;

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      if(i != j){
        double scalar = (data[i*COLS+i] != 0) ? data[j*COLS +i]/data[i*COLS+i]: 0;
        for(int k = 0; k < COLS; ++k){
          if(data[j*COLS + k] == 0){ //whats this?

          }
          data[j*COLS+k] = data[j*COLS+k] - scalar*data[i*COLS+k];
        }
      }
    }
  }

  for(int i = 0; i < ROWS; ++i){
    double scalar = data[i*COLS + i];
    if (scalar){
      for(int j = 0; j < COLS; ++j){
        data[i*COLS + j] /= scalar;
      }
    }
  }

  for(int i = 0; i < ROWS; ++i){
    int track = 0;
      for(int j = 0; j < COLS; ++j){
        if(data[i*COLS + j] == 0) {
          ++track;
        }

        if(track == COLS) swap(data,ROWS, COLS, i, ROWS-1);

      }
    }
}

bool is_invertible(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  struct Matrix *c = clone(x);
  GE(c);

  bool result = rank(c) == ROWS;
  free(c);
  return result;
}

double det2(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));
  assert(get_rows(x) == 2 && get_columns(x) == 2);

  const double *data = get_data(x);

  double result = data[0]*data[3] - data[1]*data[2];
  return result;
}

struct Matrix *inverse(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));
  assert(is_invertible(x));

  const int ROWS = get_rows(x);
  const double *o_data = get_data(x);
  struct Matrix *re;;

  if(ROWS == 2){
    const double det = det2(x); //det won't be zero
    double *data = malloc(sizeof(data)*ROWS*ROWS);
    data[0] = (1/det)*o_data[3];
    data[1] = (1/det)*-1*o_data[1];
    data[2] = (1/det)*-1*o_data[2];
    data[3] = (1/det)*o_data[0];

    re = matrix_create(ROWS, ROWS, data);
  } else {
    double *i_data = malloc(sizeof(double)*ROWS*ROWS);

    for(int i = 0; i < ROWS; ++i){
        for(int j = 0; j < ROWS; ++j){
          i_data[i*ROWS + i] = 1;
        }
    }
    struct Matrix *I = matrix_create(ROWS, ROWS, i_data);

    struct Matrix *aug = merger(x, I);
    GE(aug);
    matrix_destroy(I);
    double *aug_data = get_data(aug);

    for(int i = 0; i < ROWS; ++i){
      double scalar = aug_data[i*2*ROWS + i];
      for(int j = 0; j < 2*ROWS; ++j){
        aug_data[i*2*ROWS + j] /= scalar;
      }
    }

    double *data = malloc(sizeof(data)*ROWS*ROWS);
    for(int i = 0; i < ROWS; ++i){
      for(int j = 0, k = ROWS; j < ROWS; ++j, ++k){
        data[i*ROWS + j] = aug_data[i*2*ROWS + k];
      }
    }

    re = matrix_create(ROWS, ROWS, data);
  }
  return re;
}


double *partial_clone(const double *x, const int ROWS, int row, int col){
  assert(x);
  assert((0 <= row) && (row < ROWS));
  assert((0 <= col) && (col < ROWS));

  double *re = malloc(sizeof(double)*(ROWS-1)*(ROWS-1));
  int track = 0;

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < ROWS; ++j){
      if((i != row) && (j != col)){
      re[track] = x[i*ROWS + j];
      ++track;
    }
  }
}
  return re;
}

double cofactor(const struct Matrix *x, const int row, const int col){
  assert(x);
  assert(is_matrix_valid(x));
  assert(get_rows(x) == get_columns(x));
  assert(get_rows(x) >= 2);

  const int ROWS = get_rows(x);
  double re = 0;

  if(ROWS == 2){
    int exp = pow(-1, row+col);
    re = pow(-1, row+col)*det2(x);
  }else{
    double *partial = get_data(x);
    partial = partial_clone(partial, ROWS, row, col);
    struct Matrix *y = matrix_create(ROWS-1, ROWS-1, partial);
    re += pow(-1, row+col)*det(y);
    matrix_destroy(y);
  }

  return re;
}

double det(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));
  assert(get_rows(x) == get_columns(x));

  const int ROWS = get_rows(x);
  double *partial = get_data(x);
  double re;

  if(ROWS == 1){
    return partial[0];
  } else if (ROWS == 2){
    re += det2(x);
  } else {
    for(int i = 0; i < ROWS; ++i){
      double scalar = get_elem(x,i,0);
      double factor = cofactor(x,i,0);
      re += scalar*factor;
    }
  }

  return re;
}

static int *lead_free(const struct Matrix *x, bool y){
  assert(x);
  assert(is_matrix_valid(x));

  int track  = 0;
  int num = 0;
  const int ROWS = get_rows(x);

  for(int i = 0; i < ROWS; ++i){
    if(y){
      if(get_elem(x,i,i) != 0){
        ++num;
      }
    } else{
      if(get_elem(x,i,i) == 0){
        ++num;
      }
    }
  }

  int *indices = malloc(sizeof(int)*num+1);

  for(int i = 0; i < ROWS; ++i){
    if(y){
      if(get_elem(x,i,i) != 0){
        indices[track] = i;
        ++track;
      }
    }else{
      if(get_elem(x,i,i) == 0){
        indices[track] = i;
        ++track;
      }
    }
  }
  indices[-1] = num; //length
  return indices;
}

struct Matrix *col(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  const int COL = get_columns(x);

  struct Matrix *cl = clone(x);
  GE(cl);

  int *indices = lead_free(cl,1);
  const int LENGTH = indices[-1];

  double *data = malloc(sizeof(double)*LENGTH*ROWS);

  for(int i = 0; i < LENGTH; ++i){
    for(int j = 0; j < ROWS; ++j){
      data[j*LENGTH + i] = get_elem(x,j,indices[i]);
    }
  }

  struct Matrix *re = matrix_create(ROWS, LENGTH, data);
  return re;
}

struct Matrix *null(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int COLS = get_columns(x);

  struct Matrix *cl = clone(x);
  GE(cl);

  int *indices = lead_free(cl,0);
  const int LENGTH = indices[-1];
  int track = 0;

  double *data = malloc(sizeof(double)*LENGTH*COLS);

  for(int i = 0; i < LENGTH; ++i){
    for(int j = 0; j < COLS; ++j){
      if(j == indices[track]){
        data[j*LENGTH + i] = 1;
      }else{
        data[j*LENGTH + i] = -1*get_elem(cl, j,indices[track]);
      }
    }
    ++track;
  }

  struct Matrix *re = matrix_create(COLS, LENGTH, data);
  return re;


}
struct Matrix *row(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  //const int ROWS = get_rows(x);
  const int COLS = get_columns(x);

  struct Matrix *cl = clone(x);
  GE(cl);

  int *indices = lead_free(cl,1);
  const int LENGTH = indices[-1];
  int track = 0;

  double *data = malloc(sizeof(double)*LENGTH*COLS);

  for(int i = 0; i < LENGTH; ++i){
    for(int j = 0; j < COLS; ++j){
      data[j*LENGTH + i] = get_elem(cl,indices[track],j);
    }
    ++track;
  }

  struct Matrix *re = matrix_create(COLS, LENGTH, data);
  return re;
}

struct Matrix *null_T(struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  struct Matrix *cl = tranpose(x);
  cl = null(cl);

  return cl;
}

static double *get_vector(const struct Matrix *x, const int col){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  const int COLS = get_columns(x);

  double *vector = malloc(sizeof(double)*ROWS);
  for(int i = 0; i < ROWS; ++i){
    vector[i] = get_elem(x,i,col);
  }

  return vector;

}

static double cof(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(y));

  double dot = dot_product(x, y);
  double mag = magnitude(y)*magnitude(y);

  return dot/mag;
}

struct Matrix *b_matrix(const struct Matrix *x, const struct Matrix *o_basis){
  assert(x);
  assert(o_basis);
  assert(is_matrix_valid(x));
  assert(is_matrix_valid(o_basis));

  const int ROWS = get_rows(o_basis);
  const int COLS = get_columns(o_basis);
  const int SIZE = sizeof_data(x);
  const double *x_data = get_data(x);

  double *cor_data = malloc(sizeof(double)*COLS);

  for(int i = 0; i < COLS; ++i){
    double *v_data = get_vector(o_basis, i);

    struct Matrix *vector = matrix_create(ROWS,1,v_data);

    cor_data[i] = cof(x, vector);
    matrix_destroy(vector);
  }

  struct Matrix *re = matrix_create(COLS, 1, cor_data);
  return re;
}

// struct Martix *QR(const struct Matrix *x){
//   assert(x);
//   assert(is_matrix_valid(x));
//
//   const int ROWS = get_rows(x);
//   const int COLS = get_columns(x);
//
//   struct Matrix *q = GSP(x);
//   double *r_data = malloc(sizoeof(double)*ROWS*COLS);
//
//   for(int i = 0; i < ROWS; ++i){
//     for(int j = 0; j < COLS; ++j){
//       if(j < i){
//         double *x_data = get_vector(x,j);
//         double *q_data = get_vector(q, i);
//         struct Matrix *x_vector = matrix_create(ROWS, 1, x_data);
//         struct Matrix *q_vector = matrix_create(ROWS, 1, q_data);
//         r_data[i*COLS + j] = dot_product)=(x_vector, q_vector);
//         matrix_destroy(x_data);
//         matrix_destroy(q_data);
//       }else{
//         r_data[i*COLS + j] = 0;
//       }
//     }
//   }
//
//   struct Matrix *R = matrix_create(ROWS, COLS, r_data);
//   return R;
// }

struct Matrix *least_sqaures(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(is_matrix_valid(x));
  assert(y);
  assert(is_matrix_valid(y));

  struct Matrix *x_tranpose = tranpose(x);
  struct Matrix *product = mmulti(x_tranpose, x);
  struct Matrix *inverse_product = inverse(product);
  struct Matrix *y_product = mmulti(x_tranpose, y);

  struct Matrix *re = mmulti(inverse_product, y_product);

  return re;
}
