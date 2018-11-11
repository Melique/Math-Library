#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include "matrix_ADT.h"

struct Matrix{
  int rows;
  int cols;
  double *data;
};

struct Matrix *matrix_create(const int m, const int n, double *r){
  assert(m > 0);
  assert(n > 0);

  struct Matrix *x = (struct Matrix *) malloc(sizeof(struct Matrix));

  x->rows = m;
  x->cols = n;
  x->data = r;

  return x;
}

void matrix_destroy(struct Matrix *x){
  assert(x);
  assert(x->data);

  //free(x->data); Was never being allocated
  free(x);
  x = NULL;
}

bool is_matrix_valid(const struct Matrix *x){
  return (x && (x->rows > 0) && (x->cols > 0) && x->data);
}

int get_rows(const struct Matrix *x){
  assert(is_matrix_valid(x));
  return x->rows;
}

int get_columns(const struct Matrix *x){
  assert(is_matrix_valid(x));
  return x->cols;
}

double *get_data(const struct Matrix *x){
  assert(is_matrix_valid(x));

  return x->data;
}

double get_elem(const struct Matrix *x, const int m, const int n){
  assert(x);
  assert(m >= 0);
  assert(n >=  0);
  assert(x->data);

  return x->data[m * x->cols + n];
}

int sizeof_data(const struct Matrix *x){
  assert(is_matrix_valid(x));

  return sizeof(x->data)/sizeof(x->data[0]);
}

void print_matrix(const struct Matrix *x){
  assert(is_matrix_valid(x));

  const int ROWS = x->rows;
  const int COLS = x->cols;

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      printf("%f ", x->data[i *COLS + j]);
    }
    printf("\n");
  }

}

struct Matrix *tranpose(const struct Matrix *x){
  assert(is_matrix_valid(x));

  const int ROWS = x->cols;
  const int COLS = x->rows;

  struct Matrix *trans_x = matrix_create(ROWS, COLS, x->data);

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      trans_x->data[i * COLS + j] = x->data[j * ROWS + i];
    }
  }

  return trans_x;
}

struct Matrix *add(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(x->rows == y->rows);
  assert(x->cols == y->cols);
  assert(x->rows > 0);
  assert(x->cols > 0);

  struct Matrix *re = matrix_create(x->rows, x->cols, x->data);

  const int rows = re->rows;
  const int cols = re->cols;

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      re->data[i * cols + j] = x->data[i * cols + j] + y->data[i * cols + j];
    }
  }

  return re;
}

struct Matrix *smulti(const struct Matrix *x, const double C){
  assert(is_matrix_valid(x));

  struct Matrix *re = matrix_create(x->rows, x->cols, x->data);
  const int rows = re->rows;
  const int cols = re->cols;

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      re->data[i * cols + j] = C*x->data[i * cols + j];
    }
  }

  return re;
}


struct Matrix *mmulti(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(y);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->cols == y->rows);

  const int ROWS = x->rows;
  const int COLS = y->cols;
  struct Matrix *product = matrix_create(ROWS, COLS, x->data);
  const int length = x->cols;

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      for(int k = 0; k < length; ++k){
        product->data[i * COLS + j] += x->data[i * COLS + k]*y->data[k*COLS + j];
      }
    }
  }
  return product;
}

struct Matrix *smERO(const struct Matrix *x, const int row_num, const double C) {
  assert(is_matrix_valid(x));
  assert(0 < row_num <= x->rows);
  assert(C);

  const int ROWS = x->rows;
  const int COLS = x->cols;
  struct Matrix *re = matrix_create(ROWS, COLS, x->data);

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      if(i == row_num){
        re->data[i * COLS + j] = C*x->data[i * COLS + j];
      }
      else{
        re->data[i * COLS + j] = C*x->data[i * COLS + j];
      }
    }
  }

  return re;
}

struct Matrix *aERO(const struct Matrix *x, const int l, const int m, const double C){
  assert(is_matrix_valid(x));
  assert(0 < l <= x->rows);
  assert(0 < m <= x->rows);

  const int ROWS = x->rows;
  const int COLS = x->cols;
  struct Matrix *re = matrix_create(ROWS, COLS, x->data);

  double *row_l = malloc(sizeof(double)*COLS);

  for(int i = 0; i < COLS; ++i){
    row_l[i] = x->data[(l-1)*COLS + i] + C*x->data[(m-1)*COLS + i];

  }

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      if(i == (l-1)){
        re->data[i * COLS + j] = row_l[j];
      }
      else{
        re->data[i * COLS + j] = x->data[i * COLS + j];
      }
    }
  }
  free(row_l);
  return re;
}

struct Matrix *swapERO(const struct Matrix *x, const int l, const int m){
  assert(is_matrix_valid(x));
  assert(0 < l <= x->rows);
  assert(0 < m <= x->rows);

  const int ROWS = x->rows;
  const int COLS = x->cols;

  double *row_l = malloc(sizeof(double)*COLS);
  double *row_m = malloc(sizeof(double)*COLS);

  struct Matrix *re = matrix_create(ROWS, COLS, x->data);

  for(int i = 0; i < COLS; ++i){
    row_l[i] = x->data[l * COLS + i];
    row_m[i] = x->data[m * COLS + i];
  }

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      if(i == l){
        re->data[i * COLS + j] = row_l[j];
      }
      else if(i == m){
        re->data[i * COLS + j] = row_m[j];
      }
      else{
        re->data[i * COLS + j] = x->data[i * COLS + j];
      }
    }
  }

  free(row_l);
  free(row_m);

  return re;
}
