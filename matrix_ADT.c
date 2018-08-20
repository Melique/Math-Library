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

// typedef struct Matrix *Matrix;

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
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  free(x->data);
  free(x);
  x = NULL;
}

int get_row(const struct Matrix *x){
  assert(x);
  return x->rows;
}

int get_column(const struct Matrix *x){
  assert(x);
  return x->cols;
}

void print_matrix(const struct Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  const int ROWS = x->rows;
  const int COLS = x->cols;

  // for(int i = 0; i < ROWS*COLS; ++i){
  //   if(i%COLS == 0)
  //     printf("\n");
  //
  //   printf("%f ", x->data[i * COLS + (i%COLS)]);
  // }

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      printf("%f ", x->data[i *COLS + j]);
    }
    printf("\n");
  }

}

double get_elem(const struct Matrix *x, const int m, const int n){
  assert(x);
  assert(m > 0);
  assert(n > 0);
  assert(x->data);

  return x->data[m * x->cols + n];
}

struct Matrix *tranpose(const struct Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

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

struct Matrix *add_matrix(const struct Matrix *x, const struct Matrix *y){
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

struct Matrix *scalar_multiply(const struct Matrix *x, const double C){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

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


struct Matrix *matrix_multiply(const struct Matrix *x, const struct Matrix *y){
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
