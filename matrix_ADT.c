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
  int size;
};

struct Matrix *matrix_create(const int m, const int n, double *r){
  assert(m > 0);
  assert(n > 0);

  struct Matrix *x = (struct Matrix *) malloc(sizeof(struct Matrix));

  x->rows = m;
  x->cols = n;
  x->data = r;
  x->size = m*n;

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
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  return 1;
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

struct Matrix *clone(const struct Matrix *x){
  assert(x);
  assert(is_matrix_valid(x));

  const int ROWS = get_rows(x);
  const int COLS = get_columns(x);
  const double *DATA = get_data(x);
  double *r = malloc(sizeof(double)*ROWS*COLS);

  for(int i = 0; i < ROWS*COLS; ++i){
    r[i] = DATA[i];
  }

  struct Matrix *re = matrix_create(ROWS, COLS, r);

  return re;
}

int sizeof_data(const struct Matrix *x){
  assert(is_matrix_valid(x));

  return x->size;
}

void data_change(struct Matrix *x, double *r){
  assert(x);
  assert(r);

  x->data = r;

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

  const int ROWS = x->rows;
  const int COLS = x->cols;
  const int N = ROWS*COLS;

  double *data_t = malloc(sizeof(double)*N);

  for(int i = 0; i < N; ++i){
    int a = i/ROWS;
    int b = i%ROWS;

    data_t[i] = x->data[b*COLS + a];
  }

  struct Matrix *re = matrix_create(COLS, ROWS, data_t);

  return re;
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
  double *data = malloc(sizeof(double)*ROWS*COLS);
  const int length = x->cols;

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      for(int k = 0; k < length; ++k){
        data[i*COLS +j] += get_elem(x,i,k) * get_elem(y,k,j);
      }
    }
  }

  struct Matrix *product = matrix_create(ROWS, COLS, data);
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

bool equality(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(is_matrix_valid(x));
  assert(y);
  assert(is_matrix_valid(y));

  const int X_ROWS = get_rows(x);
  const int X_COLS = get_columns(x);
  const int Y_ROWS = get_rows(y);
  const int Y_COLS = get_columns(y);

  const double *X_DATA = get_data(x);
  const double *Y_DATA = get_data(y);


  if (X_ROWS != Y_ROWS || X_COLS != Y_COLS)return 0;

  for(int i = 0; i < X_ROWS; ++i){
    for(int j = 0; j < X_COLS; ++j){
      if(X_DATA[i*X_COLS + j] != Y_DATA[i*X_COLS + j]) return 0;
    }
  }

  return 1;
}

struct Matrix *merger(const struct Matrix *x, const struct Matrix *y){
  assert(x);
  assert(is_matrix_valid(x));
  assert(y);
  assert(is_matrix_valid(x));
  assert(get_rows(x) == get_columns(x));
  assert(get_rows(x) == get_rows(y));

  const int ROWS = get_rows(x);
  const double *data_x = get_data(x);
  const double *data_y = get_data(y);
  double *data = malloc(sizeof(double)*2*ROWS*ROWS);

  for(int i = 0; i < 2*ROWS; ++i){
    for(int j = 0; j < ROWS; ++j){
      if(i < ROWS){
        data[i*2*ROWS + j] = data_x[i*ROWS + j];
      } else {
        data[j*2*ROWS + i] = data_y[(i%ROWS)*ROWS + j];
      }
    }
  }

  struct Matrix *re = matrix_create(ROWS, 2*ROWS, data);
  return re;
}
