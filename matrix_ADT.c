#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <string.h>

struct Matrix{
  int rows;
  int cols;
  double *data;
};

typedef struct Matrix Matrix;

Matrix *matrix_create(const int m, const int n){
  assert(m > 0);
  assert(n > 0);

  Matrix *x = malloc(sizeof(Matrix));

  x->rows = m;
  x->cols = n;
  x->data = malloc(sizeof(double)*m*n);

  return x;
}

void matrix_destroy(Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  free(x->data);
  free(x);
  x = NULL;
}

void print_matrix(const Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  const int ROWS = x->rows;
  const int COLS = x->cols;

  for(int i = 0; i < ROWS*COLS; ++i){
    int track;

    if(track == COLS){
      track = 0;
      printf("\n");
    }

    printf("%f\n", x->data[i]);

  }
  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      printf("%f ", x->data[i *COLS + j]);
    }
    printf("\n");
  }
}

Matrix *tranpose(const Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  const int ROWS = x->cols;
  const int COLS = x->rows;

  Matrix *trans_x = matrix_create(ROWS, COLS);


  if((ROWS == 0) || (COLS == 0)){
    trans_x->data = x->data;
    return trans_x;
  }

  for(int i = 0; i < ROWS; ++i){
    for(int j = 0; j < COLS; ++j){
      trans_x->data[i * COLS + j] = x->data[j * COLS + i];
    }
  }

  return trans_x;
}

Matrix *add_matrix(const Matrix *x, const Matrix *y){
  assert(x);
  assert(y);
  assert(x->rows == y->rows);
  assert(x->cols == y->cols);
  assert(x->rows > 0);
  assert(x->cols > 0);

  Matrix *re = matrix_create(x->rows, x->cols);

  const int rows = re->rows;
  const int cols = re->cols;

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      re->data[i * cols + j] = x->data[i * cols + j] + y->data[i * cols + j];
    }
  }

  return re;
}

Matrix *scalar_multiply(const Matrix *x, const double C){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  Matrix *re = matrix_create(x->rows, x->cols);
  const int rows = re->rows;
  const int cols = re->cols;

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      re->data[i * cols + j] = C*x->data[i * cols + j];
    }
  }

  return re;
}

Matrix *matrix_multiply(const Matrix *x, const Matrix *y){
  assert(x);
  assert(y);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->cols == y->rows);

  const int rows = x->rows;
  const int cols = y->cols;
  Matrix *re = matrix_create(rows, cols);

  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      re->data[i * cols + j] += ((x->data[i * cols + j]) * (y->data[j *rows+i]));
    }
  }

  return re;
}

int main(){
  double a[6] = {3,4,-5,1,0,2};
  double b[3] = {2,-1,6};
  double c[2] = {1,2};
  double d[4] = {1,2,3,4};
  double e[6] = {1,2,3,4,5,6};
  double f[9] = {1,3,1,-4,0,2,5,9,-3};

  Matrix *dole = matrix_create(2,3);
  dole->data = a;
  Matrix *chris = matrix_create(3,1);
  chris->data = b;
  Matrix *pavi = matrix_create(2,1);
  pavi->data = c;
  Matrix *kabi = matrix_create(2,2);
  kabi->data = d;
  Matrix *kyara = matrix_create(3,2);
  kyara->data = e;
  Matrix *rachel = matrix_create(3,3);
  rachel->data = f;

  Matrix *dole_three = scalar_multiply(dole, 3);
  Matrix *chris_T = tranpose(chris);
  Matrix *pavi_T = tranpose(pavi);
  Matrix *kabi_T = tranpose(kabi);
  Matrix *kyara_T = tranpose(kyara);
  Matrix *dole_chris = matrix_multiply(dole, chris);

  printf("%d %d\n", dole_chris->rows, dole_chris->cols);
  print_matrix(chris);
  printf("\n");
  print_matrix(chris_T);
  printf("\n");
  print_matrix(pavi);
  printf("\n");
  print_matrix(pavi_T);




  // matrix_destroy(dole);
  // matrix_destroy(chris);
  // matrix_destroy(pavi);
  // matrix_destroy(kabi);
  // matrix_destroy(kyara);
  // matrix_destroy(rachel);
  // matrix_destroy(dole_chris);
  // matrix_destroy(dole_three);
  // matrix_destroy(chris_two);
  // matrix_destroy(pavi_T);
  // matrix_destroy(kabi_T);
  // matrix_destroy(kyara_T);
  // matrix_destroy(rachel_T);


}
