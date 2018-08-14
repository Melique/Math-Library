#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

struct Matrix{
  int rows;
  int cols;
  double *data;
}

typedef struct Matrix Matrix;

Matrix * matrix_create(const int m, const int n){
  assert(m > 0);
  assert(n > 0);

  Matrix *x = malloc(sizeof(Matrix));

  x->rows = m;
  x->cols = n;
  x->data = malloc(sizeof(double)*m*n);

  return x;
}

void matrix_destroy(Matrix x){
  assert(x);
  assert(x->data);

  free(x->data);
  free(x);
  x = NULL;
}

double dot_product(const double *x, const double *y, const int length){
  assert(x);
  assert(y);
  assert(length > 0);

  double re = 0;

  for(int i=0; i < length; ++i){
    re += (x[i]*y[i]);
  }
  return re;
}

double magnitude(const double *x, const int length){
  assert(x);
  assert(length > 0);

  double dot = dot_product(x, x, length);
  return sqrt(dot);
}

bool is_unit_vector(const double *x, const int length){
  assert(x);
  assert(length > 0);

  if(magnitude(x, length) == 1.0)
    return 1;

  return 0;
}

bool is_orthogonal(const double *x, const double *y, const int length){
  assert(x);
  assert(y);
  assert(length > 0);

  if(!(dot_product(x, y, length)))
    return 1;

  return 0;
}

double *cross_product(const double *x, const double *y){
  assert(x);
  assert(y);

  double *cross = malloc(sizeof(double)*3);
  cross[0] = (x[1]*y[2]) - (x[2]*y[1]);
  cross[1] = (x[2]*y[0]) - (x[0]*y[2]);
  cross[2] = (x[0]*y[1]) - (x[1]*y[0]);

  return cross;
}

double *projection(const double *x, const double *y, const int length){
  assert(u);
  assert(v);
  assert(length > 0);

  double *re = malloc(sizeof(double)*length);
  double dot = dot_product(x, y, length);
  double mag = pow(magnitude(y, length), 2);

  for(int i =0; i < length; ++i){
    re[i] = (dot/mag)*v[i];
  }

  return re;
}

double *perpendicular(const double *x, const double *y, const int length){
  assert(u);
  assert(v);
  assert(length > 0);

  double *re = malloc(sizeof(double)*length);

  double *proj = projection(x, y, length);

  for(int i=0; i < length; ++i){
    re[i] = x[i] - proj[i];
  }

  free(proj);
  return re;
}

bool is_RREF(const Matrix *x){
  assert(x);
  assert(x->rows > 0);
  assert(x->cols > 0);
  assert(x->data);

  int pivot[2] = {0,0};

  for(int i = 0; i < x->rows; ++i){

    //Find pivot
    for(int j =0; j < x->cols; ++j){
      //Non-zero and not 1
      if(x->data[i * x->cols + j]){
        if(x->data[i * x->cols + j] != 1){
          printf("\n%s", "Rule 1");
          return 0;
        }
        else{
          pivot[0] = i;
          pivot[1] = j;
        }
      }
    }

    for(int j = 0; j < x->cols; ++j){
      //Non-zero and not itself
      if(x->data[i * x->cols + pivot[1]] && (i != pivot[0])){
        printf("\n%s", "Rule 2");
        return 0;
      }

    }

    int track;
    int track_below;

    for(int j = 0; j < (x->cols); ++j){
      track += x->data[i * x->cols + j];
      track_below += x->data[(i+1) * x->cols + j];
    }

    if((track == 0) && track_below != 0){
      printf("\n%s", "Rule 3a");
      return 0;
    }

    else if(pivot[1] < pivot[0]){
      printf("\n%s", "Rule 3b");
      return 0;
    }
  }

    return 1;
  }

int main(){
  double x[2] = {2,3};
  double y[2] = {3,4};
  double z[3] = {1,3,-1};
  double v[3] = {0,2,1};

  double *chris = perpendicular(x,y,2);
  double *dole = perpendicular(v,z,3);

  Matrix *helen = matrix_construct(2,3);
  double test[12] = {1,0,7,0,0,0,0,0,0,1,0,1};
  helen->data = test;
  bool kelsey = is_RREF(helen);
  printf("\n%s", kelsey ? "true" : "false");

}
