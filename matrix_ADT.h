#include <stdio.h>
#include <stdlib.h>

// matrix_ADT.h [INTERFACE]

// This interface contains the matrix struct and its opertations

#ifndef MATRIX_H_
#define MATRIX_H_

struct Matrix;
// typedef struct Matrix *Matrix;

// matrix_create(m, n) creates an m x n martix
// effects: allocates memory (client must call matrix_destroy)
// requires:  m, n are postivie integers
struct Matrix *matrix_create(const int m, const int n, double *r);

int get_row(const struct Matrix *x);

int get_column(const struct Matrix *x);

// matrix_destroy(x) removes memory
// effects: x is no longer valid
void matrix_destroy(struct Matrix *x);

// print_matrix(x) outputs the x as a matrix
// effects: output
// requires: x is valid
// time: O(n^2) || O(n)
void print_matrix(const struct Matrix *x);

// tranpose(x) returns the tranpose of x
// effects: allocates memory (client must call matrix_destroy)
// requires: x is valid
struct Matrix *tranpose(const struct Matrix *x);

// add_matrix(x, y) returns the sum of x and y
// effects: allocates memory (client must call matrix_destroy)
// requires: x, y are valid
struct Matrix *add_matrix(const struct Matrix *x, const struct Matrix *y);

// scalar_multiply(x, C) returns scalar multiplication of C and x
// effects: allocates memory (client must call matrix_destroy)
// requires: x is valid
struct Matrix *scalar_multiply(const struct Matrix *x, const double C);

// matrix_multiply(x, y) returns matrix multiplication of x and y
// effects: allocates memory (client must call matrix_destroy)
// requires: x, y are valid
struct Matrix *matrix_multiply(const struct Matrix *x, const struct Matrix *y);

#endif
