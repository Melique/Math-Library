#include <stdio.h>
#include <stdlib.h>

// matrix_ADT.h [INTERFACE]

// This interface contains the matrix struct and its opertations

#ifndef MATRIX_H_
#define MATRIX_H_

struct Matrix;

// matrix_create(m, n) creates an m x n martix
// effects: allocates memory (client must call matrix_destroy)
// requires:  m, n are postivie integers
struct Matrix *matrix_create(const int m, const int n, double *r);

// matrix_destroy(x) removes memory
// effects: x is no longer valid
void matrix_destroy(struct Matrix *x);

// is_matrix_valid(x) returns true of matrix is valid, false otherwise
bool is_matrix_valid(const struct Matrix *x);

// get_row(x) returns the number of rows in x
int get_rows(const struct Matrix *x);

// get_row(x) returns the number of rows in x
int get_columns(const struct Matrix *x);

// get_row(x) returns the data in x
double *get_data(const struct Matrix *x);

// get_elem(x, m, n) returns the element at postion m,n in x
double get_elem(const struct Matrix *x, const int m, const int n);

// clone(x) returns a copy of matrix x
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *clone(const struct Matrix *x);

// sizeof_data(x) returns the total number of elements
int sizeof_data(const struct Matrix *x);

// data_change(x, r) changes x->data to point a r
// effects: original data is lost
void data_change(struct Matrix *x, double *r);

// print_matrix(x) outputs the x as a matrix
// time: O(mn))
void print_matrix(const struct Matrix *x);

// tranpose(x) returns the tranpose of x
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *tranpose(const struct Matrix *x);

// add_matrix(x, y) returns the sum of x and y
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *add(const struct Matrix *x, const struct Matrix *y);

// scalar_multiply(x, C) returns scalar multiplication of C and x
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *smulti(const struct Matrix *x, const double C);

// matrix_multiply(x, y) returns matrix multiplication of x and y
// effects: allocates memory (client must call matrix_destroy)
// time: O(mnl)
struct Matrix *mmulti(const struct Matrix *x, const struct Matrix *y);

// scalar_multiply_row(x, row_num, C) multiplies row_num in x by C and returns
//  a new matrix
// effects: allocates memory (client must call matrix_destroy)
// requires: x is valid, row_num is a valid row, and C is an non-zero integer
// time: O(mn)
struct Matrix *smERO(const struct Matrix *x, const int row_num, const double C);

// addition_row(x, l, m, C) returns a new matrix of row_l + C*row_m
// effects: allocates memory (client must call matrix_destroy)
// requires: C is an integer
// time: O(mn)
struct Matrix *aERO(const struct Matrix *x, const int l, const int m, const double C);

// swap(x, l, m) swaps row l and row m and returns a new matrix
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *swapERO(const struct Matrix *x, const int l, const int m);

// equality(x,y) returns 1 if x and y are equal. 0 otherwise
// time: O(mn)
bool equality(const struct Matrix *x, const struct Matrix *y);
 // merger(x,y) combines x and y by extending the rows. The new matrix will be
//  m x 2*n
// effects: allocates memory (client must call matrix_destroy)
// time: O(mn)
struct Matrix *merger(const struct Matrix *x, const struct Matrix *y);

#endif
