// linear_algebra.h [INTERFACE]

// This module includes algorithms and computational task done in Linear Algebra 1
// (MATH 136) offered at the Univeristy of Waterloo

#ifndef LINALG_H_
#define LINALG_H_


// dot_product(x, y, length) returns the dot product of x and y
// requires: x and y are valid and length > 0
// time: O(n)
double dot_product(const struct Matrix *x, const struct Matrix *y);

// magnitude(x, length) returns the vector length of x
// requires: x is valid and length > 0
// time: O(n)
double magnitude(const struct Matrix *x);

// is_unit_vector(x, length) returns true if x is a unit vector, false otherwise
// requires: x is valid and length > 0
// time: O(n)
bool is_unit_vector(const struct Matrix *x);

// is_orthogonal(x, y, length) returns true if x and y are orthogonal, false otherwise
// requires: x and y are valid and length > 0
// time: O(n)
bool is_orthogonal(const struct Matrix *x, const struct Matrix *y);

// cross_product(x, y) returns the cross product of x and y
// effects: allocates memory (client must free memory)
// requires: x, y are in R^3 and valid
struct Matrix *cross_product(const struct Matrix *x, const struct Matrix *y);

// projection(x, y, length) returns the projection of x onto y
// effects: allocates memory (client must free memory)
// requires: x and y are valid and length > 0
// time: O(n)
struct Matrix *projection(const struct Matrix *x, const struct Matrix *y);

// perpendicular(x, y, length) returns the perpendicular of x onto y
// effects: allocates memory (client must free memory)
// requires: x and y are valid and length > 0
// time: O(n)
struct Matrix *perpendicular(const struct Matrix *x, const struct Matrix *y);

void swap(double *x, const int ROWS, const int COLS, int l, int m);

int arj_max(double *x, int a, int b, int ROWS, int COLS);

// rank(x) returns the rank of x
// requires: x be a valid matrix. x be in REF
// time: O(n)
int rank(struct Matrix *x);

// GE(x) returns -1 if the system is incosistent
//       returns 0 is the system is consistent with a unqiue solution
//       returns 1 if the system is consistent with infinite solutions
// effects: x is changed
// requires: x be a valid matrix
// time:  O(4mn)
void GE(struct Matrix *x);

// is_invertible(x) returns non-zero integer if x is invertible. 0 otherwise
// requires: x is valid
// time: O(mn)
bool is_invertible(const struct Matrix *x);

// deter2(x) returns the determinant of a 2x2 matrix x
// requires: x is valid
double deter2(const struct Matrix *x);

// deter(x) returns the determinant of a
// requires: x is valid
// time:
double deter(const struct Matrix *x);

// cofactor(x, i, j) returns the cofactor at row i or column j
// requires: x is valid.
// effects: new matrix is created
// time:
float cofactor(const struct Matrix *x, const int row, const int col);
#endif
