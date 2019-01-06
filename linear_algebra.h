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

// inverse(x) returns the inverse of x
// requires: x is valid and invertible
// effects: new matrix created
// time:
void *inverse(const struct Matrix *x);

// deter2(x) returns the determinant of a 2x2 matrix x
// requires: x is valid
double det2(const struct Matrix *x);

// deter(x) returns the determinant of a
// requires: x is valid
// time:
double det(const struct Matrix *x);

// cofactor(x, i, j) returns the cofactor at row i or column j
// requires: x is valid.
// effects: new matrix is created
// time:
double cofactor(const struct Matrix *x, const int row, const int col);

// col(x) returns the column space of x
// requires: x is valid
// effects: new matrix allocated
struct Matrix *col(struct Matrix *x);

// null(x) returns the null space of x
// requires: x is valid
// effects: new matrix allocated
struct Matrix *null(const struct Matrix *x);

// row(x) returns the row space of x
// requires: x is valid
// effects: new matrix allocated
struct Matrix *row(const struct Matrix *x);

// null_T(x) returns the left null space of x
// requires: x is valid
// effects: new matrix allocated
struct Matrix *null_T(const struct Matrix *x);

// b_matrix(x, o_basis) returns the B Matrix of x
// requires: x, o_basis are valid and o_basis is orthogonal basis
// effects: new matrix allocated
// note: uses standard inner product (dot product).
//       the orthongoal basis should be read as size of of one vector by
//       the number of vectors
//       Ex. An orthogonal basis basis with 3 elements for M2X2(R). The matrix is
//       4x3
struct Matrix *b_matrix(const struct Matrix *x, const struct Matrix *y);

// GSP(x) returns an orthgonal basis
// requires: x is valid
// effetcs: new matrix is allocated
// note: standard inner product (dot product)
//       If the zero vector appears simply remove it 
struct Matrix *GSP(const struct Matrix *x);

// GSP(x) returns a QR decompostion of x
// requires: x is valid. dim Col(x) = n
// effetcs: new matrix is allocated
struct Matrix *QR(const struct Matrix *x);

// least_squares(x, y) applies the method of least squares to find a polynomial of
//  best fit
// requires: x and y are valid
struct Matrix *least_sqaures(const struct Matrix *x, const struct Matrix *y);

#endif
