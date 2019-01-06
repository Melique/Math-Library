// linear_algebra.h [INTERFACE]

// This module includes algorithms and computational task done in Linear Algebra
// Note: All vectors are column vectors

#ifndef LINALG_H_
#define LINALG_H_


// dot_product(x, y) returns the dot product of x and y
// time: O(n)
double dot_product(const struct Matrix *x, const struct Matrix *y);

// magnitude(x) returns the vector length of x
// time: O(n)
double magnitude(const struct Matrix *x);

// is_unit_vector(x) returns true if x is a unit vector, false otherwise
// time: O(n)
bool is_unit_vector(const struct Matrix *x);

// is_orthogonal(x, y) returns true if x and y are orthogonal, false otherwise
// time: O(n)
bool is_orthogonal(const struct Matrix *x, const struct Matrix *y);

// cross_product(x, y) returns the cross product of x and y
// requires: x, y are in R^3
// effects: allocates memory (client must free memory)
struct Matrix *cross_product(const struct Matrix *x, const struct Matrix *y);

// projection(x, y) returns the projection of x onto y
// effects: allocates memory (client must free memory)
// time: O(n)
struct Matrix *projection(const struct Matrix *x, const struct Matrix *y);

// perpendicular(x, y) returns the perpendicular of x onto y
// effects: allocates memory (client must free memory)
// requires: x and y are valid and length > 0
// time: O(n)
struct Matrix *perpendicular(const struct Matrix *x, const struct Matrix *y);

// rank(x) returns the rank of x
// requires: x be in REF
// time: O(n)
int rank(struct Matrix *x);

// GE(x) changes x in the RREF x
// time:  O(mn)
void GE(struct Matrix *x);

// is_invertible(x) returns non-zero integer if x is invertible. 0 otherwise
// time: O(mn)
bool is_invertible(const struct Matrix *x);

// deter2(x) returns the determinant of a 2x2 matrix x
double det2(const struct Matrix *x);

// inverse(x) returns the inverse of x
// requires: x is invertible
// effects: allocates memory (client must free memory)
// time: O(mn)
void *inverse(const struct Matrix *x);

// deter(x) returns the determinant of a
// time:
double det(const struct Matrix *x);

// cofactor(x, i, j) returns the cofactor at row i or column j
// effects: allocates memory (client must free memory)
// time:
double cofactor(const struct Matrix *x, const int row, const int col);

// col(x) returns the column space of x
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *col(struct Matrix *x);

// null(x) returns the null space of x
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *null(const struct Matrix *x);

// row(x) returns the row space of x
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *row(const struct Matrix *x);

// null_T(x) returns the left null space of x
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *null_T(const struct Matrix *x);

// normalizer(x) normalizes x
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *normalizer(const struct Matrix *x);

// b_matrix(x, o_basis) returns the B Matrix of x
// requires:  o_basis is orthogonal basis
// effects: allocates memory (client must free memory)
// note: uses standard inner product (dot product).
//       the orthongoal basis should be read as size of of one vector by
//       the number of vectors
//       Ex. An orthogonal basis basis with 3 elements for M2X2(R). The matrix is
//       4x3
struct Matrix *b_matrix(const struct Matrix *x, const struct Matrix *y);

// GSP(x) returns an orthgonal basis
// requires: x is valid
// effects: allocates memory (client must free memory)
// note: standard inner product (dot product)
//       If the zero vector appears simply remove it
// time: O(mn)
struct Matrix *GSP(const struct Matrix *x);

// QR(x) returns a QR decompostion of x
// requires: dim Col(x) = n
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *QR(const struct Matrix *x);

// least_squares(x, y) applies the method of least squares to find a polynomial of
//  best fit
// effects: allocates memory (client must free memory)
// time: O(mn)
struct Matrix *least_sqaures(const struct Matrix *x, const struct Matrix *y);

#endif
