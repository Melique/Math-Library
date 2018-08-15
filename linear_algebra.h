// linear_algebra.h [INTERFACE]

// This module includes all algorithms and computational task done in Linear Algebra 1
// (MATH 136) offered at the Univeristy of Waterloo

// dot_product(x, y, length) returns the dot product of x and y
// requires: x and y are valid and length > 0
// time: O(n)
double dot_product(const double *x, const double *y, const int length);

// magnitude(x, length) returns the vector length of x
// requires: x is valid and length > 0
// time: O(n)
double magnitude(const double *x, const int length);

// is_unit_vector(x, length) returns true if x is a unit vector, false otherwise
// requires: x is valid and length > 0
// time: O(n)
bool is_unit_vector(const double *x, const int length);

// is_orthogonal(x, y, length) returns true if x and y are orthogonal, false otherwise
// requires: x and y are valid and length > 0
// time: O(n)
bool is_orthogonal(const double *x, const double *y, const int length);

// cross_product(x, y) returns the cross product of x and y
// effects: allocates memory (client must free memory)
// requires: x, y are in R^3 and valid
double *cross_product(const double *x, const double *y);

// projection(x, y, length) returns the projection of x onto y
// effects: allocates memory (client must free memory)
// requires: x and y are valid and length > 0
// time: O(n)
double *projection(const double *u, const double *v, const int length)

// perpendicular(x, y, length) returns the perpendicular of x onto y
// effects: allocates memory (client must free memory)
// requires: x and y are valid and length > 0
// time: O(n)
double *perpendicular(const double *x, const double *y, const int length);

// is_RREF(x) returns true if x is in RREF, false otherwise
// requires: x is valid
// time: O(n^2)
bool is_RREF(const Matrix *x);
