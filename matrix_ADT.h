// matrix_ADT.h [INTERFACE]

// This interface contains the matrix struct and its opertations

struct Matrix;
typedef struct Matrix Matrix;

// matrix_create(m, n) creates an m x n martix
// effects: allocates memory (client must call matrix_destroy)
// requires:  m, n are postivie integers
Matrix *matrix_create(const int m, const int n);

// matrix_destroy(x) removes memory
// effects: x is no longer valid
void matrix_destroy(Matrix x);

// print_matrix(x) outputs the x as a matrix
// effects: output
// requires: x is valid
// time: O(n^2) || O(n)
void print_matrix(const Matrix *x);

// tranpose(x) returns the tranpose of x
// effects: allocates memory (client must call matrix_destroy)
// requires: x is valid
Matrix *tranpose(const Matrix *x);

// add_matrix(x, y) returns the sum of x and y
// effects: allocates memory (client must call matrix_destroy)
// requires: x, y are valid
Matrix *add_matrix(const Matrix *x, const Matrix *y);

// scalar_multiply(x, C) returns scalar multiplication of C and x
// effects: allocates memory (client must call matrix_destroy)
// requires: x is valid
Matrix *scalar_multiply(const Matrix *x, const double C);

// matrix_multiply(x, y) returns matrix multiplication of x and y
// effects: allocates memory (client must call matrix_destroy)
// requires: x, y are valid
Matrix *matrix_multiply(const Matrix *x, const Matrix *y);
