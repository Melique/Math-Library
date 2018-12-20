// algebra.h [INTERFACE]

// This module includes algorithms and computational task done in Algebra
// (MATH 135) offered at the University of Waterloo

#ifndef ALG_H_
#define ALG_H_

// is_prime(n) determines if n is a prime number
// requires: n > 1
// time: O(√n)
bool is_prime(const int n);

// division_alg(a, b) returns the quotient and remainder when a is divided by b
// requires: b > 0 and caller to free the malloc memory
// time: O(1)
int* division_alg(const int a, const int b);

// EE(a, b) returns the gcd(a,b) by the Euclidean Algorithm
// time: O(n)
int EE(const int a, const int b);

// EEA(a, b) returns x, y, d=gcd(a,b), and q[uotient] of the Linear Diophantine
// equation ax + by = c by the Extended Euclidean Algorithm
// requires: a,b > 0, the caller frees memory
// effects: a, b could swap
// time: O(n)
int* EEA(int a, int b);

// PF(a) returns the Prime Factorization of a
// requires: a > 1
// effects: memory is dyanmically allocated, caller frees memory
// time: O(m√n)
int** PF(const int a);

// RSA_setUp(p, q) returns the private and public key
// requires: p and q are prime, caller has to free memory
// time: O(n)
int** RSA_setUp(const int p, const int q);

#endif
