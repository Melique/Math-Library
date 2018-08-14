// calculus.h [INTERFACE]

// This module includes all algorithms taught in MATH 137 and 138 offered at the
// University of Waterloo

// heron(alpha,a, n) returns the square root of alpha by the Heron Algorithm
//  where a is the starting value and n is the number of iterations
// requires: alpha and a be real numbers > 0 and n be non-negative
//  integer > 0
// time: O(n)
double heron(const double alpha, double a, int n);

// bisection(epi, a, b, fp) returns an approximation of a root of fp by the bisection
//  method
// requires: fp is a continuous function, EPI > 0, a < b, and fp(a)fp(b) < 0
// time: O(logn)
double bisection(double a, double b, double (*fp)(double));

// newton(EPI, x, fp, fprimep) uses Newton's method to find approximations such as
//  roots of functions and square roots
// requires: x is close to point c such that fp(c)=0 (IVT) and fp is differentiable
//  a x
double newton(double x, const double (*fp)(double), const double (*fprimep)(double));
