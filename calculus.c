#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#define E 2.71828

double heron(const double alpha, double a, int n){
  assert(alpha > 0);
  assert(a > 0);
  assert(n > 0);

  if(n == 0){
    return a;
  }
  else{
      a = 0.5*(a + (alpha/a));
      return heron(alpha, a, --n);
  }
}

double f(double x){
  return (x*x -3);
}

double g(double x){
  return (pow(E, -1*x)*(3.2*sin(x)-0.5*cos(x)));
}

double h(double x){
  return (pow(E,x));
}

double hprime(double x){
  return 2*x;
}

double bisection(const double epi, double a, double b, double (*fp)(double)){
  assert(epi > 0);
  assert(a < b);
  assert(fp(a)*fp(b) < 0);

  double d = (a + b)/2;
  if(fp(d) == 0.0)
    return d;
  else if((b-a) < epi)
    return d;
  else{
    if(fp(d)*fp(a) < 0)
      return bisection(epi, a, d, fp);
    else
      return bisection(epi, d, b, fp);
    }
}

double newton(double x, const double (*fp)(double), const double (*fprimep)(double)){
  double next = x - (fp(x)/fprimep(x));

  if ((next - x) == 0.0)
    return next;
  else if (fp(x) == 0.0)
    return next;
  else
    return newton(next, fp, fprimep);
}

double num_dev(double (*fp)(double), const double x){
  const double DELTA = 1.0e8;
  double x1 = x - DELTA;
  double x2 = x + DELTA;
  double y1 = fp(x1);
  double y2 = fp(x2);
  return (y2 - y1)/(x2 - x1);
}

int first_dev_test()
int main(){
  double test = heron(17, 4, 10);
  printf("%.10lf\n", test);

  double (*fp)(double) = f;
  double (*sin_)(double) = sin;
  double (*gp)(double) = g;
  double (*hp)(double) = h;
  double (*hprimep)(double) = hprime;

  test = num_dev(fp, 3);
  printf("%.10lf\n", test);

  test = num_dev(gp, -2);
  printf("%.10lf\n", test);

  test = num_dev(hp, 0.5);
  printf("%.10lf\n", test);

  // test = bisection(0.001, 3, 4, gp);
  // printf("%.10lf\n", test);
  //
  // test = bisection(0.00001, 1, 99, sin_);
  // printf("%.10lf\n", test);
  //
  // test = newton(0.000000001, 1, hp, hprimep);
  // printf("%.10lf\n", test);
}
