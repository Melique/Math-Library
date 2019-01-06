#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

bool is_prime(const int n){
  assert(n > 1);

  const int ROOT_N = sqrt(n);
  int track = 0;

  for(int i =2 ; i <= ROOT_N; ++i){
    if(n % i == 0) {
      track = 1;
      break;
    }
  }
  if(track == 0) return 1;
  return 0;
}


int *division_alg(const int a, const int b){
  assert(b > 0);

  int *re = malloc(sizeof(int) * 2);
  int q = a/b;
  int r = a-b*q;

  if (r < 0){
    r += b;
    q = (a-r)/b;
  }

  re[0] = q;
  re[1] = r;
  return re;
}

int EE(const int a, const int b) {

  if(b == 0) return abs(a);
  else if (b == a) return abs(a) ;
  else{
    int *re = division_alg(a,b);
    return EE(b,re[1]);
  }
}

int* EEA(int a, int b){
  assert(a > 0);
  assert(b > 0);

  if(a < b){
    int temp = a;
    a = b;
    b = temp;
  }

  int track_2[4] = {1,0,a,0};
  int track_1[4] = {0,1,b,0};
  int current[4];

  do {
    current[3] = track_2[2]/track_1[2];
    current[0] = track_2[0]-current[3]*track_1[0];
    current[1] = track_2[1]-current[3]*track_1[1];
    current[2] = track_2[2]-current[3]*track_1[2];
    memcpy(track_2, track_1, sizeof(int)*4);
    memcpy(track_1, current, sizeof(int)*4);
  }while(current[2] != 0);


  if(a < 0) track_2[0] *= -1;
  if(b < 0) track_2[1] *= -1;

  const int length = 4;
  int *re = malloc(sizeof(int)*length);
  for(int i = 0; i < length; ++i){
    re[i] = track_2[i];
  }

  return re;
}

int** PF(const int a){
  assert(a > 1);

  int maxlen = 1;
  int len = 1;
  int *primes = malloc(maxlen * sizeof(int));
  int *occur = malloc(maxlen * sizeof(int));

  int **re = malloc(3 * sizeof(int*));

  if(is_prime(a)){
    primes[0] = a;
    occur[0] = 1;
    re[0] = primes;
    re[1] = occur;
    re[2] = &len;
    return re;
  }

  int track = 0;
  const int root_a = sqrt(a);

  for(int i =2; i < root_a; ++i){
    if(is_prime(i) && (a % i == 0)) {

      if(len == maxlen){
        maxlen *= 2;
        primes = realloc(primes, maxlen * sizeof(int));
        occur = realloc(occur, maxlen * sizeof(int));
      }

      primes[track] = i;
      int temp = a;
      int j = 0;

      while(temp % i == 0){
        occur[track] = ++j;
        temp /= i;
      }
      ++track;
      ++len;
    }
  }

  re[0] = primes;
  re[1] = occur;
  len = (len == 1) ? len: --len;
  re[2] = &len;

  return re;
}
