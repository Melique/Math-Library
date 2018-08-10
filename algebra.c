#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>

//Todo:
// -destroy functions
// -cosnt
// -documentation

bool is_prime(const int n){
  assert(n > 1);

  const int root_n = sqrt(n);
  int track = 0;

  for(int i =2 ; i <= root_n; ++i){
    if(n % i == 0) {
      track = 1;
      break;
    }
  }
  if(track == 0) return 1;
}


int* division_alg(const int a, const int b){
  assert(b > 0);

  int *re = malloc(sizeof(int)*2);
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
    EE(b,re[1]);
  }
}

int* EEA(int a, int b){
  assert(a > 0);
  assert(b > 0);

  if(a < b){
    int temp = a;
    a = b;
    b= temp;
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

// long SM(double a, const double e, const double n){
//   int i = 2;
//   const double ori = a;
//   printf("ori: %f\n", ori);
//
//   while(i < e){
//     a = (long) pow(a, 2) % (long) n;
//     printf("a: %f\n", a);
//     i *= 2;
//   }
//
//   i /= 2;
//
//   long rest = (long) pow(ori, (e-i)) % (long) n;
//   printf("rest: %ld\n", rest);
//   long re = (long) (a*rest) % (long) n;
//   return re;
// }
//
int ** RSA_setUp(const int p, const int q){
  assert(is_prime(p));
  assert(is_prime(q));

  const int n = p*q;
  const int r = (p-1)*(q-1);
  int e;

  for(int i = 2; i < r; ++i){
    if(EE(i,r) == 1){
      e=i;
      break;
    }
  }

  const int *re = EEA(r,e);
  int *public = malloc(sizeof(int)*2);
  int *private = malloc(sizeof(int)*2);

  const int d = (re[1] > 0) ? re[1]: re[1]+r;

  public[0] = e;
  public[1] = n;
  private[0] = d;
  private[1] = n;

  int **test = malloc(sizeof(int*)*2);
  test[0] = public;
  test[1] = private;

  return test;
}
//
// long RSA_encrypt(int m, int e, int n){
//   assert(0 <= m);
//   assert(m < n);
//
//   double re = pow((double) m, (double) e);
//   printf("TEST 2: %d\n", re);
//   const long C = (int) re % n;
//   return C;
// }
//
// long RSA_decrypt(const int C, int d, int n){
//   double re = pow((double) C, (double) d);
//   printf("TEST 1: %d\n", re);
//   const long int R = (int) re % n;
//   return R;
// }


//DESTROY FUNCIONS
int main(){
  // int *test = EEA(1386,322);
  // printf("x=%d, y=%d, d=%d, q=%d\n", test[0],test[1],test[2],test[3]);
  // free(test);
  // int *neg = EEA(231,660);
  // printf("x=%d, y=%d, d=%d, q=%d\n", neg[0],neg[1],neg[2],neg[3]);
  // free(neg);
  // int *temp = divisionAlg(24750,2);
  // printf("%d, %d\n", temp[0],temp[1]);
  // printf("%d\n", is_prime(8));
  // int **dole = PF(434511);
  // int length = *(dole[2]);
  //
  // for(int i = 0; i < length; ++i){
  //   printf("%d: %d \n", dole[0][i], dole[1][i]);
  // }

  // int **jacob = RSA_setUp(7, 19);
  // printf("\npublic: %d,%d\n", jacob[0][0], jacob[0][1]);
  // printf("private: %d,%d\n", jacob[1][0], jacob[1][1]);
  // long C = RSA_encrypt(6, 5, 133);
  // printf("C: %ld\n", C);
  // long R = RSA_decrypt(62, 65, 133);
  // printf("R: %ld\n", R);

  printf("%ld\n", SM(30, 19, 391));


}
