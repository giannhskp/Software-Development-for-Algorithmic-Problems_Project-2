#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

// these constants used at normal distribution
#define SIGMA 1.00
#define MI 0.00

double uniform_distribution(int rangeLow, int rangeHigh) {
    double myRand = rand()/(1.0 + RAND_MAX);
    int range = rangeHigh - rangeLow;
    double myRand_scaled = (myRand * range) + rangeLow;
    return myRand_scaled;
}

double rand_gen() {
   // return a uniformly distributed random value
   return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}
double normalRandom() {
   // return a normally distributed random value
   double v1=rand_gen();
   double v2=rand_gen();
   return (cos(2*3.14*v2)*sqrt(-2.*log(v1)))*SIGMA+MI;
}

double dot_product(double *v, double *u,int d){
  // calculate the dot product between the given vectors and return it
  double result = 0.0;
  for (int i = 0; i < d; i++)
      result += v[i]*u[i];
  return result;
}

int mod_Int_Int(int a, int b){
  int r = a % b;
  int result = (r < 0) ? (r + b) : r;
  return result;
}

long long int mod_LLI_UI(long long int a, unsigned int b){
  long long  int r = a % b;
  long long int result = (r < 0) ? (r + b) : r;
  return result;
}

long long int mod_LLI_I(long long int a,int b){
  long long  int r = a % b;
  long long int result = (r < 0) ? (r + b) : r;
  return result;
}
