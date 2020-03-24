#ifndef __SSE_H__
#define __SSE_H__

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sys/time.h>
#include "globals.h"

const float eps                = 1e-7;       // single precision epsilon
const float inv4PI             = 0.25/M_PI;  // Laplace kernel coefficient

template<typename T>
class vec3 {
public:
  T x;
  T y;
  T z;
};

template<typename T>
class vec4 {
public:
  T x;
  T y;
  T z;
  T w;
};

extern vec3<float> *bodyAccel;
extern vec4<float> *bodyPos;
extern vec3<float> *bodyVel;
extern double tic,t[9];
double get_time(void);
 void direct(int numParticles);
 void direct_seq(int numParticles);
 double get_nearwall_potential(float x, float y);

 void getEFromElectrons(vec3<float> &bodyAccel_, double x, double y, double z,  int n);



#endif // __FMM_H__
