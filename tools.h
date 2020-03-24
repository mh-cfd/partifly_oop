#ifndef DISTRIBUTEDRANDOM_H
#define DISTRIBUTEDRANDOM_H

#include<xmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <random>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/
#include  <math.h>
#include <sys/time.h>
#include <iostream>
#include <vector>
#include <thread>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#define W_WIDTH 1400
#define W_HEIGHT 900
#define RNX 3000
#define RNY 2000

#define THREADNUM 1

double get_time(void);
int strCompare(const void* a, const void* b);

class vec3 {
public:
    float x, y, z;
};

class vec4 {
public:
    float x, y, z, rho, r;
    int count;
    int partNumber;
};

class MultythreadRandom
{
public:
    MultythreadRandom();
    float get(int threadIdx);
private:
    void init();
    drand48_data randBuffers[32];
    unsigned int seeds[32];
    uint32_t xr[32],yr[32],zr[32],wr[32];
};


class DistributedRandom
{
public:
    DistributedRandom();
    ~DistributedRandom(){ delete m_rand; }
    double get(int threadIdx);
private:
    void dist();
    void integ();
    void interpol();
    typedef struct  {
        double x;
        double y;
    }vec2;
    vec2 m_pdf_x[RNX];
    vec2 m_int_y[RNY];
    vec2 m_int_x[RNX];
    MultythreadRandom *m_rand;
};

#endif // DISTRIBUTEDRANDOM_H
