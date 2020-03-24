#include "sse_sum.h"
#include "globals.h"

// summation
#define AX      "%xmm0"
#define AY      "%xmm1"
#define AZ      "%xmm2"
#define PHI     "%xmm3"
// j particle
#define XJ      "%xmm4"
#define YJ      "%xmm5"
#define ZJ      "%xmm6"
#define MJ      "%xmm7"
// temporary
#define RINV    "%xmm8"
#define X2      "%xmm9"
#define Y2      "%xmm10"
#define Z2      "%xmm11"
// fixed i particle
#define XI      "%xmm12"
#define YI      "%xmm13"
#define ZI      "%xmm14"
#define R2      "%xmm15"

#define XORPS(a, b) asm("xorps "  a  ","  b );
#define LOADPS(mem, reg) asm("movaps %0, %" reg:: "m" (mem));
#define STORPS(reg, mem) asm("movaps %" reg " , %0" :: "m"(mem));
#define MOVAPS(src, dst) asm("movaps " src "," dst);
#define MOVQ(src, dst) asm("movq " src "," dst);
#define BCAST0(reg) asm("shufps $0x00, " reg ","  reg);
#define BCAST1(reg) asm("shufps $0x55, " reg ","  reg);
#define BCAST2(reg) asm("shufps $0xaa, " reg ","  reg);
#define BCAST3(reg) asm("shufps $0xff, " reg ","  reg);
#define MULPS(src, dst) asm("mulps " src "," dst);
#define ADDPS(src, dst) asm("addps " src ","  dst);
#define SUBPS(src, dst) asm("subps "  src "," dst);
#define RSQRTPS(src, dst) asm("rsqrtps " src "," dst);
#define MOVHLPS(src, dst) asm("movhlps " src "," dst);
#define DEBUGPS(reg)

#define ALIGN16 __attribute__ ((aligned(16)))

#define NJMAX (1<<24)

typedef double v2df __attribute__ ((vector_size(16)));
typedef float  v4sf __attribute__ ((vector_size(16)));
typedef int    v4si __attribute__ ((vector_size(16)));
typedef short  v8hi __attribute__ ((vector_size(16)));

typedef struct iptdata{
  float x[4];
  float y[4];
  float z[4];
  float eps2[4]; // not used in this implementation
} Ipdata ALIGN16;
typedef struct fodata{
  float ax[4];
  float ay[4];
  float az[4];
  float phi[4];
} Fodata ALIGN16;
typedef struct jpdata{
  float x, y, z, m;
} Jpdata ALIGN16;

/*void p2p_kernel(Ipdata *ipdata, Fodata *fodata, Jpdata *jpdata, int nj){
  int j;
  assert(((unsigned long)jpdata & 15) == 0);
  assert(((unsigned long)ipdata & 15) == 0);
  assert(((unsigned long)fodata & 15) == 0);

  XORPS(AX, AX);               // AX = 0
  XORPS(AY, AY);               // AY = 0
  XORPS(AZ, AZ);               // AZ = 0
  XORPS(PHI, PHI);             // PHI = 0

  LOADPS(*ipdata->x, XI);      // XI = *ipdata->x
  LOADPS(*ipdata->y, YI);      // YI = *ipdata->y
  LOADPS(*ipdata->z, ZI);      // ZI = *ipdata->z
  LOADPS(*ipdata->eps2, R2);   // R2 = *ipdata->eps2

  LOADPS(jpdata[0], MJ);       // MJ = *jpdata->x,y,z,m
  MOVAPS(MJ, X2);              // X2 = MJ
  MOVAPS(MJ, Y2);              // Y2 = MJ
  MOVAPS(MJ, Z2);              // Z2 = MJ

  BCAST0(X2);                  // X2 = *jpdata->x
  BCAST1(Y2);                  // Y2 = *jpdata->y
  BCAST2(Z2);                  // Z2 = *jpdata->z
  BCAST3(MJ);                  // MJ = *jpdata->m

  SUBPS(XI, X2);               // X2 = X2 - XI
  SUBPS(YI, Y2);               // Y2 = Y2 - YI
  SUBPS(ZI, Z2);               // Z2 = Z2 - ZI

  MOVAPS(X2, XJ);              // XJ = X2
  MOVAPS(Y2, YJ);              // YJ = Y2
  MOVAPS(Z2, ZJ);              // ZJ = Z2

  MULPS(X2, X2);               // X2 = X2 * X2
  MULPS(Y2, Y2);               // Y2 = Y2 * Y2
  MULPS(Z2, Z2);               // Z2 = Z2 * Z2

  ADDPS(X2, R2);               // R2 = R2 + X2
  ADDPS(Y2, R2);               // R2 = R2 + Y2
  ADDPS(Z2, R2);               // R2 = R2 + Z2

  LOADPS(jpdata[1], X2);       // X2 = *jpdata->x,y,z,m
  MOVAPS(X2, Y2);              // Y2 = X2
  MOVAPS(X2, Z2);              // Z2 = X2
  for(j=0;j<nj;j++){
    RSQRTPS(R2, RINV);         // RINV = rsqrt(R2)
    jpdata++;
    LOADPS(*ipdata->eps2, R2); // R2 = *ipdata->eps2
    BCAST0(X2);                // X2 = *jpdata->x
    BCAST1(Y2);                // Y2 = *jpdata->y
    BCAST2(Z2);                // Z2 = *jpdata->z
    SUBPS(XI, X2);             // X2 = X2 - XI
    SUBPS(YI, Y2);             // Y2 = Y2 - YI
    SUBPS(ZI, Z2);             // Z2 = Z2 - ZI

    MULPS(RINV, MJ);           // MJ = MJ * RINV
    SUBPS(MJ, PHI);            // PHI = PHI - MJ
    MULPS(RINV, RINV);         // RINV = RINV * RINV
    MULPS(MJ, RINV);           // RINV = MJ * RINV
    LOADPS(jpdata[0], MJ);     // MJ = *jpdata->x,y,z,m
    BCAST3(MJ);                // MJ = *jpdata->m

    MULPS(RINV, XJ);           // XJ = XJ * RINV
    ADDPS(XJ, AX);             // AX = AX + XJ
    MOVAPS(X2, XJ);            // XJ = X2
    MULPS(X2, X2);             // X2 = X2 * X2
    ADDPS(X2, R2);             // R2 = R2 + X2
    LOADPS(jpdata[1], X2);     // X2 = *jpdata->x,y,z,m

    MULPS(RINV, YJ);           // YJ = YJ * RINV
    ADDPS(YJ, AY);             // AY = AY + YJ
    MOVAPS(Y2, YJ);            // YJ = Y2
    MULPS(Y2, Y2);             // Y2 = Y2 * Y2
    ADDPS(Y2, R2);             // R2 = R2 + Y2
    MOVAPS(X2, Y2);            // Y2 = X2

    MULPS(RINV, ZJ);           // ZJ = ZJ * RINV
    ADDPS(ZJ, AZ);             // AZ = AZ + ZJ
    MOVAPS(Z2, ZJ);            // ZJ = Z2
    MULPS(Z2, Z2);             // Z2 = Z2 * Z2
    ADDPS(Z2, R2);             // R2 = R2 + Z2
    MOVAPS(X2, Z2);            // Z2 = X2
  }
  STORPS(AX, *fodata->ax);     // AX = *fodata->ax
  STORPS(AY, *fodata->ay);     // AY = *fodata->ay
  STORPS(AZ, *fodata->az);     // AZ = *fodata->az
  STORPS(PHI, *fodata->phi);   // PHI = *fodata->phi
}*/

/*static inline void v4sf_transpose(
    v4sf *d0, v4sf *d1, v4sf *d2, v4sf *d3,
    v4sf  s0, v4sf  s1, v4sf  s2, v4sf  s3)
{
  *d0 = __builtin_ia32_unpcklps(
        __builtin_ia32_unpcklps(s0, s2),
        __builtin_ia32_unpcklps(s1, s3));
  *d1 = __builtin_ia32_unpckhps(
        __builtin_ia32_unpcklps(s0, s2),
        __builtin_ia32_unpcklps(s1, s3));
  *d2 = __builtin_ia32_unpcklps(
        __builtin_ia32_unpckhps(s0, s2),
        __builtin_ia32_unpckhps(s1, s3));
  *d3 = __builtin_ia32_unpckhps(
        __builtin_ia32_unpckhps(s0, s2),
        __builtin_ia32_unpckhps(s1, s3));
}*/

/*static inline void v3sf_store_sp(v4sf vec, float *d0, float *d1, float *d2){
  *d0 = __builtin_ia32_vec_ext_v4sf(vec, 0);
  *d1 = __builtin_ia32_vec_ext_v4sf(vec, 1);
  *d2 = __builtin_ia32_vec_ext_v4sf(vec, 2);
}
*/




vec3<float> *bodyAccel;
vec3<float> *bodyVel;
vec4<float> *bodyPos;
double get_time(void) {                          // a simple timer
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);
  return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}



// direct summation kernel
void direct_wall(vec3<float> *bodyAccel_, vec4<float> *bodyPos_, int n) {
  int i,j;
  vec3<float> dist;
  float invDist,invDistCube;
  for( i=0; i<n; i++ ) 
  {
    vec3<float> ai = {0.0, 0.0, 0.0};
    for( j=0; j<N_X-1; j++ )
    {
      dist.x = bodyPos_[i].x-(w_x0+ dx*j);
      dist.y = bodyPos_[i].y-w_y0;
      dist.z = 0.0;//bodyPos_[i].z-wallP[j].z;
      invDist = 1.0/sqrtf(dist.x*dist.x+dist.y*dist.y+eps*eps);
      invDistCube =10.0* BoundaryLayerGauss[j]*invDist*invDist*invDist;
      ai.x -= dist.x*invDistCube;
      ai.y -= dist.y*invDistCube;
      //ai.z -= dist.z*invDistCube;
    }
    bodyAccel_[i].x += inv4PI*ai.x;
    bodyAccel_[i].y += inv4PI*ai.y;
    //bodyAccel_[i].z = 0.0;//inv4PI*ai.z;
  }

}


double get_nearwall_potential(float x, float y)
{
    int i,j;
    float Dx,Dy;
    double phi=0.0;
    float invDist,invDistW,invDistW0;
    Dy = y-w_y0;
          
      for( j=0; j<N_X; j++ )
      {
        Dx = x-(w_x0+ dx*j);
        invDistW = BoundaryLayerGauss[j]/sqrtf(Dx*Dx+Dy*Dy+dx*dx);
        invDistW0 = BoundaryLayerGauss[j]/sqrtf(Dx*Dx+dx*dx);
        phi -= invDistW-invDistW0;
      }
        return -phi*inv4PI;
}

// direct summation kernel
/*void direct(int n) {
  int ii,i,offset;
  Ipdata iptcl;
  Fodata fout;
  Jpdata *jptcl;
  jptcl = (Jpdata *) malloc(sizeof(Jpdata)*NJMAX);
  for( i=0; i<n; i++ )
  {
    *(v4sf *)(jptcl+i) = (v4sf) {bodyPos[i].x,bodyPos[i].y,bodyPos[i].z,bodyPos[i].w};
  }
  for( ii=0; ii<n/4; ii++ )
  {
    offset = 4*ii;
    for(i=0;i<4;i++)
    {
      iptcl.x[i] = bodyPos[offset+i].x;
      iptcl.y[i] = bodyPos[offset+i].y;
      iptcl.z[i] = bodyPos[offset+i].z;
      iptcl.eps2[i] = eps*eps;
    }
    p2p_kernel(&iptcl, &fout, jptcl, n);
    v4sf ax = *(v4sf *)(fout.ax);
    v4sf ay = *(v4sf *)(fout.ay);
    v4sf az = *(v4sf *)(fout.az);
    v4sf phi= -*(v4sf *)(fout.phi);
    v4sf f0, f1, f2, f3;
    v4sf_transpose(&f0, &f1, &f2, &f3, ax, ay, az, phi);
    v3sf_store_sp(f0, &bodyAccel[offset+0].x, &bodyAccel[offset+0].y, &bodyAccel[offset+0].z);
    v3sf_store_sp(f1, &bodyAccel[offset+1].x, &bodyAccel[offset+1].y, &bodyAccel[offset+1].z);
    v3sf_store_sp(f2, &bodyAccel[offset+2].x, &bodyAccel[offset+2].y, &bodyAccel[offset+2].z);
    v3sf_store_sp(f3, &bodyAccel[offset+3].x, &bodyAccel[offset+3].y, &bodyAccel[offset+3].z);
    for(i=0;i<4;i++)
    {
      bodyAccel[offset+i].x *= inv4PI;
      bodyAccel[offset+i].y *= inv4PI;
      bodyAccel[offset+i].z *= inv4PI;
    }
  }
  free(jptcl);

  //direct_wall(bodyAccel,bodyPos,n);
}*/

void getEFromElectrons(vec3<float> &bodyAccel_, double x, double y, double z,  int n) {
    vec3<float> dist;
    float invDist,invDistCube;
      vec3<float> ai = {0.0, 0.0, 0.0};
      for( int j=0; j<n; j++ ){
        dist.x = bodyPos[j].x-x;
        dist.y = bodyPos[j].y-y;
        dist.z = bodyPos[j].z-z;
        invDist = 1.0/sqrtf(dist.x*dist.x+dist.y*dist.y+dist.z*dist.z+1e-14);
        invDistCube = bodyPos[j].w*invDist*invDist*invDist;
        ai.x -= dist.x*invDistCube;
        ai.y -= dist.y*invDistCube;
        ai.z -= dist.z*invDistCube;
      }

      bodyAccel_.x = inv4PI*ai.x/eps0;
      bodyAccel_.y = inv4PI*ai.y/eps0;
      bodyAccel_.z = inv4PI*ai.z/eps0;
}

// direct summation kernel
void direct_seq(int n) {
  int i,j;
  vec3<float> dist;
  float invDist,invDistCube;
  for( i=0; i<n; i++ ) {
    vec3<float> ai = {0.0, 0.0, 0.0};
    for( j=0; j<n; j++ ){
      dist.x = bodyPos[i].x-bodyPos[j].x;
      dist.y = bodyPos[i].y-bodyPos[j].y;
      dist.z = bodyPos[i].z-bodyPos[j].z;
      invDist = 1.0/sqrtf(dist.x*dist.x+dist.y*dist.y+dist.z*dist.z+1e-18/*+eps*/);
      invDistCube = bodyPos[j].w*invDist*invDist*invDist;
      ai.x -= dist.x*invDistCube;
      ai.y -= dist.y*invDistCube;
      ai.z -= dist.z*invDistCube;
    }
   /* for( j=0; j<N_X-1; j++ ){
      dist.x = bodyPos_[i].x-(dx*j-0.5);
      dist.y = bodyPos_[i].y+0.15;
      dist.z = 0.0;//bodyPos_[i].z-wallP[j].z;
      invDist = 1.0/sqrtf(dist.x*dist.x+dist.y*dist.y+eps);
      invDistCube = BoundaryLayer[j]*invDist*invDist*invDist;
      ai.x -= dist.x*invDistCube;
      ai.y -= dist.y*invDistCube;
      //ai.z -= dist.z*invDistCube;
    }*/


    bodyAccel[i].x = inv4PI*ai.x/eps0;
    bodyAccel[i].y = inv4PI*ai.y/eps0;
    bodyAccel[i].z = inv4PI*ai.z/eps0;
  }

//    direct_wall(bodyAccel_, bodyPos_, n);
}


