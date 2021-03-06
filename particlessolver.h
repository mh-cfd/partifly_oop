#ifndef PARTICLESSOLVER_H
#define PARTICLESSOLVER_H
#include "tools.h"
#include "fieldsloader.h"

#define MAXPARTICLES 31250
#define NXDepos 30
#define NYDepos 30

class ParticlesSolver
{
public:
    ParticlesSolver();
    void sweep();
    void particlesStep(int threadIdx, int steps);
    void saveParticlesInFile(double time);

    double m_dtParticles;
    double m_tParticles[THREADNUM];
    double m_dtGrid;
    double m_tGrid[THREADNUM];
    double m_dtSave;
    double m_tPrevSave;
    int m_counter;
    int m_numParticles[THREADNUM];
    int m_currFile[THREADNUM];
    vec3* m_bodyAccel[THREADNUM];
    vec4* m_bodyPos[THREADNUM];
    vec3* m_bodyVel[THREADNUM];
    vec3* m_bodyTurbVel[THREADNUM];
    MultythreadRandom* m_multyRand;
    DistributedRandom* m_gaussRand;
    FieldsLoader* m_fieldsLoader;
    int m_NX, m_NY, m_NZ;
    int* m_dims;

    double* m_x;
    double* m_y;
    double* m_z;
    double**** m_u;
    double**** m_v;
    double**** m_w;
    double**** m_uSig;
    double**** m_vSig;
    double**** m_wSig;
    double m_deposPart[THREADNUM][NXDepos][NYDepos];

    double x_save[MAXPARTICLES*2][1000];
    double y_save[MAXPARTICLES*2][1000];
    double z_save[MAXPARTICLES*2][1000];
    double u_save[MAXPARTICLES*2][1000];
    double v_save[MAXPARTICLES*2][1000];
    double w_save[MAXPARTICLES*2][1000];
    int totalNum_save = 0;
    int num_save[MAXPARTICLES*2];

private:
    vec4 initPos(int threadIdx);
    vec3 initVel(int threadIdx);
    void createRandomParticles(int threadIdx);
    void deleteParticle(int particlesIdx, int threadIdx);
    vec3 interpolateValue(int fileIdx, const vec4&bodyP, double**** u_, double**** v_, double**** w_);
    void calcTurbulentVel(int fileIdx, int threadIdx, const vec4&bodyP, int particleIdx, vec3&em_uSig);
};

#endif // PARTICLESSOLVER_H
