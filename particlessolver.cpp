#include "particlessolver.h"

ParticlesSolver::ParticlesSolver()
{
    m_dtParticles = 5.0;
    m_dtGrid = 60.;
    m_dtSave = 0.0005;
    m_tPrevSave = 0.0;
    for(int i = 0; i < THREADNUM; i++)
    {
        m_tParticles[i] = 0.0;
        m_tGrid[i] = 0.0;
        m_currFile[i] = 1;
        m_numParticles[i] = 0;
        m_bodyAccel[i] = new vec3[MAXPARTICLES];
        m_bodyVel[i] = new vec3[MAXPARTICLES];
        m_bodyPos[i] = new vec4[MAXPARTICLES];
        m_bodyTurbVel[i] = new vec3[MAXPARTICLES];
    }

    for( int i = 0; i < MAXPARTICLES; i++ )
        num_save[i] = 0;

    m_multyRand = new MultythreadRandom();
    m_gaussRand = new DistributedRandom();
    m_fieldsLoader = new FieldsLoader();
    m_fieldsLoader->loadVTK();

    m_dims = m_fieldsLoader->m_dims;
    m_x = m_fieldsLoader->m_x;
    m_y = m_fieldsLoader->m_y;
    m_z = m_fieldsLoader->m_z;
    m_u = m_fieldsLoader->m_u;
    m_v = m_fieldsLoader->m_v;
    m_w = m_fieldsLoader->m_w;
    m_uSig = m_fieldsLoader->m_uSig;
    m_vSig = m_fieldsLoader->m_vSig;
    m_wSig = m_fieldsLoader->m_wSig;
    m_NX = m_dims[0];
    m_NY = m_dims[1];
    m_NZ = m_dims[2];
    m_counter=0;
}

void ParticlesSolver::particlesStep(int threadIdx, int steps)
{
    int i, k;
    for(k = 0; k < steps; k++)
    {
        m_tParticles[threadIdx] += m_dtParticles;
        if(m_tParticles[threadIdx] > m_tGrid[threadIdx])
        {
            m_tGrid[threadIdx] += m_dtGrid;
            m_currFile[threadIdx]++;
        }

        createRandomParticles(threadIdx);

        for( i=0; i<m_numParticles[threadIdx]; i++ )
        {
            double a = (m_tParticles[threadIdx] - (m_tGrid[threadIdx] - m_dtGrid))/(m_dtGrid);
            vec3 eu, euPrev, euSig, euSigPrev,euSigInter;
            eu = interpolateValue(m_currFile[threadIdx] , m_bodyPos[threadIdx][i], m_u, m_v, m_w);
            euPrev = interpolateValue(m_currFile[threadIdx] - 1, m_bodyPos[threadIdx][i], m_u, m_v, m_w);
            euSig = interpolateValue(m_currFile[threadIdx] , m_bodyPos[threadIdx][i], m_uSig, m_vSig, m_wSig);
            euSigPrev = interpolateValue(m_currFile[threadIdx] - 1, m_bodyPos[threadIdx][i], m_uSig, m_vSig, m_wSig);
            euSigInter.x = euSig.x * a + euSigPrev.x * (1.0-a);
            euSigInter.y = euSig.y * a + euSigPrev.y * (1.0-a);
            euSigInter.z = euSig.z * a + euSigPrev.z * (1.0-a);
            calcTurbulentVel(m_currFile[threadIdx] , threadIdx, m_bodyPos[threadIdx][i], i, euSigInter);

            vec3 V;
            V.x = (eu.x * a + euPrev.x * (1.0-a));
            V.y = (eu.y * a + euPrev.y * (1.0-a));
            V.z = (eu.z * a + euPrev.z * (1.0-a));

            double rhoAir = 0.999;

            m_bodyVel[threadIdx][i].x = 10000*V.x + 0.01 * m_bodyTurbVel[threadIdx][i].x;
            m_bodyVel[threadIdx][i].y = 10000*V.y + 0.01 * m_bodyTurbVel[threadIdx][i].y;
            m_bodyVel[threadIdx][i].z = 10000*V.z + 0.01 * m_bodyTurbVel[threadIdx][i].z - 0.01 * 9.8 * m_dtParticles * (1.0 - rhoAir / m_bodyPos[threadIdx]->rho);
        }

        std::vector<int> numToDel;
        for( i=0; i<m_numParticles[threadIdx]; i++ )
        {

            m_bodyPos[threadIdx][i].x +=  m_dtParticles * m_bodyVel[threadIdx][i].x;
            m_bodyPos[threadIdx][i].y +=  m_dtParticles * m_bodyVel[threadIdx][i].y;
            m_bodyPos[threadIdx][i].z +=  m_dtParticles * m_bodyVel[threadIdx][i].z;
            if((m_bodyPos[threadIdx][i].x > m_x[m_NX-1]) || (m_bodyPos[threadIdx][i].x < m_x[0])
                    || (m_bodyPos[threadIdx][i].y > m_y[m_NY-1]) || (m_bodyPos[threadIdx][i].y < m_y[0])
                    || (m_bodyPos[threadIdx][i].z > m_z[m_NZ-1]))
                numToDel.push_back(i);
            else if(m_bodyPos[threadIdx][i].z < m_z[0])
            {
                numToDel.push_back(i);
                int xIdx =  m_bodyPos[threadIdx][i].x / ((m_x[m_NX - 1] - m_x[0]) / NXDepos);
                int yIdx =  m_bodyPos[threadIdx][i].y / ((m_y[m_NY - 1] - m_y[0]) / NYDepos);
                m_deposPart[threadIdx][xIdx][yIdx] += m_bodyPos[threadIdx][i].count * m_bodyPos[threadIdx][i].rho * (4.0 / 3.0) * M_PI * m_bodyPos[threadIdx][i].r * m_bodyPos[threadIdx][i].r * m_bodyPos[threadIdx][i].r;
            }
        }
        for( i=0; i<numToDel.size(); i++)
            deleteParticle(numToDel[i], threadIdx);
    }
}

void ParticlesSolver::sweep()
{
    double time1=get_time();
    std::vector <std::thread> th_vec;
    for (int i = 0; i < THREADNUM; ++i)
    {
        th_vec.push_back(std::thread(&ParticlesSolver::particlesStep, this , i, 30));
    }
    for (int i = 0; i < THREADNUM; ++i)
    {
        th_vec.at(i).join();
    }
    double time2=get_time();

    printf("t2=%e \n", time2-time1);

    m_counter++;
    if (m_counter>1)
    {
        for(int i=0; i<m_numParticles[0]; i++ )
        {
            x_save[i][num_save[i]]=m_bodyPos[0][i].x;
            y_save[i][num_save[i]]=m_bodyPos[0][i].y;
            z_save[i][num_save[i]]=m_bodyPos[0][i].z;
            u_save[i][num_save[i]]=m_bodyVel[0][i].x;
            v_save[i][num_save[i]]=m_bodyVel[0][i].y;
            w_save[i][num_save[i]]=m_bodyVel[0][i].z;
            num_save[i]++;
            if (num_save[i]>1000) num_save[i]=0;
        }
        m_counter=0;
    }
    /*if(fabs(tParticles - tPrevSave - dtSave) < 0.000001)
    {
      tPrevSave += dtSave;
    save_particles_in_file(tParticles);
    }*/
}

vec3 ParticlesSolver::interpolateValue(int fileIdx, const vec4&bodyP, double**** u_, double**** v_, double**** w_)
{
    int i_,j_,k_;
    double a,b,c;
    i_ = (int)((m_NX-1) * (bodyP.x - m_x[0]) / (m_x[m_NX-1] - m_x[0]));
    j_ = (int)((m_NY-1) * (bodyP.y - m_y[0]) / (m_y[m_NY-1] - m_y[0]));
    k_ = 0;
    while ((m_z[k_] < bodyP.z) && (k_ < m_NZ-2))
        k_++;

    a = (bodyP.x - m_x[i_]) / (m_x[i_ + 1] - m_x[i_]);
    b = (bodyP.y - m_y[j_]) / (m_y[j_ + 1] - m_y[j_]);
    c = (bodyP.z - m_z[k_]) / (m_z[k_ + 1] - m_z[k_]);

    i_ = fmax(1, i_);
    i_ = fmin(m_NX - 2, i_);
    j_ = fmax(1, j_);
    j_ = fmin(m_NY - 2, j_);
    k_ = fmax(1, k_);
    k_ = fmin(m_NZ - 2, k_);

    a = fmax(0.0, a);
    a = fmin(1.0, a);
    b = fmax(0.0, b);
    b = fmin(1.0, b);
    c = fmax(0.0, c);
    c = fmin(1.0, c);

    vec3 res;
    res.x = (1.0 - c) * ((1.0 - b) * ((1.0 - a) * u_[fileIdx][i_][j_][k_]         + (a) * u_[fileIdx][i_ + 1][j_][k_])
            + (b)     * ((1.0 - a)              * u_[fileIdx][i_][j_ + 1][k_]     + (a) * u_[fileIdx][i_ + 1][j_ + 1][k_]))
            + (c)     * ((1.0 - b) * ((1.0 - a) * u_[fileIdx][i_][j_][k_ + 1]     + (a) * u_[fileIdx][i_ + 1][j_][k_ + 1])
            + (b)     * ((1.0 - a)              * u_[fileIdx][i_][j_ + 1][k_ + 1] + (a) * u_[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    res.y = (1.0 - c) * ((1.0 - b) * ((1.0 - a) * v_[fileIdx][i_][j_][k_]         + (a) * v_[fileIdx][i_ + 1][j_][k_])
            + (b)     * ((1.0 - a)              * v_[fileIdx][i_][j_ + 1][k_]     + (a) * v_[fileIdx][i_ + 1][j_ + 1][k_]))
            + (c)     * ((1.0 - b) * ((1.0 - a) * v_[fileIdx][i_][j_][k_ + 1]     + (a) * v_[fileIdx][i_ + 1][j_][k_ + 1])
            + (b)     * ((1.0 - a)              * v_[fileIdx][i_][j_ + 1][k_ + 1] + (a) * v_[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    res.z = (1.0 - c) * ((1.0 - b) * ((1.0 - a) * w_[fileIdx][i_][j_][k_]         + (a) * w_[fileIdx][i_ + 1][j_][k_])
            + (b)                  * ((1.0 - a) * w_[fileIdx][i_][j_ + 1][k_]     + (a) * w_[fileIdx][i_ + 1][j_ + 1][k_]))
            + (c)     * ((1.0 - b) * ((1.0 - a) * w_[fileIdx][i_][j_][k_ + 1]     + (a) * w_[fileIdx][i_ + 1][j_][k_ + 1])
            + (b)                  * ((1.0 - a) * w_[fileIdx][i_][j_ + 1][k_ + 1] + (a) * w_[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    return res;
}

void ParticlesSolver::calcTurbulentVel(int fileIdx, int threadIdx,const vec4&bodyP, int particleIdx, vec3&em_uSig)
{
    int i_,j_,k_;
    double a,b,c;

    i_ = (int)((m_NX-1) * (bodyP.x - m_x[0]) / (m_x[m_NX-1] - m_x[0]));
    j_ = (int)((m_NY-1) * (bodyP.y - m_y[0]) / (m_y[m_NY-1] - m_y[0]));
    k_ = 0;
    while ((m_z[k_] < bodyP.z) && (k_ < m_NZ-2))
        k_++;

    a = (bodyP.x - m_x[i_]) / (m_x[i_ + 1] - m_x[i_]);
    b = (bodyP.y - m_y[j_]) / (m_y[j_ + 1] - m_y[j_]);
    c = (bodyP.z - m_z[k_]) / (m_z[k_ + 1] - m_z[k_]);

    i_ = fmax(1, i_);
    i_ = fmin(m_NX - 2, i_);
    j_ = fmax(1, j_);
    j_ = fmin(m_NY - 2, j_);
    k_ = fmax(1, k_);
    k_ = fmin(m_NZ - 2, k_);

    a = fmax(0.0, a);
    a = fmin(1.0, a);
    b = fmax(0.0, b);
    b = fmin(1.0, b);
    c = fmax(0.0, c);
    c = fmin(1.0, c);

    double a_ = 1.0;
    double uxp =  (1.0 - c) * ((1.0 - b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_][k_]         + (a_) * m_uSig[fileIdx][i_ + 1][j_][k_])
            +                        (b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_ + 1][k_]     + (a_) * m_uSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +           (c) * ((1.0 - b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_][k_ + 1]     + (a_) * m_uSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                        (b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_ + 1][k_ + 1] + (a_) * m_uSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    a_ = 0.0;
    double ux =  (1.0 - c) * ((1.0 - b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_][k_]         + (a_) * m_uSig[fileIdx][i_ + 1][j_][k_])
            +                       (b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_ + 1][k_]     + (a_) * m_uSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +          (c) * ((1.0 - b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_][k_ + 1]     + (a_) * m_uSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                       (b) * ((1.0 - a_) * m_uSig[fileIdx][i_][j_ + 1][k_ + 1] + (a_) * m_uSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    double b_ = 1.0;
    double uyp =  (1.0 - c) * ((1.0 - b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_][k_]         + (a) * m_vSig[fileIdx][i_ + 1][j_][k_])
            +                        (b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_ + 1][k_]     + (a) * m_vSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +           (c) * ((1.0 - b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_][k_ + 1]     + (a) * m_vSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                        (b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_ + 1][k_ + 1] + (a) * m_vSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    b_ = 0.0;
    double uy =  (1.0 - c) * ((1.0 - b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_][k_]         + (a) * m_vSig[fileIdx][i_ + 1][j_][k_])
            +                       (b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_ + 1][k_]     + (a) * m_vSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +          (c) * ((1.0 - b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_][k_ + 1]     + (a) * m_vSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                       (b_) * ((1.0 - a) * m_vSig[fileIdx][i_][j_ + 1][k_ + 1] + (a) * m_vSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));


    double c_ = 1.0;
    double uzp =  (1.0 - c_) * ((1.0 - b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_][k_]         + (a) * m_wSig[fileIdx][i_ + 1][j_][k_])
            +                         (b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_ + 1][k_]     + (a) * m_wSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +           (c_) * ((1.0 - b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_][k_ + 1]     + (a) * m_wSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                         (b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_ + 1][k_ + 1] + (a) * m_wSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    c_ = 0.0;
    double uz =  (1.0 - c_) * ((1.0 - b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_][k_]         + (a) * m_wSig[fileIdx][i_ + 1][j_][k_])
            +                        (b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_ + 1][k_]     + (a) * m_wSig[fileIdx][i_ + 1][j_ + 1][k_]))
            +          (c_) * ((1.0 - b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_][k_ + 1]     + (a) * m_wSig[fileIdx][i_ + 1][j_][k_ + 1])
            +                        (b) * ((1.0 - a) * m_wSig[fileIdx][i_][j_ + 1][k_ + 1] + (a) * m_wSig[fileIdx][i_ + 1][j_ + 1][k_ + 1]));

    double C0, eps, k, tau;
    C0 = 0.6;
    eps = 1.0;
    k = 1.0;
    tau = 0.01*(0.5 + 0.75 * C0) * eps / k;

    m_bodyTurbVel[threadIdx][particleIdx].x += -m_bodyTurbVel[threadIdx][particleIdx].x * m_dtParticles * tau + em_uSig.x * sqrt(2.0 * m_dtParticles * tau) * 1000.0 * m_gaussRand->get((threadIdx))   + m_dtParticles * (uxp - ux) / (m_x[i_ + 1] - m_x[i_]);
    m_bodyTurbVel[threadIdx][particleIdx].y += -m_bodyTurbVel[threadIdx][particleIdx].y * m_dtParticles * tau + em_uSig.y * sqrt(2.0 * m_dtParticles * tau) * 1000.0 * m_gaussRand->get((threadIdx))   + m_dtParticles * (uyp - uy) / (m_y[j_ + 1] - m_y[j_]);
    m_bodyTurbVel[threadIdx][particleIdx].z += -m_bodyTurbVel[threadIdx][particleIdx].z * m_dtParticles * tau + em_uSig.z * sqrt(2.0 * m_dtParticles * tau) * 1000.0 * m_gaussRand->get((threadIdx))   + m_dtParticles * (uzp - uz) / (m_z[k_ + 1] - m_z[k_]);
}


vec4 ParticlesSolver::initPos(int threadIdx)
{
    vec4 res;
    res.x = m_x[0] + m_multyRand->get(threadIdx) * (m_x[m_NX-1] - m_x[0]); //0.5*(m_x[0]+m_x[m_NX-1]) + m_multyRand->get(threadIdx) * 0.16*(m_x[m_NX-1] - m_x[0]);//m_x[0] + m_multyRand->get(threadIdx) * (m_x[m_NX-1] - m_x[0]);
    res.y = m_y[0] + m_multyRand->get(threadIdx) * (m_y[m_NY-1] - m_y[0]); //0.5*(m_y[0]+m_y[m_NY-1]) + m_multyRand->get(threadIdx) * 0.16*(m_y[m_NY-1] - m_y[0]);//m_y[0] + m_multyRand->get(threadIdx) * (m_y[m_NY-1] - m_y[0]);
    res.z = m_z[0] + m_multyRand->get(threadIdx) * 0.1 * (m_z[m_NZ-1] - m_z[0]);//0.25*(m_z[0]+m_z[m_NZ-1]) + m_multyRand->get(threadIdx) * 0.16 * (m_z[m_NZ-1] - m_z[0]);//m_z[0] + m_multyRand->get(threadIdx) * 0.1 * (m_z[m_NZ-1] - m_z[0]);
    res.count = 1;
    res.rho = 1.0;
    res.r = 0.01;
    return res;
}

vec3 ParticlesSolver::initVel(int threadIdx)
{
    vec3 vel;
    vel.x = m_multyRand->get(threadIdx)*0.1-0.05;
    vel.y = m_multyRand->get(threadIdx)*0.1-0.05;
    vel.z = m_multyRand->get(threadIdx)*0.1-0.05;
    return vel;
}

void ParticlesSolver::createRandomParticles(int threadIdx)
{
    int curNum = m_numParticles[threadIdx];
    if  (m_multyRand->get(threadIdx)>0.9)
    {
        int numToAdd = std::min(MAXPARTICLES, MAXPARTICLES - m_numParticles[threadIdx]);//std::min(int(rand()*1.0/RAND_MAX * 300.0), maxParticles - numParticles-1);
        m_numParticles[threadIdx] += numToAdd;
        for (int i = curNum; i < curNum + numToAdd; ++i)
        {
            vec4 pos=initPos(threadIdx);
            m_bodyPos[threadIdx][i].x=pos.x;
            m_bodyPos[threadIdx][i].y=pos.y;
            m_bodyPos[threadIdx][i].z=pos.z;
            m_bodyPos[threadIdx][i].rho=pos.rho;
            m_bodyPos[threadIdx][i].r=pos.r;
            m_bodyPos[threadIdx][i].count=pos.count;
            m_bodyPos[threadIdx][i].partNumber = totalNum_save;
            totalNum_save++;
            m_bodyVel[threadIdx][i] = interpolateValue(m_currFile[threadIdx], m_bodyPos[threadIdx][i], m_fieldsLoader->m_u, m_fieldsLoader->m_v, m_fieldsLoader->m_w);
            m_bodyTurbVel[threadIdx][i].x=0.0;
            m_bodyTurbVel[threadIdx][i].y=0.0;
            m_bodyTurbVel[threadIdx][i].z=0.0;
            m_bodyAccel[threadIdx][i].x=0.0;
            m_bodyAccel[threadIdx][i].y=0.0;
            m_bodyAccel[threadIdx][i].z=0.0;


        }
        if (threadIdx==0)
        {
            for (int i=curNum;i<m_numParticles[threadIdx];i++)
                num_save[i]=0;
        }
    }

}

void ParticlesSolver::deleteParticle(int particlesIdx, int threadIdx)
{
    //printf("delete\n");
    m_bodyPos[threadIdx][particlesIdx] = m_bodyPos[threadIdx][m_numParticles[threadIdx]-1];
    m_bodyVel[threadIdx][particlesIdx] = m_bodyVel[threadIdx][m_numParticles[threadIdx]-1];
    m_bodyTurbVel[threadIdx][particlesIdx] = m_bodyTurbVel[threadIdx][m_numParticles[threadIdx]-1];
    m_bodyAccel[threadIdx][particlesIdx] = m_bodyAccel[threadIdx][m_numParticles[threadIdx]-1];
    if(threadIdx == 0)
    {
        for(int i = 0; i< num_save[m_numParticles[threadIdx]-1]; i++)
        {
            x_save[particlesIdx][i]=x_save[m_numParticles[threadIdx]-1][i];
            y_save[particlesIdx][i]=y_save[m_numParticles[threadIdx]-1][i];
            z_save[particlesIdx][i]=z_save[m_numParticles[threadIdx]-1][i];
            u_save[particlesIdx][i]=u_save[m_numParticles[threadIdx]-1][i];
            v_save[particlesIdx][i]=v_save[m_numParticles[threadIdx]-1][i];
            w_save[particlesIdx][i]=w_save[m_numParticles[threadIdx]-1][i];
        }
        num_save[particlesIdx] = num_save[m_numParticles[threadIdx]-1];
    }
    m_numParticles[threadIdx]-= 1;
}

void ParticlesSolver::saveParticlesInFile(double time)
{
    char filename[64];
    sprintf(filename, "particles%2f.txt", time);
    FILE *file_data=fopen(filename,"w");
    for (int k = 0; k < THREADNUM; k++)
    {
        for (int i = 0; i < m_numParticles[k]; i++)
        {
            fprintf(file_data,"%d %lf %lf %lf\n",m_bodyPos[k][i].partNumber, m_bodyPos[k][i].x, m_bodyPos[k][i].y, m_bodyPos[k][i].z);
        }
    }
    fclose(file_data);
}
