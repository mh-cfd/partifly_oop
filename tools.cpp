#include "tools.h"

double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}



int strCompare(const void * a, const void * b)
{
  return strcmp(*(char**)a,*(char**)b);
}
MultythreadRandom::MultythreadRandom()
{
    init();
}

void MultythreadRandom::init()
{
    for (int i=0;i<THREADNUM;i++)
    {
        xr[i] = 123456789+i;
        yr[i] = 362436069+2*i;
        zr[i] = 521288629+3*i;
        wr[i] = 88675123+5*i;
        seeds[i]=i*10;
        srand48_r(time(NULL), &(randBuffers[i]));
    }
}

float MultythreadRandom::get(int threadIdx) {
    /* uint32_t t;
    t = xr[i] ^ (xr[i] << 11);
    xr[i] = yr[i]; yr[i] = zr[i]; zr[i] = wr[i];
    wr[i] = wr[i] ^ (wr[i] >> 19) ^ (t ^ (t >> 8));
    return wr[i]*0.5/0x7FFFFFFF;*/
    // return std::uniform_real_distribution<double> zeroToOne(0.0, 1.0);
    // return rand_r(&(seeds[i]))*1.0/RAND_MAX;
    double res;
    drand48_r(&randBuffers[threadIdx], &res);
    return res;//drand48(); // jrand48()
}

DistributedRandom::DistributedRandom()
{
    dist ();
    integ();
    interpol ();
    m_rand = new MultythreadRandom();
}

void DistributedRandom::dist ()
{
    for (int i=0;i<RNX;i++)
    {
        m_pdf_x[i].x=6*i*1.0/RNX-3;
        m_pdf_x[i].y=exp(-(m_pdf_x[i].x)*(m_pdf_x[i].x));
    }
}

void DistributedRandom::integ (){

    m_int_x[0].x=m_pdf_x[0].x;
    m_int_x[0].y=0.0;
    for (int i=1;i<RNX;i++){

        double dx=(m_pdf_x[i].x-m_pdf_x[i-1].x)/RNX;
        m_int_x[i].x=m_pdf_x[i].x;
        m_int_x[i].y=m_int_x[i-1].y+dx*0.5*(m_pdf_x[i].y+m_pdf_x[i-1].y);
    }
    double max=m_int_y[0].y;
    for (int i=1;i<RNX;i++){
        if(max<m_int_x[i].y)
            max=m_int_x[i].y;
    }
    for (int i=0;i<RNX;i++)
        m_int_x[i].y/=max;
}

void DistributedRandom::interpol ()
{
    for (int i=0;i<RNY;i++){
        m_int_y[i].y=i*1.0/(RNY-1);
        int j=0;
        while (m_int_y[i].y>m_int_x[j].y)
        {
            j++;
        }
        double alpha=(-m_int_y[i].y+m_int_x[j].y)/(m_int_x[j].y-m_int_x[j-1].y);
        m_int_y[i].x=m_int_x[j-1].x*(alpha)+m_int_x[j].x*(1.0-alpha);
    }
    m_int_y[0].x=m_int_x[0].x;
    m_int_y[RNY-1].x=m_int_x[RNX-1].x;
}

double DistributedRandom::get(int threadIdx)
{
    double rand_num = m_rand->get(threadIdx);
    double dy=1.0/RNY;
    double alpha=rand_num/dy;
    int n=(int)alpha;
    alpha-=n;
    return m_int_y[n].x*(1.0-alpha)+m_int_y[n+1].x*alpha;
}
