#include "tools.h"
#include "particlessolver.h"

#define VIEW_U 0
#define VIEW_V 1
#define VIEW_W 2

int itn=0;
int clear_w=1.0;
float rx=0;
float ry=0;
int mx0,my0;
int rotate=0;
float rx0=0;
float ry0=0;
double d_temp;
double mouse_x,mouse_y;
double r_2g=0.0;
int kCur= 0;
double sc=1.0;
int view=VIEW_U;
int redr=0;
double ck=2.0;
double cv=0.001;
bool clearc=true;
void display(void);
void init();

void sweep();
void sweepInit();

ParticlesSolver* partSolver;
int NX, NY, NZ;
double* x;
double* y;
double* z;

void display(void)
{

    if (redr==1)
    {
        for (int i=0;i<1;i++)
            sweep();
    }

    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();

    glTranslatef(x[NX/2],z[NZ/2],y[NY/2]);
    glScalef(sc,sc,sc);
    glRotatef(-90,1,0,0);
    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);
    //glRotatef(90,1,0,0);
    glTranslatef(-x[NX/2],-z[NZ/2],-y[NY/2]);

    glColor3f(1,1,1);
    double l_2;

    for (int i=0;i<NX-1;i++)
    {
        /*glBegin(GL_TRIANGLE_STRIP);


        for (int j=0;j<NY;j++)
        {
            if (view==VIEW_U)
                l_2=ck*(u[i][j][kCur]);
            if (view==VIEW_V)
                l_2=ck*(v[i][j][kCur]);
            if (view==VIEW_W)
                l_2=ck*(w[i][j][kCur]);
            glColor3f(l_2,l_2,-l_2);
            glVertex3f(x[i],y[j],z[kCur]);

            if (view==VIEW_U)
                l_2=ck*(u[i+1][j][kCur]);
            if (view==VIEW_V)
                l_2=ck*(v[i+1][j][kCur]);
            if (view==VIEW_W)
                l_2=ck*(w[i+1][j][kCur]);

            glColor3f(l_2,l_2,-l_2);
            glVertex3f(x[i+1],y[j],z[kCur]);
        }

        glEnd();*/
    }

    double scale =150;


    double dxDepos = (x[NX - 1] - x[0]) / (NXDepos-1);
    double dyDepos = (y[NY - 1] - y[0]) / (NYDepos-1);
    double maxSumMass = - 1e-10;
    double sumMass[NXDepos-1][NXDepos-1];
    for( int i=0; i < NXDepos-1; i++ )
    {
        for( int j=0; j < NYDepos-1; j++ )
        {
            sumMass[i][j] = 0;
            for( int k=0; k < THREADNUM; k++ )
            {
                sumMass[i][j] += partSolver->m_deposPart[k][i][j];
            }
            if(sumMass[i][j] > maxSumMass)
                maxSumMass = sumMass[i][j];
        }
    }

    glBegin(GL_QUADS);
    for( int i=0; i < NXDepos-1; i++ )
    {
        for( int j=0; j < NYDepos-1; j++ )
        {

            glColor3f(sumMass[i][j]  / maxSumMass ,0,0);
            glVertex3f(scale*i * dxDepos ,  scale*z[0] - 0.1, scale*j * dyDepos);
            glVertex3f(scale*i * dxDepos, scale*z[0] - 0.1 , scale*(j+1) * dyDepos);
            glVertex3f(scale*(i+1) * dxDepos , scale*z[0] - 0.1, scale*(j+1) * dyDepos);
            glVertex3f(scale*(i+1) * dxDepos , scale*z[0] - 0.1, scale*j * dyDepos);
        }
    }
    glEnd();

    glEnable(GL_BLEND);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(2.5);

    double vel_scale=0.0;
    for( int i=0; i<partSolver->m_numParticles[0]; i++ ) {
        for( int j=1; j<partSolver->num_save[i]; j++ )
        {
            vel_scale+=(partSolver->u_save[i][j]*partSolver->u_save[i][j]+partSolver->v_save[i][j]*partSolver->v_save[i][j]+partSolver->w_save[i][j]*partSolver->w_save[i][j]);
        }
        //vel_scale/=num_save[i];
    }
    //printf("numPart=%d\n", numParticles[0]);
    vel_scale=sqrt(vel_scale/partSolver->m_numParticles[0]+0.00001);


    for( int i=0; i<partSolver->m_numParticles[0]; i++ )
    {
        glBegin(GL_POINTS);//(GL_LINE_STRIP);
        for( int j=0; j<partSolver->num_save[i]; j++ )
        {
            glColor3f(ck*fabs(partSolver->u_save[i][j])/vel_scale,ck*fabs(partSolver->v_save[i][j])/vel_scale,ck*fabs(partSolver->w_save[i][j])/vel_scale);
            glVertex3f(scale*partSolver->x_save[i][j],scale*partSolver->z_save[i][j],scale*partSolver->y_save[i][j]);
        }
        glEnd();
    }


    glLineWidth(2);

    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_LOOP);
    glVertex3f(scale*x[0],scale*z[0],scale*y[0]);
    glVertex3f(scale*x[NX-1],scale*z[0],scale*y[0]);
    glVertex3f(scale*x[NX-1],scale*z[0],scale*y[NY-1]);
    glVertex3f(scale*x[0],scale*z[0],scale*y[NY-1]);
    glEnd();

    glColor3f(1,1,1);
    glBegin(GL_LINE_LOOP);
    glVertex3f(scale*x[0],scale*z[NZ-1],scale*y[0]);
    glVertex3f(scale*x[NX-1],scale*z[NZ-1],scale*y[0]);
    glVertex3f(scale*x[NX-1],scale*z[NZ-1],scale*y[NY-1]);
    glVertex3f(scale*x[0],scale*z[NZ-1],scale*y[NY-1]);
    glEnd();

    glPointSize(3.0);
    glLineWidth(1.0);
    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();

}

void m_m(int x,int y) //mouse move
{
    if (rotate==1)
    {
        rx=rx0+0.05*(x-mx0);
        ry=ry0+0.05*(y-my0);
    }
    glutPostRedisplay();
}

void m_d(int button, int state,int x, int y)  //mouse down
{
    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;
    }

    mouse_x=(1.0*x)/W_WIDTH;
    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;
    glutPostRedisplay();
}


void kb(unsigned char key, int x, int y)
{
    if (key==']')
    {
        if (kCur<NZ)
            kCur++;
    }

    if (key=='[')
    {
        if (kCur>0)
            kCur--;
    }

    if (key=='q')
    {
        sc*=1.1;
    }

    if (key=='e')
    {
        sc/=1.1;
    }

    if (key=='.')
    {
        ck*=1.1;
    }

    if (key==',')
    {
        ck/=1.1;
    }

    if (key=='1')
    {
        view=VIEW_U;
        printf("viewing U \n");
    }

    if (key=='2')
    {
        view=VIEW_V;
        printf("viewing V \n");
    }

    if (key=='2')
    {
        view=VIEW_W;
        printf("viewing W \n");
    }

    if (key=='m')
    {


    }

    if (key=='n')
    {
        clearc=!clearc;
    }

    if (key=='s')
    {
        sweep();
    }

    if (key==' ')
    {
        redr=!redr;
    }


    glutPostRedisplay();
}

void sweep_init()
{
    partSolver = new ParticlesSolver();
    NX = partSolver->m_NX;
    NY = partSolver->m_NY;
    NZ = partSolver->m_NZ;
    x = partSolver->m_x;
    y = partSolver->m_y;
    z = partSolver->m_z;
}

void sweep()
{
    partSolver->sweep();
}

void init()
{
    sweep_init();
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(x[0]-1000, x[NX-1]+1000,  z[0]-1000, z[NZ-1]+1000,y[0]-100000, y[NY-1]+100000);
    glMatrixMode (GL_MODELVIEW);
}

int main(int argc, char** argv)
{

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT,W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
}
