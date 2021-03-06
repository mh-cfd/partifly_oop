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
double sc=0.01;
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

double view_x=14.507903;
double view_y=8.300000;
double view_z=24.2;

double o_x=0.0;
double o_y=0.0;
double o_z=0.0;

void display(void)
{   
    int i,j,k,l,ii,jj;
    double l_2,tx,ty,tx0,ty0,vx,vy,v0x,v0y;


    double orient_x=0.0;
    double orient_y=0.0;
    double orient_z=5.0;



    glEnable(GL_DEPTH_TEST);
   glClearColor(0.0, 0.0, 0.0, 0.0);
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

   glLoadIdentity();

    glLineWidth(1);
    glScalef(sc,sc,sc);
   o_x=cos(-ry/20)*cos(-rx/20);
   o_y=cos(-ry/20)*sin(-rx/20);
   o_z=sin(-ry/20);


   //printf("rx=%f ry=%f view_x=%f view_y=%f view_z=%f angle=%f \n",rx,ry,view_x,view_y,view_z,angle);

   orient_x=view_x+o_x;
   orient_y=view_y+o_y;
   orient_z=view_z+o_z;

   gluLookAt(view_x,view_y,view_z,orient_x,orient_y,orient_z,0,0,1);

    //////////////////////////
    if (redr==1)
    {
        for (int i=0;i<1;i++)
            sweep();
    }


    glColor3f(1,1,1);
    //double l_2;




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
            glVertex3f(i * dxDepos ,    j * dyDepos, z[0] - 0.1);
            glVertex3f(i * dxDepos,(j+1) * dyDepos, z[0] - 0.1 );
            glVertex3f((i+1) * dxDepos ,(j+1) * dyDepos, z[0] - 0.1);
            glVertex3f((i+1) * dxDepos ,j * dyDepos, z[0] - 0.1);
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

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glLineWidth(2.5);

    for( int i=0; i<partSolver->m_numParticles[0]; i++ )
    {
        glBegin/*(GL_POINTS);*/(GL_LINE_STRIP);
        for( int j=0; j<partSolver->num_save[i]; j++ )
        {
            glColor3f(ck*fabs(partSolver->u_save[i][j])/vel_scale,ck*fabs(partSolver->v_save[i][j])/vel_scale,ck*fabs(partSolver->w_save[i][j])/vel_scale);
            glVertex3f(partSolver->x_save[i][j],partSolver->y_save[i][j],partSolver->z_save[i][j]);
        }
        glEnd();
    }


    glLineWidth(2);

    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_LOOP);
    glVertex3f(x[0], y[0],  z[0]);
    glVertex3f(x[NX-1],y[0],z[0]);
    glVertex3f(x[NX-1],y[NY-1],z[0]);
    glVertex3f(x[0],   y[NY-1],z[0]);
    glEnd();

    glColor3f(1,1,1);
    glBegin(GL_LINE_LOOP);
    glVertex3f(x[0], y[0],  z[NZ-1]);
    glVertex3f(x[NX-1],y[0],z[NZ-1]);
    glVertex3f(x[NX-1],y[NY-1],z[NZ-1]);
    glVertex3f(x[0], y[NY-1],  z[NZ-1]);
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




    if (key=='w')
      {
        view_x+=(o_x)*100.1;  //sf*=1.1;
        view_y+=(o_y)*100.1;
        view_z+=(o_z)*100.1;
      }

    if (key=='s')
      {

        view_x-=(o_x)*100.1;  //sf*=1.1;
        view_y-=(o_y)*100.1;
        view_z-=(o_z)*100.1;


      }


    if (key=='q')
      {

         view_z+=30.1;

      }

    if (key=='e')
      {

         view_z-=30.1;



      }



    if (key=='a')
      {
        double l2=sqrt(o_y*o_y+o_x*o_x);

        view_y+=(o_x)*100.1/l2;  //sf*=1.1;


        view_x+=-(o_y)*100.1/l2;
      }

    if (key=='d')
      {
        double    l2=sqrt(o_y*o_y+o_x*o_x);
        view_x+=(o_y)*100.1/l2;  //sf*=1.1;

        view_y+=-(o_x)*100.1/l2;


      }





/*    if (key=='m')
    {


    }

    if (key=='n')
    {
        clearc=!clearc;
    }

    if (key=='s')
    {
        sweep();
    }*/

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

    double x0=x[0];
    double y0=y[0];
    double z0=z[0];


    double x1=x[NX-1];
    double y1=y[NY-1];
    double z1=z[NZ-1];



    view_x=x0-0.1*(x1+x0);
    view_y=0.5*(y0+y1);
    view_z=0.5*(z0+z1);

   // gluLookAt(view_x,view_y,view_z,orient_x,orient_y,orient_z,0,1,0);

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

    /* set clear color to black */

    gluPerspective(45.0f, W_WIDTH*1.0/W_HEIGHT, 0.1f, 250.0f);

    // glOrtho(-15.0, 15.0, -15.0*(W_WIDTH*1.0/W_HEIGHT),15.0*(W_WIDTH*1.0/W_HEIGHT), -400.0, 400.0);
    // glFrustum(-15.0, 15.0, -15.0*(W_WIDTH*1.0/W_HEIGHT),15.0*(W_WIDTH*1.0/W_HEIGHT), 0.1, 20.0);
    glMatrixMode (GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);


//    glOrtho(x[0]-1000, x[NX-1]+1000,  z[0]-1000, z[NZ-1]+1000,y[0]-100000, y[NY-1]+100000);
  //  glMatrixMode (GL_MODELVIEW);
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
