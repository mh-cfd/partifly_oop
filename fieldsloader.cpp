#include "fieldsloader.h"

FieldsLoader::FieldsLoader()
{
    m_readFileNum = 0;
    std::string  nameFile,nameForOpen;
    std::string nameDir("/mnt/storage/MicroHH/RFFI_ATMOSPH/av3D/mean/vtk/");//"/mnt/storage/MicroHH/RFFI_ATMOSPH/av3D/mean/vtk/"///"/mnt/storage/MicroHH/boris/rb/vtk3D/"
    DIR *mydir = opendir(nameDir.data());
    if(mydir == NULL) {
        perror("opendir");
        //return -1;
    }
    struct dirent *entry;
    while ((entry = readdir(mydir)))
    {
        int len = strlen (entry->d_name);
        if (len >= 4) {
            if (strcmp (".vtk", entry->d_name + len - 4) == 0)
            {
                nameFile=entry->d_name;
                nameForOpen = nameDir+ nameFile;
                char * str = (char*)malloc(sizeof(char)*strlen((char*)(nameForOpen.data())));
                strcpy(str,(char*)(nameForOpen.data()));
                m_fileNames.push_back(str);
                printf("FilesName= %s\n",str);
            }
        }
    }
    char *sort_char_array[m_fileNames.size()];
    for (int i = 0 ;i < m_fileNames.size();i++)
        sort_char_array[i] = m_fileNames[i];
    qsort (sort_char_array, m_fileNames.size(), sizeof(sizeof(char**)), strCompare);
    for (int i =0 ;i<m_fileNames.size();i++)
    {
        m_fileNames[i] = sort_char_array[i];
        printf("FilesName= %s\n",m_fileNames[i]);
    }
    closedir(mydir);
}

FieldsLoader::~FieldsLoader()
{
    for (int i =0 ;i<m_fileNames.size();i++)
    {
        free(m_fileNames[i]);
    }
    delete m_x;
    delete m_y;
    delete m_z;
    for (int i=0;i<m_dims[0];i++)
    {
        for (int j=0;j<m_dims[1];j++)
        {
            delete m_u[i][j];
            delete m_v[i][j];
            delete m_w[i][j];
            delete m_uPrev[i][j];
            delete m_vPrev[i][j];
            delete m_wPrev[i][j];
        }
        delete m_u[i];
        delete m_v[i];
        delete m_w[i];
        delete m_uPrev[i] ;
        delete m_vPrev[i] ;
        delete m_wPrev[i] ;
    }
    delete m_u;
    delete m_v;
    delete m_w;
    delete m_uPrev;
    delete m_vPrev;
    delete m_wPrev;
    delete m_dims;
}

void FieldsLoader::loadVTK()
{
    printf("Load VTK file: %s\n", m_fileNames[m_readFileNum]);
    if(!m_isCreated)
    {
        m_reader = vtkStructuredGridReader::New();
    }
    else
    {
        for (int i=0;i<m_dims[0];i++)
        {
            for (int j=0;j<m_dims[1];j++)
            {
                for (int k=0;k<m_dims[2];k++)
                {
                    m_uPrev[i][j][k]=m_u[i][j][k];
                    m_vPrev[i][j][k]=m_v[i][j][k];
                    m_wPrev[i][j][k]=m_w[i][j][k];
                    m_uSigPrev[i][j][k]=m_uSig[i][j][k];
                    m_vSigPrev[i][j][k]=m_vSig[i][j][k];
                    m_wSigPrev[i][j][k]=m_wSig[i][j][k];
                }
            }
        }
    }

    vtkDataArray *vel;
    vtkDataArray *vel2;

    m_reader->ReadAllScalarsOn();
    m_reader->ReadAllVectorsOn();
    m_reader->SetFileName(m_fileNames[m_readFileNum]);

    m_reader->Update();
    vel = m_reader->GetOutput()->GetPointData()->GetArray("Vel_mean");
    vel2 = m_reader->GetOutput()->GetPointData()->GetArray("Vel_2");
    m_structured_grid = m_reader->GetOutput();
    m_dims = m_structured_grid->GetDimensions();

    if(!m_isCreated)
    {
        m_x=new double [m_dims[0]];
        m_y=new double [m_dims[1]];
        m_z=new double [m_dims[2]];

        m_u=new double **[m_dims[0]];
        m_v=new double **[m_dims[0]];
        m_w=new double **[m_dims[0]];
        m_uSig = new double **[m_dims[0]];
        m_vSig = new double **[m_dims[0]];
        m_wSig = new double **[m_dims[0]];
        m_uPrev=new double **[m_dims[0]];
        m_vPrev=new double **[m_dims[0]];
        m_wPrev=new double **[m_dims[0]];
        m_uSigPrev=new double **[m_dims[0]];
        m_vSigPrev=new double **[m_dims[0]];
        m_wSigPrev=new double **[m_dims[0]];

        for (int i=0;i<m_dims[0];i++)
        {
            m_u[i]=new double *[m_dims[1]];
            m_v[i]=new double *[m_dims[1]];
            m_w[i]=new double *[m_dims[1]];
            m_uSig[i]=new double *[m_dims[1]];
            m_vSig[i]=new double *[m_dims[1]];
            m_wSig[i]=new double *[m_dims[1]];
            m_uPrev[i]=new double *[m_dims[1]];
            m_vPrev[i]=new double *[m_dims[1]];
            m_wPrev[i]=new double *[m_dims[1]];
            m_uSigPrev[i]=new double *[m_dims[1]];
            m_vSigPrev[i]=new double *[m_dims[1]];
            m_wSigPrev[i]=new double *[m_dims[1]];

            for (int j=0;j<m_dims[1];j++)
            {
                m_u[i][j]=new double [m_dims[2]];
                m_v[i][j]=new double [m_dims[2]];
                m_w[i][j]=new double [m_dims[2]];
                m_uSig[i][j]=new double [m_dims[2]];
                m_vSig[i][j]=new double [m_dims[2]];
                m_wSig[i][j]=new double [m_dims[2]];
                m_uPrev[i][j]=new double [m_dims[2]];
                m_vPrev[i][j]=new double [m_dims[2]];
                m_wPrev[i][j]=new double [m_dims[2]];
                m_uSigPrev[i][j]=new double [m_dims[2]];
                m_vSigPrev[i][j]=new double [m_dims[2]];
                m_wSigPrev[i][j]=new double [m_dims[2]];
                for (int k=0;k<m_dims[2];k++)
                {
                    m_u[i][j][k]=0.0;
                    m_v[i][j][k]=0.0;
                    m_w[i][j][k]=0.0;
                    m_uSig[i][j][k]=0.0;
                    m_vSig[i][j][k]=0.0;
                    m_wSig[i][j][k]=0.0;
                    m_uPrev[i][j][k]=0.0;
                    m_vPrev[i][j][k]=0.0;
                    m_wPrev[i][j][k]=0.0;
                    m_uSigPrev[i][j][k]=0.0;
                    m_vSigPrev[i][j][k]=0.0;
                    m_wSigPrev[i][j][k]=0.0;
                }
            }
        }
        m_isCreated=true;
    }


    m_reader->CloseVTKFile();
    double p[3];

    for (int k = 0; k < m_dims[2]; k++)
    {
        for (int j = 0; j < m_dims[1]; j++)
        {
            for (int i = 0; i < m_dims[0]; i++)
            {
                m_structured_grid->GetPoint(i, j, k, p);
                m_x[i]=p[0];
            }
            m_y[j]=p[1];
        }
        m_z[k]=p[2];
    }


    for (int k = 0; k < m_dims[2]; k++)
    {
        for (int j = 0; j < m_dims[1]; j++)
        {
            for (int i = 0; i < m_dims[0]; i++)
            {
                int ijk=i+j*m_dims[0]+k*m_dims[1]*m_dims[1];
                m_u[i][j][k]=vel->GetComponent(ijk, 0);
                m_v[i][j][k]=vel->GetComponent(ijk, 1);
                m_w[i][j][k]=vel->GetComponent(ijk, 2);
                m_uSig[i][j][k]=sqrt(fabs(vel2->GetComponent(ijk, 0)));
                m_vSig[i][j][k]=sqrt(fabs(vel2->GetComponent(ijk, 1)));
                m_wSig[i][j][k]=sqrt(fabs(vel2->GetComponent(ijk, 2)));
            }
        }
    }
    m_readFileNum++;
}

void FieldsLoader::saveVTK(char *name)
{
    vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
    vtkStructuredGrid *structured_grid_new = vtkStructuredGrid::New();
    structured_grid_new->SetDimensions(m_dims[0], m_dims[1], m_dims[2]);
    vtkPoints *nodes = vtkPoints::New();
    nodes->Allocate(m_dims[0]*m_dims[1]*m_dims[2]);
    vtkDoubleArray* velArray = vtkDoubleArray::New();
    velArray->SetNumberOfComponents(3);
    velArray->SetName("Velocity");

    for (int k = 0; k < m_dims[2]; k++)
    {
        for (int j = 0; j < m_dims[1]; j++)
        {
            for (int i = 0; i < m_dims[0]; i++)
            {
                nodes->InsertNextPoint(m_x[i], m_y[j], m_z[k]);
                velArray->InsertNextTuple3(m_u[i][j][k], m_v[i][j][k], m_w[i][j][k]);
            }
        }
    }

    structured_grid_new->SetPoints(nodes);
    structured_grid_new->GetPointData()->AddArray(velArray);


#if (VTK_MAJOR_VERSION >=6)
    writer->SetInputData(structured_grid_new);
#else
    writer->SetInput(structured_grid_new);
#endif

    writer->SetFileName(name);
    writer->SetFileTypeToBinary();
    writer->Write();

    structured_grid_new->Delete();
    writer->Delete();
    nodes->Delete();
    velArray->Delete();
}

