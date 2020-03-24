#ifndef FIELDSLOADER_H
#define FIELDSLOADER_H

#include "tools.h"

class FieldsLoader
{
public:
    FieldsLoader();
    ~FieldsLoader();

    void loadVTK();
    void saveVTK(char *name);

    vtkStructuredGridReader* m_reader;
    vtkStructuredGrid* m_structured_grid;
    bool m_isCreated=false;
    int* m_dims;
    double* m_x;
    double* m_y;
    double* m_z;
    double*** m_u;
    double*** m_v;
    double*** m_w;
    double*** m_uSig;
    double*** m_vSig;
    double*** m_wSig;
    double*** m_uPrev;
    double*** m_vPrev;
    double*** m_wPrev;
    double*** m_uSigPrev;
    double*** m_vSigPrev;
    double*** m_wSigPrev;
    int m_readFileNum;
    std::vector <char*>  m_fileNames;
};

#endif // FIELDSLOADER_H
