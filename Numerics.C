#include <iostream>
#include <cstdlib>
#include "tensor.H"//must include a .H file that is defined in the Foam namespace
#include "Numerics.H"

namespace Foam
{ 
Numerics::Numerics(void)//try to remove void later!
{
   Info <<"Minghan is in Constructor"<<endl;
   
    m_A_ij              = new double* [3];
    m_newA_ij           = new double* [3];
    m_Eig_Vec           = new double* [3];
    m_New_Eig_Vec       = new double* [3];
    Corners             = new double* [3];
    m_Eig_Val           = new double [3];
    Barycentric_Coord   = new double [2];
    New_Coord           = new double [2];
    for (unsigned short i = 0; i < 3; i++)
    {
        m_A_ij[i]           = new double [3];
        m_newA_ij[i]        = new double [3];
        m_Eig_Vec[i]        = new double [3];
        m_New_Eig_Vec[i]    = new double [3];
        Corners[i]          = new double [2];//i=3 means 3 row pointers and [2] means actual two elements each row
        m_Eig_Val[i]        = 0; //initializing 
    }
    /* define barycentric traingle corner points */
    Corners[0][0] = 1.0;
    Corners[0][1] = 0.0;
    Corners[1][0] = 0.0;
    Corners[1][1] = 0.0;
    Corners[2][0] = 0.5;
    Corners[2][1] = 0.866025;
}

Numerics::~Numerics(void)
{
    Info <<"Minghan is in Destructor"<<endl;
    for (unsigned short i = 0; i < 3; i++)
    {
        delete [] m_A_ij[i];
        delete [] m_newA_ij[i];
        delete [] m_Eig_Vec[i];
        delete [] m_New_Eig_Vec[i];
        delete [] Corners[i];
    }
        delete [] m_A_ij;
        delete [] m_newA_ij;
        delete [] m_Eig_Vec;
        delete [] m_New_Eig_Vec;
        delete [] Corners;
        delete [] m_Eig_Val;
        delete [] Barycentric_Coord;
        delete [] New_Coord;
}

void Numerics::EigenRecomposition(double** A_ij, double** Eig_Vec, double* Eig_Val, unsigned short n)
{
    unsigned short i, j, k;//Note for loop can recognize i, j, k inside the scope
    double** tmp = new double* [n];
    double** deltaN = new double* [n];
    //double small = 10e-10;

    for(i = 0; i < n; i++)
    {
        tmp[i] = new double [n];
        deltaN[i] = new double [n];
    }

    for(i = 0; i < n; i++)//I matrix
    {
        for(j = 0; j < n; j++)
        {
            if (i == j) deltaN[i][j] = 1.0;
            else deltaN[i][j] = 0.0;
            //std::cout << "deltaN["<<i<<"]["<<j<<"]= "<<deltaN[i][j]<<std::endl;
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            tmp[i][j] = 0.0;//initializing every entry with 0
            for (k = 0; k < n; k++)
            {
                tmp[i][j] += Eig_Vec[i][k] * Eig_Val[k] * deltaN[k][j];
                //std::cout << "tmp["<<i<<"]["<<j<<"]= "<<tmp[i][j]<<std::endl;
                //std::cout << "Eig_Vec_init["<<i<<"]["<<k<<"]= "<<Eig_Vec[i][k]<<std::endl;
                //std::cout << "Eig_Val_init["<<k<<"]= "<<Eig_Val[k]<<std::endl;
                //std::cout << "DeltaN_init["<<k<<"]["<<j<<"]= "<<deltaN[i][k]<<std::endl;

            }
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A_ij[i][j] = 0.0;
            for (k = 0; k < n; k++)
            {
                A_ij[i][j] += tmp[i][k] * Eig_Vec[j][k];
            }
        }
        
    }

    for (i = 0; i < n; i++)
    {
        delete [] tmp[i];
        delete [] deltaN[i];
    }
    delete [] tmp;
    delete [] deltaN;
}


}//End namespace foam
