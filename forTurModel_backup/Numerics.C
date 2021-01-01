#include <iostream>
#include <cstdlib>
#include "tensor.H"//must include a .H file that is defined in the Foam namespace
#include "Numerics.H"

namespace Foam
{ 
Numerics::Numerics(void)//try to remove void later!
{

    unsigned short i,j;

    /* --- Initializing variables for the UQ methodology --- */
    m_delta = NULL;
    m_delta3 = NULL;

    m_delta = new double* [nDim];
    for (i = 0; i < nDim; i++)
    {
        m_delta[i] = new double [nDim];
    }

    for (i = 0; i < nDim; i++)
    {
        for (j = 0; j < nDim; j++)
        {
            if (i == j) m_delta[i][j] = 1.0;
            else m_delta[i][j] = 0.0;
        }
    }

    m_delta3 = new double* [3];
    for (i = 0; i < 3; i++)
    {
        m_delta3[i] = new double [3];
    }

    for (i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            if (i == j) m_delta3[i][j] = 1.0;
            else m_delta3[i][j] = 0.0;
        }
    }

    m_MeanReynoldsStress    = new double* [3];
    m_MeanPerturbedRSM      = new double* [3];
    m_PerturbedStrainRate   = new double* [nDim];
    m_B_ij                  = new double* [3];
    m_newB_ij               = new double* [3];
    m_Eig_Vec               = new double* [3];
    m_New_Eig_Vec           = new double* [3];
    m_Corners               = new double* [3];
    m_Eig_Val               = new double [3];
    m_Barycentric_Coord     = new double [2];
    m_New_Coord             = new double [2];
    for (i = 0; i < 3; i++)
    {
        m_MeanReynoldsStress[i]     = new double [3];
        m_MeanPerturbedRSM[i]       = new double [3];
        m_PerturbedStrainRate[i]    = new double [nDim];
        m_B_ij[i]                   = new double [3];
        m_newB_ij[i]                = new double [3];
        m_Eig_Vec[i]                = new double [3];
        m_New_Eig_Vec[i]            = new double [3];
        m_Corners[i]                = new double [2];//i=3 means 3 row pointers and [2] means actual two elements each row
        m_Eig_Val[i]                = 0; //initializing 
        
        if(i < 2)
        {
            m_Barycentric_Coord[i]  = 0;
            m_New_Coord[i]          = 0;
        }
    }

    for (i = 0; i < 3; i++)
	{  
		for (j = 0; j < 3; j++)
		{
			m_MeanReynoldsStress[i][j] 		= 0.0;
			m_MeanPerturbedRSM[i][j]		= 0.0;
            m_PerturbedStrainRate[i][j]     = 0.0;
			m_B_ij[i][i]                    = 0.0;
            m_newB_ij[i][j]					= 0.0;
			m_Eig_Vec[i][j] 		   		= 0.0;
			m_New_Eig_Vec[i][j] 			= 0.0; 
		}
		
	}

    /* define barycentric traingle corner points */
    m_Corners[0][0] = 1.0;
    m_Corners[0][1] = 0.0;
    m_Corners[1][0] = 0.0;
    m_Corners[1][1] = 0.0;
    m_Corners[2][0] = 0.5;
    m_Corners[2][1] = 0.866025;


}

Numerics::~Numerics(void)
{

    if (m_delta != NULL)
    {
        for (unsigned short i = 0; i < nDim; i++)
        {
            delete [] m_delta[i];
        }
        delete [] m_delta;
    }
    
    if (m_delta3 != NULL)
    {
        for (unsigned short i = 0; i < 3; i++)
        {
            delete [] m_delta3[i];
        }
        delete [] m_delta3;
    }

    for (unsigned short i = 0; i < 3; i++)
    {
        delete [] m_MeanReynoldsStress[i];
        delete [] m_MeanPerturbedRSM[i];
        delete [] m_PerturbedStrainRate[i];
        delete [] m_B_ij[i];
        delete [] m_newB_ij[i];
        delete [] m_Eig_Vec[i];
        delete [] m_New_Eig_Vec[i];
        delete [] m_Corners[i];
    }
        delete [] m_MeanReynoldsStress;
        delete [] m_MeanPerturbedRSM;
        delete [] m_PerturbedStrainRate;
        delete [] m_B_ij;
        delete [] m_newB_ij;
        delete [] m_Eig_Vec;
        delete [] m_New_Eig_Vec;
        delete [] m_Corners;
        delete [] m_Eig_Val;
        delete [] m_Barycentric_Coord;
        delete [] m_New_Coord;
}

void Numerics::EigenRecomposition(double** m_B_ij, double** m_Eig_Vec, double* m_Eig_Val, unsigned short n)
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
            //Info<< "deltaN["<<i<<"]["<<j<<"]= "<<deltaN[i][j]<<endl;
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            tmp[i][j] = 0.0;//initializing every entry with 0
            for (k = 0; k < n; k++)
            {
                tmp[i][j] += m_Eig_Vec[i][k] * m_Eig_Val[k] * deltaN[k][j];
                //Info<< "tmp["<<i<<"]["<<j<<"]= "<<tmp[i][j]<<endl;
                //Info<< "Eig_Vec_init["<<i<<"]["<<k<<"]= "<<Eig_Vec[i][k]<<endl;
                //Info<< "Eig_Val_init["<<k<<"]= "<<Eig_Val[k]<<endl;
                //Info<< "DeltaN_init["<<k<<"]["<<j<<"]= "<<deltaN[i][k]<<endl;

            }
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            m_B_ij[i][j] = 0.0;
            for (k = 0; k < n; k++)
            {
                m_B_ij[i][j] += tmp[i][k] * m_Eig_Vec[j][k];
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
