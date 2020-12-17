#pragma once

#include <iostream>


namespace Foam{

class Numerics 
{
     /*The following declarations are defined in the "Numerics.H" (base) however are actually used in the subclasses*/
      protected:
        double **m_MeanReynoldsStress;
        double **m_MeanPerturbedRSM;
        double **m_delta3;

        double **m_B_ij, **m_Eig_Vec, **m_newB_ij, **m_New_Eig_Vec, **m_Corners; //2D array
        double *m_Eig_Val, *m_Barycentric_Coord, *m_New_Coord; //1D array
 
      protected:
        /*constructor*/
        Numerics(void); //try without specifying the type void later!
        

        /*destructor*/
        virtual ~Numerics(void);

        /*member functions (note m_Aij, m_Eig_Vec and m_Eig_Val could be any names!)*/
        static void EigenRecomposition(double** B_ij, double** Eig_Vec, double* Eig_Val, unsigned short n);

};

}//End namespace foam



