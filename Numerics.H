#pragma once

#include <iostream>


namespace Foam{

class Numerics 
{
     /*The following declarations are defined in the "Numerics.H" (base) however are actually used in the subclasses*/
      protected:
        double **m_MeanReynoldsStress, **m_MeanPerturbedRSM, **m_PerturbedStrainRate; //corresponding to <uiuj>, <uiuj>* and <Sij>* (* means perturbation)
       
        double **m_delta;   /*!< \brief Identity matrix. */
        double **m_delta3;  /*!< \brief 3 row Identity matrix. */

      

        double **m_B_ij, **m_Eig_Vec, **m_newB_ij, **m_New_Eig_Vec, **m_Corners; //2D array
        double *m_Eig_Val, *m_Barycentric_Coord, *m_New_Coord; //1D array
 
        const unsigned short nDim = 3;//nDim could be 2 or 3. However for research it should always be 3

      protected:
      
        /*constructor*/
        Numerics(void); //try without specifying the type void later!
 

        /*destructor*/
        virtual ~Numerics(void);

        /*member functions (note m_Aij, m_Eig_Vec and m_Eig_Val could be any names!)*/
        static void EigenRecomposition(double** m_B_ij, double** m_Eig_Vec, double* m_Eig_Val, unsigned short n);

};

}//End namespace foam



