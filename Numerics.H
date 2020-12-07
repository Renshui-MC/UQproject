#pragma once

#include <iostream>


namespace Foam{

class Numerics 
{
    protected:
        /*constructor*/
        Numerics(void);

        /*Constructor*/
        virtual ~Numerics(void);

        virtual void InitLog();

       // void EigenRecomposition(double **A_ij, double **Eig_Vec, double *Eig_Val, unsigned short n);

};

}//End namespace foam



