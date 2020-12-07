#include <iostream>
#include "tensor.H"//must include a .H file that is defined in the Foam namespace
#include "Numerics.H"

namespace Foam
{ 
Numerics::Numerics(void)
{
   Info <<"Minghan is in Constructor"<<endl;
}

Numerics::~Numerics(void)
{
    Info <<"Minghan is in Destructor"<<endl;
}

void Numerics::InitLog()
{
    Info <<"Minghan is inside log with namespace"<<endl;
}
/*void Numerics::EigenRecomposition(double** A_ij, double** Eig_Vec, double* Eig_Val, unsigned short n)
{
    unsigned short i,j,k;
    double** tmp = new double* [n];
    double** deltaN = new double* [n];

}*/

/*void InitLog()
{
	std::cout << "Yes I am in after dozens of tests!!!!" << std::endl;
}*/
}//End namespace foam
