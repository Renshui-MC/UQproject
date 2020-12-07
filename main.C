#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include "tensor.H"
#include "MyUQ.H"

using namespace Foam;

int main()
{
	
	
	UQ eigenspace;


	eigenspace.EigenvalueSort();
	eigenspace.Eigenvector();
	eigenspace.PerturbedAij();
	
	


	
	 int intArray[7] = {5, 3, 32, -1, 1, 104, 53};
	 std::sort(intArray, intArray+7);
	 	for (int value : intArray)
		std::cout << value << std::endl;
  


	tensor T(1, -3, 3, 3, -5, 3, 6, -6, 4);
	
	Info << "T.inv()= " << T.inv() << endl;
	Info << "inv(T)= " << inv(T) << endl;
	
	{
			Info << "rows: " << nl;
			for (direction i = 0; i < 3; ++i)
		{
			Info << " (" << i << ") = " << T.row(i) << nl;
		}
	}
{
	Info << "cols:" << nl;
	for (direction i=0; i < 3; ++i)
	{
		Info<< " (" << i << ") = " << T.col(i) << nl;
	}
	Info<< "col<0> = " << T.col<0>() << nl;
	Info<< "col<1> = " << T.col<1>() << nl;
	Info<< "col<2> = " << T.col<2>() << nl;

	T.col<1>({0, 2, 4});
	Info<< "replaced col<1> = " << T.col<1>() << nl;
	Info<< "tensor " << T << nl;

	T.row<2>(Zero);
	Info<< "replaced row<2> = " << T.row<2>() << nl;
	Info<< "tensor " << T << nl;
	

	
}
Info<< nl;
	
	//Info<< "tensor" << T << endl;
    //vector e = eigenValues(T);
    //Info<< "eigenvalues " << e << endl;

	//tensor ev = eigenVectors(T);
    //Info<< "eigenvectors " << ev << endl;

	std::cin.get();

}

