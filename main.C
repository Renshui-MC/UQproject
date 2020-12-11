#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include "tensor.H"
#include "MyUQ.H"//note only include the subclass (superclass) then everything is included, i.e. the parent class

using namespace Foam;

int main()
{
	
	UQ uncertainty;


	uncertainty.EigenSpace();
	
	//uncertainty.PerturbedAij();
	
	


	
	 int intArray[7] = {5, 3, 32, -1, 1, 104, 53};
	 std::sort(intArray, intArray+7);
	 	for (int value : intArray)
		std::cout << value << std::endl;
  


	tensor T(1, -3, 3, 3, -5, 3, 6, -6, 4);
	Info <<"T= "<< T << endl;
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
	
    tensor t6(1,0,-4,0,5,4,-4,4,3);
    vector e = eigenValues(t6);
    Info<< "eigenvalues " << e << endl;

    tensor ev = eigenVectors(t6);
    Info<< "eigenvectors " << ev << endl;
    Info<< "xx()= "<<ev.xx()<<" xy()= "<<ev.xy()<<" xz()= "<<ev.xz()<<"\n";
    Info<< "yx()= "<<ev.yx()<<" yy()=  "<<ev.yy()<<" yz()= "<<ev.yz()<<"\n";
	Info<< "zx()= "<<ev.zx()<<" zy()=  "<<ev.zy()<<" zz()= "<<ev.zz()<<endl;

    Info<< "Check eigenvectors "
        << (eigenVectors(t6, e) & t6) << " " 
        << (e.x() * eigenVectors(t6, e)).x()
        << (e.x() * eigenVectors(t6, e)).y()
        << (e.x() * eigenVectors(t6, e)).z()
        << endl;    

        Info<< symm(t6) <<endl;
        
    
}
Info<< nl;
	
	//Info<< "tensor" << T << endl;
    //vector e = eigenValues(T);
    //Info<< "eigenvalues " << e << endl;

	//tensor ev = eigenVectors(T);
    //Info<< "eigenvectors " << ev << endl;

	std::cin.get();

}

