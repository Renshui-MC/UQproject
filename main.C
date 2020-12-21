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
	scalar k = 0.0090;
	scalar nut = 0.5;
	symmTensor Aij_OF;
	symmTensor uiuj_OF;
	symmTensor Sij_OF;

	symmTensor Aij[5];
	symmTensor Rij[5];
	symmTensor perturbed_Aij[5];
	symmTensor perturbed_uiuj[5];
	symmTensor perturbed_Sij[5];

	unsigned short i, x, y;
	/* newAij could be any name and passed by reference by calling "EigenSpace" function
	   Note newAij must be exactly same as defined in the Numerics.H otherwise an linking
	   error will be generated */
	double** newAij 	= new double* [3];
	double** newuiuj 	= new double* [3];
	double** newSij		= new double* [3];
	for (i = 0; i < 3; i++)
	{
		newAij[i] 	= new double [3];
		newuiuj[i]	= new double [3];
		newSij[i]   = new double [3];
	} 


	Aij[0] = {1, 0, -4, 5, 4, 3};//to read in fields
	Rij[0] = {0, 0, 0, 0, 0, 0};


	//read in Reynold stress field
	Info<<"*********** Start reading anisotropy tensor in each cell **************\n\n";
	for (i = 1; i < 5; i++)
	{
		Aij[i] = {0.6592491, 0.04143975, 0, -0.3259158, 0, -0.333333};
		Rij[i] = {0.0179, 7.4577e-4, 0, 1.3349e-4, 0, 0};
		//Aij[i] = {1, 1, -1, 1, 2, 0, -1, 0, 5};
		//Aij[i] = {3, -2, 4, -2, 6, 2, 4, 2, 3};
	}



	for (i = 0; i < 5; i++)
	{
	 UQ uncertainty(k, nut, Rij[i], i);
	 
	 uncertainty.EigenSpace(Aij_OF, uiuj_OF, Sij_OF, 
	 						newAij, newuiuj, newSij,
	 						Aij[i], i);

	 Info<<"************** celli= "<<i<<" inside MAIN! **************"<<endl;
	 	for (x = 0; x < 3; x++)
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"newAij"<<"["<< x <<"]"<<"["<<y<<"]= "<< newAij[x][y]<<endl; 
			}
		}

		for (x = 0; x < 3; x++)
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"newuiuj"<<"["<< x <<"]"<<"["<<y<<"]= "<< newuiuj[x][y]<<endl; 
			}
		}

		for (x = 0; x < 3; x++)
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"newSij"<<"["<< x <<"]"<<"["<<y<<"]= "<< newSij[x][y]<<endl; 
			}
		}


		perturbed_Aij[i] 	= Aij_OF;
		perturbed_uiuj[i] 	= uiuj_OF;
		perturbed_Sij[i]    = Sij_OF;
		Info<<"perturbed_Aij= "<<perturbed_Aij[i]<<endl;
		Info<<"perturbed_uiuj= "<<perturbed_uiuj[i]<<endl;
		Info<<"perturbed_Sij= "<<perturbed_Sij[i]<<endl;

	 

	}
	
	Info<<"************ UQ loop has ended! Check the store perturbed Aij ************"<<endl;
	Info<<"perturbed_Aij[0]= "<<perturbed_Aij[0]<<"\n";
	Info<<"perturbed_Aij[1]= "<<perturbed_Aij[1]<<"\n";
	Info<<"perturbed_Aij[1]= "<<perturbed_Aij[2]<<"\n";
	Info<<"perturbed_Aij[1]= "<<perturbed_Aij[3]<<"\n";
	Info<<"perturbed_Aij[2]= "<<perturbed_Aij[4]<<"\n";
	Info<<"************ End MAIN *****************************************"<<endl;
	
	for (i = 0; i < 3; i++)
	{
		delete [] newAij[i];
		delete [] newuiuj[i];
	}
	delete [] newAij;
	delete [] newuiuj;

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

        Info<<"symmetric tensor t6= "<< symm(t6) <<endl;
        Info<<"two symm t6= "<< twoSymm(t6) <<endl;
		Info<<"dev(t6)= "<<dev(t6)<<endl;
		Info<<"dev2(t6)= "<<dev2(t6)<<endl;
}
Info<< nl;



	std::cin.get();

}

