#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include "tensor.H"
#include "triad.H"
#include "MyUQ.H"


//using namespace Foam;//using directive namespace has the scope that encloses entire program

namespace Foam{//specifies the Foam namespace for restricted scopes

UQ::UQ()
{
		
	Rij[0] = {1, -3, 3, 3, -5, 3, 6, -6, 4};//to read in fields


	//read in Reynold stress field
	for ( int i = 1; i < 5; i++)
	{
		Rij[i] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
		Info<<"i= "<<i-1<<" Rij= "<<Rij[i-1]<<Rij[i-1].row(1)<<nl;
		if(i==4)
		{
			Info<<"i= "<<i<<" Rij= "<<Rij[i]<<Rij[i].row(1)<<nl;
		}
	}
}

//**********Member functions************//
void UQ::EigenvalueSort ()
{
	for (int i = 0; i < 5; i++)
	{
		e[i] = eigenValues(Rij[i]);//Note here "eigenValues" is the intrinsic function
								   //offered by OpenFoam 
		e_sort[i][0] = e[i].x() ;
		e_sort[i][1] = e[i].y() ;
		e_sort[i][2] = e[i].z() ;
		std::sort(e_sort[i], e_sort[i]+3, std::greater<double>());
		Info<< "i= "<<i<<" eigenvalues= " << e[i]<< endl;	
	}

	for (int i = 0; i < 5; i++)
	{
		for (double value : e_sort[i])//i stands for number of rows
		Info<<"non-decreasing eigenvalue= "<< value << endl;
		
		/*for (int j = 0; j < 3; j++)
		{
			Info<<" e_sort= "<< e_sort[i][j] <<endl;
			
		}*/	
	}
}
	
void UQ::Eigenvector()
{
	for (int i = 0; i < 5; i++)
	{
		ev[i] = eigenVectors(Rij[i]);
		ev_sole[i][0] = ev[i].xx();
		ev_sole[i][1] = ev[i].xy();
		ev_sole[i][2] = ev[i].xz();
		ev_sole[i][3] = ev[i].yx();
		ev_sole[i][4] = ev[i].yy();
		ev_sole[i][5] = ev[i].yz();
		ev_sole[i][6] = ev[i].zx();
		ev_sole[i][7] = ev[i].zy();
		ev_sole[i][8] = ev[i].zz();	
		Info<< " eigenvectors= " << ev[i] <<endl;
	}

		for (int i = 0; i < 5; i++)
		{
		
			for (double value : ev_sole[i])
			Info<<"i= "<<i<<" Individual eigenvectors= "<<value<<endl;
			
		}

}

void UQ::PerturbedAij()
{
	/* declare metrics within member functions means 
	1. you must use it (not only assign values)
	2. the lifetime is within this scope
	Here array variables are defined instead of normal variables, 
	I will write a second version using normal variables compatible
	with the way of looping in OF */
		double c1c[5];
		double c2c[5];
		double c3c[5];

	/*for loop below must read in the baseline eigenvalues and 
	calculate the corresponding c1c c2c and c3c */
	for (int i = 0; i < 5; i++) 
	{
		c1c[i] = e_sort[i][0] - e_sort[i][1];	
		c2c[i] = 2.0*(e_sort[i][1] - e_sort[i][2]);
		c3c[i] = 3.0*e_sort[i][2] + 1.0;	
	}

	for (double value : c1c)//i stands for number of rows
		Info<<"c1c= "<< value << endl;
	for (double value : c2c)//i stands for number of rows
		Info<<"c2c= "<< value << endl;
	for (double value : c3c)//i stands for number of rows
		Info<<"c3c= "<< value << endl;
	
/* define barycentric triangle corner points */
/* Note The corners is a 2D array containing corner points */
/* later maybe try to put it into 1D array */
  double Corners[3][2];
  Corners[0][0] = 1.0; //x1
  Corners[0][1] = 0.0; //y1
  Corners[1][0] = 0.0; //x2
  Corners[1][1] = 0.0; //y2
  Corners[2][0] = 0.5; //x3
  Corners[2][1] = 0.866025; //y3

  /* define barycentric cordinates for baseline prediction */
  double Barycentric_Coord[5][2];

  for (int i = 0; i < 5; i++ )
  {
	Barycentric_Coord[i][0] = Corners[0][0] * c1c[i] + Corners[1][0] * c2c[i] + Corners[2][0] * c3c[i];
	Barycentric_Coord[i][1] = Corners[0][1] * c1c[i] + Corners[1][1] * c2c[i] + Corners[2][1] * c3c[i];
	  
	  for (double value : Barycentric_Coord[i])//i stands for number of rows
		Info<<"Baseline BaryTri coordinates["<<i<<"]= "<< value << endl;
  }

 /* perturbation at 1c, 2c, or 3c state */
	//1c state
	double New_Coord[2];//Only one of these three cases is required therefore [2] for x and y
						//later may change it to if statements see SU2

	//New_Coord[0] = Corners[0][0];
	//New_Coord[1] = Corners[0][1];
	
	//2c state
  	New_Coord[0] = Corners[1][0];
	New_Coord[1] = Corners[1][1];
	
	//3c state
	//New_Coord[0] = Corners[2][0];
	//New_Coord[1] = Corners[2][1];

	/* Calculate perturbed barycentric coordinates:
	note here the baseline Barycentric_coord is replaced with 
	perturbed barycentric_coord */
	
	double uq_delta_b = 0.5; //deltaB

  for (int i = 0; i < 5; i++)
  {
	Barycentric_Coord[i][0] = Barycentric_Coord[i][0] + (uq_delta_b)*(New_Coord[0] - Barycentric_Coord[i][0]);

	Barycentric_Coord[i][1] = Barycentric_Coord[i][1] + (uq_delta_b)*(New_Coord[1] - Barycentric_Coord[i][1]);
	  
	  for (double value : Barycentric_Coord[i])//i stands for number of rows
		Info<<"new BaryTri coordinates["<<i<<"]= "<< value << endl;
  }

  /* rebuild c1c, c2c, and c3c based on new barycentric coordinates */
  for (int i = 0; i < 5; i++)
  {
	
	c3c[i] = Barycentric_Coord[i][1] / Corners[2][1];
	c1c[i] = Barycentric_Coord[i][0] - Corners[2][0]*c3c[i];
	c2c[i] = 1 - c1c[i] - c3c[i];

	e_sort[i][0] = c1c[i] + e_sort[i][1];
	e_sort[i][1] = 0.5*c2c[i] + e_sort[i][2];
	e_sort[i][2] = (c3c[i] - 1.0)/3.0;   
  	 
	for (double value : e_sort[i])//i stands for number of rows
		Info<<"perturbed eigenvalue["<<i<<"]= "<< value << endl;  
  }

  {

		for (double value : c1c)//i stands for number of rows
			Info<<"perturbed c1c= "<< value << endl;
		for (double value : c2c)//i stands for number of rows
			Info<<"perturbed c2c= "<< value << endl;
		for (double value : c3c)//i stands for number of rows
			Info<<"perturbed c3c= "<< value << endl; 

  }

  InitLog();

}

}//End namespace foam

