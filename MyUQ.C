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
		

	unsigned short i,j;


	for (i = 0; i < 3; i++)
	{  
		m_Eig_Val[i] = 0.0;
		for (j = 0; j < 3; j++)
		{
			m_Eig_Vec[i][j] 	= 0.0;
			m_New_Eig_Vec[i][j] = 0.0; 
		}
	}
		

/* The following block of code will disappear in the final implementation
   and will be replaced with calculated anisotropy tensor */
/*************************** Start ******************************************/   
	Rij[0] = {1,0,-4,0,5,4,-4,4,3};//to read in fields

	//read in Reynold stress field
	Info<<"*********** Start reading anisotropy tensor in each cell **************\n\n";
	for ( int i = 1; i < 5; i++)
	{
		Rij[i] = {0.6592491, 0.04143975, 0, 0.04143974, -0.3259158, 0, 0, 0, -0.333333};
		//Rij[i] = {1, 1, -1, 1, 2, 0, -1, 0, 5};
		//Rij[i] = {3, -2, 4, -2, 6, 2, 4, 2, 3};
		Info<<"i= "<<i-1<<" Rij= "<<Rij[i-1]<<Rij[i-1].row(1)<<nl;
		if(i==4)
		{
			Info<<"i= "<<i<<" Rij= "<<Rij[i]<<Rij[i].row(1)<<"\n\n";

		}
	}
	Info<<"*********** End ********************************************************\n\n"<<endl;
/*************************** End ******************************************/

}

//**********Member functions************//
void UQ::EigenSpace ()
{
	tensor ev_inv[5];//An temporary variable declared here to help store the inverse of 
					 //eigenvectors, because eigenvectors (3 components) corresponding
					 //to each eigenvalue (3 in total) need to be arranged in column manner
					 //so that "EigenRecomposition" can work properly
	unsigned short i, x, y;

	for (i = 0; i < 5; i++)
	{
		
		e[i] = eigenValues(Rij[i]);//Note here "eigenValues" is the intrinsic function
								   //offered by OpenFoam 
		ev[i] = eigenVectors(Rij[i]);
		ev_inv[i] = inv(ev[i]); //find the inverse of eigenvectors will set each set of eigenvectors
								//(for each eigenvalue) in column manner which is directly used to
								//calculate recomposed anisotropy tensor 

		m_Eig_Val[0] = e[i].x();
		m_Eig_Val[1] = e[i].y();
		m_Eig_Val[2] = e[i].z();

		ev_sole[i][0] = ev_inv[i].xx();
		ev_sole[i][1] = ev_inv[i].xy();
		ev_sole[i][2] = ev_inv[i].xz();
		ev_sole[i][3] = ev_inv[i].yx();
		ev_sole[i][4] = ev_inv[i].yy();
		ev_sole[i][5] = ev_inv[i].yz();
		ev_sole[i][6] = ev_inv[i].zx();
		ev_sole[i][7] = ev_inv[i].zy();
		ev_sole[i][8] = ev_inv[i].zz();	

		Info<<"****************"<<" i= "<<i<< " Setting eigenvectors in square matrix ********************\n\n";
		for (x = 0; x < 3; x++)
		{
			int count = x*3;
			for (y = 0; y < 3; y++)
			{
				m_Eig_Vec[x][y] 	= ev_sole[i][y+count];
				m_New_Eig_Vec[x][y] = ev_sole[i][y+count];
				Info << "m_Eig_Vec["<<x<<"]["<<y<<"]= "<<m_Eig_Vec[x][y] <<endl;   
			}
		}

		std::sort(m_Eig_Val, m_Eig_Val+3);//keep the original eigenvalue matrix as default (increasing order) to 
										  //be compatible with the calculated eigenvalues, i.e. eigenValues(en_inv[i])

		//std::sort(m_Eig_Val, m_Eig_Val+3, std::greater<double>());
		Info<< "i= "<<i<<" m_Eig_Val[0]= " << m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;	
		Info<<"******************************** End **********************************************\n\n"<<endl;
	
		/* The following block of code verifies the methodology that was implemented to recompose the original matrix from 
		   its eigenvectors and eigenvalues (from CS point of view since variables m_A_ij, m_Eig_Vec, m_Eig_Val, etc are stored on heap
		   they have the lifetime till the program comes to the end of the global scope, keep that in mind and ENSURE your perturbed aij 
		   matrices are passed to the functions that need these perturbed matrices, e.g. pk and pw ) */

		/*Info<<"****************"<<" i_before= "<<i<< " Start calculating Recomposed anisotropy tensor ********************\n\n";
		EigenRecomposition(m_A_ij, m_Eig_Vec, m_Eig_Val, 3);

		for (x = 0; x < 3; x++)
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"m_A_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_A_ij[x][y]<<endl; 
			}
		}
		Info<<"******************************** End **********************************************\n\n"<<endl;*/



		Info<<"****************"<<" i_before= "<<i<< " Start calculating Perturbed anisotropy tensor ********************\n\n";

		PerturbedAij();

		Info<<"******************************** End **********************************************\n\n"<<endl;

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
		double c1c;
		double c2c;
		double c3c;

		unsigned short x, y;
	//std::sort(m_Eig_Val, m_Eig_Val+3, std::greater<double>()); //This reorders eigenvalues in 
																 //a non-increasing order. However
																 //not useful just for reference.
											

	/*for loop below must read in the baseline eigenvalues and 
	calculate the corresponding c1c c2c and c3c */

		c1c = m_Eig_Val[2] - m_Eig_Val[1];	
		c2c = 2.0*(m_Eig_Val[1] - m_Eig_Val[0]);
		c3c = 3.0*m_Eig_Val[0] + 1.0;	
	
	Info<< "**********Start unperturbed (baseline) eigenvalues**********\n\n";  
	Info<< "c1c= "<<c1c<<" c2c= "<<c2c<<" c3c= "<<c3c<<endl;
	Info<< "m_Eig_Val[0]= "<<m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;
	Info<< "**********End unperturbed (baseline) eigenvalues**********\n\n"<<endl; 

	
/* define barycentric triangle corner points */
/* Note The corners is a 2D array containing corner points */
/* later maybe try to put it into 1D array */
 
  Corners[0][0] = 1.0; //x1
  Corners[0][1] = 0.0; //y1
  Corners[1][0] = 0.0; //x2
  Corners[1][1] = 0.0; //y2
  Corners[2][0] = 0.5; //x3
  Corners[2][1] = 0.866025; //y3

  /* define barycentric cordinates for baseline prediction */
	Barycentric_Coord[0] = Corners[0][0] * c1c + Corners[1][0] * c2c + Corners[2][0] * c3c;
	Barycentric_Coord[1] = Corners[0][1] * c1c + Corners[1][1] * c2c + Corners[2][1] * c3c;
	Info<<"Bary_Coord[0]= "<<Barycentric_Coord[0]<<" Bary_Coord[1]= "<<Barycentric_Coord[1]<<endl;  

 /* perturbation at 1c, 2c, or 3c state: Only one of these three cases is required 
 therefore [2] for x and y later may change it to if statements see SU2 */
	
	//1c state
	//New_Coord[0] = Corners[0][0];
	//New_Coord[1] = Corners[0][1];
	
	//2c state
  	New_Coord[0] = Corners[1][0];
	New_Coord[1] = Corners[1][1];
	Info<<"New_Coord[0]= "<<New_Coord[0]<<" New_Coord[1]= "<<New_Coord[1]<<endl;  
	
	//3c state
	//New_Coord[0] = Corners[2][0];
	//New_Coord[1] = Corners[2][1];

	/* Calculate perturbed barycentric coordinates:
	note here the baseline Barycentric_coord is replaced with 
	perturbed barycentric_coord */
	
	double uq_delta_b = 0.5; //deltaB is manually defined in the code
							 //later try to Get in the "controdict"

	Barycentric_Coord[0] = Barycentric_Coord[0] + (uq_delta_b)*(New_Coord[0] - Barycentric_Coord[0]);

	Barycentric_Coord[1] = Barycentric_Coord[1] + (uq_delta_b)*(New_Coord[1] - Barycentric_Coord[1]);
	  
	Info<<"Bary_Coord_perturb[0]= "<<Barycentric_Coord[0]<<" Bary_Coord_perturb[1]= "<<Barycentric_Coord[1]<<endl;  

  

  /* rebuild c1c, c2c, and c3c based on new (perturbed) barycentric coordinates */	
	c3c = Barycentric_Coord[1] / Corners[2][1];
	c1c = Barycentric_Coord[0] - Corners[2][0]*c3c;
	c2c = 1 - c1c - c3c;  

	m_Eig_Val[0] = (c3c - 1.0)/3.0;   
	m_Eig_Val[1] = 0.5*c2c + m_Eig_Val[0];
	m_Eig_Val[2] = c1c + m_Eig_Val[1];
  	 
	Info<< "**********Start perturbed eigenvalues**********\n\n";  
	Info<< "c1c= "<<c1c<<" c2c= "<<c2c<<" c3c= "<<c3c<<endl;
	Info<< "m_Eig_Val[0]= "<<m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;
	Info<< "**********End perturbed eigenvalues**********\n\n"<<endl; 
	


	EigenRecomposition(m_newA_ij, m_New_Eig_Vec, m_Eig_Val, 3);
	
	for (x = 0; x < 3; x++)
	{
		for (y = 0; y < 3; y++)
		{
			//Info<<"m_A_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_A_ij[x][y]<<endl; 
			//Info<<"m_NewEig_Vec"<<"["<< x <<"]"<<"["<<y<<"]= "<<m_New_Eig_Vec[x][y]<<endl;
			Info<<"m_newA_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_newA_ij[x][y]<<endl; 
		}
	}

}

}//End namespace foam

