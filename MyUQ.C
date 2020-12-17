#include <iostream>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include "tensor.H"
#include "triad.H"
#include "MyUQ.H"


//using namespace Foam;//using directive namespace has the scope that encloses entire program
/* Note for variables prefixed with "_" means the variable is passed (always by reference) into
   the UQ translation unit from OpenFOAM environment (already calculated in OpenFOAM and always 
   need to be converted into the compatible type with that used in the UQ methodology, e.g. restore _Rij
   into m_MeanReynoldsStress[x][y] that will then be used to calculate perturbed Reynolds stresses) 
   
   prefix meaning:

   "-"    form OpenFoam
   "tmp"  temporary
   "m"    member variables */

namespace Foam{//specifies the Foam namespace for restricted scopes

UQ::UQ(scalar& _k, symmTensor& _Rij, unsigned short _celli)
{
	unsigned short i,j;

	Info<<"C*********************** UQ methodology start here ***************************C\n";
	Info<<"C		**************************************                        C\n";                           
	Info<<"C			***************************                           C\n";
	Info<<"C				*********                                     C\n";
	Info<<"C**************** "<<"celli= "<<_celli<<" Setting eigenvectors in square matrix *************C\n\n";

	Info<<"turbulence kinetic energy= "<< _k <<endl;
	for (i = 0; i < 6; i++)
	{
		Info<<"Rij["<<i<<"]= "<<_Rij[i]<<endl;
	}
	
	for (i = 0; i < 3; i++)
	{  
		m_Eig_Val[i] = 0.0;
		for (j = 0; j < 3; j++)
		{
			m_MeanReynoldsStress[i][j] 	= 0.0;
			m_MeanPerturbedRSM[i][j]	= 0.0;
			m_Eig_Vec[i][j] 		   	= 0.0;
			m_New_Eig_Vec[i][j] 		= 0.0; 
			m_newB_ij[i][j]				= 0.0;
		}
	}

 	
	m_turb_ke = _k;//turb_ke is not declared within function, e.g. UQ::QU(scalar& m_turb_ke, ..), therefore must be declared in MyUQ.H file
	m_MeanReynoldsStress[0][0] = _Rij[0];		
 	m_MeanReynoldsStress[0][1] = _Rij[1];
 	m_MeanReynoldsStress[0][2] = _Rij[2];
 	m_MeanReynoldsStress[1][0] = _Rij[1];
 	m_MeanReynoldsStress[1][1] = _Rij[3];
 	m_MeanReynoldsStress[1][2] = _Rij[4];
 	m_MeanReynoldsStress[2][0] = _Rij[2];
 	m_MeanReynoldsStress[2][1] = _Rij[4];
 	m_MeanReynoldsStress[2][2] = _Rij[5];
	
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			Info<<"m_MeanReynoldsStress["<<i<<"]["<<j<<"]"<<m_MeanReynoldsStress[i][j]<<endl;
	


/* The following block of code will disappear in the final implementation
   and will be replaced with calculated anisotropy tensor */
/*************************** Start ******************************************/   
//	Aij[0] = {1,0,-4,0,5,4,-4,4,3};//to read in fields


//	//read in Reynold stress field
//	Info<<"*********** Start reading anisotropy tensor in each cell **************\n\n";
//	for ( int i = 1; i < 5; i++)
//	{
//		Aij[i] = {0.6592491, 0.04143975, 0, 0.04143974, -0.3259158, 0, 0, 0, -0.333333};
//		//Aij[i] = {1, 1, -1, 1, 2, 0, -1, 0, 5};
//		//Aij[i] = {3, -2, 4, -2, 6, 2, 4, 2, 3};
//		Info<<"i= "<<i-1<<" Aij= "<<Aij[i-1]<<Aij[i-1].row(1)<<nl;
//		if(i==4)
//		{
//			Info<<"i= "<<i<<" Aij= "<<Aij[i]<<=Aij[i].row(1)<<"\n\n";

//		}
//	}
//	Info<<"*********** End ********************************************************\n\n"<<endl;
/*************************** End ******************************************/

}

//**********Member functions************//
void UQ::EigenSpace (symmTensor& m_newBij_OF, symmTensor& m_newuiuj_OF, double** m_newB_ij, double** m_MeanPerturbedRSM, //pass to main()
					 symmTensor& _Bij, unsigned short& _celli)//passed from main()
{
	tensor tmp_ev_inv;//An temporary variable declared here to help store the inverse of 
					 //eigenvectors, because eigenvectors (3 components) corresponding
					 //to each eigenvalue (3 in total) need to be arranged in column manner
					 //so that "EigenRecomposition" can work properly
	unsigned short x, y;

	
		
		m_e = eigenValues(_Bij);//Note here "eigenValues" is the intrinsic function
								   //offered by OpenFoam 
		m_ev = eigenVectors(_Bij);
		tmp_ev_inv = inv(m_ev); //find the inverse of eigenvectors will set each set of eigenvectors
								//(for each eigenvalue) in column manner which is directly used to
								//calculate recomposed anisotropy tensor 

		m_Eig_Val[0] = m_e.x();
		m_Eig_Val[1] = m_e.y();
		m_Eig_Val[2] = m_e.z();

		m_ev_sole[0] = tmp_ev_inv.xx();
		m_ev_sole[1] = tmp_ev_inv.xy();
		m_ev_sole[2] = tmp_ev_inv.xz();
		m_ev_sole[3] = tmp_ev_inv.yx();
		m_ev_sole[4] = tmp_ev_inv.yy();
		m_ev_sole[5] = tmp_ev_inv.yz();
		m_ev_sole[6] = tmp_ev_inv.zx();
		m_ev_sole[7] = tmp_ev_inv.zy();
		m_ev_sole[8] = tmp_ev_inv.zz();	

	
		for (x = 0; x < 3; x++)
		{
			int count = x*3;
			for (y = 0; y < 3; y++)
			{
				m_Eig_Vec[x][y] 	= m_ev_sole[y+count];
				m_New_Eig_Vec[x][y] = m_ev_sole[y+count];
				Info << "m_Eig_Vec["<<x<<"]["<<y<<"]= "<<m_Eig_Vec[x][y] <<endl;   
			}
		}

		std::sort(m_Eig_Val, m_Eig_Val+3);//keep the original eigenvalue matrix as default (increasing order) to 
										  //be compatible with the calculated eigenvalues, i.e. eigenValues(en_inv[i])

		//std::sort(m_Eig_Val, m_Eig_Val+3, std::greater<double>());
		Info<< "celli= "<<_celli<<" m_Eig_Val[0]= " << m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;	
		Info<<"******************************** End **********************************************\n\n"<<endl;
	
		/* The following block of code verifies the methodology that was implemented to recompose the original matrix from 
		   its eigenvectors and eigenvalues (from CS point of view since variables m_B_ij, m_Eig_Vec, m_Eig_Val, etc are stored on heap
		   they have the lifetime till the program comes to the end of the global scope, keep that in mind and ENSURE your perturbed aij 
		   matrices are passed to the functions that need these perturbed matrices, e.g. pk and pw ) */

		/*Info<<"****************"<<" i_before= "<<i<< " Start calculating Recomposed anisotropy tensor ********************\n\n";
		EigenRecomposition(m_B_ij, m_Eig_Vec, m_Eig_Val, 3);

		for (x = 0; x < 3; x++)
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"m_B_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_B_ij[x][y]<<endl; 
			}
		}
		Info<<"******************************** End **********************************************\n\n"<<endl;*/



		Info<<"**************** "<<"celli= "<<_celli<< " Start calculating Perturbed anisotropy tensor ********************\n";
		Info<<"**************** ********************************************************************* *****\n\n";

		PerturbedAij(m_newB_ij, m_MeanPerturbedRSM);//perturbed anisotropy matrix at an individual cell is returned 
		m_newBij_OF = {m_newB_ij[0][0], m_newB_ij[0][1], m_newB_ij[0][2], 
				  m_newB_ij[1][1], m_newB_ij[1][2], m_newB_ij[2][2]};//Aij_OF is a "symmTensor" type variable which is a 1D array with size 
				  													 //1 in this particular case. It stores the calculated perturbed anisotropy
																	 //matrix at an individual cell and return to where the "EigenSpace" function
																	 //is being called (main function) through passing by reference
																	 
		m_newuiuj_OF = {m_MeanPerturbedRSM[0][0], m_MeanPerturbedRSM[0][1], m_MeanPerturbedRSM[0][2],
						m_MeanPerturbedRSM[1][1], m_MeanPerturbedRSM[1][2], m_MeanPerturbedRSM[2][2]};
		for (x = 0; x < 3; x++)//this for loop is created to test if the perturbed anisotropy matrix is returned correctly
		{
			for (y = 0; y < 3; y++)
			{
				Info<<"m_newB_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_newB_ij[x][y]<<endl; 
			}
		}
		Info<<"******************************** End **********************************************\n";
		Info<<"******************************** End **********************************************\n\n"<<endl;
}

void UQ::PerturbedAij(double** m_newB_ij, double** m_MeanPerturbedRSM)
{
	/* declare metrics within member functions means 
	1. you must use it (not only assign values)
	2. the lifetime is within this scope
	Here array variables are defined instead of normal variables, 
	I will write a second version using normal variables compatible
	with the way of looping in OF */
		double tmp_c1c;
		double tmp_c2c;
		double tmp_c3c;

		unsigned short x, y;
	//std::sort(m_Eig_Val, m_Eig_Val+3, std::greater<double>()); //This reorders eigenvalues in 
																 //a non-increasing order. However
																 //not useful just for reference.
											

	/*for loop below must read in the baseline eigenvalues and 
	calculate the corresponding c1c c2c and c3c */

		tmp_c1c = m_Eig_Val[2] - m_Eig_Val[1];	
		tmp_c2c = 2.0*(m_Eig_Val[1] - m_Eig_Val[0]);
		tmp_c3c = 3.0*m_Eig_Val[0] + 1.0;	
	
	Info<< "********** Start unperturbed (baseline) eigenvalues **********\n\n";  
	Info<< "c1c= "<<tmp_c1c<<" c2c= "<<tmp_c2c<<" c3c= "<<tmp_c3c<<endl;
	Info<< "m_Eig_Val[0]= "<<m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;
	Info<< "********** End unperturbed (baseline) eigenvalues **********\n\n"<<endl; 

	
/* define barycentric triangle corner points */
/* Note The corners is a 2D array containing corner points */
/* later maybe try to put it into 1D array */
 
  m_Corners[0][0] = 1.0; //x1
  m_Corners[0][1] = 0.0; //y1
  m_Corners[1][0] = 0.0; //x2
  m_Corners[1][1] = 0.0; //y2
  m_Corners[2][0] = 0.5; //x3
  m_Corners[2][1] = 0.866025; //y3

  /* define barycentric coordinates for baseline prediction */
	m_Barycentric_Coord[0] = m_Corners[0][0] * tmp_c1c + m_Corners[1][0] * tmp_c2c + m_Corners[2][0] * tmp_c3c;
	m_Barycentric_Coord[1] = m_Corners[0][1] * tmp_c1c + m_Corners[1][1] * tmp_c2c + m_Corners[2][1] * tmp_c3c;
	Info<<"Bary_Coord[0]= "<<m_Barycentric_Coord[0]<<" Bary_Coord[1]= "<<m_Barycentric_Coord[1]<<endl;  

 /* perturbation at 1c, 2c, or 3c state: Only one of these three cases is required 
 therefore [2] for x and y later may change it to if statements see SU2 */
	
	//1c state
	//m_New_Coord[0] = m_Corners[0][0];
	//m_New_Coord[1] = m_Corners[0][1];
	
	//2c state
  	m_New_Coord[0] = m_Corners[1][0];
	m_New_Coord[1] = m_Corners[1][1];
	Info<<"New_Coord[0]= "<<m_New_Coord[0]<<" New_Coord[1]= "<<m_New_Coord[1]<<endl;  
	
	//3c state
	//m_New_Coord[0] = m_Corners[2][0];
	//m_New_Coord[1] = m_Corners[2][1];

	/* Calculate perturbed barycentric coordinates:
	note here the baseline Barycentric_coord is replaced with 
	perturbed barycentric_coord */
	
	double uq_delta_b = 0.5; //deltaB is manually defined in the code
							 //later try to Get in the "controdict"

	m_Barycentric_Coord[0] = m_Barycentric_Coord[0] + (uq_delta_b)*(m_New_Coord[0] - m_Barycentric_Coord[0]);

	m_Barycentric_Coord[1] = m_Barycentric_Coord[1] + (uq_delta_b)*(m_New_Coord[1] - m_Barycentric_Coord[1]);
	  
	Info<<"Bary_Coord_perturb[0]= "<<m_Barycentric_Coord[0]<<" Bary_Coord_perturb[1]= "<<m_Barycentric_Coord[1]<<endl;  

  

  /* rebuild c1c, c2c, and c3c based on new (perturbed) barycentric coordinates */	
	tmp_c3c = m_Barycentric_Coord[1] / m_Corners[2][1];
	tmp_c1c = m_Barycentric_Coord[0] - m_Corners[2][0]*tmp_c3c;
	tmp_c2c = 1 - tmp_c1c - tmp_c3c;  

	m_Eig_Val[0] = (tmp_c3c - 1.0)/3.0;   
	m_Eig_Val[1] = 0.5*tmp_c2c + m_Eig_Val[0];
	m_Eig_Val[2] = tmp_c1c + m_Eig_Val[1];
  	 
	Info<< "********** Start perturbed eigenvalues **********\n\n";  
	Info<< "c1c= "<<tmp_c1c<<" c2c= "<<tmp_c2c<<" c3c= "<<tmp_c3c<<endl;
	Info<< "m_Eig_Val[0]= "<<m_Eig_Val[0]<<" m_Eig_Val[1]= "<<m_Eig_Val[1]<<" m_Eig_Val[2]= "<<m_Eig_Val[2]<<endl;
	Info<< "********** End perturbed eigenvalues **********\n\n"<<endl; 
	


	EigenRecomposition(m_newB_ij, m_New_Eig_Vec, m_Eig_Val, 3);
	
	for (x = 0; x < 3; x++)
	{
		for (y = 0; y < 3; y++)
		{
			//Info<<"m_B_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_B_ij[x][y]<<endl; 
			//Info<<"m_NewEig_Vec"<<"["<< x <<"]"<<"["<<y<<"]= "<<m_New_Eig_Vec[x][y]<<endl;
			Info<<"m_newB_ij"<<"["<< x <<"]"<<"["<<y<<"]= "<< m_newB_ij[x][y]<<endl; 
		}
	}

	double uq_urlx = 0.1;

	for (x = 0; x < 3; x++)
	{
		for (y = 0; y < 3; y++)
		{
			m_MeanPerturbedRSM[x][y] = 2.0 * m_turb_ke * (m_newB_ij[x][y] + 1.0/3.0 * m_delta3[x][y]);
			m_MeanPerturbedRSM[x][y] = m_MeanReynoldsStress[x][y] + uq_urlx*(m_MeanPerturbedRSM[x][y] - m_MeanReynoldsStress[x][y]);
			
			Info<<"m_MeanPerturbedRSM["<<x<<"]["<<y<<"]= "<<m_MeanPerturbedRSM[x][y]<<
			" m_MeanReynoldsStress["<<x<<"]["<<y<<"]= "<<m_MeanReynoldsStress[x][y]<<endl;
			
		}
	}
	
}



}//End namespace foam

