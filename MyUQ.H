#pragma once

#include <iostream>
#include "Numerics.H"

namespace Foam
{
class UQ : public Numerics
{
	
	protected:
		/***** The following variables are declared and defined in main() and will compile fine without 
		 errors and calculations inside the MyUQ.C translation unit can proceed thoroughly through passing by reference */

		//scalar _k;
		//symmTensor _Rij; 
		//symmTensor _Bij; 
		//symmTensor m_newBij_OF;
		//symmTensor m_newuiuj_OF;
		/****************** End ********************/

		scalar m_turb_ke; //turbulence kinetic energy passed from main(), which could be perturbed
		scalar m_nut; //kinematic eddy viscosity passed from main()
		vector m_e; 
		tensor m_ev;
		
		const unsigned short n = 3;
		


	public:
//*********** Constructor ***************//
		UQ (scalar& _k, scalar& _nut, symmTensor& _Rij, label& _celli);//already declared in main and passed into MyUQ.C by reference
																  //do not declare them again (refer to the above comments)

//********* Member functions ************//
		virtual void EigenSpace (symmTensor& m_newBij_OF, symmTensor& m_newuiuj_OF, symmTensor& m_newSij_OF,
								 double** m_newB_ij, double** m_MeanPerturbedRSM, double** m_PerturbedStrainRate,
								 symmTensor& _Bij, label& _celli);


		virtual void PerturbedBij (double** m_newB_ij, double** m_MeanPerturbedRSM, double** m_PerturbedStrainRate); //return the perturbed normalized anisotropy matrix

		virtual void SetPerturbedStrainMag (double** m_MeanPerturbedRSM, double** m_PerturbedStrainRate);
		
};

}//End namespace foam
