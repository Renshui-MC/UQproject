#pragma once

#include <iostream>
#include "Numerics.H"

namespace Foam
{
class UQ : public Numerics
{
	
	protected:
		scalar m_turb_ke;
		scalar _k;
		/***** Even though the following two variables are declared and defined in main() and will compile fine without 
		 errors and calculations inside the MyUQ.C translation unit can proceed thoroughly through passing by reference 
		 under certain conditions, the BIG problem is once being passed into the MyUQ translation unit there is no room 
		 (memory) to hold them which must be assigned to any defined variables. In short, always allocate variables to
		 memory so that you can use them when you want */
		symmTensor _Rij; 
		symmTensor _Bij; 
		/****************** End ********************/

		symmTensor m_newBij_OF;
		symmTensor m_newuiuj_OF;
		vector m_e; 
		tensor m_ev;
		double m_e_sort[3];//the first [] is rwo and the second [] is col. **later might change to double** e_sort = new double* [];**
		double m_ev_sole[9];


	public:
//*********** Constructor ***************//
		UQ (scalar& _k, symmTensor& _Rij, unsigned short _celli);

//********* Member functions ************//
		virtual void EigenSpace (symmTensor& m_newBij_OF, symmTensor& m_newuiuj_OF, double** m_newB_ij, double** m_MeanPerturbedRSM,
								 symmTensor& _Bij, unsigned short& _celli);


		virtual void PerturbedAij (double** m_newB_ij, double** m_MeanPerturbedRSM);
	
};

}//End namespace foam
