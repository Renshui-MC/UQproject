#pragma once

#include <iostream>
#include "Numerics.H"

namespace Foam
{
class UQ : public Numerics
{
	
	protected:
		symmTensor Aij; //Note 5 means 5 rows for this particular example
		symmTensor Aij_OF;
		vector e; 
		tensor ev;
		double e_sort[3];//the first [] is rwo and the second [] is col. **later might change to double** e_sort = new double* [];**
		double ev_sole[9];


	public:
//*********** Constructor ***************//
		UQ ();

//********* Member functions ************//
		virtual void EigenSpace (symmTensor& Aij_OF, double** m_newA_ij, symmTensor Aij, unsigned short& celli);
		

		virtual void PerturbedAij (double** m_newA_ij);
	
};

}//End namespace foam
