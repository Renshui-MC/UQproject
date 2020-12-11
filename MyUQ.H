#pragma once

#include <iostream>
#include "Numerics.H"

namespace Foam
{
class UQ : public Numerics
{
	
	protected:
		tensor Rij[5]; //Note 5 means 5 rows for this particular example
		vector e[5]; 
		tensor ev[5];
		double e_sort[5][3];//the first [] is rwo and the second [] is col. **later might change to double** e_sort = new double* [];**
		double ev_sole[5][9];


	public:
//*********** Constructor ***************//
		UQ ();

//********* Member functions ************//
		void EigenSpace ();
		

		void PerturbedAij ();
	
};

}//End namespace foam
