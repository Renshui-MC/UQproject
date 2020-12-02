#pragma once

#include <iostream>


using namespace Foam;


class UQ
{
	protected:
		tensor Rij[5]; //Note 5 means 5 rows for this particular example
		vector e[5]; 
		tensor ev[5];
		double e_sort[5][3];//the first [] is rwo and the second [] is col. **later might change to double** e_sort = new double* [];**
	

	
	public:
//*********** Constructor ***************//
		UQ ();

//********* Member functions ************//
		void EigenvalueSort ();
		
		void Eigenvector ();

		void PerturbedAij ();
		void InitLog();
		
};