/*
 * Copyright (C) 2004-2018 David Bernstein <david.h.bernstein@gmail.com>
 *
 * This file is part of OPBE.
 *
 * OPBE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OPBE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPBE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _mkproblem_h_
#define _mkproblem_h_

#include "problem.h"
#include "realmatrix.h"
#include "hermitepolynomial.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class MKProblem : public Problem {
	 public:
        MKProblem(void);
		~MKProblem(void) { };
        
		// copy constructor
		MKProblem(const MKProblem &sol);
		
        // run
        void Run(const std::string &fileName);
		
		// test
		void Test(const std::string &fileName);
						
	private:
		// run
		void Run(void);
		
		// input
		void ReadInputFile(const std::string &fileName);
				
		// Initialization
		void Initialize(void);
		
		// initial data
		void ComputeInitialData(void);
		
		// density
		void SetInitialDensities(void);
		
		// volterra equation
		void UpdateVolterraCoefficients(long timeStep);
		void UpdateVolterraFAverage(long timeStep);
		
		// averaging
		double GaussianRandomVariable(double mean, double sigma) const;
		void SetBigS(void);
		
		// IO 
		void WriteVolterraFFile(void);

		// member data
	private:
		// monte carlo
		long mNumMonteCarloRuns;
		long mRunCount;
		
		// volterra equation
		short mFiniteRankSize;
		Matrix<double> mVolterraF0;
		double mBigS;
	};



	inline MKProblem::MKProblem()
	{
		mNumMonteCarloRuns = 0;
		mFiniteRankSize = 1;
		mRunCount = 0;
		
		return;
	} 
}

#endif // _mkproblem_h_	
