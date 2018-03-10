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

#ifndef _averagingproblem_h_
#define _averagingproblem_h_

#include "problem.h"
#include "realmatrix.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class AveragingProblem : public Problem {
	 public:
        AveragingProblem(void);
		~AveragingProblem(void) { };
        
		// copy constructor
		AveragingProblem(const AveragingProblem &sol);
		
        // run
        void Run(const std::string &fileName);
		
		// testing
		void Test(const std::string &fileName);
				
	private:
		// input
		void ReadInputFile(const std::string &fileName);
		
		// run
		void RunOneUnresolved(void);
		void RunFixedSpacingOneUnresolved(void);

		// Initialization
		void Initialize(long numSystems = 1);
		void SetInitialDataOneUnresolved(void);
		void CheckDensity(void) const;
		
		// quadrature
		void SetQuadratureType(std::string quadratureType);
		
		// averaging
		void AverageOneUnresolved(void);
		
		// member data
	private:
		// quadrature 
		short mInitialQuadratureLevel;
		double mQuadratureTolerance;
		QuadratureType mQuadratureType;
		short mQuadratureNumGridPoints;
		
		// averages
		Matrix<double> mAverage;
	};



	inline AveragingProblem::AveragingProblem()
	{
		mInitialQuadratureLevel = 0;
		mQuadratureTolerance = -1.0;
		mQuadratureType = NO_QUADRATURE_TYPE;
		mQuadratureNumGridPoints = 0;
		
		return;
	} 
}

#endif // _averagingproblem_h_	
