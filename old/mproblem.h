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

#ifndef _mproblem_h_
#define _mproblem_h_

#include "problem.h"
#include "realmatrix.h"
#include "hermitepolynomial.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class MProblem : public Problem {
	 public:
        MProblem(void);
		~MProblem(void) { };
        
		// copy constructor
		MProblem(const MProblem &sol);
		
        // run
        void Run(const std::string &fileName);
						
	private:
		// input
		void ReadInputFile(const std::string &fileName);
				
		// Initialization
		void Initialize(void);
		
		// initial data
		
		// volterra equation
		void ComputeVolterraCoefficients(long timeStep);
		void ComputeVolterraF(long timeStep);
		void SetVolterraTimeSchedule(long numSteps);
		
		// averaging
		double ComputeAverage(const Array<Array<double> > &f) const;
		double ComputeAverage(const Array<double> &f, ModeType modeType) const;
		
		// IO 
		void WriteVolterraFToFile(std::string fileName) const;

		// member data
	private:
		// selected resolved and unresolved modes
		long mResolvedMode;
		
		// initial data
		Array<double> mResolvedID;
		
		// volterra equation
		Array<double> mVolterraT;
		Array<double> mVolterraF;
		double mSecondModeSpacing;
		
		// finite rank expansion
		HermitePolynomial mHermitePolynomial;
		short mFiniteRankSize;
		Array<double> mFiniteRankCoefficient;
	};



	inline MProblem::MProblem()
	{
		mResolvedMode = 0;
		mSecondModeSpacing = -1.0;
		mFiniteRankSize = 0;
		
		return;
	} 
}

#endif // _mproblem_h_	
