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

#ifndef _sproblem_h_
#define _sproblem_h_

#include "problem.h"
#include "realmatrix.h"
#include "hermitepolynomial.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class SProblem : public Problem {
	 public:
        SProblem(void);
		~SProblem(void) { };
        
		// copy constructor
		SProblem(const SProblem &sol);
		
        // run
        void Run(const std::string &fileName);
		
		// testing
		void Test(const std::string &fileName);
				
	private:
		// input
		void ReadInputFile(const std::string &fileName);
		
		// modes
		void GetMode(Array<Array<double> > &u, long modeIndex) const;
		
		// density
		void CheckDensities(void);
		
		// Initialization
		void Initialize(void);
		
		// initial data
		void SetInitialDataAndDensity(void);
		void SetInitialDataAndDensity(ModeType modeType, short numPts);
		double Spacing(ModeType modeType) const;
		
		// volterra equation
		void ComputeVolterraCoefficients(void);
		void ComputeLiouvillian(const Array<Array<double> > &f, Array<Array<double> > &liouvillian) const;
		void ComputeProjection(const Array<Array<double> > &f, Array<double> &projection) const;
		
		// averaging
		double ComputeAverage(const Array<Array<double> > &f) const;
		double ComputeAverage(const Array<double> &f, ModeType modeType) const;

		// member data
	private:
		// selected resolved and unresolved modes
		long mResolvedMode;
		long mUnresolvedMode;
		
		// initial data
		Array<double> mResolvedID, mUnresolvedID;
		Array<double> mDensityResolved, mDensityUnresolved;
		
		// volterra equation
		long mNumVolterraTimeSteps;
		Array<Array<double> > mVolterraF;
		
		// finite rank expansion
		HermitePolynomial mHermitePolynomial;
		short mFiniteRankSize;
		Array<double> mFiniteRankCoefficient;
	};



	inline SProblem::SProblem()
	{
		mResolvedMode = 0;
		mUnresolvedMode = 0;
		
		mNumVolterraTimeSteps = 0;
		mFiniteRankSize = 0;
		
		return;
	} 
	
	
	
	inline double SProblem::Spacing(ModeType modeType) const
	{
		switch (modeType) {
		case RESOLVED_MODE:
			return mResolvedID[1] - mResolvedID[0];
		
		case UNRESOLVED_MODE:
			return mUnresolvedID[1] - mUnresolvedID[0];
		
		default:
			ThrowException("SProblem::Spacing : bad mode type");
			return -1.0;
		}
	}
}

#endif // _sproblem_h_	
