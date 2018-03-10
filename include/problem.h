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

#ifndef _problem_h_
#define _problem_h_

#include "utility.h"
#include "opbeconst.h"
#include "system.h"
#include "runcontrol.h"
#include "density.h"
#include "modeindex.h"
#include "opbeparameter.h"

#include <string>
#include <iostream>

#include <gsl/gsl_rng.h>

namespace NAMESPACE {
	class Problem {
	 public:
        Problem(void);
		~Problem(void) { };
        
		// copy constructor
		Problem(const Problem &sol);
				
		// problem type
		ProblemType GetProblemType(const std::string &fileName);
		
	protected:
		// Initialization
		virtual void ReadInputFile(const std::string &fileName);
		void ReadInitialConditions(const std::string &fileName);
		void Reset(void);
		void InitializeRandomNumberGenerator(void);
		
		// evolution
		void Evolve(double t1);
		
		// IO
		void WriteOutput(void);
		void WriteModes(void);
		void WriteEnergy(void);
		void WriteMoments(void);
		void WriteTModelRatio(void);
		void PrintCurrentTime(void) const;
		
		// member data
	protected:
		// run control parameters
		RunControl mRunControl;
		
		// time
		double mCurrentTime;
		
		// state
		ProblemState mState;
		
		// parameters
		OPBEParameter mOPBEParameter;
		
		// modes
		long mNumModes;
		long mNumResolvedModes;
		long mNumUnresolvedModes;
		
		// mode index function
		ModeIndex mModeIndex;
		
		// systems
		Array<System> mSystem;
		
		// initial density
		Array<Density> mInitialDensity;
		
		// fixed initial conditions (for all Systems)
		Array<double> mInitialCondition;
		
		// random number generator
		gsl_rng *mpGSLRandomNumberGenerator;
	};



	inline Problem::Problem()
	{
		mNumModes = -1;
		mNumResolvedModes = -1;
		mNumUnresolvedModes = -1;

		mCurrentTime = 0.0;
		mState = NO_PROBLEM_STATE;
		
		mpGSLRandomNumberGenerator = NULL;
		
		return;
	} 
}

#endif // _problem_h_	
