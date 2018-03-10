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

#ifndef _system_h_
#define _system_h_

#include "array.h"
#include "opbeenums.h"
#include "runcontrol.h"
#include "modeindex.h"
#include "opbeparameter.h"

#include <fstream>

namespace NAMESPACE {
	class System {
	 public:
        System(void);
		~System(void) { };
        
		// copy constructor
		System(const System &sol);
		
        // run
		void InitializeSolver(void);
		void Evolve(double t1);
		void SetRunControl(const RunControl *pRunControl);
		void SetModeIndex(const ModeIndex *pModeIndex);
		void CleanUpSolver(void);
		
		// time
		void SetCurrentTime(double t);
		
		// modes
		void SetNumModes(long numModes);
		double GetMode(long modeIndex) const;
		double U(long i) const;
		double U(long i, long j, long k) const;
		void RHS(Array<double> &rhs) const;
		void RatioTModel(Array<double> &ratio) const;
		
		// initial conditions
		void SetInitialCondition(long modeIndex, double value);
		void SetInitialConditions(const Array<double> &ic);
		void SetToInitialCondition(void);
		double InitialCondition(long modeIndex) const;
		void PrintInitialConditions(void) const;

		// parameters
		void SetOPBEParameter(OPBEParameter *pParam);
		
		// properties of System
		double Energy(long modeIndex = -1) const;
		void ComputeMoments(Array<double> &moment, short maxMoment) const;
		double SqrNorm(long modeIndex = -1) const;
		double Norm(long modeIndex = -1) const;
		double ResolvedNoise(long i) const;
		double ResolvedNoise(long i, long j, long k) const;
		
	private:
		void GSLEvolve(double t1);
		
		// member data
	protected:
		// mode values
		Array<double> mMode;
		Array<double> mInitialCondition;
		
		// run controller
		RunControl *mpRunControl;
		
		// parameters
		OPBEParameter *mpOPBEParameter;
		
		// mode index
		ModeIndex *mpModeIndex;
		
		// time
		double mCurrentTime;
	};



	inline System::System()
	{
		mCurrentTime = 0.0;
		mpRunControl = NULL;
		mpModeIndex = NULL;
		mpOPBEParameter = NULL;
		
		return;
	} 
	
	
	
	inline void System::SetRunControl(const RunControl *pRunControl)
	{
		mpRunControl = (RunControl*) pRunControl;
		return;
	}
	
	
	
	inline void System::SetModeIndex(const ModeIndex *pModeIndex)
	{
		mpModeIndex = (ModeIndex*) pModeIndex;
		return;
	}
	
	
	
	inline void System::SetOPBEParameter(OPBEParameter *pParam)
	{
		mpOPBEParameter = pParam;
		return;
	}

	
	
	inline void System::SetCurrentTime(double t)
	{
		mCurrentTime = t;
		return;
	}
	
	
	
	inline void System::SetInitialCondition(long modeIndex, double value)
	{	
		mInitialCondition[modeIndex] = value;
		return;
	}
	
	
	
	inline double System::InitialCondition(long modeIndex) const
	{
		return mInitialCondition[modeIndex];
	}
	
	
	
	inline double System::GetMode(long modeIndex) const
	{
		return mMode[modeIndex];
	}
	
	
	
	inline double System::U(long i) const
	{
		return mMode[(*mpModeIndex)(i)];
	}
	
	
	
	inline double System::U(long i, long j, long k) const
	{
		return mMode[(*mpModeIndex)(i, j, k)];
	}
}

#endif // _system_h_	
