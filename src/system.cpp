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

#include "system.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace NAMESPACE;
using namespace std;

System::System(const System &sol)
{
	ThrowException("System : copy constructor not implemented");
	return;
}



void System::Evolve(double t1)
{
	GSLEvolve(t1);
    return;
}



void System::InitializeSolver()
{
	if (mpRunControl->State() != SYSTEM_INITIALIZE)
		ThrowException("System::InitializeSolver : RunControl does not indicate initialization");
		
	GSLEvolve(0.0);
	
	return;
}



void System::CleanUpSolver()
{
	if (mpRunControl->State() != SYSTEM_STOP)
		ThrowException("System::CleanUpSolver : RunControl does not indicate initialization");
	
	GSLEvolve(0.0);
}



void System::SetInitialConditions(const Array<double> &ic)
{
	if (ic.Size() != mInitialCondition.Size())
		ThrowException("System::SetInitialConditions : initial data array wrong size");
	
	mInitialCondition = ic;
	
	return;
}



void System::SetToInitialCondition()
{
	mMode = mInitialCondition;
	return;
}

	
	
void System::SetNumModes(long numModes)
{
	if (numModes <= 0)
		ThrowException("System::SetNumTotalAndResolvedModes : total number of modes must be positive");
		
	mMode.SetSize(numModes);
	mInitialCondition.SetSize(numModes);
		
	// all modes are initialized to zero
	for (long i = 0; i < numModes; ++i)
		mMode[i] = 0.0;
		
	return;
}



double System::ResolvedNoise(long index) const
{
	// for Burgers equation
	
	// make sure index is in resolved range
	long m = index + 1;
	if (!mpModeIndex->InResolvedRange(m))
		ThrowException("System::ResolvedNoise : index not in resolved range");
		
	long numResolved = mpModeIndex->NumResolvedModes();
	long numModes = mpModeIndex->NumModes();
	
	double sum = 0.0;
	for (long kp = numResolved - m + 1; kp <= numModes - m; ++kp)
		sum += U(kp) * U(m + kp);
	
	return 0.5 * m * sum;
}



double System::ResolvedNoise(long i, long j, long k) const
{
	// for NS
	return 0.0;
}



double System::Energy(long modeIndex) const
{
	return ONE_OVER_2PI * SqrNorm(modeIndex);
}



void System::ComputeMoments(Array<double> &moment, short maxMoment) const
{
	moment.SetSize(maxMoment + 1);
	for (short m = 0; m <= maxMoment; ++m)
		moment[m] = 0.0;
	
	moment[0] = mMode.Size();
	
	for (long k = 0; k < mMode.Size(); ++k) {
		double power = mMode[k];
		
		for (short m = 1; m <= maxMoment; ++m) {
			moment[m] += power;
			power *= mMode[k];
		}
	}
	
	// divide by PI
	for (short m = 0; m <= maxMoment; ++m)
		moment[m] *= 2.0 * ONE_OVER_2PI;
		
	
	return;
}



double System::Norm(long modeIndex) const
{
	return sqrt(SqrNorm(modeIndex));
}



double System::SqrNorm(long modeIndex) const
{
	long n = (modeIndex == -1) ? mMode.Size() : modeIndex;
	
	double sum = 0.0;
	for (long k = 0; k < n; ++k) 
		sum += mMode[k] * mMode[k];
	
	return sum;
}



void System::PrintInitialConditions() const
{
	for (long i = 0; i < mMode.Size(); ++i)
		cout << i << " " << mInitialCondition[i] << endl;
	
	
	return;
}