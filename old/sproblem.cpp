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

#include "sproblem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"
#include "hermitepolynomial.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace NAMESPACE;
using namespace std;

SProblem::SProblem(const SProblem &sol)
{
	ThrowException("SProblem : copy constructor not implemented");
	return;
}



void SProblem::Run(const std::string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array and set initial conditions
	Initialize();
	
	// run
	mRunControl.SetState(SYSTEM_RUN);
		
	double volterraDt = (mRunControl.EndTime() - mRunControl.StartTime()) / mNumVolterraTimeSteps;
	
	ComputeVolterraCoefficients();
	
	for (long n = 0; n < mNumVolterraTimeSteps; ++n) {
		Evolve(mCurrentTime + volterraDt);
		ComputeVolterraCoefficients();
	}
	
	mRunControl.SetState(SYSTEM_STOP);
	mSystem[0].CleanUpSolver();

    return;
}



void SProblem::ComputeVolterraCoefficients()
{
	
	return;
}



double SProblem::ComputeAverage(const Array<Array<double> > &f) const
{
	Array<double> f1(f.Size());
	
	for (long i = 0; i < f.Size(); ++i)
		f1[i] = ComputeAverage(f[i], RESOLVED_MODE); 
	
	return ComputeAverage(f1, UNRESOLVED_MODE);
}



double SProblem::ComputeAverage(const Array<double> &f, ModeType modeType) const
{
	// extended trapezoidal rule
	const double *pRho = (modeType == RESOLVED_MODE) ? mDensityResolved.Begin() : mDensityUnresolved.Begin();
	
	double sum= 0.5 * pRho[0] * f[0];
	long i;
	for (i = 1; i < f.Size() - 1; ++i)
		sum += pRho[i] * f[i];
		
	sum += 0.5 * pRho[i] * f[i];
				
	return sum * Spacing(modeType);
}



void SProblem::ComputeLiouvillian(const Array<Array<double> > &f, Array<Array<double> > &liouvillian) const
{
	Array<Array<double> > u;
	GetMode(u, 1);

	return;
}



void SProblem::ComputeProjection(const Array<Array<double> > &f, Array<double> &projection) const
{


	return;
}



void SProblem::GetMode(Array<Array<double> > &u, long modeIndex) const
{	
	long numResolvedID = mResolvedID.Size();
	long numUnresolvedID = mUnresolvedID.Size();
	
	u.SetSize(numResolvedID);
	for (short i = 0; i < numResolvedID; ++i) {
		u[i].SetSize(numUnresolvedID);
		for (short j = 0; j < numUnresolvedID; ++j) {
			long index = i * numResolvedID + j;
			u[i][j] = mSystem[index].GetMode(index);
		}
	}
	

	return;
}



void SProblem::Initialize()
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
	
	SetInitialDataAndDensity();
	
	// set current time
	mCurrentTime = mRunControl.StartTime();
	
	long numResolvedID = mResolvedID.Size();
	long numUnresolvedID = mUnresolvedID.Size();
	
	mSystem.SetSize(numResolvedID * numUnresolvedID);
	for (short i = 0; i < mResolvedID.Size(); ++i) {
		for (short j = 0; j < mUnresolvedID.Size(); ++j) {			
			long index = i * numResolvedID + j;
			
			mSystem[index].SetRunControl(&mRunControl);
			mSystem[i].SetModeIndex(&mModeIndex);
			mSystem[index].SetNumModes(mNumModes);
			mSystem[index].SetOPBEParameter(&mOPBEParameter);
			mSystem[index].SetCurrentTime(mRunControl.StartTime());
			
			mSystem[index].SetInitialCondition(mResolvedMode, mResolvedID[i]);
			mSystem[index].SetInitialCondition(mUnresolvedMode, mUnresolvedID[j]);
		}
	}
	
	// initialize all systems
	mSystem[0].InitializeSolver();
	
	
	return;
}



void SProblem::SetInitialDataAndDensity()
{
	short numResolvedID = 10, numUnresolvedID = 10;
	
	SetInitialDataAndDensity(RESOLVED_MODE, numResolvedID);
	SetInitialDataAndDensity(UNRESOLVED_MODE, numUnresolvedID);
	
	return;
}



void SProblem::SetInitialDataAndDensity(ModeType modeType, short numPts)
{
	Density *pDensity;
	if (modeType == RESOLVED_MODE)
		pDensity = &mInitialDensity[0];
	else
		pDensity = &mInitialDensity[1];
	
	double dMin, dMax;
	//HermitePolynomial hermite;
	//hermite.DomainBounds(DEFAULT_HERMITE_INTEGRATION_LIMIT, dMin, dMax);
	
	double spacing = (dMax - dMin) / (numPts - 1);
	
	double *pX, *pData;
	if (modeType == RESOLVED_MODE) {
		mResolvedID.SetSize(numPts);
		pX = mResolvedID.Begin();
		
		mDensityResolved.SetSize(numPts);
		pData = mDensityResolved.Begin();
	}
	else {
		mUnresolvedID.SetSize(numPts);
		pX = mUnresolvedID.Begin();
		
		mDensityUnresolved.SetSize(numPts);
		pData = mDensityUnresolved.Begin();
	}

	for (short i = 0; i < numPts; ++i) {
		pX[i] = dMin + i * spacing;
		pData[i] = pDensity->Value(pX[i]);
	}


	return;
}



void SProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	// find resolved and uresolved modes
	Parser parser(fileName);
	
	Array<double> parameter(3);
	if (parser.FindBracedFloats("resolvedmode={", parameter) == false)
		ThrowException("SProblem::ReadInputFile : file " + fileName + " does not contain resolved mode");

	mResolvedMode = mModeIndex(parameter);
	
	if (parser.FindBracedFloats("unresolvedmode={", parameter) == false)
		ThrowException("SProblem::ReadInputFile : file " + fileName + " does not contain unresolved mode");

	mUnresolvedMode = mModeIndex(parameter);
	
	CheckDensities();
	
	// volterra equation
	if (parser.FindInteger("numberofvolterratimesteps=", mNumVolterraTimeSteps) == false)
		ThrowException("SProblem::ReadInputFile : didn't find number of Volterra time steps");
	
	if (mNumVolterraTimeSteps < 1)
		ThrowException("SProblem::ReadInputFile : number of Volterra time steps less than 1");
		
	long finiteRankSize;
	if (parser.FindInteger("finiterankexpansionsize=", finiteRankSize) == false)
		ThrowException("SProblem::ReadInputFile : didn't find number of finite rank size");
	
	if (finiteRankSize < 1)
		ThrowException("SProblem::ReadInputFile : finite rank size less than 1");
		
	mFiniteRankSize = (short) finiteRankSize;
	
	return;
}



void SProblem::CheckDensities()
{
	// makes sure that resolved mode is always first in mInitialDensity
	
	Array<short> k;
	mInitialDensity[0].GetModeIndices(k);

	if (k[0] != mResolvedMode) 
		swap(mInitialDensity[0], mInitialDensity[1]);
		
	return;
}



void SProblem::Test(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array and set initial conditions
	Initialize();
	
	Array<double> f;
	f = mDensityResolved;
	
	for (short i = 0; i < f.Size(); ++i)
		f[i] = 1.0;
		
	double fBar = ComputeAverage(f, RESOLVED_MODE);
	
	cout << fBar << endl;
	
	return;
}