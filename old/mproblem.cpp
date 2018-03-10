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

#include "mproblem.h"
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

MProblem::MProblem(const MProblem &sol)
{
	ThrowException("MProblem : copy constructor not implemented");
	return;
}



void MProblem::Run(const std::string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array and set initial conditions
	Initialize();
	//SetToInitialConditions();
	
	// run
	mRunControl.SetState(SYSTEM_RUN);
				
	for (long n = 0; n < mVolterraT.Size(); ++n) {
		Evolve(mVolterraT[n]);
		ComputeVolterraCoefficients(n);
		PrintCurrentTime();
	}
		
	mRunControl.SetState(SYSTEM_STOP);
	mSystem[0].CleanUpSolver();

	WriteVolterraFToFile("f.dat");
	
	
    return;
}



void MProblem::ComputeVolterraCoefficients(long timeStep)
{
	ComputeVolterraF(timeStep);
	
	return;
}



void MProblem::ComputeVolterraF(long timeStep)
{
	static Array<double> f1, f;
	f1.SetSize(mSystem.Size());
	f.SetSize(mResolvedID.Size());
	
	long numResolvedID = mResolvedID.Size();
	
	double epsilon = mOPBEParameter.ViscosityCoefficient();
	
	// noise 
	for (long i = 0; i < mSystem.Size(); ++i)
		f1[i] = mSystem[i].ResolvedNoise(1);
		
	// integrand
	double dx = mResolvedID[1] - mResolvedID[0];
	
	long i;
	double x1, df1dx1, term1, df1dx2, term2;
	for (i = 1; i < numResolvedID - 1; ++i) {
		x1 = mResolvedID[i];
		
		df1dx1 = (f1[i + 1] - f1[i - 1]) / (2.0 * dx);
		term1 = - epsilon * x1 * df1dx1;
		
		df1dx2 = (f1[2 * numResolvedID + i] - f1[numResolvedID + i]) / (2.0 * mSecondModeSpacing);
		term2 = -0.5 * x1 * x1 * df1dx2;
		
		f[i] = term1 + term2;
	}
	
	i = 0;
	x1 = mResolvedID[i];
	df1dx1 = (f1[i + 1] - f1[i]) / (dx);
	term1 = -epsilon * x1 * df1dx1;
	df1dx2 = (f1[2 * numResolvedID + i] - f1[numResolvedID + i]) / (2.0 * mSecondModeSpacing);
	term2 = -0.5 * x1 * x1 * df1dx2;
	f[i] = term1 + term2;
	
	i = numResolvedID - 1;
	x1 = mResolvedID[i];
	df1dx1 = (f1[i] - f1[i - 1]) / (dx);
	term1 = -epsilon * x1 * df1dx1;
	df1dx2 = (f1[2 * numResolvedID + i] - f1[numResolvedID + i]) / (2.0 * mSecondModeSpacing);
	term2 = -0.5 * x1 * x1 * df1dx2;
	f[i] = term1 + term2;
	
	mVolterraF[timeStep] = mHermitePolynomial.InnerProduct(f, mResolvedID, 0);
	
	
	return;
}



double MProblem::ComputeAverage(const Array<Array<double> > &f) const
{
	Array<double> f1(f.Size());
	
	for (long i = 0; i < f.Size(); ++i)
		f1[i] = ComputeAverage(f[i], RESOLVED_MODE); 
	
	return ComputeAverage(f1, UNRESOLVED_MODE);
}



double MProblem::ComputeAverage(const Array<double> &f, ModeType modeType) const
{
/*
	// extended trapezoidal rule
	const double *pRho = (modeType == RESOLVED_MODE) ? mDensityResolved.Begin() : mDensityUnresolved.Begin();
	
	double sum= 0.5 * pRho[0] * f[0];
	long i;
	for (i = 1; i < f.Size() - 1; ++i)
		sum += pRho[i] * f[i];
		
	sum += 0.5 * pRho[i] * f[i];
				
	return sum * Spacing(modeType);
	*/
	return 0.0;
}



void MProblem::Initialize()
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
	
	mHermitePolynomial.Initialize(mInitialDensity[0], mFiniteRankSize);
	
	mResolvedID = mHermitePolynomial.MakeGrid();
	
	// Volterra coefficients
	mVolterraF.SetSize(mVolterraT.Size());
		
	// set current time
	mCurrentTime = mRunControl.StartTime();
	
	long numResolvedID = mResolvedID.Size();
	
	mSystem.SetSize(3 * numResolvedID);
	
	for (short i = 0; i < mSystem.Size(); ++i) {
		mSystem[i].SetRunControl(&mRunControl);
		mSystem[i].SetModeIndex(&mModeIndex);
		mSystem[i].SetNumModes(mNumModes);
		mSystem[i].SetOPBEParameter(&mOPBEParameter);
		mSystem[i].SetCurrentTime(mRunControl.StartTime());
	}
	
	// set initial conditions
	for (short i = 0; i < mResolvedID.Size(); ++i) {
		long i1 = numResolvedID + i;
		long i2 = 2 * numResolvedID + i;
		
		mSystem[i].SetInitialCondition(mResolvedMode, mResolvedID[i]);
		mSystem[i1].SetInitialCondition(mResolvedMode, mResolvedID[i]);
		mSystem[i2].SetInitialCondition(mResolvedMode, mResolvedID[i]);
		
		mSystem[i1].SetInitialCondition(1, -mSecondModeSpacing);
		mSystem[i2].SetInitialCondition(1,  mSecondModeSpacing);
	}
	
	
	// initialize all systems
	mSystem[0].InitializeSolver();
	
	
	return;
}



void MProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	// find resolved and uresolved modes
	Parser parser(fileName);
	
	Array<double> parameter(3);
	if (parser.FindBracedFloats("resolvedmode={", parameter) == false)
		ThrowException("MProblem::ReadInputFile : file " + fileName + " does not contain resolved mode");

	mResolvedMode = mModeIndex(parameter);
	
	mNumResolvedModes = 1;
	mNumUnresolvedModes = mNumModes - 1;
	mModeIndex.SetNumResolvedAndUnresolvedModes(mNumResolvedModes, mNumUnresolvedModes);
		
	// volterra equation
	long numVolterraTimeSteps;
	if (parser.FindInteger("numberofvolterratimesteps=", numVolterraTimeSteps) == false)
		ThrowException("MProblem::ReadInputFile : didn't find number of Volterra time steps");
	
	if (numVolterraTimeSteps < 1)
		ThrowException("MProblem::ReadInputFile : number of Volterra time steps less than 1");
	
	SetVolterraTimeSchedule(numVolterraTimeSteps);
	
	long finiteRankSize;
	if (parser.FindInteger("finiterankexpansionsize=", finiteRankSize) == false)
		ThrowException("MProblem::ReadInputFile : didn't find number of finite rank size");
	
	if (finiteRankSize < 1)
		ThrowException("MProblem::ReadInputFile : finite rank size less than 1");
		
	mFiniteRankSize = (short) finiteRankSize;
	
	if (parser.FindFloat("secondmodespacing=", mSecondModeSpacing) == false)
		ThrowException("MProblem::ReadInputFile : didn't find second mode spacing");
	
	
	return;
}



void MProblem::SetVolterraTimeSchedule(long numSteps)
{
	mVolterraT.SetSize(numSteps + 1);
	
	double dT = (mRunControl.EndTime() - mRunControl.StartTime()) / numSteps;
	
	for (long i = 0; i < numSteps + 1; ++i)
		mVolterraT[i] = mRunControl.StartTime() + i * dT;
	
	return;
}



void MProblem::WriteVolterraFToFile(string fileName) const
{
	fileName = mRunControl.OutputDirectory() + fileName;
	ofstream file;
	OpenOutputFile(fileName, file);
	
	for (long n = 0; n < mVolterraT.Size(); ++n) {
		file << mVolterraT[n] << " ";
		file << mVolterraF[n] << endl;
	}
	
	
	return;
}