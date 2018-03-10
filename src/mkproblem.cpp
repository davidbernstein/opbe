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

#include "mkproblem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"
#include "hermitepolynomial.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include <gsl/gsl_randist.h>

using namespace NAMESPACE;
using namespace std;

MKProblem::MKProblem(const MKProblem &sol)
{
	ThrowException("MKProblem : copy constructor not implemented");
	return;
}



void MKProblem::Run(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array, etc.
	Initialize();
	
	Clock clock;
	if (mRunControl.RunClockOn())
		clock.Start();

	for (mRunCount = 1; mRunCount <= mNumMonteCarloRuns; ++mRunCount) {
		Run();
		mRunControl.PrintRunCount(mRunCount);
	}
	
	clock.StopAndPrintTime();
	
	mSystem[0].CleanUpSolver();

	WriteVolterraFFile();
	
	
    return;
}



void MKProblem::Run()
{
	ComputeInitialData();
	Reset();
	
	SetBigS();
	
	mRunControl.SetState(SYSTEM_RUN);
	
	for (long i = 0; i < mRunControl.NumOutputTimes(); ++i) {
		Evolve(mRunControl.OutputTime(i));
		UpdateVolterraCoefficients(i);
	}		
	
		
	mRunControl.SetState(SYSTEM_STOP);
	
    return;
}



void MKProblem::UpdateVolterraCoefficients(long timeStep)
{
	UpdateVolterraFAverage(timeStep);

	return;
}



void MKProblem::UpdateVolterraFAverage(long timeStep)
{
	static double fa, fb;
	
	if (mRunCount == 1) {
		fa = 0.0;
		fb = 1.0;
	}
	else {
		fa = 1.0 / mRunCount;
		fb = 1.0 - fa;
	}
	
	for (short i = 0; i < mNumResolvedModes; ++i) {
		double newValue = -mSystem[0].ResolvedNoise(i) * mBigS;
		//double oldValue = mVolterraF0(timeStep, i);
		mVolterraF0(timeStep, i) = fb * newValue + fa * mVolterraF0(timeStep, i);
	}
	
	
	return;
}



void MKProblem::SetBigS() 
{
	Array<double> rhs;
	
	double sumKSquared = mNumModes * (2.0 * mNumModes * mNumModes + 3.0 * mNumModes + 1) / 6.0;
	
	double divR = -sumKSquared * mOPBEParameter.ViscosityCoefficient();
	for (long k = 2; k <= mNumModes; k = k + 2) 
		divR += 0.5 * k * mSystem[0].U(k);
	
	mSystem[0].RHS(rhs);
	double sum = 0.0;
	for (long i = 0; i < mNumModes; ++i) {
		double mean = mInitialDensity[i].GetParameter(0);
		double sigma = mInitialDensity[i].GetParameter(1);
		
		double densityTerm = -(mSystem[0].InitialCondition(i) - mean) / (sigma * sigma);
		sum += rhs[i] * densityTerm;
	}
	
	mBigS = divR + sum;
	
	return;
}
	


void MKProblem::ComputeInitialData()
{
	for (long i = 0; i < mNumModes; ++i) {
		if (mInitialDensity[i].GetModeIndex() != i)
			ThrowException("MKProblem::ComputeInitialData : bad initial density array");
		
		double mean = mInitialDensity[i].GetParameter(0);
		double sigma = mInitialDensity[i].GetParameter(1);
		
		double x = GaussianRandomVariable(mean, sigma);
		
		mSystem[0].SetInitialCondition(i, x);
	}
	
	
	return;
}



void MKProblem::Initialize()
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
	
	InitializeRandomNumberGenerator();
	
	SetInitialDensities();
				
	// Volterra coefficients
	mVolterraF0.SetSize(mRunControl.NumOutputTimes(), mNumResolvedModes);
		
	// set current time
	mCurrentTime = mRunControl.StartTime();
		
	mSystem.SetSize(1);
	
	mSystem[0].SetRunControl(&mRunControl);
	mSystem[0].SetModeIndex(&mModeIndex);
	mSystem[0].SetNumModes(mNumModes);
	mSystem[0].SetOPBEParameter(&mOPBEParameter);
	mSystem[0].SetCurrentTime(mRunControl.StartTime());
	
	// initialize all systems
	mSystem[0].InitializeSolver();
	
	
	return;
}



double MKProblem::GaussianRandomVariable(double mean, double sigma) const
{
	return gsl_ran_gaussian(mpGSLRandomNumberGenerator, sigma) + mean;
}



void MKProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	// find resolved and uresolved modes
	Parser parser(fileName);
		
	// random number generator
	string rngName;
	if (parser.FindString("gslrandomnumbergenerator=", rngName))
		mRunControl.SetGSLRandomNumberGeneratorName(rngName);
	
	long randomSeed;
	if (parser.FindInteger("randomseed=", randomSeed))
		mRunControl.SetRandomSeed(abs(randomSeed));
	
	// number of runs
	if (parser.FindInteger("numberofruns=", mNumMonteCarloRuns) == false)
		ThrowException("MKProblem::ReadInputFile : didn't find number of monte carlo runs");
	
	// print run count
	long increment;
	if (parser.FindInteger("printruncountincrement=", increment))
		mRunControl.SetPrintRunCountIncrement(increment);
		
		
	return;
}



void MKProblem::SetInitialDensities()
{
	// incoming initial density should already have some data
	if (mInitialDensity.Size() != 2)
		ThrowException("MKProblem::SetInitialDensities : density information missing");
	
	Array<Density> dum(mNumModes - 2);
	mInitialDensity += dum;
	
	for (long i = 2; i < mNumModes; ++i) {
		mInitialDensity[i] = mInitialDensity[1];
		mInitialDensity[i].SetModeIndex(i);
	}
	
	
	return;
}



void MKProblem::WriteVolterraFFile() 
{
	ofstream& fileStream = mRunControl.GetOutputStream(VOLTERRA_F0_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;
	
	for (long n = 0; n < mRunControl.NumOutputTimes(); ++n) {
		fileStream << mRunControl.OutputTime(n) << " ";
		
		for (short i = 0; i < mNumResolvedModes; ++i) {
			fileStream << mVolterraF0(n, i);
			
			if (i != mNumResolvedModes - 1)
				fileStream << " ";
		}
		
		fileStream << endl;
	}
	
	
	return;
}



void MKProblem::Test(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array and set initial conditions
	Initialize();
	
	for (short i = 0; i < 1000; ++i)
		cout << GaussianRandomVariable(1.0, 0.1) << endl;
	
    return;
}
