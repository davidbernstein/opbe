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

#include "averagingproblem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace NAMESPACE;
using namespace std;

AveragingProblem::AveragingProblem(const AveragingProblem &sol)
{
	ThrowException("AveragingProblem : copy constructor not implemented");
	return;
}



void AveragingProblem::Run(const std::string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// run
	if (mNumUnresolvedModes == 1) {
		RunOneUnresolved();
		return;
	}
	
	ThrowException("AveragingProblem::Run : number of unresolved modes not supported");
	
	
    return;
}



void AveragingProblem::RunOneUnresolved()
{	
	switch (mQuadratureType) {
	case ADAPTIVE_QUADRATURE:
		ThrowException("AveragingProblem::Run : adaptive integration not implemented");
		break;
	
	case MONTE_CARLO_QUADRATURE:
		ThrowException("AveragingProblem::Run : Monte Carlo integration not implemented");
		break;
		
	case FIXED_SPACING_QUADRATURE:
		RunFixedSpacingOneUnresolved();
		break;
	
	default:
		ThrowException("AveragingProblem::Run : quadrature type not set");
	}
	
	
	return;
}



void AveragingProblem::RunFixedSpacingOneUnresolved()
{
	double xMin, xMax;
	mInitialDensity[0].DomainBounds(DEFAULT_QUADRATURE_TOLERANCE, xMin, xMax);

	// set initial conditions
	Initialize(mQuadratureNumGridPoints);
	double dX = (xMax - xMin) / (mQuadratureNumGridPoints - 1.0);
	for (long i = 0; i < mQuadratureNumGridPoints; ++i) {
		double x = xMin + i * dX;
		mSystem[i].SetInitialCondition(mNumModes - 1, x);
	}
	
	Initialize();
	mState = PROBLEM_START;
	AverageOneUnresolved();
	
	mRunControl.SetState(SYSTEM_RUN);
	for (long i = 0; i < mRunControl.NumOutputTimes(); ++i) {
		Evolve(mRunControl.OutputTime(i));
		AverageOneUnresolved();
	}
	
	mState = PROBLEM_DONE;
	
	return;
}



void AveragingProblem::AverageOneUnresolved()
{
/*
	static ofstream modeFile, averageFile;
	
	WriteModesToStream(modeFile);
	
	averageFile << mCurrentTime << " ";
	double dX = mSystem[1].InitialCondition(mNumModes - 1) - mSystem[0].InitialCondition(mNumModes - 1);
	for (long i = 0; i < mNumModes - 1; ++i) {
		double average = 0.0;
		
		// endpoints
		long j = 0;
		average += 0.5 * mSystem[j].GetMode(i) * mInitialDensity[0].Value(mSystem[j].InitialCondition(mNumModes - 1));
		
		j = mSystem.Size() - 1;
		average += 0.5 * mSystem[j].GetMode(i) * mInitialDensity[0].Value(mSystem[j].InitialCondition(mNumModes - 1));
		
		for (long j = 1; j < mSystem.Size() - 1; ++j) 
			average += mSystem[j].GetMode(i) * mInitialDensity[0].Value(mSystem[j].InitialCondition(mNumModes - 1));
		
		average *= dX;
		
		averageFile << average << " ";
	}
	averageFile << endl;
	*/
	return;
}



void AveragingProblem::Initialize(long numSystems)
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
	
	mAverage.SetSize(mRunControl.NumOutputTimes(), mNumResolvedModes);
	
	// set current time
	mCurrentTime = mRunControl.StartTime();
	
	mSystem.SetSize(numSystems);
		
	for (long i = 0; i < numSystems; ++i) {
		mSystem[i].SetRunControl(&mRunControl);
		mSystem[i].SetModeIndex(&mModeIndex);
		mSystem[i].SetNumModes(mNumModes);
		mSystem[i].SetOPBEParameter(&mOPBEParameter);
		mSystem[i].SetCurrentTime(mRunControl.StartTime());
		mSystem[i].SetInitialConditions(mInitialCondition);
	}
	// initialize solver for all systems
	mSystem[0].InitializeSolver();
	
	
	return;
}



void AveragingProblem::SetInitialDataOneUnresolved()
{


	return;
}



void AveragingProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	// find fixed initial conditions
	Problem::ReadInitialConditions(fileName);
	
	// find resolved and uresolved modes
	Parser parser(fileName);

	// quadrature type
	string quadratureType;
	if (parser.FindString("quadraturetype=", quadratureType) == false) 
		ThrowException("AveragingProblem::ReadInputFile : quadrature type not found in file " + fileName);
	
	SetQuadratureType(quadratureType);
		
	// initial quadrature level
	long initialLevel;
	if (parser.FindInteger("initialquadraturelevel=", initialLevel) == false) {
		mInitialQuadratureLevel = DEFAULT_INITIAL_QUADRATURE_LEVEL;
	}
	else {
		if (initialLevel <= 0)
			ThrowException("AveragingProblem::ReadInputFile : non-positive initial quadrature level");

		mInitialQuadratureLevel = initialLevel;
	}
	
	// quadrature tolerance
	mQuadratureTolerance = DEFAULT_QUADRATURE_TOLERANCE;
	
	// grid size for fixed spacing quadrature
	long gridSize;
	if (mQuadratureType == FIXED_SPACING_QUADRATURE) {
		if (parser.FindInteger("numberofquadraturepoints=", gridSize) == false) 
			mQuadratureNumGridPoints = DEFAULT_QUADRATURE_NUM_POINTS;
		else 
			mQuadratureNumGridPoints = gridSize;
	}
	
	// check to make sure all required densities have been set
	CheckDensity();
	
	
	return;
}



void AveragingProblem::CheckDensity() const
{
	ThrowException("AveragingProblem::CheckDensity : not implemented");
	
	return;
}



void AveragingProblem::SetQuadratureType(string quadratureType)
{
	if (quadratureType == "adaptive") {
		mQuadratureType = ADAPTIVE_QUADRATURE;
		return;
	}
	
	if (quadratureType == "montecarlo") {
		mQuadratureType = MONTE_CARLO_QUADRATURE;
		return;
	}
	
	if (quadratureType == "fixedspacing") {
		mQuadratureType = FIXED_SPACING_QUADRATURE;
		return;
	}

	ThrowException("AveragingProblem::SetQuadratureType : no type corresponds to string " + quadratureType);
	
	return;
}



void AveragingProblem::Test(const string &fileName)
{
	
	
	return;
}