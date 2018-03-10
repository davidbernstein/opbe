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

#include "deltaproblem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include <gsl/gsl_randist.h>

using namespace NAMESPACE;
using namespace std;

DeltaProblem::DeltaProblem(const DeltaProblem &sol)
{
	ThrowException("DeltaProblem : copy constructor not implemented");
	return;
}



void DeltaProblem::Run(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// allocate system array, etc.
	Initialize();
	
	Clock clock;
	if (mRunControl.RunClockOn())
		clock.Start();

	Run();
	
	clock.StopAndPrintTime();
	
	mSystem[0].CleanUpSolver();
	
	
    return;
}



void DeltaProblem::Run()
{
	Reset();
		
	mRunControl.SetState(SYSTEM_RUN);
	
	for (long i = 0; i < mRunControl.NumOutputTimes(); ++i) {
		Evolve(mRunControl.OutputTime(i));		
		WriteOutput();
		WriteVolterraFFile();
	}		
	
		
	mRunControl.SetState(SYSTEM_STOP);
	
    return;
}



void DeltaProblem::Initialize()
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
		
	// set current time
	mCurrentTime = mRunControl.StartTime();
		
	mSystem.SetSize(3);
	
	for (short i = 0; i < 3; ++i) {
		mSystem[i].SetRunControl(&mRunControl);
		mSystem[i].SetModeIndex(&mModeIndex);
		mSystem[i].SetNumModes(mNumModes);
		mSystem[i].SetOPBEParameter(&mOPBEParameter);
		mSystem[i].SetCurrentTime(mRunControl.StartTime());
	}
	
	// set initial conditions
	for (short i = 0; i < 3; ++i) 
		mSystem[i].SetInitialConditions(mInitialCondition);
	
	double x0 = mSystem[1].InitialCondition(0);
	mSystem[1].SetInitialCondition(0, x0 + mDeltaX1);
	mSystem[2].SetInitialCondition(1, mDeltaX2);
	
	// initialize all systems
	mSystem[0].InitializeSolver();
	

	return;
}



void DeltaProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	Problem::ReadInitialConditions(fileName);
	
	// find delta x1 and delta x2
	Parser parser(fileName);
	
	if (parser.FindFloat("deltax1=", mDeltaX1) == false)
		ThrowException("DeltaProblem::ReadInputFile : didn't find delta x1");
	
	if (parser.FindFloat("deltax2=", mDeltaX2) == false)
		ThrowException("DeltaProblem::ReadInputFile : didn't find delta x2");
		
		
	return;
}



void DeltaProblem::WriteVolterraFFile() 
{
	ofstream& fileStream = mRunControl.GetOutputStream(VOLTERRA_F0_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;
	
	fileStream << mCurrentTime<< " ";
		
	for (short i = 0; i < mNumResolvedModes; ++i) {
		fileStream << VolterraF0(i);
			
		if (i != mNumResolvedModes - 1)
			fileStream << " ";
	}
		
	fileStream << endl;
	
	
	return;
}



double DeltaProblem::VolterraF0(short modeIndex) const
{
	double f0 = mSystem[0].ResolvedNoise(modeIndex);
	double f1 = mSystem[1].ResolvedNoise(modeIndex);
	double f2 = mSystem[2].ResolvedNoise(modeIndex);
	
	double diff1 = (f1 - f0) / mDeltaX1;
	double diff2 = (f2 - f0) / mDeltaX2;
	
	double epsilon = mOPBEParameter.ViscosityCoefficient();
	double a1 = mSystem[0].InitialCondition(0);
	
	return -epsilon * a1 * diff1 - 0.5 * a1 * a1 * diff2;
}