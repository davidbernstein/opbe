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

#include "fixedicproblem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace NAMESPACE;
using namespace std;

FixedICProblem::FixedICProblem(const FixedICProblem &sol)
{
	ThrowException("FixedICProblem : copy constructor not implemented");
	return;
}



void FixedICProblem::Run(const string &fileName)
{
	// read input file
    ReadInputFile(fileName);
	
	// set initial conditions
	Initialize();
	Reset();
	
	mState = PROBLEM_START;
	mRunControl.SetState(SYSTEM_RUN);

	Clock clock;
	if (mRunControl.RunClockOn()) {
		clock.SetPrintMode(PRINT_SECONDS);
		clock.Start();
	}
		
	for (long i = 0; i < mRunControl.NumOutputTimes(); ++i) {
		Evolve(mRunControl.OutputTime(i));
		WriteOutput();
	}
	
	clock.StopAndPrintTime();
	
	mState = PROBLEM_DONE;
	
    return;
}



void FixedICProblem::Initialize()
{	
	mRunControl.SetState(SYSTEM_INITIALIZE);
		
	// set current time
	mCurrentTime = mRunControl.StartTime();
	
	mSystem.SetSize(1);

	mSystem[0].SetRunControl(&mRunControl);
	mSystem[0].SetModeIndex(&mModeIndex);
	mSystem[0].SetNumModes(mNumModes);
	mSystem[0].SetOPBEParameter(&mOPBEParameter);	
	mSystem[0].SetCurrentTime(mRunControl.StartTime());
	mSystem[0].SetInitialConditions(mInitialCondition);
	
	// initialize solver for all systems
	mSystem[0].InitializeSolver();
	
	
	return;
}



void FixedICProblem::ReadInputFile(const string &fileName)
{
	// read file
	Problem::ReadInputFile(fileName);
	
	// find fixed initial conditions
	Problem::ReadInitialConditions(fileName);
	
	return;
}
