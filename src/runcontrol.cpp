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

#include "runcontrol.h"

#include <iostream>
#include <vector>

using namespace NAMESPACE;
using namespace std;



void RunControl::SetSystemType(string systemType)
{
	if (mSystemType != NO_SYSTEM_TYPE)
		ThrowException("Solution::SetSystemType : system type already set");
		
	if (systemType == "burgersequation") {
		mSystemType = BURGERS_EQUATION;
		return;
	}
	
	if (systemType == "navierstokes") {
		mSystemType = NAVIER_STOKES;
		return;
	}
	
	ThrowException("Solution::SetSystemType : invalid system type");
	
	return;
}



void RunControl::MakeOutputSchedule(double timeStep)
{
	switch (mOutputScheduleMode) {
	case OUTPUT_SCHEDULE_LINEAR:
		MakeOutputScheduleLinear(timeStep);
		break;
	
	case OUTPUT_SCHEDULE_LOGARITHMIC:
		MakeOutputScheduleLogarithmic(timeStep);
		break;
	
	default:
		ThrowException("RunControl::MakeOutputSchedule : output schedule type not set");
		break;
	}

	return;
}



void RunControl::MakeOutputScheduleLinear(double timeStep)
{
	// mOutputSchedule includes the end points mStartTime and mEndTime
	
	if (timeStep <= 0.0)
		ThrowException("RunControl::MakeOutputScheduleLinear : non-positive output time step");
	
	vector<double> t;

	double t1 = mStartTime;
	t.push_back(t1);
	
	while (t1 < mEndTime) {
		t1 += timeStep;
		
		if (mEndTime - t1 < 0.5 * timeStep)
			t1 = mEndTime;
		
		t.push_back(t1);
	}
	
	mOutputSchedule = t;

	return;
}



void RunControl::MakeOutputScheduleLogarithmic(double timeStep)
{
	// timestep is interpreted as a multiplier in this routine
	if (timeStep <= 1.0)
		ThrowException("RunControl::MakeOutputScheduleLogarithmic : timeStep is multiplier, but is less than 1");
		
	double startTime = mStartTime;
	double endTime = mEndTime;
	
	mStartTime = log10(DEFAULT_FIRST_OUTPUT_TIME);
	mEndTime = log10(endTime - startTime);
	
	MakeOutputScheduleLinear(log10(timeStep));
	
	mStartTime = startTime;
	mEndTime = endTime;
	
	for (long i = 0; i < mOutputSchedule.Size(); ++i) 
		mOutputSchedule[i] = mStartTime + Power(10.0, mOutputSchedule[i]);

	return;
}



void RunControl::SetOutputScheduleMode(string outputScheduleMode)
{
	if (outputScheduleMode == "linear") {
		mOutputScheduleMode = OUTPUT_SCHEDULE_LINEAR;
		return;
	}
	
	if (outputScheduleMode == "logarithmic") {
		mOutputScheduleMode = OUTPUT_SCHEDULE_LOGARITHMIC;
		return;
	}
	
	ThrowException("RunControl::SetOutputScheduleMode : bad input string = " + outputScheduleMode);
	return;
}

	
	
void RunControl::PrintRunCount(long runCount) const
{
	if (mPrintRunCountIncrement > 0) {
		if (runCount % mPrintRunCountIncrement == 0)
			std::cout << "Run count = " << runCount << endl;
	}
	
		
	return;
}