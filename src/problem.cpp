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

#include "problem.h"
#include "parser.h"
#include "opbeconst.h"
#include "constants.h"
#include "clock.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace NAMESPACE;
using namespace std;

Problem::Problem(const Problem &sol)
{
	ThrowException("Problem : copy constructor not implemented");
	return;
}



void Problem::Evolve(double t1)
{
	// evolve all systems from current time to t1
	for (long i = 0; i < mSystem.Size(); ++i) 
		mSystem[i].Evolve(t1);
		
	mCurrentTime = t1;
	mState = PROBLEM_RUNNING;
	    
		
    return;
}



void Problem::Reset()
{
	// set current time
	mCurrentTime = mRunControl.StartTime();
	
	for (short i = 0; i < mSystem.Size(); ++i) {
		mSystem[i].SetCurrentTime(mRunControl.StartTime());
		mSystem[i].SetToInitialCondition();
	}

	return;
}



ProblemType Problem::GetProblemType(const std::string &fileName)
{
	Parser parser(fileName);
	
	string problemName;
	if (parser.FindString("problemtype=", problemName) == false)
			ThrowException("Problem::GetProblemType : file " + fileName + " does not contain problem type");

	if (problemName == "memorykernelproblem")
		return MK_PROBLEM;
		
	if (problemName == "averagingproblem")
		return AVERAGING_PROBLEM;
	
	if (problemName == "fixedicproblem")
		return FIXED_IC_PROBLEM;
		
	if (problemName == "deltaproblem")
		return DELTA_PROBLEM;
		
	ThrowException("Problem::GetProblemType : undefined problem type");
	
	
	return NO_PROBLEM_TYPE;
}



void Problem::ReadInputFile(const string &fileName)
{	
	Parser parser(fileName);
	
	mRunControl.SetInputDirectory(fileName);

	// system type
	string systemType;
	if (parser.FindString("system=", systemType))
		mRunControl.SetSystemType(systemType);
	else
		ThrowException("Problem::ReadInputFile : system type not found in file " + fileName);
		
	// reynolds number
	double reynoldsNumber;
	if (parser.FindFloat("reynoldsnumber=", reynoldsNumber) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain Reynolds number");
		
	mOPBEParameter.SetReynoldsNumber(reynoldsNumber);
	
	// start and end times
	double startTime, endTime;
	if (parser.FindFloat("starttime=", startTime) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain start time");
		
	if (parser.FindFloat("endtime=", endTime) == false)
		ThrowException("Simulation::ReadInputFile : file " + fileName + " does not contain end time");
		
	mRunControl.SetStartAndEndTimes(startTime, endTime);
	
	// output times
	string outputScheduleMode;
	if (parser.FindString("outputtimemode=", outputScheduleMode))
		mRunControl.SetOutputScheduleMode(outputScheduleMode);
	else
		mRunControl.SetOutputScheduleMode(OUTPUT_SCHEDULE_LINEAR);
		
	double outputTimeStep;
	if (parser.FindFloat("outputtimestep=", outputTimeStep)) 
		mRunControl.MakeOutputSchedule(outputTimeStep);
	else
		mRunControl.TurnOffOutput();
		
	// number of modes
	Array<double> parameter(3);
	if (parser.FindBracedFloats("numberofmodes={", parameter) == false)
		ThrowException("SProblem::ReadInputFile : file " + fileName + " does not contain number of modes");

	mModeIndex.Set(parameter);
	mNumModes = mModeIndex.Max();
	
	long numResolved = 0;
	if (parser.FindInteger("numberofresolvedmodes=", numResolved) == false) {
		mNumUnresolvedModes = 0;
		mNumResolvedModes = mNumModes;
	}
	else {
		if (numResolved <= 0)
			ThrowException("AveragingProblem::ReadInputFile : non-positive number of unresolved modes");
		
		mNumUnresolvedModes = mNumModes - numResolved;
		mNumResolvedModes = numResolved;
	}
	
	mModeIndex.SetNumResolvedAndUnresolvedModes(mNumResolvedModes, mNumUnresolvedModes);
	
	// t-model
	string dum;
	if (parser.FindString("t-model=on", dum)) 
		mRunControl.TurnOnTModel();
	
	// gsl solver
	string solverName;
	parser.FindString("gslsolver=", solverName);
	mRunControl.SetGSLSolverName(solverName);
	
	double odeError;
	if (parser.FindFloat("gslrelativeerror=", odeError))
		mRunControl.SetLocalRelativeError(odeError);
	else
		mRunControl.SetLocalRelativeError(DEFAULT_LOCAL_RELATIVE_ERROR);
	
	if (parser.FindFloat("gslabsoluteerror=", odeError))
		mRunControl.SetLocalAbsoluteError(odeError);
	else
		mRunControl.SetLocalAbsoluteError(DEFAULT_LOCAL_ABSOLUTE_ERROR);
			
	// initial densities
	list<Density> densityList;
	Density density;
	string prefix = "densitytype";
	for (short i = 1; i <= mNumModes; ++i) {
		string n = ConvertIntegerToString(i);
		string name;
		if (parser.FindString(prefix + n + "=", name)) {
			density.SetType(name);
			parameter.SetSize(density.NumParametersToSpecify());
			if (parser.FindBracedFloats("density" + n + "={", parameter) == false)
				ThrowException("Simulation::ReadInputFile : didn't find all parameters for density " + n);
			
			density.SetParameters(parameter);
			densityList.push_back(density);
		}
	}   
	
	mInitialDensity = densityList;
	
	// output directory
	string outputName;
	if (parser.FindFileName("outputdirectory=", outputName) == false) {
		mRunControl.SetOutputDirectoryToInputDirectory();
	}
	else {
		mRunControl.SetOutputDirectory(outputName);
	}
	cout << "Output directory is " + mRunControl.OutputDirectory() << endl;
	
	// open output file streams
	if (parser.FindFileName("modefile=", outputName)) 
		mRunControl.OpenOutputStream(MODE_OUTPUT_STREAM, outputName);
	
	if (parser.FindFileName("energyfile=", outputName)) 
		mRunControl.OpenOutputStream(ENERGY_OUTPUT_STREAM, outputName);
		
	if (parser.FindFileName("momentsfile=", outputName)) 
		mRunControl.OpenOutputStream(MOMENTS_OUTPUT_STREAM, outputName);
			
	if (parser.FindFileName("tmodelratiofile=", outputName)) 
		mRunControl.OpenOutputStream(TMODEL_RATIO_OUTPUT_STREAM, outputName);
	
	if (parser.FindFileName("volterraffile=", outputName))
		mRunControl.OpenOutputStream(VOLTERRA_F0_OUTPUT_STREAM, outputName);
		
	// clock
	if (parser.FindString("runclock=on", dum))
		mRunControl.TurnOnRunClock();
		
	// print run time
	if (parser.FindString("printruntime=on", dum))
		mRunControl.TurnOnRunClock();
	
	// print time
	if (parser.FindString("printoutputtime=on", dum))
		mRunControl.TurnOnPrintOutputTime();
		
		
	return;
}



void Problem::ReadInitialConditions(const std::string &fileName)
{
	mInitialCondition.SetSize(mNumModes);
	for (long i = 0; i < mNumModes; ++i)
		mInitialCondition[i] = 0.0;
		
	Parser parser(fileName);
	
	// first look for initial conditions in file fileName
	bool found = false;
	Array<double> parameter(4);
	for (long i = 0; i < mNumModes; ++i) {
		string n = ConvertIntegerToString(i);
		if (parser.FindBracedFloats("initialcondition" + n + "={", parameter)) {
			mInitialCondition[mModeIndex(parameter)] = parameter[3];
			cout << parameter[0] << " " << parameter[1] << " " << parameter[2] << " " << parameter[3] << endl;
			found = true;
		}
	}

	if (found)
		return;

	// else look for initial conditions in another file (specified in file fileName)
	string icFileName;
	if (parser.FindFileName("initialconditionsfile=", icFileName) == false) {
		cout << "Didn't find initial conditions in " + fileName << endl;
	}
	else {
		// icFileName is assume to be in the same directory as fileName
		icFileName = mRunControl.InputDirectory() + icFileName;
		
		ifstream file;
		OpenInputFile(icFileName, file);
		
		for (long i = 0; i < mNumModes; ++i) {
			if (file.eof() == false) {
				double x;
				file >> x;
				if (file.eof() == false)
					mInitialCondition[i] = x;
			}
			else {
				mInitialCondition[i] = 0.0;
			}
		}
	}
	
	return;
}



void Problem::InitializeRandomNumberGenerator()
{	
	string rngName = mRunControl.RandomNumberGeneratorName();
	if (rngName == "taus")
		mpGSLRandomNumberGenerator = gsl_rng_alloc(gsl_rng_taus);
					
	if (mpGSLRandomNumberGenerator == NULL)
		ThrowException("Problem : random number generator initialization failed for type " + rngName);
	
	gsl_rng_set(mpGSLRandomNumberGenerator, mRunControl.RandomSeed());
	
	
	return;
}



void Problem::WriteOutput()
{
	if (mRunControl.PrintOutputTime())
		cout << "time " << mCurrentTime << endl;
		
	WriteModes();
	WriteEnergy();
	WriteMoments();
	WriteTModelRatio();
	
	return;
}



void Problem::WriteModes()
{
	ofstream& fileStream = mRunControl.GetOutputStream(MODE_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;

	for (long s = 0; s < mSystem.Size(); ++s) {
		fileStream << mCurrentTime << " ";
		
		for (long i = 0; i < mNumModes; ++i) {
			fileStream << setprecision(10) << mSystem[s].GetMode(i);
			if (i != mNumModes - 1)
				fileStream << " ";
		}
	}
	
	fileStream << endl;
	
	return;
}



void Problem::WriteEnergy()
{
	// currently writes only the energy of mSystem[0] to the output stream
	
	ofstream& fileStream = mRunControl.GetOutputStream(ENERGY_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;
		 
	fileStream << mCurrentTime << " ";
	
	fileStream << mSystem[0].Energy() << " " << mSystem[0].Energy(mNumResolvedModes) << endl;
	

	return;
}



void Problem::WriteMoments()
{
	// writes only moments of mSystem[0] to file
	ofstream& fileStream = mRunControl.GetOutputStream(MOMENTS_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;

	fileStream << mCurrentTime << " ";
	
	Array<double> moment;
	mSystem[0].ComputeMoments(moment, DEFAULT_NUM_MOMENTS);
	
	for (short m = 1; m <= DEFAULT_NUM_MOMENTS; ++m) {
		fileStream << moment[m]; " ";
	
		if (m != DEFAULT_NUM_MOMENTS)
			fileStream << " ";
	}
	
	fileStream << endl;

	return;
}



void Problem::WriteTModelRatio()
{
	ofstream& fileStream = mRunControl.GetOutputStream(TMODEL_RATIO_OUTPUT_STREAM);

	if (fileStream.is_open() == false)
		return;

	Array<double> ratio;
	mSystem[0].RatioTModel(ratio);
	
	fileStream << mCurrentTime << " ";
		
	for (long i = 0; i < mNumModes; ++i) {
		fileStream << ratio[i];
		if (i != mNumModes - 1)
			fileStream << " ";
	}
	
	fileStream << endl;
	
	return;
}



void Problem::PrintCurrentTime() const
{
	cout << "t = " << mCurrentTime << endl;
	
	return;
}