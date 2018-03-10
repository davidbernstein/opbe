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

#ifndef _runcontrol_h_
#define _runcontrol_h_

#include "array.h"
#include "namespace.h"
#include "opbeconst.h"
#include "opbeenums.h"
#include "utility.h"

namespace NAMESPACE {
	class RunControl {
	 public:
        RunControl(void);
		~RunControl(void) { };
		
		// system type
		void SetSystemType(std::string systemType);
		SystemType GetSystemType(void) const;
		
		// state
		void SetState(SystemState state);
		SystemState State(void) const;
		
		// time
		void SetStartAndEndTimes(double startTime, double endTime);
		double StartTime(void) const;
		double EndTime(void) const;
		void TurnOnPrintOutputTime(void);
		bool PrintOutputTime(void);
		
		// output
		long NumOutputTimes(void) const;
		double OutputTime(long n) const;
		void MakeOutputSchedule(double timeStep);
		void TurnOffOutput(void);
		void SetOutputScheduleMode(std::string outputScheduleMode);
		void SetOutputScheduleMode(OutputScheduleMode outputScheduleMode);
		
		// ensemble switches
		void SetPrintRunCountIncrement(long increment);
		void PrintRunCount(long runCount) const;
		
		// t-model
		void TurnOnTModel(void);
		void TurnOffTModel(void);
		bool TModelOn(void) const;
		
		// clock
		void TurnOnRunClock(void);
		bool RunClockOn(void) const;
		
		// random number generator
		void SetGSLRandomNumberGeneratorName(std::string rngName);
		std::string RandomNumberGeneratorName(void) const;
		void SetRandomSeed(unsigned long int seed);
		unsigned long int RandomSeed(void) const;
		
		// gsl solver name
		void SetGSLSolverName(std::string solverName);
		std::string SolverName(void) const;
		
        // scalar error controls
        void SetLocalRelativeError(double value);
		void SetLocalAbsoluteError(double value);
		double GetLocalRelativeError(void) const;
		double GetLocalAbsoluteError(void) const;
	
		// directories
		void SetOutputDirectory(std::string directoryName);
		void SetOutputDirectoryToInputDirectory(void);
		std::string OutputDirectory(void) const;
		void SetInputDirectory(std::string inputFile);
		std::string InputDirectory(void) const;
		
		// output file streams
		void OpenOutputStream(OutputFileStreamType streamType, std::string fileName);
		std::ofstream& GetOutputStream(OutputFileStreamType streamType);
	
	protected:
		void MakeOutputScheduleLinear(double timeStep);
		void MakeOutputScheduleLogarithmic(double timeStep);
		
		// member data
	protected:
		// system type and state
		SystemType mSystemType;
		SystemState mState;
		
		// time
		double mStartTime;
		double mEndTime;
		bool mPrintOutputTime;
		
		// output
		OutputScheduleMode mOutputScheduleMode;
		Array<double> mOutputSchedule;
		
		// t-model
		bool mTModelOn;
		
		// clock
		bool mRunClockOn;
		
		// ensemble
		long mPrintRunCountIncrement;
		
		// random number generator
		std::string mGSLRandomNumberGeneratorName;
		unsigned long int mRandomSeed;
		
		// solver type
		std::string mGSLSolverName;
		
		// scalar error controls
		double mLocalRelativeError;
		double mLocalAbsoluteError;
		
		// directories 
		std::string mInputDirectory;
		std::string mOutputDirectory;
		
		// output file streams
		Array<std::ofstream> mOutputStream;
	};



	inline RunControl::RunControl()
	{
		mSystemType = NO_SYSTEM_TYPE;
		mState = NO_SYSTEM_STATE;
		
		mStartTime = 0.0;
		mEndTime = 0.0;
		
		mOutputScheduleMode = NO_OUTPUT_SCHEDULE_MODE;
		
		mTModelOn = false;
		
		mRunClockOn = false;
		mPrintOutputTime = false;
		
		mPrintRunCountIncrement = 0;
		
		mLocalRelativeError = -1.0;
		mLocalAbsoluteError = -1.0;
		
		mGSLSolverName = DEFAULT_GSL_SOLVER;
		
		mGSLRandomNumberGeneratorName = DEFAULT_GSL_RANDOM_NUMBER_GENERATOR;
		mRandomSeed = DEFAULT_RANDOM_SEED;
		
		mInputDirectory = "";
		mOutputDirectory = DEFAULT_OUTPUT_DIRECTORY;
		
		mOutputStream.SetSize(END_OUTPUT_STREAM);
		
		
		return;
	} 
	
	
	
	inline void RunControl::SetState(SystemState state)
	{
		mState = state;
		return;
	}
	
	
	
	inline SystemState RunControl::State() const
	{
		return mState;
	}
	
	
	
	inline SystemType RunControl::GetSystemType() const
	{
		return mSystemType;
	}
	
	
	
	inline void RunControl::SetStartAndEndTimes(double startTime, double endTime)
	{
		if (endTime < startTime)
			ThrowException("RunControl:SetStartAndEndTimes : end time is less than start time");
		
		mStartTime = startTime;
		mEndTime = endTime;
		
		return;
	}



	inline void RunControl::SetLocalRelativeError(double value)
	{
		if (value < 0.0)
			ThrowException("SetRelativeLocalError : attempt to set negative relative error");
		
		mLocalRelativeError = value;
		return;
	}
	
	
	
	inline void RunControl::SetLocalAbsoluteError(double value)
	{
		if (value < 0.0)
			ThrowException("SetRelativeLocalError : attempt to set negative absolute error");
		
		mLocalAbsoluteError = value;
		return;
	}
	
	
	
	inline double RunControl::GetLocalRelativeError() const
	{
		return mLocalRelativeError;
	}
	
	
	
	inline double RunControl::GetLocalAbsoluteError() const
	{
		return mLocalAbsoluteError;
	}
	
	
	
	inline void RunControl::TurnOnTModel()
	{
		mTModelOn = true;
		return;
	}
	
	
	
	inline void RunControl::TurnOffTModel()
	{
		mTModelOn = false;
		return;
	}

	
	
	
	inline bool RunControl::TModelOn() const
	{
		return mTModelOn;
	}
	
	
	
	inline double RunControl::StartTime() const
	{
		return mStartTime;
	}
	
	
	
	inline double RunControl::EndTime() const
	{
		return mEndTime;
	}
	


	inline void RunControl::TurnOnRunClock()
	{
		mRunClockOn = true;
	}
	
	
	
	inline bool RunControl::RunClockOn() const
	{
		return mRunClockOn;
	}
	
	
	
	inline void RunControl::SetGSLRandomNumberGeneratorName(std::string rngName)
	{
		mGSLRandomNumberGeneratorName = rngName;
		return;
	}
	
	
	
	inline std::string RunControl::RandomNumberGeneratorName(void) const
	{
		return mGSLRandomNumberGeneratorName;
	}
	
	
	
	inline void RunControl::SetRandomSeed(unsigned long int seed)
	{
		mRandomSeed = seed;
		return;
	}
	
	
	inline unsigned long int RunControl::RandomSeed() const
	{
		return mRandomSeed;
	}
	
	
	
	inline void RunControl::SetGSLSolverName(std::string solverName)
	{
		mGSLSolverName = solverName;
		return;
	}
	
	
	
	inline std::string RunControl::SolverName() const
	{
		return mGSLSolverName;
	}
	
	
	
	inline long RunControl::NumOutputTimes() const
	{
		return mOutputSchedule.Size();
	}
	
	
	
	inline double RunControl::OutputTime(long n) const
	{
		return mOutputSchedule[n];
	}
	
	
		
	
	inline void RunControl::SetOutputDirectory(std::string directoryName)
	{
		mOutputDirectory = directoryName;
		return;
	}
	
	
	
	inline void RunControl::SetOutputDirectoryToInputDirectory()
	{
		mOutputDirectory = mInputDirectory;
		return;
	}
	
	
	
	inline std::string RunControl::OutputDirectory() const
	{
		return mOutputDirectory;
	}
	
	
	
	inline void RunControl::SetInputDirectory(std::string inputFile)
	{
		mInputDirectory = ExtractDirectoryName(inputFile);
		return;
	}
	
	
	
	inline std::string RunControl::InputDirectory() const
	{
		return mInputDirectory;
	}
	
	
	
	inline void RunControl::TurnOffOutput()
	{
		mOutputSchedule.Erase();
		return;
	}
	
	
	
	inline void RunControl::SetOutputScheduleMode(OutputScheduleMode outputScheduleMode)
	{
		mOutputScheduleMode = outputScheduleMode;
		return;
	}
	
	
	
	inline bool RunControl::PrintOutputTime()
	{
		return mPrintOutputTime;
	}
	
	
	
	inline void RunControl::TurnOnPrintOutputTime()
	{
		mPrintOutputTime = true;
		return;
	}
	
	
	
	inline void RunControl::SetPrintRunCountIncrement(long increment)
	{
		if (increment < 1)
			mPrintRunCountIncrement = 0;
		else
			mPrintRunCountIncrement = increment;
			
		return;
	}
	
	

	inline void RunControl::OpenOutputStream(OutputFileStreamType streamType, std::string fileName)
	{
		OpenOutputFile(mOutputDirectory + fileName, mOutputStream[streamType]);
		return;
	}

	
	
	inline std::ofstream& RunControl::GetOutputStream(OutputFileStreamType streamType)
	{
		return mOutputStream[streamType];
	}
}

#endif // _runcontrol_h_	
