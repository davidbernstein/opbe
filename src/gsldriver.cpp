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
#include "modeindex.h"

#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

using namespace NAMESPACE;
using namespace std;

// struct for passing parameters to TimeDerivative
struct gsl_parameters {
	double mEpsilon;
	bool mTModelOn;
	long mNumModes;
	SystemType mSystemType;
};

// global variables for this file
ModeIndex modeIndex;
double *y;

int TimeDerivative(double t, const double u[], double uDot[], void *params);
int BurgersEquation(double t, const double u[], double uDot[], gsl_parameters *pParams);
int NavierStokes(double t, const double u[], double uDot[], gsl_parameters *pParams);

inline double U(long i, const double u[]) {return u[modeIndex(i)];}
inline double U(long i, long j, long k, const double u[]) {return u[modeIndex(i, j, k)];}

void System::GSLEvolve(double t1)
{
	// note: static variables here are shared among all instances of System
	static gsl_parameters params;
	static const gsl_odeiv_step_type *pStepType = NULL;
	static gsl_odeiv_step *pStepAlloc = NULL;
	static gsl_odeiv_control *pStepControl = NULL;
	static gsl_odeiv_evolve *pEvolve = NULL;
	static Array<double> yArray;
	static double h;
		
	// initialize if first call
	if (mpRunControl->State() == SYSTEM_INITIALIZE) {
		params.mEpsilon = mpOPBEParameter->ViscosityCoefficient();
		params.mTModelOn = mpRunControl->TModelOn();
		params.mNumModes = mMode.Size();
		params.mSystemType = mpRunControl->GetSystemType();
		
		// solver
		string solverName = mpRunControl->SolverName();
		if (solverName == "rk2")
			pStepType = gsl_odeiv_step_rk2;
			
		if (solverName == "rk4")
			pStepType = gsl_odeiv_step_rk4;
		
		if (solverName == "rkf45")
			pStepType = gsl_odeiv_step_rkf45;
		
		if (solverName == "rkck")
			pStepType = gsl_odeiv_step_rkck;
		
		if (solverName == "rk8pd")
			pStepType = gsl_odeiv_step_rk8pd;
		
		if (solverName == "rk2imp")
			pStepType = gsl_odeiv_step_rk2imp;
		
		if (solverName == "rk4imp")
			pStepType = gsl_odeiv_step_rk4imp;
		
		if (solverName == "bsimp")
			pStepType = gsl_odeiv_step_bsimp;
		
		if (solverName == "gear1")
			pStepType = gsl_odeiv_step_gear1;
			
		if (solverName == "gear2")
			pStepType = gsl_odeiv_step_gear2;
			
		if (pStepType == NULL)
			ThrowException("System::GSLEvolve : invalid gsl solver name: " + solverName);
		else
			pStepAlloc = gsl_odeiv_step_alloc(pStepType, mMode.Size());
	
		// step control
		double localAbsoluteError = mpRunControl->GetLocalAbsoluteError();
		double localRelativeError = mpRunControl->GetLocalRelativeError();
		pStepControl = gsl_odeiv_control_y_new(localAbsoluteError, localRelativeError);
	
		// evolver
		pEvolve = gsl_odeiv_evolve_alloc(mMode.Size());
		
		h = DEFAULT_TIME_STEP;
		
		modeIndex = *mpModeIndex;
					
		return;
	}
	
	// clean up if done
	if (mpRunControl->State() == SYSTEM_STOP) {
		pStepType = NULL;
		
		gsl_odeiv_evolve_free(pEvolve);
		pEvolve = NULL;
		
		gsl_odeiv_control_free(pStepControl);
		pStepControl = NULL;
		
		gsl_odeiv_step_free(pStepAlloc);
		pStepAlloc = NULL;
		
		return;
	}
	
	//gsl_odeiv_system system = {TimeDerivative, Jacobian, mNumModes, &params};
	gsl_odeiv_system system = {TimeDerivative, NULL, params.mNumModes, &params};
		
	double t = mCurrentTime;
	//double *y = mMode.Begin();
	y = mMode.Begin();
	
	// call ode solver
	while (t < t1) {
		int status = gsl_odeiv_evolve_apply(pEvolve, pStepControl, pStepAlloc, &system, &t, t1, &h, y);
	
		if (status != GSL_SUCCESS) 
			ThrowException("System::GSLEvolve : gsl step unsuccessful, gsl_status = " + status);
	}
	
	// update current time
	mCurrentTime = t;
	
	return;
}

	
	
int TimeDerivative(double t, const double u[], double uDot[], void *params)
{
	gsl_parameters *pParams = (gsl_parameters *) params;
	
	switch (pParams->mSystemType) {
	case BURGERS_EQUATION:
		return BurgersEquation(t, u, uDot, pParams);
		break;
	
	case NAVIER_STOKES:
		return NavierStokes(t, u, uDot, pParams);
		break;
	
	default:
		ThrowException("TimeDerivative : System type not set");
		return GSL_FAILURE;
		break;
	}
}




int BurgersEquation(double t, const double u[], double uDot[], gsl_parameters *pParams)
{
	static Array<double> work;
	
	double epsilon = pParams->mEpsilon;
	long numModes = pParams->mNumModes;
	bool tModelOn = pParams->mTModelOn;
			
	double sum1 = 0.0, sum2 = 0.0;

	for (long k = 1; k <= numModes; ++k) {
		double viscosityTerm = -epsilon * k * k * U(k, u);

		sum1 = 0.0;
		for (long kp = 1; kp <= numModes - k; ++kp) 
			sum1 += U(kp, u) * U(kp + k, u);
		
		sum2 = 0.0;
		for (long kp = 1; kp <= k - 1; ++kp)
			sum2 += kp * U(kp, u) * U(k - kp, u);
		
		uDot[modeIndex(k)] = viscosityTerm + 0.5 * (k * sum1 - sum2);
	}
	
	if (tModelOn) {
		if (work.Empty())
			work.SetSize(numModes + 1);
		
		for (long mpp = 1; mpp <= numModes; ++mpp) {
			work[mpp] = 0.0;
			for (long mp = mpp; mp <= numModes; ++mp) 
				work[mpp] += (mpp + numModes - mp) * U(mp, u) * U(mpp + numModes - mp, u);
		}
		
		for (long m = 1; m <= numModes; ++m) {
			sum1 = 0.0;
			for (long mpp = 1; mpp <= m; ++mpp) 
				sum1 += work[mpp] * U(mpp + numModes - m, u);
			
			uDot[modeIndex(m)] += -0.25 * t * m * sum1;
		}
	}
	
	
	return GSL_SUCCESS;
}



int NavierStokes(double t, const double u[], double uDot[], gsl_parameters *pParams)
{


	return GSL_SUCCESS;
}



void System::RHS(Array<double> &rhs) const
{
	gsl_parameters params;
	params.mEpsilon = mpOPBEParameter->ViscosityCoefficient();
	params.mTModelOn = mpRunControl->TModelOn();
	params.mNumModes = mMode.Size();
	params.mSystemType = mpRunControl->GetSystemType();
		
	rhs.SetSize(mMode.Size());
	
	TimeDerivative(mCurrentTime, mMode.Begin(), rhs.Begin(), &params);
	
	return;
}



void System::RatioTModel(Array<double> &ratio) const
{
	gsl_parameters params;
	params.mEpsilon = mpOPBEParameter->ViscosityCoefficient();
	params.mTModelOn = mpRunControl->TModelOn();
	params.mNumModes = mMode.Size();
	params.mSystemType = mpRunControl->GetSystemType();
		
	ratio.SetSize(mMode.Size());
	Array<double> rhsOff(mMode.Size());
	Array<double> rhsOn(mMode.Size());
	
	params.mTModelOn = false;
	TimeDerivative(mCurrentTime, mMode.Begin(), rhsOff.Begin(), &params);
	
	params.mTModelOn = true;
	TimeDerivative(mCurrentTime, mMode.Begin(), rhsOn.Begin(), &params);
	
	for (long k = 0; k < mMode.Size(); ++k) {
		ratio[k] = rhsOn[k];
	/*
		if (rhsOff[k] != 0.0)
			ratio[k] = (rhsOn[k] - rhsOff[k]) / rhsOff[k];
		else
			ratio[k] = 0.0;
		*/
	}
		
		
	return;
}

