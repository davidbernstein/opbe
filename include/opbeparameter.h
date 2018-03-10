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

#ifndef _opbeparameter_h_
#define _opbeparameter_h_

#include "array.h"
#include "namespace.h"
#include "utility.h"
#include "constants.h"

namespace NAMESPACE {
	class OPBEParameter {
	 public:
        OPBEParameter(void);
		~OPBEParameter(void) { };
		
		// Reynolds number and viscosity coefficient
		void SetReynoldsNumber(double re);
		void SetViscosityCoefficient(double epsilon);
		double ReynoldsNumber(void) const;
		double ViscosityCoefficient(void) const;
		
		// t-model parameters
		void SetTModelNumUnresolved(long numUnresolved);
		void SetZeroSigmaLimitFlagOn(void);
		void SetAllSigmasEqual(double sigma);
		void SetSigma(const Array<double> &sigma);

		// member data
	protected:
		// viscosity coefficient
		double mViscosityCoefficient;
		
		// t-model parameters
		long mTModelNumUnresolved;
		bool mZeroSigmaLimitOn;
		Array<double> mSigma;
		bool mAllSigmasEqual;
	};



	inline OPBEParameter::OPBEParameter()
	{
		mViscosityCoefficient = -1.0;
				
		mTModelNumUnresolved = -1;
		mZeroSigmaLimitOn = false;
		mAllSigmasEqual = false;
		
		
		return;
	} 
	
	
	
	inline double OPBEParameter::ViscosityCoefficient() const
	{
		return mViscosityCoefficient;
	}
	
	
		
	inline double OPBEParameter::ReynoldsNumber() const
	{
		// should really return a flag for infinite Reynolds number
		if (mViscosityCoefficient > 0.0)
			return PI / mViscosityCoefficient;
		else
			return -1.0;
	}
	
	
	
	inline void OPBEParameter::SetViscosityCoefficient(double epsilon)
	{
		if (epsilon < 0.0)
			ThrowException("Problem:SetViscosityCoefficient : negative viscosity coefficient");
		
		mViscosityCoefficient = epsilon;
			
		return;
	}
			
	
	
	inline void OPBEParameter::SetTModelNumUnresolved(long numUnresolved)
	{
		if (numUnresolved < 1)
			ThrowException("");
		
		mTModelNumUnresolved = numUnresolved;
		return;
	}
	
	
	
	inline void OPBEParameter::SetZeroSigmaLimitFlagOn(void)
	{
		mZeroSigmaLimitOn = true;
		return;
	}
	
	
	
	inline void OPBEParameter::SetAllSigmasEqual(double sigma)
	{
		if (sigma <= 0.0)
			ThrowException("SetAllSigmasEqual : sigma must be positive");
			
		mAllSigmasEqual = true;
		mSigma.SetSize(1);
		mSigma[0] = sigma;
		
		return;
	}
	
	
	
	inline void OPBEParameter::SetSigma(const Array<double> &sigma)
	{
		if (sigma.Empty())
			ThrowException("OPBEParameter::SetSigma : array empty");
		
		mSigma = sigma;
		
		// check
		for (long i = 0; i < mSigma.Size(); ++i) {
			if (mSigma[i] <= 0.0)
				ThrowException("OPBEParameter::SetSigma : one or more values non-positive");
		}
		
		return;
	}
}

#endif // _opbeparameter_h_	
