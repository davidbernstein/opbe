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

#ifndef _mode_h_
#define _mode_h_

#include "opbeconst.h"
#include "density.h"
#include "wavenumber.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class Mode {
	 public:
        Mode(void);
		~Mode(void) { };
        
		// copy constructor
		Mode(const Mode &sol);
		
		// value
		void SetValue(double value);
		double Value(void) const;
		
		// initial condition
		void SetInitialCondition(double value);
		void SetToInitialCondition(void);
		
		// wavenumber
		void SetWaveNumbers(short kx, short ky = 0, short kz = 0);
		
		// density
		void SetDensity(const Density *pDensity);
		Density* GetDensity(void) const;
				
	private:
		// member data
	protected:
		// current value
		double mValue;
		
		// initial condition
		double mInitialValue;
		
		// pointer to density
		Density *mpInitialDensity;
		
		// wavenumber(s)
		WaveNumber mWaveNumber;
	};



	inline Mode::Mode()
	{
		mValue = 0.0;
		mInitialValue = 0.0;
		
		mpInitialDensity = NULL;
		
		return;
	} 
	
	
	
	inline void Mode::SetValue(double value)
	{
		mValue = value;
		return;
	}
	
	
	
	inline double Mode::Value() const
	{
		return mValue;
	}
		
	
	
	inline void Mode::SetWaveNumbers(short kx, short ky, short kz)
	{
		mWaveNumber.Set(kx, ky, kz);
		return;
	}
	
			
	
	inline void Mode::SetDensity(const Density *pDensity)
	{
		mpInitialDensity = (Density*) pDensity;
		return;
	}
	
	
	
	inline Density* Mode::GetDensity() const
	{
		return mpInitialDensity;
	}
	
	
	
	inline void Mode::SetInitialCondition(double value)
	{
		mInitialValue = value;
		return;
	}
	
	
	
	inline void Mode::SetToInitialCondition()
	{
		mValue = mInitialValue;
		return;
	}
}

#endif // _mode_h_	
