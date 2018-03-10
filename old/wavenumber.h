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

#ifndef _wavenumber_h_
#define _wavenumber_h_

#include "namespace.h"
#include "array.h"
#include <math.h>

namespace NAMESPACE {
	class WaveNumber {
	 public:
        WaveNumber(void);
		~WaveNumber(void) { };
		
		// wavenumber
		void Set(short kx, short ky = 0, short kz = 0);
		void Set(const Array<double> &a);
		void Get(short &kx, short &ky, short &kz) const;
		short X(void) const;
		short Y(void) const;
		short Z(void) const;
		
		// utility
		double Norm(void) const;
		double SqrNorm(void) const;
		void GetNormalized(Array<double> &normalized) const;
		bool IsZero(void) const;
				
	private:
		// wavenumber(s)
		short mpWaveNumber[3];
	};



	inline WaveNumber::WaveNumber()
	{
		mpWaveNumber[0] = 0;
		mpWaveNumber[1] = 0;
		mpWaveNumber[2] = 0;
				
		return;
	} 
	
	
	
	inline void WaveNumber::Set(short kx, short ky, short kz)
	{
		mpWaveNumber[0] = kx;
		mpWaveNumber[1] = ky;
		mpWaveNumber[2] = kz;
				
		return;
	}
	
	
	
	inline void WaveNumber::Get(short &kx, short &ky, short &kz) const
	{
		kx = mpWaveNumber[0];
		ky = mpWaveNumber[1];
		kz = mpWaveNumber[2];
	}
	
	
	
	inline short WaveNumber::X() const
	{
		return mpWaveNumber[0];
	}
	
	
	
	inline short WaveNumber::Y() const
	{
		return mpWaveNumber[1];
	}
	
	
	
	inline short WaveNumber::Z() const
	{
		return mpWaveNumber[2];
	}
	
	
	
	inline double WaveNumber::Norm() const
	{
		return sqrt(SqrNorm());
	}
	
	
	
	inline double WaveNumber::SqrNorm() const
	{
		return mpWaveNumber[0] * mpWaveNumber[0] + mpWaveNumber[1] * mpWaveNumber[1] + mpWaveNumber[2] * mpWaveNumber[2];
	}
	
	
	
	inline void WaveNumber::GetNormalized(Array<double> &normalized) const
	{
		if (IsZero())
			ThrowException("WaveNumber::GetNormalized : wavenumber is zero");
			
		normalized.SetSize(3);
		
		double factor = 1.0 / Norm();
		
		normalized[0] = mpWaveNumber[0] * factor;
		normalized[1] = mpWaveNumber[1] * factor;
		normalized[2] = mpWaveNumber[2] * factor;
		
		return;
	}
	
	
	
	inline bool WaveNumber::IsZero() const
	{
		return (mpWaveNumber[0] == 0) & (mpWaveNumber[1] == 0) & (mpWaveNumber[2] == 0);
	}
}

#endif // _wavenumber_h_	
