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

#ifndef _density_h_
#define _density_h_

#include "opbeenums.h"
#include "opbeconst.h"

namespace NAMESPACE {
	class Density {
	public:
		// Constructor
		Density(DensityType type = NO_DENSITY_TYPE);

		// Destructor
		~Density(void) { };

		// value
		double Value(double x) const;
		
		// mode index
		short GetModeIndex(void) const;
		void SetModeIndex(short index);
		
		// domain
		void DomainBounds(double epsilon, double &min, double &max) const;
		void OneDGaussianDomain(double epsilon, double &min, double &max) const;
		
		// type
		void SetType(DensityType type);
		void SetType(const std::string &name);
        DensityType Type(void) const;
        DensityType GetDensityType(const std::string &name);
		
		// parameters
		short NumParametersToSpecify(void);
		void SetParameters(const Array<double> &input);
		double GetParameter(short index) const;
		
		// I/O
		void Print(void) const;
		
	private:
		double OneDGaussian(double x) const;
		short GetOneDGaussianModeIndex(void) const;
		void SetOneDGaussianModeIndex(short index);
		
	private:
		// type of function
		DensityType mType;
		
		// parameters
		Array<double> mParameter;
	};
	
	
	
	inline Density::Density(DensityType type)
	{
		SetType(type);
		return;
	}
	
	
	
	inline void Density::SetType(DensityType type)
	{
		if (type == NO_DENSITY_TYPE) {
			mParameter.SetSize(0);
			mType = NO_DENSITY_TYPE;
			return;
		}
				
		if (type == ONE_DIMENSIONAL_GAUSSIAN) {
			mParameter.SetSize(3);
			mType = ONE_DIMENSIONAL_GAUSSIAN;
			return;
		}
		
		ThrowException("Density::SetType : bad type");
		
		
		return;
	}
	
	
	
	inline short Density::NumParametersToSpecify()
	{
		switch (mType) {
		case NO_DENSITY_TYPE: 
			return 0;
			break;
				
		case ONE_DIMENSIONAL_GAUSSIAN:
			return 3;
			break;
		
		default:
			ThrowException("Density::NumParameters : bad type");
			return -1;
		}
	}

    
    
    inline DensityType Density::Type() const
    {
        return mType;
    }
	
	
	
	inline double Density::GetParameter(short index) const
	{
		return mParameter[index];
	}
}


#endif // _density_h_

