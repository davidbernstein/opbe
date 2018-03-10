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

#include <math.h>
#include <list>

#include "density.h"
#include <iostream>

using namespace NAMESPACE;
using namespace std;

double Density::Value(double x) const
{
	switch (mType) {
	case NO_DENSITY_TYPE:
		ThrowException("Density::Value : unspecified function type");
		break;
			
	case ONE_DIMENSIONAL_GAUSSIAN:
		return OneDGaussian(x);
		break;
        
	default:
		ThrowException("Density::Value : unspecified function type");
		break;
	}
	
	return -1.0;
}



double Density::OneDGaussian(double x) const
{
	// Gaussian is of the form
	// 
	// (2 Pi s^2)^{-1/2} exp(-(x-m)/ (2 s^2)) 
	// m   = mean
	// s   = standard deviation
		
	// mParameter[0] : m
	// mParameter[1] : s
	// mParameter[2] : mode index
	
	if (mParameter[1] <= 0.0)
		ThrowException("Density::Gaussian : Variance is non-positive");
		
	double arg = (x - mParameter[0]) / mParameter[1];
	double expTerm = exp(-0.5 * arg * arg);
	
	return ONE_OVER_SQRT_2PI * expTerm / mParameter[1];
}



void Density::DomainBounds(double epsilon, double &min, double &max) const
{
	switch (mType) {
	case NO_DENSITY_TYPE:
		ThrowException("Density::Value : unspecified function type");
		break;
			
	case ONE_DIMENSIONAL_GAUSSIAN:
		return OneDGaussianDomain(epsilon, min, max);
		break;
        
	default:
		ThrowException("Density::DomainBounds : unspecified function type");
		break;
	}
	
	return;
}



void Density::OneDGaussianDomain(double epsilon, double &min, double &max) const
{
	double term = log(epsilon * mParameter[1] / ONE_OVER_SQRT_2PI);
	if (term >= 0.0)
		ThrowException("Density::GaussianDomain : weight minimum not small enough");
		
	term = mParameter[1] * sqrt(-2.0 * term);
	
	min = mParameter[0] - term;
	max = mParameter[0] + term;
	
	return;
}



short Density::GetModeIndex() const
{
	switch (mType) {
	case NO_DENSITY_TYPE:
		ThrowException("Density::GetModeIndex : unspecified function type");
		break;
			
	case ONE_DIMENSIONAL_GAUSSIAN:
		return GetOneDGaussianModeIndex();
		break;
        
	default:
		ThrowException("Density::DomainBounds : unspecified function type");
		break;
	}
	
	return -1;
}



void Density::SetModeIndex(short index) 
{
	switch (mType) {
	case NO_DENSITY_TYPE:
		ThrowException("Density::SetModeIndex : unspecified function type");
		break;
			
	case ONE_DIMENSIONAL_GAUSSIAN:
		return SetOneDGaussianModeIndex(index);
		break;
        
	default:
		ThrowException("Density::DomainBounds : unspecified function type");
		break;
	}
	
	return;
}



short Density::GetOneDGaussianModeIndex() const
{
	return (short) NearestInteger(mParameter[2]);
}



void Density::SetOneDGaussianModeIndex(short index) 
{
	if (index < 0)
		ThrowException("Density::SetOneDGaussianModeIndex : negative index");
		
	mParameter[2] = index;
	
	return;
}



void Density::SetParameters(const Array<double> &input)
{
	if (input.Size() != NumParametersToSpecify())
		ThrowException("Density::SetParameters : input is wrong size");
		
	mParameter = input;
	
	return;
}



void Density::SetType(const string &name)
{
	SetType(GetDensityType(name));
    return;
}



DensityType Density::GetDensityType(const string &name)
{
	string tmp = name;
	LowerCase(tmp);
    
	if (tmp == "1dgaussian")
		return ONE_DIMENSIONAL_GAUSSIAN;
    
	ThrowException("GetDensityType : bad input string : " + name);
    
	return NO_DENSITY_TYPE;
}



void Density::Print() const
{
	cout << mType << " ";
	
	for (short i = 0; i < mParameter.Size(); ++i)
		cout << mParameter[i] << " ";
	
	cout << endl;
	
	return;
}