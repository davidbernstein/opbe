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

#ifndef _hermitepolynomial_h_
#define _hermitepolynomial_h_

#include "namespace.h"
#include "utility.h"
#include "opbeconst.h"
#include "density.h"

// these hermite polynomials are orthonormal on (-infinity, +infinity) with respect to the Gaussian weight
// (1 / sqrt(2 pi) * sigma) exp(-(x - x_0)^2 / (2 sigma^2)), where x_0 is mMean below

namespace NAMESPACE {
	class HermitePolynomial {
	 public:
        HermitePolynomial(void);
		~HermitePolynomial(void) { };
		
		// value
		double Evaluate(double x, short n)  const;
		
		// initiliaze
		void Initialize(double mean, double variance, short maxPower);
		void Initialize(const Density &density, short maxPower);
		
		// domain
		Array<double> MakeGrid(double tolerance = 0.0) const;
		
		// inner product
		double InnerProduct(const Array<double> &f, const Array<double> &grid, short n) const;
		
	private:
		void ComputeNormalizationFactors(short numFactors);
		double EvaluateStandard(double x, short n) const;
		double InnerProduct(short n1, short n2, const Array<double> &grid) const;
		double OrthogonalityError(const Array<double> &grid) const;
		void MakeGrid(Array<double> &grid, double xMin, double xMax, short numPoints) const;
		void GetInitialGrid(Array<double> &grid) const;
		void DecreaseGridSpacing(Array<double> &grid) const;
		void IncreaseGridRange(Array<double> &grid) const;
		double DecreaseGridSpacing(Array<double> &grid, double tolerance) const;
		
		// member data
	private:
		// density
		Density mDensity;
		
		// convenient to store these separately
		double mSigma;
		double mMean;
		
		// static data common to all instances
		static Array<double> mNormalizationFactor;
		static bool mNormalizationFactorsReady;
	};



	inline HermitePolynomial::HermitePolynomial()
	{
		mSigma = 0.0;
		mMean = 0.0;
		
		mNormalizationFactorsReady = false;
		
		return;
	} 
}

#endif // _hermitepolynomial_h_	
