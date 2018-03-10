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

#include "hermitepolynomial.h"
#include "constants.h"

#include <iostream>
#include <iomanip>

using namespace NAMESPACE;
using namespace std;

// initialize static member data
Array<double> HermitePolynomial::mNormalizationFactor(0);
bool HermitePolynomial::mNormalizationFactorsReady = false;

double HermitePolynomial::Evaluate(double x, short n) const
{
	if (mDensity.Type() == NO_DENSITY_TYPE)
		ThrowException("HermitePolynomial::Evaluate : not initialized");
		
	double y = (x - mMean) / (SQRT_TWO * mSigma);
	return EvaluateStandard(y, n) * mNormalizationFactor[n];
}



double HermitePolynomial::EvaluateStandard(double x, short n) const
{
	// the standard hermite polynomials are the ones normal w.r.t exp(-x^2),
	// i.e., 1, 2x, 4x^2 - 2, 8x^3 - 12x, etc.
	
	if (n < 0)
		ThrowException("HermitePolynomial::Evaluate : n < 0");
		
	if (n == 0)
		return 1.0;
	
	if (n == 1)
		return 2.0 * x;
	
	double him1 = 2.0 * x;
	double him2 = 1.0;
	double hi;
	for (short i = 2; i <= n; ++i) {
		hi = 2.0 * x * him1 - 2.0 * (i - 1.0) * him2;
		him2 = him1;
		him1 = hi;
	}
	
	return hi;
}



void HermitePolynomial::ComputeNormalizationFactors(short numFactors)
{	
	if (mNormalizationFactorsReady)
		return;
		
	if (numFactors <= 0)
		ThrowException("HermitePolynomial::ComputeNormalizationFactors : number of factors non-positive");
	
	mNormalizationFactor.SetSize(numFactors);
	
	mNormalizationFactor[0] = 1.0;
	
	for (short i = 1; i < numFactors; ++i) {
		double factor = 1.0;
		for (short n = i; n > 0; --n) {
			factor *= sqrt(2.0 * n);
		}
		mNormalizationFactor[i] = 1.0 / factor;
		//cout << i << " " << mNormalizationFactor[i] << endl;
	}
	
	return;
}



void HermitePolynomial::Initialize(double mean, double variance, short maxPower)
{
	mDensity.SetType(ONE_DIMENSIONAL_GAUSSIAN);
	
	mMean = mean;
	mSigma = variance;
	
	Array<double> parameter(3);
	parameter[0] = mean;
	parameter[1] = variance;
	parameter[2] = 0.0;
	
	mDensity.SetParameters(parameter);
	
	ComputeNormalizationFactors(maxPower);
	
	return;
}



void HermitePolynomial::Initialize(const Density &density, short maxPower)
{
	mDensity = density;
	
	mMean = mDensity.GetParameter(0);
	mSigma = mDensity.GetParameter(1);
	
	ComputeNormalizationFactors(maxPower);
	
	return;
}



Array<double> HermitePolynomial::MakeGrid(double tolerance) const
{
	Array<double> grid;
	
	if (tolerance <= 0.0) {
		cout << "HermitePolynomial::ComputeGrid : tolerance is not positive, using default value of ";
		cout << DEFAULT_HERMITE_INNER_PRODUCT_TOLERANCE << endl;
		tolerance = DEFAULT_HERMITE_INNER_PRODUCT_TOLERANCE;
	}
	
	GetInitialGrid(grid);
	
	double error = DecreaseGridSpacing(grid, tolerance);
	while (error > tolerance) {
		IncreaseGridRange(grid); 
		//cout << grid[0] << " " << grid[grid.Size() - 1] << " " << grid.Size() << endl;
		error = DecreaseGridSpacing(grid, tolerance);
	}
	
	return grid;
}



double HermitePolynomial::DecreaseGridSpacing(Array<double> &grid, double tolerance) const
{
	double errorPrev = OrthogonalityError(grid);
	DecreaseGridSpacing(grid);
	double errorNew = OrthogonalityError(grid);
	double fractionalChange = fabs(errorNew - errorPrev) / errorPrev;
	
	while (fractionalChange > tolerance) {
		errorPrev = errorNew;
		DecreaseGridSpacing(grid);
		errorNew = OrthogonalityError(grid);
		fractionalChange = fabs(errorNew - errorPrev) / errorPrev;
	}


	return errorNew;
}



double HermitePolynomial::OrthogonalityError(const Array<double> &grid) const
{
	double error, errorMax = -1.0;
	
	short nMax = mNormalizationFactor.Size();
	
	for (short i = 0; i < nMax; ++i) {
		for (short j = 0; j <= i; ++j) {
			double innerProduct = InnerProduct(i, j, grid);
			error = (i == j) ? fabs(1.0 - innerProduct) : fabs(innerProduct);
			errorMax = max(error, errorMax);
		}
	}

	return errorMax;
}



double HermitePolynomial::InnerProduct(short n1, short n2, const Array<double> &grid) const
{
	short i = 0;
	double x = grid[i];
	double sum = 0.5 * Evaluate(x, n1) * Evaluate(x, n2) * mDensity.Value(x);
	
	for (i = 1; i < grid.Size() - 1; ++i) {
		double x = grid[i];
		sum += Evaluate(x, n1) * Evaluate(x, n2) * mDensity.Value(x);
	}
	
	x = grid[i];
	sum += 0.5 * Evaluate(x, n1) * Evaluate(x, n2) * mDensity.Value(x);
	

	return sum * (grid[1] - grid[0]);
}



double HermitePolynomial::InnerProduct(const Array<double> &f, const Array<double> &grid, short n) const
{
	// the grid is assume to be evenly spaced and f is assumed to be evaluated at the points grid[i]
	if (f.Size() != grid.Size())
		ThrowException("HermitePolynomial::InnerProduct : f and grid have different sizes");
		
	short i = 0;
	double x = grid[i];
	double sum = 0.5 * Evaluate(x, n) * f[i] * mDensity.Value(x);
	
	for (i = 1; i < grid.Size() - 1; ++i) {
		x = grid[i];
		sum += Evaluate(x, n) * f[i] * mDensity.Value(x);
	}
	
	i = grid.Size() - 1;
	x = grid[i];
	sum += 0.5 * Evaluate(x, n) * f[i] * mDensity.Value(x);
	

	return sum * (grid[1] - grid[0]);

}



void HermitePolynomial::MakeGrid(Array<double> &grid, double xMin, double xMax, short numPoints) const
{
	grid.SetSize(numPoints);
	double dX = (xMax - xMin) / (numPoints - 1);
	
	for (short i = 0; i < numPoints; ++i) 
		grid[i] = xMin + i * dX;
	
	return;
}



void HermitePolynomial::GetInitialGrid(Array<double> &grid) const
{
	MakeGrid(grid, mMean - mSigma, mMean + mSigma, DEFAULT_QUADRATURE_NUM_POINTS);
	
	return;
}



void HermitePolynomial::DecreaseGridSpacing(Array<double> &grid) const
{
	short numPoints = NearestInteger(HERMITE_INTEGRATION_MULTIPLIER * grid.Size());
	
	MakeGrid(grid, grid[0], grid[grid.Size() - 1], numPoints);

	return;
}



void HermitePolynomial::IncreaseGridRange(Array<double> &grid) const
{
	double xMin = grid[0];
	double xMax = grid[grid.Size() - 1];
	
	double range = xMax - xMin;
	double center = 0.5 * (xMax + xMin);
	
	double a = center - 0.5 * range * HERMITE_INTEGRATION_MULTIPLIER;
	double b = center + 0.5 * range * HERMITE_INTEGRATION_MULTIPLIER;
	
	MakeGrid(grid, a, b, grid.Size());
	
	
	return;
}