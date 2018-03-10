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

#ifndef _opbeconst_h_
#define _opbeconst_h_

#include "constants.h"

namespace NAMESPACE {
	// quadrature constants
	const double DEFAULT_QUADRATURE_TOLERANCE = 1.0e-4;
	const short DEFAULT_INTEGRATION_INITIAL_LEVEL = 4;
	const double DEFAULT_HERMITE_INNER_PRODUCT_TOLERANCE = 1.0E-3;
	const short DEFAULT_INITIAL_QUADRATURE_LEVEL = 4;
	const short DEFAULT_QUADRATURE_NUM_POINTS = 5;
	const double HERMITE_INTEGRATION_MULTIPLIER = 1.1;
	
	// gsl
	const std::string DEFAULT_GSL_SOLVER = "rk8pd";
	const double DEFAULT_TIME_STEP = 1.0e-3;
	const double DEFAULT_LOCAL_RELATIVE_ERROR = 0.0;
	const double DEFAULT_LOCAL_ABSOLUTE_ERROR = 1.0e-6;
	const std::string DEFAULT_GSL_RANDOM_NUMBER_GENERATOR = "taus";
	const unsigned long int DEFAULT_RANDOM_SEED = 0;
	
	// numerical constants
	const double ONE_OVER_SQRT_2PI = 1.0 / sqrt(2.0 * PI);
	const double ONE_OVER_2PI = 1.0 / (2.0 * PI);
	const double SQRT_TWO_OVER_PI = sqrt(2.0 / PI);
	
	// default output directory is root
	const std::string DEFAULT_OUTPUT_DIRECTORY = "/";
	
	// output 
	const double DEFAULT_FIRST_OUTPUT_TIME = 0.1;
	const short DEFAULT_NUM_MOMENTS = 10;
	
	// hermite polynomials
	const short DEFAULT_FINITE_RANK_SIZE = 10;
	const double HERMITE_DOMAIN_STEP = 0.1;
}

#endif // _opbeconst_h_
