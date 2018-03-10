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

#ifndef _opbeenums_h_
#define _opbeenums_h_

namespace OPBurgersEquation {
	enum SystemState{NO_SYSTEM_STATE, SYSTEM_INITIALIZE, SYSTEM_RUN, SYSTEM_STOP};
	
	enum ProblemState{NO_PROBLEM_STATE, PROBLEM_START, PROBLEM_RUNNING, PROBLEM_DONE};
	
	enum DensityType{NO_DENSITY_TYPE, ONE_DIMENSIONAL_GAUSSIAN};
	
	enum SystemType{NO_SYSTEM_TYPE, BURGERS_EQUATION, NAVIER_STOKES};
	
	enum ModeType{RESOLVED_MODE, UNRESOLVED_MODE};
	
	enum ProblemType{NO_PROBLEM_TYPE, AVERAGING_PROBLEM, MK_PROBLEM, DELTA_PROBLEM, FIXED_IC_PROBLEM};
	
	enum QuadratureType{NO_QUADRATURE_TYPE, 
						FIXED_SPACING_QUADRATURE, 
						ADAPTIVE_QUADRATURE, 
						MONTE_CARLO_QUADRATURE};
	
	enum OutputScheduleMode{NO_OUTPUT_SCHEDULE_MODE, 
							OUTPUT_SCHEDULE_LINEAR, 
							OUTPUT_SCHEDULE_LOGARITHMIC};
							
	enum OutputFileStreamType{MODE_OUTPUT_STREAM, 
							  ENERGY_OUTPUT_STREAM, 
							  MOMENTS_OUTPUT_STREAM,
							  TMODEL_RATIO_OUTPUT_STREAM,
							  VOLTERRA_F0_OUTPUT_STREAM,
							  END_OUTPUT_STREAM};
}

#endif // _opbeenums_h_
