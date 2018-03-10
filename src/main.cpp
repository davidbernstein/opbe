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

#include "problem.h"
#include "mkproblem.h"
#include "averagingproblem.h"
#include "fixedicproblem.h"
#include "deltaproblem.h"

#include <iostream>

using namespace std; 
using namespace NAMESPACE;

int main()
{	
	MKProblem mk;			
	AveragingProblem a;
	FixedICProblem f;
	DeltaProblem delta;

    try { 	
		string fileName = "/Users/dave/Documents/OP_Burgers/OP_Burgers3/Figures/OvershootFig/overshoot.txt";
        
		//string fileName = "/Users/dave/Projects/OPBE/SimulationFiles/fixedICBurgers.txt";
		Problem p;
		
		switch (p.GetProblemType(fileName)) {
		case MK_PROBLEM: 
			mk.Run(fileName);
			break;
					
		case AVERAGING_PROBLEM: 
			a.Run(fileName);
			break;
			
		case FIXED_IC_PROBLEM:
			f.Run(fileName);
			break;
			break;
			
		case DELTA_PROBLEM:
			delta.Run(fileName);
			break;

		default:
			ThrowException("bad problem type");
			break;
		}
    }
    catch (exception &standardException) {
        HandleException(standardException);
    }
    catch (string &message) {
        HandleException(message);
    }
	

    return 0;
}

