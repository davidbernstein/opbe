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

#ifndef _deltaproblem_h_
#define _deltaproblem_h_

#include "problem.h"
#include "realmatrix.h"

#include <string>
#include <iostream>

namespace NAMESPACE {
	class DeltaProblem : public Problem {
	 public:
        DeltaProblem(void);
		~DeltaProblem(void) { };
        
		// copy constructor
		DeltaProblem(const DeltaProblem &sol);
		
        // run
        void Run(const std::string &fileName);
		
		// test
		void Test(const std::string &fileName);
						
	private:
		// run
		void Run(void);
		
		// input
		void ReadInputFile(const std::string &fileName);
				
		// Initialization
		void Initialize(void);
		
		// Volterra
		double VolterraF0(short modeIndex) const;
					
		// IO 
		void WriteVolterraFFile(void);

		// member data
	private:
		double mDeltaX1;
		double mDeltaX2;
	};



	inline DeltaProblem::DeltaProblem()
	{
		mDeltaX1 = -1.0;
		mDeltaX2 = -1.0;
		
		return;
	} 
}

#endif // _deltaproblem_h_	
