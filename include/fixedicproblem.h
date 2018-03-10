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

#ifndef _fixedicproblem_h_
#define _fixedicproblem_h_

#include "problem.h"

namespace NAMESPACE {
	class FixedICProblem : public Problem {
	 public:
        FixedICProblem(void) { };
		~FixedICProblem(void) { };
        
		// copy constructor
		FixedICProblem(const FixedICProblem &sol);
		
        // run
        void Run(const std::string &fileName);
				
	private:
		// input
		void ReadInputFile(const std::string &fileName);
		
		// initialization
		void Initialize(void);
	};
}

#endif // _fixedicproblem_h_	
