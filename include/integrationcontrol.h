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

#ifndef _integrationcontrol_h_
#define _integrationcontrol_h_

#include "opbeconst.h"
#include "namespace.h"

namespace NAMESPACE {
	class IntegrationControl {
	 public:
        IntegrationControl(void);
		~IntegrationControl(void) { };
		
		// member data
	protected:
		// controls for 1D integration
		double m1DTolerance;
		double mDomainMin
		double mDomainMax;
		short mInitialLevel;
	};



	inline IntegrationControl::IntegrationControl()
	{
		m1DTolerance = DEFAULT_INTEGRATION_TOLERANCE;
		mDomainMin = 0.0;
		mDomainMax = 0.0;
		mInitialLevel = DEFAULT_INTEGRATION_INITIAL_LEVEL;
		
		return;
	} 
}

#endif // _integrationcontrol_h_	
