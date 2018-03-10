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

#include "wavenumber.h"
#include "utility.h"

using namespace NAMESPACE;
using namespace std;

namespace NAMESPACE {

void WaveNumber::Set(const Array<double> &a)
{
	if (a.Empty()) {
		mpWaveNumber[0] = 0;
		mpWaveNumber[1] = 0;
		mpWaveNumber[2] = 0;
		return;
	}
	
	if (a.Size() == 1) {
		mpWaveNumber[0] = NearestInteger(a[0]);
		mpWaveNumber[1] = 0;
		mpWaveNumber[2] = 0;
		return;
	}
	
	if (a.Size() == 2) {
		mpWaveNumber[0] = NearestInteger(a[0]);
		mpWaveNumber[1] = NearestInteger(a[1]);
		mpWaveNumber[2] = 0;
		return;
	}
	
	mpWaveNumber[0] = NearestInteger(a[0]);
	mpWaveNumber[1] = NearestInteger(a[1]);
	mpWaveNumber[2] = NearestInteger(a[2]);
	return;
}
}