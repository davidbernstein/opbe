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

#include "modeindex.h"
#include "utility.h"

using namespace NAMESPACE;
using namespace std;
	
	
void ModeIndex::Set(long iSize, long jSize, long kSize)
{
	if (iSize < 1)
		ThrowException("ModeIndex::SetISizeAndJSize : iSize should be positive");
			
	if (jSize < 0)
		ThrowException("ModeIndex::SetISizeAndJSize : jSize should be non-negative");
		
	if (kSize < 0)
		ThrowException("ModeIndex::SetISizeAndJSize : kSize should be non-negative");
			
	mISize = iSize;
	mJSize = jSize;
	mKSize = kSize;
		
	return;
}



long ModeIndex::operator()(const Array<double> &array) const
{
	long i = NearestInteger(array[0]);
	long j = NearestInteger(array[1]);
	long k = NearestInteger(array[2]);
	
	if (i <= 0)
		ThrowException("ModeIndex::operator() : i is non-positive");
		
	j = (j <= 0) ? 1 : j;
	k = (k <= 0) ? 1 : k;
	
	return (*this)(i, j, k);
}



void ModeIndex::Set(const Array<double> &parameter)
{
	Set(NearestInteger(parameter[0]), NearestInteger(parameter[1]), NearestInteger(parameter[2]));
	return;
}