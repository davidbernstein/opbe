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

#ifndef _modeindex_h_
#define _modeindex_h_

#include "array.h"

namespace NAMESPACE {
	class ModeIndex {
	 public:
        ModeIndex(void);
		~ModeIndex(void) { };
        
		// copy constructor
		ModeIndex(const ModeIndex &sol);
		
		// Set
		void Set(long iSize, long jSize, long kSize);
		void Set(const Array<double> &parameter);
		
		// number of modes
		void SetNumResolvedAndUnresolvedModes(long numResolved, long numUnresolved);
		long NumResolvedModes(void) const;
		long NumUnresolvedModes(void) const;
		long NumModes(void) const;
		bool InResolvedRange(long i) const;
		
		// indexing
		long operator()(long i) const;
		long operator()(long i, long j, long k) const;
		long operator()(const Array<double> &array) const;
		
		// max
		long Max(void);
		
		// member data
	protected:
		// modes
		long mNumResolvedModes;
		long mNumUnresolvedModes;
		
		// index ranges
		long mISize;
		long mJSize;
		long mKSize;
	};



	inline ModeIndex::ModeIndex()
	{
		mNumResolvedModes = 0;
		mNumUnresolvedModes = 0;
		
		mISize = 0;
		mJSize = 0;
		mKSize = 0;
		
		return;
	} 
	
	
    
    inline long ModeIndex::operator()(long i) const
	{
        return i - 1;
	}
	
	
	
	inline long ModeIndex::operator()(long i, long j, long k) const
	{
		--i;
		--j;
		--k;
		
		return i + mISize * (j + k * mJSize);
	}
	
	
	
	inline long ModeIndex::Max()
	{
		return mISize * (1 + mJSize * (1 + mKSize));
	}
	
	
	
	inline void ModeIndex::SetNumResolvedAndUnresolvedModes(long numResolved, long numUnresolved)
	{
		if (numResolved <= 0)
			ThrowException("OPBEParameter::SetNumResolvedAndUnresolvedModes : number of resolved modes must be positive");
			
		if (numUnresolved < 0)
			ThrowException("OPBEParameter::SetNumResolvedAndUnresolvedModes : number of unresolved modes is negative");
			
		mNumResolvedModes = numResolved;
		mNumUnresolvedModes = numUnresolved;
		
		return;
	}
	
	
	
	inline long ModeIndex::NumResolvedModes() const
	{
		return mNumResolvedModes;
	}
	
	
	
	inline long ModeIndex::NumUnresolvedModes() const
	{
		return mNumUnresolvedModes;
	}
	
	
	
	inline long ModeIndex::NumModes() const
	{
		return mNumResolvedModes + mNumUnresolvedModes;
	}
	
	
	
	inline bool ModeIndex::InResolvedRange(long i) const
	{
		return i <= mNumResolvedModes;
	}
}

#endif // _modeindex_h_	
