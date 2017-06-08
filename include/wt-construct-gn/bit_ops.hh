/*
 * Copyright (c) 2017 Tuukka Norri
 * This file is part of wt-construct-gn.
 * 
 * wt-construct-gn is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * wt-construct-gn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with wt-construct-gn.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef WT_CONSTRUCT_GN_BIT_OPS_HH
#define WT_CONSTRUCT_GN_BIT_OPS_HH

#include <climits>
#include <cstdint>
#include <sdsl/int_vector.hpp>


namespace wtcgn {
	
	enum {
		UINT32_T_BIT	= CHAR_BIT * sizeof(uint32_t),
		BIT_VECTOR_BIT	= CHAR_BIT * sizeof(sdsl::std_size_type_for_int_vector)
	};
	
	
	// t_type must be unsigned.
	template <typename t_type>
	inline t_type rol(t_type const word, t_type const n)
	{
		return ((word << n) | (word >> (CHAR_BIT * sizeof(t_type) - n)));
	}
	
	
#ifndef HAVE_PEXT
	template <typename t_type>
	inline t_type extract_bits(t_type const word, t_type mask)
	{
		// Adapted from https://github.com/fmatthew5876/stdcxx-bitops/
		t_type res(0);
		for (t_type bb(1); mask != 0; bb += bb)
		{
			if (word & mask & -mask)
				res |= bb;
			
			mask &= (mask - 1);
		}
		return res;
	}
	
#else
	inline uint32_t extract_bits(uint32_t const word, uint32_t const mask)
	{
		uint32_t dst{0};
		asm (
			"pextl %2, %1, %0"
			: "=r" (dst)
			: "r" (word), "r" (mask)
		);
		return dst;
	}


	inline uint64_t extract_bits(uint64_t const word, uint64_t const mask)
	{
		uint64_t dst{0};
		asm (
			"pextq %2, %1, %0"
			: "=r" (dst)
			: "r" (word), "r" (mask)
		);
		return dst;
	}
#endif
	
	
	inline uint32_t and_not(uint32_t const a, uint32_t const b)
	{
		uint32_t dst{0};
		asm (
			"andnl %2, %1, %0"
			: "=r" (dst)
			: "r" (a), "r" (b)
		);
		return dst;
	}
	
	
	inline uint64_t and_not(uint64_t const a, uint64_t const b)
	{
		uint64_t dst{0};
		asm (
			"andnq %2, %1, %0"
			: "=r" (dst)
			: "r" (a), "r" (b)
		);
		return dst;
	}
	
	
	inline uint16_t count_set_bits(uint16_t const word)
	{
		uint16_t count{0};
		asm (
			"popcntw %1, %0"
			: "=r" (count)
			: "r" (word)
		);
		return count;
	}

	
	inline uint32_t count_set_bits(uint32_t const word)
	{
		uint32_t count{0};
		asm (
			"popcntl %1, %0"
			: "=r" (count)
			: "r" (word)
		);
		return count;
	}


	inline uint64_t count_set_bits(uint64_t const word)
	{
		uint64_t count{0};
		asm (
			"popcntq %1, %0"
			: "=r" (count)
			: "r" (word)
		);
		return count;
	}

	
	inline uint16_t lzcnt(uint16_t const val)
	{
		uint16_t retval(0);
		__asm__ (
			"lzcntw %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}
	
	
	inline uint32_t lzcnt(uint32_t const val)
	{
		uint32_t retval(0);
		__asm__ (
			"lzcntl %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}
	

	inline uint64_t lzcnt(uint64_t const val)
	{
		uint64_t retval(0);
		__asm__ (
			"lzcntq %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}

	
	inline uint16_t tzcnt(uint16_t const val)
	{
		uint16_t retval(0);
		__asm__ (
			"tzcntw %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}
	
	
	inline uint32_t tzcnt(uint32_t const val)
	{
		uint32_t retval(0);
		__asm__ (
			"tzcntl %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}
	

	inline uint64_t tzcnt(uint64_t const val)
	{
		uint64_t retval(0);
		__asm__ (
			"tzcntq %1, %0;"
			: "=r" (retval)
			: "r" (val)
			: "cc"
		);
		
		return retval;
	}


	template <typename t_type>
	inline bool is_exact_power_of_2(t_type const val)
	{
		// From http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
		return val && !(val & (val - 1));
	}
	
	
	// t_type must be uint{16, 32, 64}_t.
	template <typename t_type>
	inline t_type gte_exp_base_2(t_type const val) 
	{
		// Check for zero.
		if (!val)
			return 1;
		
		// If the value has only one bit set, it is a power of two.
		t_type p1(lzcnt(val));
		t_type p2(tzcnt(val));
		if (p1 == CHAR_BIT * sizeof(val) - p2 - 1)
			return p1;
		
		return CHAR_BIT * sizeof(val) - p1;
	}
}

#endif
