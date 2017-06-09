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

#ifndef WT_CONSTRUCT_GN_ITEM_MASKS_HH
#define WT_CONSTRUCT_GN_ITEM_MASKS_HH

#include <sdsl/int_vector.hpp>


namespace wtcgn {
	
	template <typename t_word>
	class create_masks
	{
	public:
		enum { T_WORD_BIT = CHAR_BIT * sizeof(t_word) };
		typedef sdsl::int_vector <T_WORD_BIT> vector_type;
		
		
		create_masks() = delete;
		
		
		// Create bit masks to access the first bit of each character in a word.
		static void create_item_masks(std::size_t count, vector_type /* out */ &vec)
		{
			vector_type temp_vec(count, 0);

			for (std::size_t i(0); i < count; ++i)
			{
				t_word const initial_mask(0x1);
				t_word mask(0);
			
				for (std::size_t j(0), shift_count(T_WORD_BIT / (1 + i)); j < shift_count; ++j)
				{
					mask <<= 1 + i;
					mask |= initial_mask;
				}
			
				//mask <<= (T_WORD_BIT % (1 + i));
				temp_vec[i] = mask;
			}
			
			vec = std::move(temp_vec);
		}
		
		
		// Create bit masks to access the first tau bits of each character in a word.
		static void create_tau_masks(
			std::size_t item_bits,
			std::size_t const tau,
			vector_type const &item_masks,
			vector_type /* out */ &vec
		)
		{
			vector_type temp_vec(item_bits / tau, 0);
			
			// Create a mask with bits set for tau.
			t_word initial_mask(0);
			initial_mask = ~initial_mask;
			initial_mask >>= T_WORD_BIT - tau;
			
			std::size_t j(0);
			while (tau <= item_bits)
			{
				t_word const factor(item_masks[item_bits - 1]);
				t_word const mask(factor * initial_mask);
				temp_vec[j] = mask;
				
				++j;
				item_bits -= tau;
			}
			
			vec = std::move(temp_vec);
		}
		
		
		// Create bit masks with all possible values with tau bits repeated and zero-padded on the right.
		static void create_tau_values(
			std::size_t const tau,
			vector_type const &item_masks,
			vector_type /* out */ &vec
		)
		{
			t_word const count(1 << tau);
			vector_type temp_vec(count, 0);
			auto const mask(item_masks[tau - 1]);
			for (t_word i(0); i < count; ++i)
			{
				t_word const val(i * mask);
				temp_vec[i] = val;
			}
			
			vec = std::move(temp_vec);
		}
	};
}

#endif
