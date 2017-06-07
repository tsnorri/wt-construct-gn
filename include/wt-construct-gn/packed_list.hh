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

#ifndef WT_CONSTRUCT_GN_PACKED_LIST_HH
#define WT_CONSTRUCT_GN_PACKED_LIST_HH

#include <sdsl/int_vector.hpp>


namespace wtcgn {
	
	// Create bit masks to access the first bit of each character in a word.
	template <typename t_word>
	class item_masks
	{
	public:
		enum { T_WORD_BIT = CHAR_BIT * sizeof(t_word) };
	
	protected:
		sdsl::int_vector <T_WORD_BIT>	m_masks;
	
	public:
		std::size_t size() const { return m_masks.size(); }
	
		item_masks() = default;
		
		item_masks(std::size_t count):
			m_masks(count, 0)
		{
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
				m_masks[i] = mask;
			}
		}
	
		t_word operator[](std::size_t const i) const { return m_masks[i]; }
	};
	
	
	// Pack as many whole items of given size into each t_word.
	template <typename t_word>
	class packed_list
	{
	public:
		enum { T_WORD_BIT = CHAR_BIT * sizeof(t_word) };
		typedef sdsl::int_vector <T_WORD_BIT> int_vector;
	
	protected:
		int_vector	m_vec;
		std::size_t	m_idx{0};
		std::size_t m_item_bits{0};
		std::size_t m_items_per_word{0};
		t_word		m_item_mask{0};
	
	public:
		packed_list() = default;
	
		// Initialize with space for a given number of item_bits-sized items.
		packed_list(std::size_t const size, std::size_t item_bits):
			m_vec(1 + size / (T_WORD_BIT / item_bits), 0),
			m_item_bits(item_bits),
			m_items_per_word(T_WORD_BIT / m_item_bits),
			m_item_mask(~0)
		{
			m_item_mask >>= (T_WORD_BIT - m_item_bits);

		}
	
		// Access the underlying vector.
		int_vector const &word_vector() const { return m_vec; }
		
		// Number of m_item_bits-sized items.
		std::size_t size() const { return m_idx; }
		
		// Number of bits per item.
		std::size_t item_bits() const { return m_item_bits; }
		
		// Reserve space for item_count items.
		void reserve(std::size_t const item_count)
		{
			auto const word_count(1 + item_count / m_items_per_word);
			if (m_vec.size() < word_count)
				m_vec.resize(word_count);
		}
	
		// Append a value to the list.
		inline bool append(t_word const val)
		{
			auto const vec_idx(m_idx / m_items_per_word);
		
			auto const current_size(m_vec.size());
			if (current_size <= vec_idx)
			{
				std::cerr << "Warning: resizing." << std::endl;
				m_vec.resize(1 + vec_idx);
			}
			
			auto &ref(m_vec[vec_idx]);
			auto const shift_amt(m_item_bits * (m_idx % m_items_per_word));	// XXX remainder operation may be slow. Restricting to powers of two could be better.
			auto const shifted_val(val << shift_amt);
			ref |= shifted_val;
		
			++m_idx;
			
			// Check if a word was filled.
			return (m_idx / m_items_per_word != vec_idx);
		}
	
		inline t_word value(std::size_t idx) const
		{
			auto const vec_idx(idx / m_items_per_word);
			t_word retval(m_vec[vec_idx]);
		
			auto const shift_amt(m_item_bits * (idx % m_items_per_word));	// XXX remainder operation may be slow. Restricting to powers of two could be better.
			retval >>= shift_amt;
			retval &= m_item_mask;
		
			return retval;
		}
	};
}

#endif
