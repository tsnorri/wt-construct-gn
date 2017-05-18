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

#ifndef WT_CONSTRUCT_BIT_VECTOR_CONTAINER_HH
#define WT_CONSTRUCT_BIT_VECTOR_CONTAINER_HH

#include <sdsl/int_vector.hpp>


namespace wtcgn
{
	class bit_vector_container
	{
	protected:
		sdsl::bit_vector	*m_vector{nullptr};
		std::size_t			m_idx{0};
		
	public:
		bit_vector_container() = default;
		
		bit_vector_container(sdsl::bit_vector *vec):
			m_vector(vec)
		{
		}
		
		sdsl::bit_vector const &vector() const { return *m_vector; }
		sdsl::bit_vector &vector() { return *m_vector; }
		
		void append(sdsl::bit_vector const &vec)
		{
			auto *dst_words(m_vector->data());
			auto *src_words(vec.data());
			auto len(vec.size());
			auto const offset(m_idx / BIT_VECTOR_BIT);
			dst_words += offset;
			
			if (m_idx % BIT_VECTOR_BIT)
			{
				// Set the lowest and highest bits of each dst word.
				auto const shift_amt_1(m_idx);
				auto const shift_amt_2(BIT_VECTOR_BIT - m_idx);
				assert(BIT_VECTOR_BIT != shift_amt_1);
				assert(BIT_VECTOR_BIT != shift_amt_2);
				
				while (len >= BIT_VECTOR_BIT)
				{
					auto val_1(*src_words << shift_amt_1);
					auto val_2(*src_words >> shift_amt_2);
					++src_words;
					
					*dst_words |= val_1;
					++dst_words;
					*dst_words |= val_2;
					
					len -= BIT_VECTOR_BIT;
					m_idx += BIT_VECTOR_BIT;
				}
			}
			else
			{
				// Just copy words.
				while (len >= BIT_VECTOR_BIT)
				{
					*dst_words++ = *src_words++;
					len -= BIT_VECTOR_BIT;
					m_idx += BIT_VECTOR_BIT;
				}
			}
			
			// Set the remaining bits if needed.
			if (len)
			{
				auto const shift_amt(m_idx % BIT_VECTOR_BIT);
				auto const remaining(BIT_VECTOR_BIT - shift_amt);

				auto const val(*src_words);
				auto shifted(val);
				shifted <<= shift_amt;
				*dst_words |= shifted;

				// FIXME: the branch isn't really needed if there is enough space for one extra word.
				if (remaining < len)
				{
					auto shifted(val);
					shifted >>= remaining;
					++dst_words;
					*dst_words |= shifted;
				}
				
				m_idx += len;
			}
		}
	};
}

#endif
