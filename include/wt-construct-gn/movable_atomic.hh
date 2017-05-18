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

#ifndef WT_CONSTRUCT_GN_MOVABLE_ATOMIC_HH
#define WT_CONSTRUCT_GN_MOVABLE_ATOMIC_HH

#include <atomic>


namespace wtcgn
{

	// Atomic that may be moved before its value becomes significant.
	template <typename t_type>
	class movable_atomic
	{
	protected:
		std::atomic <t_type> m_value;
		
		
	public:
		movable_atomic(t_type const val):
			m_value(val)
		{
		}
		
		movable_atomic(movable_atomic &&other):
			m_value(other.m_value)
		{
		}
		
		movable_atomic &operator=(movable_atomic &&other) &
		{
			m_value = t_type(other.m_value);
			return *this;
		}
		
		// Debugging helpers.
		t_type operator++()
		{
			auto const retval(++m_value);
			return retval;
		}
		
		t_type operator--()
		{
			auto const retval(--m_value);
			return retval;
		}
		
		t_type operator+=(std::size_t const val)
		{
			auto const retval(m_value += val);
			return retval;
		}
	};	
}

#endif
