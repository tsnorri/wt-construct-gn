/*
 * Copyright (c) 2016 Tuukka Norri
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

#ifndef WT_CONSTRUCT_GN_TIMER_HH
#define WT_CONSTRUCT_GN_TIMER_HH

#include <chrono>


namespace wtcgn {
	
	class timer
	{
	protected:
		std::chrono::milliseconds m_start{};
		std::chrono::milliseconds m_end{};
		
	public:
		// Get the current time in milliseconds.
		static inline std::chrono::milliseconds timestamp_ms_now()
		{
			namespace chrono = std::chrono;
			return chrono::duration_cast <chrono::milliseconds>(chrono::system_clock::now().time_since_epoch());
		}
		
		timer():
			m_start(timestamp_ms_now())
		{
		}
		
		inline void start() { m_start = timestamp_ms_now(); }
		inline void stop() { m_end = timestamp_ms_now(); }
		inline std::chrono::milliseconds::rep ms_elapsed() const { return (m_end - m_start).count(); }
	};
}

#endif
