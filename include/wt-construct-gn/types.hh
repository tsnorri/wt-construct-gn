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

#ifndef WT_CONSTRUCT_GN_TYPES_HH
#define WT_CONSTRUCT_GN_TYPES_HH

#include <wt-construct-gn/construct_wt_2.hh>


namespace wtcgn {

	// Fix t_word and threads.
	template <bool t_multithreaded>
	struct construct_type_tpl
	{
		template <typename t_wt>
		using type = wtcgn::construct_wt_gn <t_wt, uint64_t, t_multithreaded>;
	};

	
	typedef sdsl::wt_pc <
		sdsl::balanced_shape,
		sdsl::bit_vector,
		sdsl::bit_vector::rank_1_type,
		sdsl::bit_vector::select_1_type,
		sdsl::bit_vector::select_0_type,
		//sdsl::byte_tree <>
		sdsl::int_tree <>
	> wt_type_sdsl;
	
	
	typedef sdsl::wt_pc <
		sdsl::balanced_shape,
		sdsl::bit_vector,
		sdsl::bit_vector::rank_1_type,
		sdsl::bit_vector::select_1_type,
		sdsl::bit_vector::select_0_type,
		//sdsl::byte_tree <>,
		sdsl::int_tree <>,
		construct_type_tpl <false>::type
	> wt_type_gn;
	
	
	typedef sdsl::wt_pc <
		sdsl::balanced_shape,
		sdsl::bit_vector,
		sdsl::bit_vector::rank_1_type,
		sdsl::bit_vector::select_1_type,
		sdsl::bit_vector::select_0_type,
		//sdsl::byte_tree <>,
		sdsl::int_tree <>,
		construct_type_tpl <true>::type
	> wt_type_gn_mt;
}

#endif
