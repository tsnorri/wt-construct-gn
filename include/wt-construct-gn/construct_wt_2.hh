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

#ifndef WT_CONSTRUCT_GN_CONSTRUCT_WT_2_HH
#define WT_CONSTRUCT_GN_CONSTRUCT_WT_2_HH

//#define INLINE_CREATION __attribute__ ((noinline))
#define INLINE_CREATION inline

#define UNROLL_LOOPS __attribute__((optimize("unroll-loops")))


#include <iostream>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/ram_fs.hpp>
#include <sdsl/util.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <string>
#include <vector>
#include <wt-construct-gn/bit_ops.hh>
#include <wt-construct-gn/bit_vector_container.hh>
#include <wt-construct-gn/create_masks.hh>
#include <wt-construct-gn/dispatch_fn.hh>
#include <wt-construct-gn/movable_atomic.hh>
#include <wt-construct-gn/packed_list.hh>


namespace wtcgn
{
	// Does not destroy objects immediately.
	template <typename t_element>
	class array
	{
	protected:
		t_element	*m_array{nullptr};
		std::size_t	m_size{0};
		std::size_t	m_idx{0};
		
	public:
		array() = default;
		
		array(t_element *array, std::size_t const size):
			m_array(array),
			m_size(size)
		{
		}
		
		std::size_t const size() { return m_idx; }
		
		t_element operator[](std::size_t const idx) const
		{
			assert(idx < m_size);
			return m_array[idx];
		}
		
		t_element &operator[](std::size_t const idx)
		{
			assert(idx < m_size);
			return m_array[idx];
		}
		
		t_element pop()
		{
			return m_array[--m_idx];
		}
		
		void push(t_element const val)
		{
			assert(m_idx < m_size);
			m_array[m_idx++] = val;
		}
		
		void clear()
		{
			m_idx = 0;
		}
	};
	
	
	template <typename t_word>
	class big_node
	{
	public:
		typedef packed_list <t_word> character_list;
	
	protected:
		character_list	m_characters{};
		t_word			m_suffix{0};
		std::size_t		m_suffix_length{0};
		std::size_t		m_node_index{0};
		bool			m_is_leaf{false};
	
	public:
		character_list &characters() { return m_characters; }
		character_list const &characters() const { return m_characters; }
		
		std::size_t node_index() const { return m_node_index; }
		void set_node_index(std::size_t const node_index) { m_node_index = node_index; }
		
		t_word suffix() const { return m_suffix; }
		void set_suffix(t_word const suffix) { m_suffix = suffix; }
	
		std::size_t suffix_length() const { return m_suffix_length; }
		void set_suffix_length(std::size_t const suffix_length) { m_suffix_length = suffix_length; }
		
		bool is_leaf() const { return m_is_leaf; }
		void set_leaf(bool val) { m_is_leaf = val; }
	};
	
	
	template <typename t_vec>
	void remove_vector_items(t_vec &vec)
	{
		auto const original_size(vec.size());
		vec.clear();
		vec.resize(original_size);
	}
	
	
	template <typename t_wt>
	struct construct_wt_gn_delegate
	{
		virtual ~construct_wt_gn_delegate() {}
		virtual std::size_t tau(t_wt const &wt) = 0;
		virtual void finish_construction() = 0;
	};
	
	
	template <typename t_word>
	inline void place_bits(
		sdsl::bit_vector &target_bv,
		t_word extracted,
		std::size_t const pos,
		std::size_t const bit_count
	)
	{
#if 1
		if (pos <= 67 && 67 <= pos + bit_count)
			std::cerr << "extracted: " << extracted << std::endl;
#endif
		
		auto const vec_idx(pos / BIT_VECTOR_BIT);
		auto const offset(pos % BIT_VECTOR_BIT);
		auto const remaining(BIT_VECTOR_BIT - offset);
		uint64_t *bv_data(target_bv.data());
		
		sdsl::std_size_type_for_int_vector src(extracted);
		src <<= offset;
		bv_data[vec_idx] |= src;
		
		if (remaining < bit_count)
		{
			sdsl::std_size_type_for_int_vector src(extracted);
			src >>= remaining;
			bv_data[vec_idx + 1] |= src;
		}
	}
	
	
	template <typename t_wt, typename t_word, bool t_multithreaded>
	class construct_wt_gn
	{
	protected:
		class creation_context;
		
	public:
		enum {
			T_WORD_BIT = CHAR_BIT * sizeof(t_word),
			TREE_NODE_BIT = CHAR_BIT * sizeof(typename t_wt::tree_strat_type::node_type)
		};
		typedef big_node <t_word>	big_node_type;
		
	protected:
		typedef sdsl::int_vector <TREE_NODE_BIT>			shortcut_vec_type;
		typedef typename t_wt::tree_strat_type::node_type	tree_node_type;
		typedef create_masks <t_word>						create_masks_type;
		
		struct node
		{
			t_word mask{0};
			tree_node_type idx{0};
			tree_node_type level{0};
			
			node() = default;
			
			node(tree_node_type const idx_, tree_node_type const level_, t_word const mask_):
				mask(mask_),
				idx(idx_),
				level(level_)
			{
			}
		};
		
	protected:
		t_wt									*m_wt{nullptr};
		construct_wt_gn_delegate <t_wt>			*m_delegate{nullptr};
		std::vector <shortcut_vec_type>			m_node_paths;
		t_word									m_tau_mask{};
		sdsl::bit_vector						m_target_bv;
		typename create_masks_type::vector_type	m_item_masks;
		typename create_masks_type::vector_type	m_tau_masks;
		typename create_masks_type::vector_type	m_tau_values;
		std::size_t								m_tau{};
		std::size_t								m_sigma_bits{};
		
	protected:
		construct_wt_gn(t_wt &wt, construct_wt_gn_delegate <t_wt> *delegate):
			m_wt(&wt),
			m_delegate(delegate),
			m_tau_mask(~0),
			m_tau(delegate->tau(wt))
		{
			m_tau_mask >>= (T_WORD_BIT - m_tau);
			print_tree(m_wt->m_tree.root(), 0);
			
			std::size_t i(0);
			for (auto const path : m_wt->m_tree.m_path)
				std::cerr << "[" << i++ << "]: " << (path & 0xFFFFFFFFFFFFFF) << " (" << (path >> 56) << ")" << std::endl;
		}
		
		
		void print_tree(typename t_wt::tree_strat_type::node_type const node, std::size_t const level)
		{
#if 1
			auto const &val(m_wt->m_tree.m_nodes[node]);
			for (std::size_t i(0); i < level; ++i)
					std::cerr << "    ";
			std::cerr << "[" << node << "] bv_pos: " << val.bv_pos << " bv_pos_rank: " << val.bv_pos_rank << " size: " << m_wt->m_tree.size(node) << std::endl;
			
			if (!m_wt->m_tree.is_leaf(node))
			{
				print_tree(m_wt->m_tree.child(node, 0), 1 + level);
				print_tree(m_wt->m_tree.child(node, 1), 1 + level);
			}
#endif
		}
		
		
		inline void extract_and_place(
			std::size_t const node_idx,
			t_word const w,
			t_word const mask,
			std::size_t const bit_count,
			std::vector <std::size_t> &bv_node_pos
		)
		{
			// extracted will contain the bits in places indicated by mask.
			t_word const extracted(extract_bits(w, mask));
			
			// Store the bits.
			place_bits(m_target_bv, extracted, bv_node_pos[node_idx], bit_count);
			bv_node_pos[node_idx] += bit_count;
		}
		
		
		template <bool t_right>
		inline void create_small_node(
			std::size_t const node_idx,
			std::size_t const level,
			t_word const w,
			t_word const mask,
			std::vector <std::size_t> &bv_node_pos,
			array <node> &node_stack
		)
		{
			auto const child_idx(m_wt->m_tree.child(node_idx, t_right));
			if (!m_wt->m_tree.is_leaf(child_idx))
			{
				node const n(child_idx, level, mask);
				node_stack.push(n);
				
				auto const bit_count(count_set_bits(mask));
				assert(bit_count);
				extract_and_place(child_idx, w, mask, bit_count, bv_node_pos);
			}
		}
		
		
		inline void create_small_nodes(
			t_word const w,
			std::vector <std::size_t> &bv_node_pos,
			array <node> &node_stack
		)
		{
			while (node_stack.size())
			{
				auto const node(node_stack.pop());
				auto const node_idx(node.idx);
				auto const level(node.level);
				auto const mask(node.mask);
				
				//std::cerr << "Node: " << node_idx << " mask: " << mask << std::endl;
				
				// Create masks for the next level.
				t_word const left_mask((~w & mask) << 1);	// Find zeros in masked positions.
				t_word const right_mask((w & mask) << 1);	// Find ones in masked positions.
				
				auto const next_level(1U + level);
				if (left_mask && next_level < m_tau)
					create_small_node <false>(node_idx, next_level, w, left_mask, bv_node_pos, node_stack);
				
				if (right_mask && next_level < m_tau)
					create_small_node <true>(node_idx, next_level, w, right_mask, bv_node_pos, node_stack);
			}
		}
		
		
		INLINE_CREATION
		void create_small_nodes(
			big_node_type const &parent_node,
			std::size_t const level,
			std::vector <std::size_t> &bv_node_pos,
			array <node> &node_stack
		)
		{
			auto &characters(parent_node.characters());
			auto const &word_vector(characters.word_vector());
			auto const parent_idx(parent_node.node_index());
			
			auto const mask_idx(m_sigma_bits - level - 1);
			t_word const mask(m_item_masks[mask_idx]);
			auto const bit_count(count_set_bits(mask));
			
			std::size_t word_idx(0);
			std::size_t const count(word_vector.size());
			if (count)
			{
				while (word_idx < count - 1)
				{
					node_stack.clear();
					node n(parent_idx, 0, mask);
					node_stack.push(n);

					auto const w(word_vector[word_idx]);
					extract_and_place(parent_idx, w, mask, bit_count, bv_node_pos);
					create_small_nodes(w, bv_node_pos, node_stack);
					++word_idx;
				}
		
				{
					node_stack.clear();
					node n(parent_idx, 0, mask);
					node_stack.push(n);

					auto const word_count(word_vector.size());
					auto const character_count(characters.size());
					auto const character_bits(characters.item_bits());
				
					t_word last_character_mask(~0);
					auto const diff(character_count / word_count * character_bits);
					last_character_mask >>= T_WORD_BIT - diff;
					last_character_mask &= mask;
					auto const bit_count(count_set_bits(mask));
				
					auto const w(word_vector[word_idx]);
					extract_and_place(parent_idx, w, mask, bit_count, bv_node_pos);
					
					create_small_nodes(w, bv_node_pos, node_stack);
				}
			}
		}
		
		
		UNROLL_LOOPS
		INLINE_CREATION
		bool create_big_node_for_word(
			t_word const w,
			t_word const tau_mask,
			t_word const tau_rep_mask,
			t_word zero_limit_mask,
			t_word const parent_suffix,
			std::size_t const non_tau,
			std::size_t const items_per_word,
			std::size_t const parent_suffix_length,
			std::size_t const child_character_bits,
			shortcut_vec_type const &paths,
			big_node_type const &parent_node,
			std::vector <big_node_type> &output_nodes
		) const
		{
			bool retval(false);
			
			// Separate the first tau bits of each character and the remaining bits.
			t_word tau_bits(extract_bits(w, tau_rep_mask));
			t_word remaining_bits(extract_bits(w, ~tau_rep_mask));
			
			std::size_t remaining_characters(items_per_word);
			while (true)
			{
				// Get the value of the lowermost tau bits.
				t_word const first_tau_val(tau_bits & tau_mask);
				
				// Get the value repeated.
				t_word const first_tau_val_rep(m_tau_values[first_tau_val]);
				
				// Use a ^ b to set zeros to the matching suffixes of length tau.
				t_word zero_mask(first_tau_val_rep ^ tau_bits);
				
				// Make sure that there are not more zeros than the maximum number of characters.
				zero_mask |= zero_limit_mask;
				
				// Count the matching bits.
				t_word const zero_count(tzcnt(zero_mask));
				
				// Use integer division to count the matching characters.
				t_word const matching_character_count(zero_count / m_tau);
				
				// Determine the number of bits in the characters.
				t_word const matching_tau_bit_count(matching_character_count * m_tau);
				t_word const matching_non_tau_bit_count(matching_character_count * non_tau);
				
				t_word extract_mask(0);
				extract_mask = ~extract_mask;
				extract_mask >>= T_WORD_BIT - matching_non_tau_bit_count;
				
				t_word const stored_bits(remaining_bits & extract_mask);
				
				// Store the bits.
				// FIXME: if stored_bits is zero, m_idx in packed_list could be incremented without storing anything.
				{
					t_word const child_suffix_part(first_tau_val << parent_suffix_length);
					t_word const child_suffix(parent_suffix | child_suffix_part);
					
					// Get the vector for the child in question.
					assert(child_suffix < output_nodes.size());
					auto &child_node(output_nodes[child_suffix]);
					if (child_node.is_leaf())
						goto update_remaining_characters;
					
					// Check that there is enough space.
					auto &child_node_characters(child_node.characters());
					if (0 == child_node_characters.size())
					{
						// Child node has not been set up; assign values.
						auto const child_suffix_length(parent_node.suffix_length() + m_tau);
				
						child_node.set_suffix(child_suffix);
						child_node.set_suffix_length(child_suffix_length);
				
						// Count the number of characters in the current node.
						assert(child_suffix < paths.size());
						auto const tree_node(paths[child_suffix]);
						auto const child_character_count(m_wt->m_tree.size(tree_node));
						//std::cerr << "Child suffix: " << child_suffix << " len: " << child_suffix_length << " tree_node: " << tree_node << " child_character_count: " << child_character_count << std::endl;
				
						child_node.set_node_index(tree_node);
				
						bool is_leaf(m_wt->m_tree.is_leaf(tree_node));
						child_node.set_leaf(is_leaf);
						if (is_leaf)
							goto update_remaining_characters;
				
						// Replace the empty list with an initialized one.
						typename big_node_type::character_list temp_list(child_character_count, child_character_bits);
						child_node_characters = std::move(temp_list);
				
						retval = true;
					}
			
					//std::cerr << "child suffix: " << child_suffix << " word vector size: " << child_node_characters.word_vector().size() << std::endl;
					child_node_characters.append(stored_bits, matching_character_count);
				}
				
				if (remaining_characters < matching_character_count)
					std::cerr << "Error" << std::endl;
				
			update_remaining_characters:
				remaining_characters -= matching_character_count;
				if (!remaining_characters)
					break;
				
				zero_limit_mask >>= matching_tau_bit_count;
				tau_bits >>= matching_tau_bit_count;
				remaining_bits >>= matching_non_tau_bit_count;
			}
			
			return retval;
		}
		
		
		UNROLL_LOOPS
		INLINE_CREATION
		bool create_big_nodes(
			big_node_type const &parent_node,
			std::vector <big_node_type> &output_nodes,
			std::vector <std::size_t> &bv_node_pos,
			array <node> &node_stack,
			std::size_t const alpha
		) const
		{
			auto const &characters(parent_node.characters());
			if (!characters.size())
				return false;

			auto const &words(characters.word_vector());
			auto const &paths(m_node_paths[alpha - 1]);
			auto const level(alpha * m_tau);
			auto const character_bits(characters.item_bits());
			auto const items_per_word(characters.items_per_word());
			auto const non_tau(character_bits - m_tau);
			auto const child_character_bits(character_bits - m_tau);
			auto const tau_rep_mask(m_tau_masks[alpha - 1]);
			auto const parent_suffix_length(parent_node.suffix_length());
			t_word const parent_suffix(parent_node.suffix());
			
			t_word zero_limit_mask(0x1);
			zero_limit_mask <<= m_tau * items_per_word;
			
			bool retval(false);
			
			t_word tau_mask(0);
			tau_mask = ~tau_mask;
			tau_mask >>= T_WORD_BIT - m_tau;
			
			auto const count(words.size());
			if (1 < count)
			{
				for (std::size_t i(0); i < count - 1; ++i)
				{
					retval |= create_big_node_for_word(
						words[i],
						tau_mask,
						tau_rep_mask,
						zero_limit_mask,
						parent_suffix,
						non_tau,
						items_per_word,
						parent_suffix_length,
						child_character_bits,
						paths,
						parent_node,
						output_nodes
					);
				}
				
				if (count)
				{
					auto const remaining_characters(characters.size() % items_per_word);
					t_word zero_limit_mask(0x1);
					zero_limit_mask <<= m_tau * remaining_characters;

					retval |= create_big_node_for_word(
						words[count - 1],
						tau_mask,
						tau_rep_mask,
						zero_limit_mask,
						parent_suffix,
						non_tau,
						remaining_characters,
						parent_suffix_length,
						child_character_bits,
						paths,
						parent_node,
						output_nodes
					);
				}
			}
			
			return retval;
		}
		
		
		// Create a list of words each of which contains as many characters as possible s.t.
		// no characters are split between words and the characters start from bit index zero (0-based).
		// append() contains a loop that should be unrolled.
		UNROLL_LOOPS
		INLINE_CREATION
		void create_big_node_initial(
			sdsl::int_vector_buffer <t_wt::tree_strat_type::int_width> &input_buf,
			std::size_t const input_size,
			big_node_type &output_node,
			std::vector <std::size_t> &bv_node_pos,
			array <node> &node_stack
		) const
		{
			auto &characters(output_node.characters());
	
			{
				assert(m_sigma_bits);
				typename big_node_type::character_list temp_list(input_size, m_sigma_bits);
				characters = std::move(temp_list);
			}
			
			auto const &word_vector(characters.word_vector());
			
			// Level zero.
			auto const mask_idx(m_sigma_bits - 1);
			t_word const mask(m_item_masks[mask_idx]);
			auto const bit_count(count_set_bits(mask));
			
			//std::cerr << "Input size: " << input_size << std::endl;
			
			std::size_t word_idx(0);
			for (std::size_t i(0); i < input_size; ++i)
			{
				auto const c(input_buf[i]);
				t_word comp(m_wt->m_tree.bit_path(c));
				comp &= 0xFFFFFFFFFFFFFF; // Set the 8 most significant bits to zero.
				//std::cerr << "Appending '" << (char) c << "', comp " << comp << std::endl;
				characters.append(comp);
			}
		}
		
		
		std::size_t find_lowest_tree_level(tree_node_type const node, std::size_t const level) const
		{
			if (m_wt->m_tree.is_leaf(node))
				return level;
			
			auto const left(m_wt->m_tree.child(node, 0));
			auto const right(m_wt->m_tree.child(node, 1));
			auto const left_res(find_lowest_tree_level(left, 1 + level));
			auto const right_res(find_lowest_tree_level(right, 1 + level));
			return std::max(left_res, right_res);
		}
		
		
		std::size_t find_lowest_tree_level() const
		{
			auto const root(m_wt->m_tree.root());
			return find_lowest_tree_level(root, 0);
		}
		
		
		// Find each valid path of length m_tau in m_wt->m_tree and map character -> node.
		void construct_shortcuts(
			shortcut_vec_type &target,
			tree_node_type const node,
			t_word const character,
			std::size_t const idx,
			std::size_t const idx_limit
		) const
		{
			if (idx == idx_limit)
			{
				target[character] = node;
			}
			else if (!m_wt->m_tree.is_leaf(node))
			{
				auto const left_child(m_wt->m_tree.child(node, 0));
				auto const right_child(m_wt->m_tree.child(node, 1));
				
				t_word mask(1 << idx);
				construct_shortcuts(target, left_child, character, 1 + idx, idx_limit);
				construct_shortcuts(target, right_child, mask | character, 1 + idx, idx_limit);
			}
		}
		
		
		// Handle each tau-th level of the Wavelet tree.
		void construct_shortcut_level(
			std::size_t const level,
			std::size_t const level_count,
			std::size_t const item_count,
			std::size_t const initial_item_count,
			shortcut_vec_type const &previous
		)
		{
			if (level < level_count)
			{
				std::size_t const idx(level * m_tau);
				std::size_t const idx_limit(idx + m_tau);
				shortcut_vec_type vec(item_count, t_wt::tree_strat_type::undef);
				for (t_word i(0), count(previous.size()); i < count; ++i)
				{
					auto const node_idx(previous[i]);
					if (t_wt::tree_strat_type::undef != node_idx)
						construct_shortcuts(vec, node_idx, i, idx, idx_limit);
				}
			
				m_node_paths.emplace_back(std::move(vec));
				
				// Tail recursion.
				construct_shortcut_level(1 + level, level_count, item_count * initial_item_count, initial_item_count, m_node_paths.back());
			}
		}
		
		
		void construct(
			sdsl::int_vector_buffer <t_wt::tree_strat_type::int_width> &input_buf,
			typename t_wt::size_type const input_size,
			std::size_t const level_nodes
		)
		{
			// Check how many characters may be packed into a word.
			m_sigma_bits = gte_exp_base_2(m_wt->m_sigma);
			assert(m_sigma_bits);
			
			// Node index stack.
			node idx_buf[level_nodes];
			array <node> node_stack(idx_buf, level_nodes);
			
			// Initialize starting positions of Wavelet tree nodes like in SDSL.
			std::vector <std::size_t> bv_node_pos(m_wt->m_tree.size(), 0);
			for (std::size_t i(0), count(bv_node_pos.size()); i < count; ++i)
				bv_node_pos[i] = m_wt->m_tree.bv_pos(i);
			
			auto const lowest_tree_level(find_lowest_tree_level());
			auto const lowest_big_node_level(0 == lowest_tree_level ? 0 : (lowest_tree_level - 1) / m_tau * m_tau);
			auto const max_big_nodes(1 << lowest_big_node_level);
			std::cerr << "Lowest tree level: " << lowest_tree_level << std::endl;
			std::cerr << "Lowest big node level: " << lowest_big_node_level << std::endl;
			
			// Create bit masks for accessing the first bits of characters.
			create_masks_type::create_item_masks(m_sigma_bits, m_item_masks);
			
			// Create bit masks for accessing the first tau bits of each character.
			create_masks_type::create_tau_masks(lowest_tree_level, m_tau, m_item_masks, m_tau_masks);
			
			// Create a vector with all possible values of tau repeated.
			if (m_tau <= lowest_tree_level)
				create_masks_type::create_tau_values(m_tau, m_item_masks, m_tau_values);
			
			// Shortcuts for m_wt->m_tree.m_nodes.
			auto const nonzero_big_node_levels(lowest_tree_level / m_tau);
			if (nonzero_big_node_levels)
			{
				auto const level_count(nonzero_big_node_levels);
				auto const initial_item_count(1 << m_tau);
				m_node_paths.reserve(level_count);
				auto const root_node(m_wt->m_tree.root());
				
				shortcut_vec_type initial_suffix(1, 0);
				construct_shortcut_level(0, level_count, initial_item_count, initial_item_count, initial_suffix);
			}
			
#if 1
			{
				std::size_t i(0);
				for (auto const &v : m_node_paths)
				{
					std::cerr << "Level " << i++ << std::endl;
					std::size_t j(0);
					for (auto const val : v)
						std::cerr << "\t" << j++ << ": " << val << std::endl;
				}
			}
#endif
			
			std::vector <big_node_type> big_nodes_1(max_big_nodes);
			std::vector <big_node_type> big_nodes_2(max_big_nodes);
			
			{
				big_node_type root;
				create_big_node_initial(input_buf, input_size, root, bv_node_pos, node_stack);
				
				create_small_nodes(root, 0, bv_node_pos, node_stack);
				big_nodes_1[0] = std::move(root);
			}
			
			for (std::size_t alpha(1); alpha * m_tau <= lowest_big_node_level; ++alpha)
			{
				for (auto const &big_node : big_nodes_1)
					create_big_nodes(big_node, big_nodes_2, bv_node_pos, node_stack, alpha);
				
				for (auto const &big_node : big_nodes_2)
				{
					if (big_node.characters().size())
						create_small_nodes(big_node, alpha * m_tau, bv_node_pos, node_stack);
				}
				
				using std::swap;
				swap(big_nodes_1, big_nodes_2);
				//std::cerr << "Removing items" << std::endl;
				remove_vector_items(big_nodes_2);
			}
			
#if 0
			std::cerr << "WT: " << std::endl;
			{
				std::size_t j(0);
				for (std::size_t i(0); i < m_target_bv.size(); ++i)
				{
					std::cerr << m_target_bv[i];
					++j;
					if (j == input_size)
					{
						j = 0;
						std::cerr << std::endl;
					}
				}
				std::cerr << std::endl;
			}
#endif
			
			m_wt->m_bv = typename t_wt::bit_vector_type(std::move(m_target_bv));
			m_wt->finish_construct();

			// Get out of WT's constructor in single-threaded mode.
			{
				auto *delegate(m_delegate);
				auto fn = [delegate](){
					delegate->finish_construction();
				};
				dispatch_async_fn(dispatch_get_main_queue(), std::move(fn));
				//fn();
			}
		}
		
		
	public:
		static void construct(
			t_wt &wt,
			sdsl::int_vector_buffer <t_wt::tree_strat_type::int_width> &input_buf,
			typename t_wt::size_type const size,
			sdsl::bit_vector &temp_bv,
			construct_wt_gn_delegate <t_wt> *delegate
		)
		{
			assert(delegate);
			
			construct_wt_gn *cwt(new construct_wt_gn(wt, delegate));
			
			// temp_bv is allocated on the stack, so swap.
			using std::swap;
			swap(cwt->m_target_bv, temp_bv);
			
			cwt->construct(input_buf, size, (2 << (cwt->m_tau + 1)) - 1);
		}
		
		construct_wt_gn() = delete;
	};
}

#endif
