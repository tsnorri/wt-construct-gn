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

#ifndef WT_CONSTRUCT_GN_CONSTRUCT_WT_HH
#define WT_CONSTRUCT_GN_CONSTRUCT_WT_HH

#include <iostream>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/ram_fs.hpp>
#include <sdsl/util.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <string>
#include <vector>
#include <wt-construct-gn/bit_ops.hh>
#include <wt-construct-gn/bit_vector_container.hh>
#include <wt-construct-gn/dispatch_fn.hh>
#include <wt-construct-gn/movable_atomic.hh>
#include <wt-construct-gn/packed_list.hh>


namespace wtcgn
{
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
	
	
	template <typename t_word>
	class small_node
	{
	public:
		enum { T_WORD_BIT = CHAR_BIT * sizeof(t_word) };
		
	protected:
		sdsl::bit_vector	m_content{};
		std::size_t			m_idx{0};
		
	public:
		small_node() = default;
		
		small_node(std::size_t const count):
			m_content(count)
		{
		}
		
		sdsl::bit_vector const &bit_vector() { return m_content; }
	
		// Append to the high end of each word.
		void append(t_word const word, std::size_t const bit_count)
		{
			assert(bit_count);
			
			// Find the position to which the bits should be placed in output_words.
			auto const word_idx(m_idx / BIT_VECTOR_BIT);
			auto const offset(m_idx % BIT_VECTOR_BIT);
			
			assert(m_content.bit_size());
			assert(word_idx <= m_content.bit_size() / BIT_VECTOR_BIT);
			
			// Store the bits.
			auto *content_data(m_content.data());
			content_data += word_idx;
			sdsl::std_size_type_for_int_vector shifted(word);
			shifted <<= offset;
			*content_data |= shifted;
			
			auto const remaining(BIT_VECTOR_BIT - offset);
			if (remaining < bit_count)
			{
				// Get the remaining bits.
				t_word shifted(word >> remaining);
				++content_data;
				*content_data |= shifted;
			}
			
			// Update the position.
			m_idx += bit_count;
		}
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
		typedef small_node <t_word>	small_node_type;
		
	protected:
		typedef sdsl::int_vector <TREE_NODE_BIT> shortcut_vec_type;
		typedef typename t_wt::tree_strat_type::node_type tree_node_type;
		
	protected:
		t_wt								*m_wt{nullptr};
		creation_context					m_ctx{};
		sdsl::int_vector <T_WORD_BIT>		m_C;
		sdsl::bit_vector					m_target_bv;
		std::vector <shortcut_vec_type>		m_node_paths;
		item_masks <t_word>					m_small_node_masks;
		t_word								m_tau_mask{};
		std::size_t							m_tau{};
		std::size_t							m_sigma_bits{};
		std::size_t							m_lowest_tree_level{};
		construct_wt_gn_delegate <t_wt>		*m_delegate{nullptr};

	protected:
		// Mutable part when creating nodes.
		class creation_context
		{
		protected:
			std::function <void(void)>						m_callback{};
			construct_wt_gn const							*m_construct{nullptr};
			bit_vector_container							m_target_vector{nullptr};
			movable_atomic <uint32_t>						m_operation_count{0};
			std::size_t										m_alpha{0};
			std::size_t										m_node_range_start{0};
			std::size_t										m_node_range_end{0};
			dispatch_ptr <dispatch_queue_t>					m_concatenation_queue;
			std::vector <dispatch_ptr <dispatch_queue_t>>	m_queues;
			std::vector <big_node_type>						m_big_nodes_1;
			std::vector <big_node_type>						m_big_nodes_2;
			std::vector <small_node_type>					m_small_nodes;
			
			
		public:
			creation_context() = default;
			

			void print_tree(typename t_wt::tree_strat_type::node_type const node, std::size_t const level)
			{
				auto const &val(m_construct->m_wt->m_tree.m_nodes[node]);
				for (std::size_t i(0); i < level; ++i)
					std::cerr << "    ";
				std::cerr << "[" << node << "] bv_pos: " << val.bv_pos << " bv_pos_rank: " << val.bv_pos_rank << " size: " << m_construct->m_wt->m_tree.size(node) << std::endl;
				
				
				if (!m_construct->m_wt->m_tree.is_leaf(node))
				{
					print_tree(m_construct->m_wt->m_tree.child(node, 0), 1 + level);
					print_tree(m_construct->m_wt->m_tree.child(node, 1), 1 + level);
				}
			}
			
			
			creation_context(construct_wt_gn const &construct, std::size_t const count, sdsl::bit_vector &target_bv, std::function <void(void)> &&callback):
				m_callback(std::move(callback)),
				m_construct(&construct),
				m_target_vector(&target_bv),
				m_concatenation_queue(dispatch_queue_create(nullptr, DISPATCH_QUEUE_SERIAL), false),
				m_queues(count),
				m_small_nodes(construct.m_wt->m_tree.m_nodes.size())
			{
				// Create the dispatch queues.
				// FIXME: Not this many queues are actually needed; one for each internal node in the last tau levels are enough
				// since each set of tau levels is processed serially.
				// FIXME: what about space complexity?
				if (t_multithreaded)
				{
					for (std::size_t i(0); i < count; ++i)
					{
						dispatch_ptr <dispatch_queue_t> queue(dispatch_queue_create(nullptr, DISPATCH_QUEUE_SERIAL), false);
						m_queues[i] = std::move(queue);
					}
				}
				
				// Big nodes.
				{
					// Find the lowest level on which big nodes are required.
					auto const sigma_bits(m_construct->m_sigma_bits);
					auto const tau(m_construct->m_tau);
					
					auto const lowest_level(m_construct->m_lowest_tree_level);
					auto const lowest_big_node_level(lowest_level / tau * tau);
					auto const max_big_nodes(1 << lowest_big_node_level);
					m_big_nodes_1.resize(max_big_nodes);
					m_big_nodes_2.resize(max_big_nodes);
				}
				
				// Small nodes.
				// Skip the last tree node.
				for (std::size_t i(0), count(construct.m_wt->m_tree.m_nodes.size() - 1); i < count; ++i)
				{
					auto const bit_count(m_construct->m_wt->m_tree.size(i));
					small_node_type temp_node(bit_count);
					m_small_nodes[i] = std::move(temp_node);
				}
				
				//print_tree(m_construct->m_wt->m_tree.root(), 0);
			}
			
			
			template <typename t_fn>
			inline void call_with_queue(t_fn &&fn, dispatch_queue_t queue)
			{
				dispatch_async_fn(queue, std::move(fn));
			}
			
			
			template <typename t_fn>
			inline void queue_concatenate(t_fn &&fn)
			{
				if (!t_multithreaded)
					fn();
				else
					call_with_queue(std::move(fn), *m_concatenation_queue);
			}

			
			template <typename t_fn>
			inline void queue_small_node_creation(t_fn &&fn, std::size_t const node_idx)
			{
				if (!t_multithreaded)
					fn();
				else
				{
					assert(node_idx < m_queues.size());
					call_with_queue(std::move(fn), *m_queues[node_idx]);
				}
			}
			
			
			template <typename t_fn>
			inline void queue_big_node_creation(t_fn &&fn)
			{
				if (!t_multithreaded)
					fn();
				else
					call_with_queue(std::move(fn), dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
			}

			
			template <typename t_fn>
			inline void queue_find_node_range(t_fn &&fn)
			{
				if (!t_multithreaded)
					fn();
				else
					call_with_queue(std::move(fn), dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0));
			}

			
			void finish_round()
			{
				bool stop(false);
				
				// Make sure that this function will not be called while still processing the previous round.
				// This works b.c. there will not be any parallel operations before some are initiated
				// from within this function.
				++m_operation_count;

				auto const tau(m_construct->m_tau);
				auto const sigma(m_construct->m_wt->m_sigma);
				auto const prev_alpha(m_alpha);
				++m_alpha;
				
				// Queue concatenating the constructed bit vectors.
				{
					auto const range_start(m_node_range_start);
					auto const range_end(m_node_range_end);
					auto const alpha(m_alpha);
					
					auto fn = [this, range_start, range_end, alpha]() mutable {
						for (std::size_t i(range_start); i < range_end; ++i)
							m_target_vector.append(m_small_nodes [i].bit_vector());
					};
					
					queue_concatenate(std::move(fn));
				}
				
				// Check if additional nodes are needed.
				if (m_alpha * tau < m_construct->m_lowest_tree_level)
				{
					// m_big_nodes_2 was prepared in construct or previous call to finish_and_check.
					using std::swap;
					swap(m_big_nodes_1, m_big_nodes_2);

					// Count the number of non-empty big nodes.
					std::size_t non_empty_count(0);
					for (auto const &big_node : m_big_nodes_1)
					{
						if (big_node.characters().size())
							++non_empty_count;
					}

					// Create the next set of big nodes if needed.
					if ((m_alpha + 1) * tau < m_construct->m_lowest_tree_level)
					{
						// Increment the operation count.
						m_operation_count += non_empty_count;
						
						// Create the new nodes.
						remove_vector_items(m_big_nodes_2);
						for (auto const &parent_node : m_big_nodes_1)
						{
							// If the node has characters, create child nodes.
							auto const characters(parent_node.characters());
							if (characters.size())
							{
								auto fn = [this, &parent_node]() mutable {
									create_big_nodes(parent_node);
								};
								queue_big_node_creation(std::move(fn));
							}
						}
					}
					
					// Create small nodes with the last set of big nodes.
					m_operation_count += non_empty_count;
					for (auto const &big_node : m_big_nodes_1)
					{
						auto const characters(big_node.characters());
						if (characters.size())
						{
							auto fn = [this, &big_node]() mutable {
								create_small_node_initial(big_node);
							};
							
							queue_small_node_creation(std::move(fn), big_node.node_index());
						}
					}
					
					// Update the node range.
					{
						++m_operation_count;
						auto fn = [this]() mutable {
							m_node_range_start = m_node_range_end;
							auto const range_end(m_construct->find_next_node_range_end(m_node_range_end));
							m_node_range_end = range_end;
							finish_and_check();
						};
						queue_find_node_range(std::move(fn));
					}
				}
				else
				{
					stop = true;
					auto fn = [this](){
						m_callback();
					};
					
					queue_concatenate(std::move(fn));
				}
				
				if (!stop)
					finish_and_check();
			}
			
			
			inline void finish_and_check()
			{
				auto const remaining_count(--m_operation_count);
				if (0 == remaining_count)
					finish_round();
			}
			
			
			template <bool t_right>
			inline void call_create_small_node(
				big_node_type const &nearest_big_node,
				t_word const w,
				t_word const mask,
				std::size_t const level,	// Actual level % m_tau.
				std::size_t const node_idx
			)
			{
				if (mask)
				{
					auto const child_idx(m_construct->m_wt->m_tree.child(node_idx, t_right));
					if (!m_construct->m_wt->m_tree.is_leaf(child_idx))
					{
						if (t_multithreaded)
							++m_operation_count;

						auto fn = [this, &nearest_big_node, w, mask, level, child_idx]() mutable {
							create_small_node <t_right>(nearest_big_node, w, mask, 1 + level, child_idx);
						};
						queue_small_node_creation(std::move(fn), child_idx);
					}
				}
			}
			
			
			template <bool t_left>
			void create_small_node(
				big_node_type const &nearest_big_node,
				t_word const w,
				t_word const mask,
				std::size_t const level,	// Actual level % m_tau.
				std::size_t const node_idx
			)
			{
				small_node_type &output_node(m_small_nodes[node_idx]);
				auto const bit_count(count_set_bits(mask));
				assert(bit_count);

				// extracted will contain the bits in places indicated by mask.
				t_word const extracted(extract_bits(w, mask));

				// Store the bits.
				output_node.append(extracted, bit_count);
			
				if (1 + level < m_construct->m_tau)
				{
					// Create masks for the next level.
					t_word const left_mask((~w & mask) << 1);	// Find zeros in masked positions.
					t_word const right_mask((w & mask) << 1);	// Find ones in masked positions.
				
					call_create_small_node <false>(nearest_big_node, w, left_mask, level, node_idx);
					call_create_small_node <true>(nearest_big_node, w, right_mask, level, node_idx);
				}

				if (t_multithreaded)
					finish_and_check();
			}
	
	
			// Create a list of words for small nodes.
			void create_small_node_initial(big_node_type const &big_node)
			{
				auto const node_idx(big_node.node_index());
				auto const tau(m_construct->m_tau);
				auto const level(m_alpha * tau);
				
				auto const mask_idx(m_construct->m_sigma_bits - level - 1);
				t_word mask(m_construct->m_small_node_masks[mask_idx]);
				auto bit_count(count_set_bits(mask));
				
				auto const &characters(big_node.characters());
				auto const &word_vector(characters.word_vector());
				auto const word_count(word_vector.size());
				auto const character_count(characters.size());
				auto const character_bits(characters.item_bits());
				
				auto &output_node(m_small_nodes[node_idx]);
				auto *output_words(output_node.bit_vector().data());
			
				for (std::size_t i(0); i < word_count; ++i)
				{
					auto const w(word_vector[i]);
					
					// Change mask for the last iteration.
					// FIXME: there's bound to be a better way to do this.
					if (i == word_count - 1)
					{
						t_word last_character_mask(~0);
						auto const diff(character_count / word_count * character_bits);
						last_character_mask >>= T_WORD_BIT - diff;
						mask &= last_character_mask;
						bit_count = count_set_bits(mask);
					}
					
					// extracted will contain the bits in places indicated by mask.
					t_word const extracted(extract_bits(w, mask));
				
					// Store the bits.
					output_node.append(extracted, bit_count);
				
					// Create masks for the next level.
					t_word const left_mask((~w & mask) << 1);	// Find zeros in masked positions.
					t_word const right_mask((w & mask) << 1);	// Find ones in masked positions.
				
					call_create_small_node <false>(big_node, w, left_mask, 0, node_idx);
					call_create_small_node <true>(big_node, w, right_mask, 0, node_idx);
				}
			
				finish_and_check();
			}
			
			
			// Similar to below.
			// append() contains a loop that should be unrolled.
			__attribute__((optimize("unroll-loops")))
			void create_big_nodes(big_node_type const &parent_node)
			{
				auto const tau(m_construct->m_tau);
				auto const tau_mask(m_construct->m_tau_mask);
				auto const &characters(parent_node.characters());
				auto const &paths(m_construct->m_node_paths[m_alpha]);
				auto const character_bits(characters.item_bits());
				auto const child_character_bits(character_bits - tau);
				auto const level(m_alpha * tau);
		
				for (std::size_t i(0), count(characters.size()); i < count; ++i)
				{
					auto const c(characters.value(i));
			
					// Get the last tau bits of the character. (WT has LSB order.)
					// Child suffix will have the bits in the correct order b.c. the parent's suffix
					// has the lowest bits.
					t_word const child_idx(c & tau_mask);
					t_word const child_suffix_part(child_idx << level);
					t_word const parent_suffix(parent_node.suffix());
					t_word const child_suffix(parent_suffix | child_suffix_part);
			
					// Get the vector for the child in question.
					auto &child_node(m_big_nodes_2[child_suffix]);
					if (child_node.is_leaf())
						continue;
			
					// Check that there is enough space.
					auto &child_node_characters(child_node.characters());
					if (0 == child_node_characters.size())
					{
						// Child node has not been set up; assign values.
						auto const child_suffix_length(parent_node.suffix_length() + tau);
						
						child_node.set_suffix(child_suffix);
						child_node.set_suffix_length(child_suffix_length);
						
						// Count the number of characters in the current node.
						auto const tree_node(paths[child_suffix]);
						auto const child_character_count(m_construct->m_wt->m_tree.size(tree_node));
						
						child_node.set_node_index(tree_node);
						
						bool is_leaf(m_construct->m_wt->m_tree.is_leaf(tree_node));
						child_node.set_leaf(is_leaf);
						if (is_leaf)
							continue;

						// Replace the empty list with an initialized one.
						typename big_node_type::character_list temp_list(child_character_count, child_character_bits);
						child_node_characters = std::move(temp_list);
					}
			
					auto const cc(c >> tau);
					child_node_characters.append(cc);
				}
				
				finish_and_check();
			}
		
		
			// Create a list of words each of which contains as many characters as possible s.t. no characters are split between words
			// and the characters start from bit index zero (0-based).
			// append() contains a loop that should be unrolled.
			__attribute__((optimize("unroll-loops")))
			void create_big_node_initial(
				sdsl::int_vector_buffer <t_wt::tree_strat_type::int_width> &input_buf,
				std::size_t const input_size,
				big_node_type &output_node
			) const
			{
				auto &characters(output_node.characters());
		
				{
					auto const sigma_bits(m_construct->m_sigma_bits);
					assert(sigma_bits);
					typename big_node_type::character_list temp_list(input_size, sigma_bits);
					characters = std::move(temp_list);
				}
		
				for (std::size_t i(0); i < input_size; ++i)
				{
					auto const c(input_buf[i]);
					t_word comp(m_construct->m_wt->m_tree.bit_path(c));
					comp &= 0xFFFFFFFFFFFFFF; // Set the 8 most significant bits to zero.
					characters.append(comp);
				}
			}
			
			
			void construct(
				sdsl::int_vector_buffer <t_wt::tree_strat_type::int_width> &input_buf,
				typename t_wt::size_type const size
			)
			{
				// Create the initial big node.
				create_big_node_initial(input_buf, size, m_big_nodes_1[0]);
			
				// Increment for the three operations below.
				m_operation_count += 3;
				
				{
					// Create small nodes for the first set of levels.
					auto fn = [this]() mutable {
						create_small_node_initial(m_big_nodes_1[0]);
					};
					
					queue_small_node_creation(std::move(fn), 0);
				}
				
				{
					// Check if additional big nodes are needed and create accordingly.
					auto const tau(m_construct->m_tau);
					if (tau < m_construct->m_lowest_tree_level)
					{
						remove_vector_items(m_big_nodes_2);
						auto fn = [this]() mutable {
							create_big_nodes(m_big_nodes_1[0]);
						};
						queue_big_node_creation(std::move(fn));
					}
				}
				
				{
					// Update the node range.
					auto fn = [this]() mutable {
						auto const range_end(m_construct->find_next_node_range_end(0));
						m_node_range_end = range_end;
						finish_and_check();
					};
					queue_find_node_range(std::move(fn));
				}
			}
		};
		
		friend class creation_context;
		
		
	protected:
		construct_wt_gn(t_wt &wt, construct_wt_gn_delegate <t_wt> *delegate):
			m_wt(&wt),
			m_tau_mask(~0),
			m_tau(delegate->tau(wt)),
			m_delegate(delegate)
		{
			m_tau_mask >>= (T_WORD_BIT - m_tau);
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
		
		
		std::size_t find_next_node_range_end(std::size_t const start) const
		{
			auto const size(m_wt->m_tree.m_nodes.size());
			std::size_t level(0);
			std::size_t max_parent(start - 1);
			std::size_t i(start);
			
			if (0 == i)
			{
				++level;
				++i;
				max_parent = 0;
				if (level == m_tau)
					return i;
			}
			
			while (i < size)
			{
				// Check the parent of the current node. If the parent
				// is greater than max_parent, increase level and set
				// max_parent to the last node of the previous level
				// (which is the node that precedes the current one).
				auto const parent(m_wt->m_tree.parent(i));
				if (max_parent < parent)
				{
					++level;
					max_parent = i - 1;
					if (level == m_tau)
						break;
				}
				++i;
			}
			
			return i;
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
			typename t_wt::size_type const size
		)
		{
			// Check how many characters may be packed into a word.
			m_sigma_bits = gte_exp_base_2(m_wt->m_sigma);
			assert(m_sigma_bits);
			
			m_lowest_tree_level = find_lowest_tree_level();
			
			// Shortcuts for m_wt->m_tree.m_nodes.
			auto const nonzero_big_node_levels(m_lowest_tree_level / m_tau);
			if (nonzero_big_node_levels)
			{
				auto const level_count(nonzero_big_node_levels);
				auto const initial_item_count(1 << m_tau);
				m_node_paths.reserve(level_count);
				auto const root_node(m_wt->m_tree.root());
				
				shortcut_vec_type initial_suffix(1, 0);
				construct_shortcut_level(0, level_count, initial_item_count, initial_item_count, initial_suffix);
			}
			
			// Create bit masks for accessing the first bits of characters.
			{
				decltype(m_small_node_masks) temp_masks(m_sigma_bits);
				m_small_node_masks = std::move(temp_masks);
			}
			
			// Node creation context.
			{
				auto callback = [this](){
					// Move the bit vector to the WT.
					m_wt->m_bv = typename t_wt::bit_vector_type(std::move(m_target_bv));
					m_wt->finish_construct();
					
					auto *delegate(m_delegate);
					
					// The instance was allocated with new.
					delete this;
					
					// Get out of WT's constructor in single-threaded mode.
					auto fn = [delegate](){
						delegate->finish_construction();
					};
					dispatch_async_fn(dispatch_get_main_queue(), std::move(fn));
				};
				
				// FIXME: not this many queues are actually needed.
				creation_context temp_ctx(*this, m_wt->m_tree.m_nodes.size(), m_target_bv, std::move(callback));
				m_ctx = std::move(temp_ctx);
			}
			
			m_ctx.construct(input_buf, size);
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
			
			cwt->construct(input_buf, size);
		}
		
		construct_wt_gn() = delete;
	};
}

#endif
