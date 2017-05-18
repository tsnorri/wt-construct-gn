/*
 * Copyright (c) 2016-2017 Tuukka Norri
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

#ifndef WT_CONSTRUCT_GN_DISPATCH_FN_HH
#define WT_CONSTRUCT_GN_DISPATCH_FN_HH

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <dispatch/dispatch.h>
#include <iostream>
#include <string>
#include <stdexcept>


namespace wtcgn { namespace detail {
	
	template <typename Fn>
	class dispatch_fn_context
	{
	public:
		typedef Fn function_type;
		
	protected:
		function_type m_fn;
		
	public:
		dispatch_fn_context(Fn &&fn):
			m_fn(std::move(fn))
		{
		}
		
		static void call_fn(void *dispatch_context)
		{
			assert(dispatch_context);
			auto *ctx(reinterpret_cast <dispatch_fn_context *>(dispatch_context));
			
			try
			{
				ctx->m_fn();
			}
			catch (std::exception const &exc)
			{
				std::cerr << "Caught exception: " << exc.what() << std::endl;
			}
			catch (...)
			{
				std::cerr << "Caught non-std::exception." << std::endl;
			}
			
			delete ctx;
		}
	};
	
	
	template <typename t_owner, void(t_owner::*t_fn)()>
	void call_member_function(void *ctx)
	{
		t_owner *owner(static_cast <t_owner *>(ctx));
		(owner->*t_fn)();
	}
}}


namespace wtcgn {
	
	// Allow passing pointer-to-member-function to dispatch_async_f without std::function.
	template <typename t_owner, void(t_owner::*t_fn)(void)>
	void dispatch_async_f(dispatch_queue_t queue, t_owner *obj)
	{
		dispatch_async_f(queue, obj, detail::call_member_function <t_owner, t_fn>);
	}

	template <typename t_owner, void(t_owner::*t_fn)()>
	void dispatch_barrier_async_f(dispatch_queue_t queue, t_owner *obj)
	{
		dispatch_barrier_async_f(queue, obj, detail::call_member_function <t_owner, t_fn>);
	}
	
	template <typename Fn>
	void dispatch_async_fn(dispatch_queue_t queue, Fn fn)
	{
		// A new expression doesn't leak memory if the object construction throws an exception.
		typedef detail::dispatch_fn_context <Fn> context_type;
		auto *ctx(new context_type(std::move(fn)));
		dispatch_async_f(queue, ctx, &context_type::call_fn);
	}
	
	template <typename Fn>
	void dispatch_barrier_async_fn(dispatch_queue_t queue, Fn fn)
	{
		typedef detail::dispatch_fn_context <Fn> context_type;
		auto *ctx(new context_type(std::move(fn)));
		dispatch_barrier_async_f(queue, ctx, &context_type::call_fn);
	}
	
	
	template <typename t_dispatch>
	class dispatch_ptr
	{
	protected:
		t_dispatch	m_ptr{};
		
	public:
		dispatch_ptr()
		{
		}
		
		dispatch_ptr(t_dispatch ptr, bool retain = false):
			m_ptr(ptr)
		{
			if (m_ptr && retain)
				dispatch_retain(m_ptr);
		}
		
		~dispatch_ptr()
		{
			if (m_ptr)
				dispatch_release(m_ptr);
		}
		
		dispatch_ptr(dispatch_ptr const &other):
			dispatch_ptr(other.m_ptr, true)
		{
		}
		
		dispatch_ptr(dispatch_ptr &&other):
			m_ptr(other.m_ptr)
		{
			other.m_ptr = nullptr;
		}
		
		bool operator==(dispatch_ptr const &other) const
		{
			return m_ptr == other.m_ptr;
		}
		
		bool operator!=(dispatch_ptr const &other) const
		{
			return !(*this == other);
		}
		
		dispatch_ptr &operator=(dispatch_ptr const &other) &
		{
			if (*this != other)
			{
				if (m_ptr)
					dispatch_release(m_ptr);
				
				m_ptr = other.m_ptr;
				dispatch_retain(m_ptr);
			}
			return *this;
		}
		
		dispatch_ptr &operator=(dispatch_ptr &&other) &
		{
			if (*this != other)
			{
				m_ptr = other.m_ptr;
				other.m_ptr = nullptr;
			}
			return *this;
		}
		
		t_dispatch operator*() { return m_ptr; }
	};
	
	
	template <typename t_dispatch>
	void swap(dispatch_ptr <t_dispatch> &lhs, dispatch_ptr <t_dispatch> &rhs)
	{
		using std::swap;
		swap(lhs.m_ptr, rhs.m_ptr);
	}
}

#endif
