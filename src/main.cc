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

#include "cmdline.h"
#include <sdsl/io.hpp>
#include <wt-construct-gn/construct_wt_2.hh>
#include <wt-construct-gn/file_handling.hh>
#include <wt-construct-gn/timer.hh>
#include <wt-construct-gn/types.hh>

#ifdef __linux__
#include <pthread_workqueue.h>
#endif


namespace {
	
	inline void stop_timer(wtcgn::timer &timer)
	{
		timer.stop();
		std::cerr << "Done. Milliseconds elapsed: " << timer.ms_elapsed() << std::endl;
	}

	
	template <typename t_wt>
	inline void write_wt_to_file(t_wt const &wt, wtcgn::file_ostream &stream)
	{
		std::cerr << "Writing the Wavelet tree to a file…" << std::endl;
		wtcgn::timer timer;
		sdsl::serialize(wt, stream);
		stop_timer(timer);
	}
	
	
	struct construct_ctx_base
	{
		virtual ~construct_ctx_base() {}
		virtual void construct() = 0;
		virtual void set_input_file_name(char const *fname) = 0;
		virtual void set_input_width(uint8_t const width) = 0;
		virtual void open_output_file(char const *fname) = 0;
		virtual wtcgn::timer &wt_construction_timer() = 0;
		virtual void set_tau(std::size_t const tau) = 0;
		virtual void set_print_wt(bool flag) = 0;
	};

	
	template <typename t_wt>
	class construct_ctx : public wtcgn::construct_wt_gn_delegate <t_wt>, public construct_ctx_base
	{
	public:
		
	protected:
		std::unique_ptr <t_wt>	m_wt;
		wtcgn::timer			m_wt_construction_timer;
		std::string				m_input_file_name;
		wtcgn::file_ostream		m_output_stream;
		std::size_t				m_tau{0};
		uint8_t					m_input_width{0};
		bool					m_output_wt{false};
		bool					m_print_wt{false};
		
	public:
		construct_ctx()
		{
		}
		
		
		virtual ~construct_ctx() {}
		
		
		void cleanup()
		{
			delete this;
		}
		
		
		virtual wtcgn::timer &wt_construction_timer() { return m_wt_construction_timer; }
		virtual void set_tau(std::size_t const tau) { m_tau = tau; }
		virtual void set_print_wt(bool const flag) { m_print_wt = flag; }
		
		
		virtual void set_input_file_name(char const *fname)
		{
			m_input_file_name = fname;
		}

		
		virtual void set_input_width(uint8_t const width)
		{
			m_input_width = width;
		}

		
		virtual void open_output_file(char const *fname)
		{
			wtcgn::open_file_for_writing(fname, m_output_stream, false);
			m_output_wt = true;
		}
		
		
		virtual void construct()
		{
			sdsl::int_vector_buffer <t_wt::alphabet_category::WIDTH>
			//text_buf(m_input_file_name, std::ios::in, 1024 * 1024, t_wt::alphabet_category::WIDTH, true);
			text_buf(m_input_file_name, std::ios::in, 1024 * 1024, m_input_width, true);

			m_wt.reset(new t_wt(text_buf, text_buf.size(), this));
		}
		
		
		virtual std::size_t tau(t_wt const &wt)
		{
			// Tau should be the largest power of two not greater than sqrt(w / log w).
			// sqrt(64 / log 64) = sqrt(64 / 6) = sqrt(10 2/3) ≈ 3.26.
			return (m_tau ?: 2);
		}
		
		
		virtual void finish_construction()
		{
			stop_timer(m_wt_construction_timer);
			
			if (m_print_wt)
			{
				std::cerr << "WT size: " << m_wt->size() << std::endl;
				if (8 == m_input_width)
				{
					for (std::size_t i(0); i < m_wt->size(); ++i)
						std::cout << (char)((*m_wt)[i]);
				}
				else
				{
					for (std::size_t i(0); i < m_wt->size(); ++i)
						std::cout << (*m_wt)[i];
				}
				std::cout << std::endl;
			}
			
			if (m_output_wt)
				write_wt_to_file(*m_wt, m_output_stream);
			
			cleanup();
			// *this is no longer valid.
			exit(EXIT_SUCCESS);
		}
	};
}



int main (int argc, char **argv)
{
	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		exit(EXIT_FAILURE);
	
	std::ios_base::sync_with_stdio(false);	// Don't use C style IO after calling cmdline_parser.
	std::cin.tie(nullptr);					// We don't require any input from the user.
	
	if (args_info.sleep_given)
	{
		std::cerr << "Waiting for " << args_info.sleep_arg << " seconds…" << std::endl;
		sleep(args_info.sleep_arg);
	}
	
	uint8_t input_width(0);
	{
		switch (args_info.input_width_arg)
		{
			case input_width_arg_8:
				input_width = 8;
				break;

			case input_width_arg_16:
				input_width = 16;
				break;

			case input_width_arg_32:
				input_width = 32;
				break;

			case input_width_arg_64:
				input_width = 64;
				break;
				
			default:
				abort();
		}
	}
	
	if (args_info.sdsl_given)
	{
		wtcgn::timer timer;

		wtcgn::file_ostream output_stream;
		if (args_info.output_file_given)
			wtcgn::open_file_for_writing(args_info.output_file_arg, output_stream, false);
		
		sdsl::int_vector_buffer <wtcgn::wt_type_sdsl::alphabet_category::WIDTH>
		//text_buf(args_info.input_file_arg, std::ios::in, 1024 * 1024, wtcgn::wt_type_sdsl::alphabet_category::WIDTH, true);
		text_buf(args_info.input_file_arg, std::ios::in, 1024 * 1024, input_width, true);

		std::cerr << "Constructing the Wavelet tree…" << std::endl;
		wtcgn::wt_type_sdsl wt(text_buf, text_buf.size());
		stop_timer(timer);
		
		if (args_info.output_file_given)
		{
			assert(output_stream.good());
			write_wt_to_file(wt, output_stream);
		}
	}
	else if (args_info.bitparallel_given)
	{
		construct_ctx_base *ctx{nullptr};
		
		if (true || args_info.no_mt_given)
			ctx = new construct_ctx <wtcgn::wt_type_gn>();
		else
		{
			// libdispatch on macOS does not need pthread_workqueue.
#ifdef __linux__
			pthread_workqueue_init_np();
#endif
			
			ctx = new construct_ctx <wtcgn::wt_type_gn_mt>();
		}
		
		// Set input file name and open output stream.
		ctx->set_input_file_name(args_info.input_file_arg);
		if (args_info.output_file_given)
			ctx->open_output_file(args_info.output_file_arg);
		
		ctx->set_input_width(input_width);
		
		if (args_info.tau_given)
		{
			auto const tau(args_info.tau_arg);
			if (! (0 < tau))
			{
				std::cerr << "Tau must be positive." << std::endl;
				exit(EXIT_FAILURE);
			}
			ctx->set_tau(tau);
		}
		
		if (args_info.print_wt_given)
			ctx->set_print_wt(true);
		
		// Construct the WT.
		std::cerr << "Constructing the Wavelet tree…" << std::endl;
		auto fn = [ctx](){
			ctx->construct();
		};
		
		wtcgn::dispatch_async_fn(dispatch_get_main_queue(), std::move(fn));
		
		// Calls pthread_exit.
		dispatch_main();
	}
	else
	{
		// Unexpected mode.
		abort();
	}
	
	// Not reached from dispatch_main.
	return EXIT_SUCCESS;
}
