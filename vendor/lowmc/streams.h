/** Sabre, an anonymous bulletin board with speedier ripostes
 *  Copyright (C) 2020  Sabre authors
 *
 *  @file    streams.h
 *  @brief   
 *
 *  @author  Ryan Henry        <ryan.henry@ucalgary.ca>
 *  @author  Adithya Vadapalli <avadapal@iu.edu>
 *  @author  Kyle Storrier     <kyle.storrier@ucalgary.ca>
 *
 *  @license GNU Public License (version 2); see LICENSE for full license text
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License along
 *    with this program; if not, write to the Free Software Foundation, Inc.,
 *    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 **/

#ifndef LOWMC_STREAMS_H__
#define LOWMC_STREAMS_H__

#include <type_traits>   // std::conditional
#include <array>         // std::array
#include <queue>         // std::queue
#include <mutex>         // std::mutex
#include <vector>        // std::vector

//#include "openssl/sha.h" // OpenSSL's SHA256 implementation

#include "transposition.h"
#include <boost/asio.hpp>

namespace lowmc
{

namespace streams
{

using boost::asio::ip::tcp;

template <typename block_t>
class input_stream
{
  public:
	input_stream() = default;
	input_stream(const input_stream &) = delete;
	~input_stream() = default;

	input_stream & operator=(const input_stream &) = delete;
	input_stream & operator=(input_stream &&) = default;

	virtual input_stream & operator>>(block_t & val) = 0; 
	//virtual input_stream & operator>>(bool & val); 
}; // class lowmc::streams::input_stream

template <typename block_t, size_t slices = 1>
class output_stream
{
  public:
	output_stream() = default;
	output_stream(const output_stream &) = delete;
	~output_stream() = default;

	output_stream & operator=(const output_stream &) = delete;
	output_stream & operator=(output_stream &&) = default;

	virtual output_stream & operator<<(const block_t & val) = 0; 
	//virtual input_stream & operator<<(bool & val);
 
}; // class lowmc::streams:: output_stream

// template <typename block_t, size_t slices = 1>
// class hashed_output_stream : public output_stream<block_t>
// {
//   public:
// 	hashed_output_stream()
// 	{
// 		for (size_t i = 0; i < slices; ++i)
// 		{
// 			if (!SHA256_Init(&contexts[i])) { /*gah!*/ }
// 		}
// 	}
// 	hashed_output_stream(const hashed_output_stream &) = delete;
// 	~hashed_output_stream() = default;

// 	hashed_output_stream & operator=(const hashed_output_stream &) = delete;
// 	hashed_output_stream & operator=(hashed_output_stream &&) = default;

// 	hashed_output_stream & operator<<(const block_t & val) override
// 	{
// 		auto transposed = transpose(val);
// 		for (size_t i = 0; i < slices; ++i)
// 		{
// 			if (!SHA256_Update(&contexts[i], &transposed[i],
// 				sizeof(transposed[i]))) { /*gah!*/ }
// 		}
// 		handle_read(val);
// 		return *this;
// 	}

// 	inline __m256i sha256_digest(size_t i = 0)
// 	{
// 		__m256i hash;
// 		// TODO: check return value; "invalidate" this stream
// 		SHA256_Final(reinterpret_cast<unsigned char*>(&hash), &contexts[i]);
// 		return hash;
// 	}

//   protected:
// 	virtual void handle_read(const block_t & val) {  }

//   private:
// 	SHA256_CTX contexts[slices];
// }; // class lowmc::hashed_output_stream

// template <typename block_t, size_t slices = 1>
// class bidirectional_stream : public input_stream<block_t>, public hashed_output_stream<block_t, slices>
// {
//   public:
// 	bidirectional_stream() = default;
// 	bidirectional_stream(const bidirectional_stream &) = delete;
// 	~bidirectional_stream() = default;

// 	bidirectional_stream & operator=(const bidirectional_stream &) = delete;
// 	bidirectional_stream & operator=(bidirectional_stream &&) = default;
// }; // class lowmc::streams:bidirectional_stream

// template <typename block_t, size_t slices = 1>
// class basic_stream : public bidirectional_stream<block_t, slices>
// {
//   public:
// 	basic_stream() = default;
// 	basic_stream(const basic_stream &) = delete;
// 	basic_stream(basic_stream &&) = default;
// 	~basic_stream() = default;

// 	inline basic_stream & operator=(const basic_stream &) = delete;
// 	inline basic_stream & operator=(basic_stream &&) = default;

// 	inline basic_stream & operator>>(block_t & val) override
// 	{
// 		while (buffer.empty()) std::this_thread::yield();
// 		const std::lock_guard<std::mutex> lock(mutex);
// 		val = buffer.front();
// 		buffer.pop();
// 		return *this;
// 	}

//   protected:
// 	inline void handle_read(const block_t & val) override
// 	{
// 		const std::lock_guard<std::mutex> lock(mutex);
// 		buffer.push(val);
// 	}

//   private:
// 	std::queue<block_t> buffer;
// 	std::mutex mutex;
// }; // class lowmc::streams::basic_stream

// template <typename block_t>
// class input_socket_stream : public input_stream<block_t>
// {
//   public:
// 	input_socket_stream() = delete;
// 	input_socket_stream(const input_socket_stream &) = delete;
// 	input_socket_stream(input_socket_stream &&) = default;
// 	input_socket_stream(tcp::socket && sin, tcp::resolver & resolver, std::string host, std::string port)
// 	    : socket(std::move(sin))
// 	{
// 		boost::asio::connect(socket, resolver.resolve({host, port}));

// 		printf("connection established\n");
// 	}
// 	~input_socket_stream() = default;

// 	inline input_socket_stream & operator=(const input_socket_stream &) = delete;
// 	inline input_socket_stream & operator=(input_socket_stream &&) = default;

// 	inline input_socket_stream & operator>>(block_t & val) override
// 	{
// 		boost::asio::read(socket, boost::asio::buffer(&val, sizeof(block_t)));
// 		printf("val written = %llu, %llu\n", val.mX[0], val.mX[1]);
// 		return *this;
// 	}
//  	inline input_socket_stream & operator>>(block<__m256i> & val)
// 	{
// 		boost::asio::read(socket, boost::asio::buffer(&val, sizeof(block<__m256i>)));
// 		printf("val written = %llu, %llu\n", val.mX[0], val.mX[1]);
// 		return *this;
// 	}

//   private:
// 	tcp::socket socket;
// }; // class lowmc::streams::input_socket_stream

// template <typename block_t>
// class output_socket_stream : public output_stream<block_t>
// {
//   public:
// 	output_socket_stream() = delete;
// 	output_socket_stream(tcp::socket && sout) : socket(std::move(sout)) { }
// 	~output_socket_stream() = default;
// 	output_socket_stream & operator=(const output_socket_stream &) = delete;
// 	inline output_socket_stream & operator=(output_socket_stream &&) = default;

// 	inline output_socket_stream & operator<<(const block_t & val) override
// 	{
// 		boost::asio::async_write(socket, boost::asio::buffer(&val, sizeof(block_t)), 
// 			[](boost::system::error_code, std::size_t){ });
// 		return *this;
// 	}

 


//   private:
// 	tcp::socket socket;
// }; // class lowmc::output_socket_stream


template <typename block_t>
class socket_stream : public input_stream<block_t>, public output_stream<block_t>
{
  public:
    socket_stream() = delete;
    socket_stream(const socket_stream &) = delete;
    socket_stream(socket_stream &&) = default;
    socket_stream(tcp::socket && sout) : socket(std::move(sout)) { }
    socket_stream(tcp::socket && sock, tcp::resolver & resolver, std::string host, std::string port)
        : socket(std::move(sock))
    {
        boost::asio::connect(socket, resolver.resolve({host, port}));
    }
    ~socket_stream() = default;

    inline socket_stream & operator=(const socket_stream &) = delete;
    inline socket_stream & operator=(socket_stream &&) = default;

    // inline socket_stream & operator>>(block_t & val) override
    // {
    //     boost::asio::read(socket, boost::asio::buffer(&val, sizeof(block_t)));
    //     return *this;
    // }




    // inline socket_stream & operator<<(const block_t & val) override
    // {
    //     //boost::asio::write(socket, boost::asio::buffer(&val, sizeof(block_t)));
    //     boost::asio::async_write(socket, boost::asio::buffer(&val, sizeof(block_t)), 
    //         [](boost::system::error_code, std::size_t){ });
    //     return *this;
    // }


	inline socket_stream & operator>>(uint8_t & val) 
    {
        boost::asio::read(socket, boost::asio::buffer(&val, sizeof(uint8_t)));
  
        return *this;
    }



    inline socket_stream & operator>>(block<__m256i> & val)
    {
        boost::asio::read(socket, boost::asio::buffer(&val, sizeof(block<__m256i>)));
        return *this;
    }

    inline socket_stream & operator>>(block<__m128i> & val)
    {
        boost::asio::read(socket, boost::asio::buffer(&val, sizeof(block<__m128i>)));
        return *this;
    }

    inline socket_stream & operator<<(const block<__m128i> & val)
    {
        boost::asio::write(socket, boost::asio::buffer(&val, sizeof(block<__m128i>)));
        // boost::asio::async_write(socket, boost::asio::buffer(&val, sizeof(block<__m256i>)), 
        //     [](boost::system::error_code, std::size_t){ });
        return *this;
    }
    
    inline socket_stream & operator<<(const block<__m256i> & val)
    {
        boost::asio::write(socket, boost::asio::buffer(&val, sizeof(block<__m256i>)));
        // boost::asio::async_write(socket, boost::asio::buffer(&val, sizeof(block<__m256i>)), 
        //     [](boost::system::error_code, std::size_t){ });
        return *this;
    }

    inline socket_stream & operator<<(const uint8_t & val)
    {
        //boost::asio::write(socket, boost::asio::buffer(&val, sizeof(block_t)));
      //  printf("bool-val = %d\n", val);
        boost::asio::async_write(socket, boost::asio::buffer(&val, sizeof(uint8_t)), 
            [](boost::system::error_code, std::size_t){ });
        return *this;
    }
  protected:
    tcp::socket socket;
};

//template <typename block_t, size_t slices = 1>
// class rewindable_stream : public bidirectional_stream<block_t, slices>
// {
//   public:
// 	rewindable_stream() : cur(0) { }
// 	rewindable_stream(const rewindable_stream &) = delete;
// 	rewindable_stream(rewindable_stream &&) = default;
// 	~rewindable_stream() = default;

// 	inline rewindable_stream & operator=(const rewindable_stream &) = delete;
// 	inline rewindable_stream & operator=(rewindable_stream &&) = default;

// 	inline rewindable_stream & rewind(size_t to = 0) { cur = to; return *this; }
// 	inline rewindable_stream & operator>>(block_t & val) override
// 	{
// 		while (cur >= buffer.size()) std::this_thread::yield();
// 		const std::lock_guard<std::mutex> lock(mutex);
// 		val = buffer[cur++];
// 		return *this;
// 	}

//   protected:
// 	inline void handle_read(const block_t & val) override
// 	{
// 		const std::lock_guard<std::mutex> lock(mutex);
// 		buffer.push_back(val);
// 	}

//   private:
// 	std::vector<block_t> buffer;
// 	std::mutex mutex;
// 	size_t cur;
// }; // class lowmc::streams::rewindable_stream

} // namespace lowmc::streams

} // namespace lowmc

#endif // LOWMC_STREAMS_H__