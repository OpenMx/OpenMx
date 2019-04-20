
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_ZLIB_HPP
#define GENFILE_ZLIB_HPP

#include <vector>
#include <stdint.h>
#include <cassert>
#include <zlib.h>
#include "zstd.h"
#include "genfile/types.hpp"

namespace genfile {

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).
	//
	// If offset is nonzero, compressed data will be written starting at position [offset].
	// The first [offset] bytes will be untouched.
	void zstd_compress(
		byte_t const* buffer,
		byte_t const* const end,
		std::vector< byte_t >* dest,
		std::size_t const offset = 0,
		int const compressionLevel = 22
	) ;

	// Compress the given data into the given destination buffer.  The destination will be resized
	// to fit the compressed data.  (Since the capacity of dest may be larger than its size,
	// to save memory you may need to copy the contents of dest elsewhere after calling
	// this function).

	template< typename T >
	void zstd_uncompress( byte_t const* begin, byte_t const* const end, std::vector< T >* dest ) {
		std::size_t const source_size = ( end - begin ) ;
		std::size_t const dest_size = dest->size() * sizeof( T ) ;
	    std::size_t const uncompressed_size = ZSTD_getDecompressedSize( reinterpret_cast< void const* >( begin ), source_size ) ;
		std::size_t const result = ZSTD_decompress(
			reinterpret_cast< void* >( &dest->operator[]( 0 ) ),
			dest_size,
			reinterpret_cast< void const* >( begin ),
			source_size
		) ;
		assert( result == uncompressed_size ) ;
		dest->resize( dest_size / sizeof( T )) ;
	}
}

#endif
