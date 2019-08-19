
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_QUERY_HPP
#define BGEN_QUERY_HPP

#include <boost/function.hpp>
#include <boost/optional.hpp>
#include <stdint.h>
#include <vector>
#include <string>

namespace genfile {
	namespace bgen {
		// Base class representing a query against a BGEN file or index
		struct Query {
		public:
			// We use std::unique_ptr to avoid using C++11 features here.
			typedef uint8_t byte_t ;
			typedef std::unique_ptr< Query > UniquePtr ;
			
			static UniquePtr create() ;

		public:
			struct GenomicRange {
				GenomicRange(): m_start(0), m_end(0) {}
				GenomicRange( GenomicRange const& other ):
					m_chromosome( other.m_chromosome ),
					m_start( other.m_start ),
					m_end( other.m_end )
				{}
					
				GenomicRange& operator=( GenomicRange const& other ) {
					m_chromosome = other.m_chromosome ;
					m_start = other.m_start ;
					m_end = other.m_end ;
					return *this ;
				}
				GenomicRange(
					std::string const& chromosome,
					uint32_t start,
					uint32_t end
				):
					m_chromosome( chromosome ),
					m_start( start ),
					m_end( end )
				{
					if( m_end < m_start ) {
						throw std::invalid_argument( "end" ) ;
					}
				}

				std::string const& chromosome() const { return m_chromosome ; }
				uint32_t const& start() const { return m_start ; }
				uint32_t const& end() const { return m_end ; }
				
			private:
				std::string m_chromosome ;
				uint32_t m_start ;
				uint32_t m_end ;
			} ;

		public:
			// Methods for building queries
			// Each method returns this object, allowing methods to be chained
			Query& include_range( GenomicRange const& range ) ;
			Query& exclude_range( GenomicRange const& range ) ;
			Query& include_rsids( std::vector< std::string > const& ids ) ;
			Query& exclude_rsids( std::vector< std::string > const& ids ) ;
			
			std::vector< std::string > const& included_rsids() const { return m_included_rsids ; }
			std::vector< std::string > const& excluded_rsids() const { return m_excluded_rsids ; }
			std::vector< GenomicRange > const& included_ranges() const { return m_included_ranges ; }
			std::vector< GenomicRange > const& excluded_ranges() const { return m_excluded_ranges ; }
			
		private:
			std::vector< std::string > m_included_rsids ;
			std::vector< std::string > m_excluded_rsids ;
			std::vector< GenomicRange > m_included_ranges ;
			std::vector< GenomicRange > m_excluded_ranges ;
		} ;
	}
}

#endif
