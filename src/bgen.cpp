
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <limits>
#include <climits>
#include <algorithm>
#include <iomanip>
#include "genfile/types.hpp"
#include "genfile/bgen/bgen.hpp"

#ifdef CHAR_BIT
#if (CHAR_BIT != 8)
#error CHAR_BIT "Sorry, this implementation assumes 8-bit bytes. It won't work on your platform"
#endif
#endif

namespace genfile {
	namespace bgen {
#if DEBUG_BGEN_FORMAT
		namespace impl {
			std::string to_hex( std::string const& str ) {
				std::ostringstream o ;
				for( std::size_t i = 0; i < str.size(); ++i ) {
					if( i % 4 == 0 )
						o << "|" ;
					o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> ( static_cast<unsigned char>( str[i] ) ) ;
				}
				return o.str() ;
			}
		}
#endif
		
		Context::Context():
			number_of_samples(0),
			number_of_variants(0),
			magic( "bgen" ),
			free_data( "" ),
			flags(0)
		{}
			
		Context::Context( Context const& other ):
			number_of_samples( other.number_of_samples ),
			number_of_variants( other.number_of_variants ),
			magic( other.magic ),
			free_data( other.free_data ),
			flags( other.flags )
		{}

		Context& Context::operator=( Context const& other ) {
			number_of_samples = other.number_of_samples ;
			number_of_variants = other.number_of_variants ;
			magic = other.magic ;
			free_data = other.free_data ;
			flags = other.flags ;
			return *this ;
		}
		
		uint32_t Context::header_size() const { return free_data.size() + 20 ; }
		
		void read_offset( std::istream& iStream, uint32_t* offset ) {
			read_little_endian_integer( iStream, offset ) ;
		}

		void write_offset( std::ostream& oStream, uint32_t const offset ) {
			write_little_endian_integer( oStream, offset ) ;
		}
		
		std::size_t read_header_block(
			std::istream& aStream,
			Context* context
		) {
			assert( context != 0 ) ;
			uint32_t
				header_size = 0,
				number_of_snp_blocks = 0,
				number_of_samples = 0,
				flags = 0 ;

			char magic[4] ;
			std::size_t fixed_data_size = 20 ;
			std::vector<char> free_data ;

			read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			read_little_endian_integer( aStream, &number_of_samples ) ;
			aStream.read( &magic[0], 4 ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( free_data.data(), free_data.size() ) ;
			read_little_endian_integer( aStream, &flags ) ;

			if(
				( magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e' || magic[3] != 'n' )
				&& ( magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0 )
			) {
				ThrowBGenError() ;
			}

			if( aStream ) {
				context->number_of_samples = number_of_samples ;
				context->number_of_variants = number_of_snp_blocks ;
				context->magic.assign( &magic[0], &magic[0] + 4 ) ;
				context->free_data.assign( free_data.begin(), free_data.end() ) ;
				context->flags = flags ;

				return( header_size ) ;
			} else {
				ThrowBGenError() ;
			}
		}

		void write_header_block(
			std::ostream& aStream,
			Context const& context
		) {
			uint32_t header_size = context.header_size() ;
			write_little_endian_integer( aStream, header_size ) ;
			write_little_endian_integer( aStream, context.number_of_variants ) ;
			write_little_endian_integer( aStream, context.number_of_samples ) ;
			aStream.write( context.magic.data(), 4 ) ;
			aStream.write( context.free_data.data(), context.free_data.size() ) ;
			write_little_endian_integer( aStream, context.flags ) ;
		}
		
		std::size_t write_sample_identifier_block(
			std::ostream& aStream,
			Context const& context,
			std::vector< std::string > const& sample_ids
		) {
			assert( sample_ids.size() == context.number_of_samples ) ;
			uint32_t block_size = 8 ;
			for( uint32_t i = 0; i < sample_ids.size(); ++i ) {
				block_size += 2 + sample_ids[i].size() ;
			}
			write_little_endian_integer( aStream, block_size ) ;
			write_little_endian_integer( aStream, context.number_of_samples ) ;
#if DEBUG_BGEN_FORMAT
			std::cerr << "genfile::bgen::write_sample_identifier_block(): writing " << sample_ids.size() << " samples...\n" ;
#endif
			for( uint32_t i = 0; i < sample_ids.size(); ++i ) {
#if DEBUG_BGEN_FORMAT
				std::cerr << "genfile::bgen::write_sample_identifier_block(): writing sample " << sample_ids[i] << ".\n" ;
#endif
				std::string const& identifier = sample_ids[i] ;
				assert( identifier.size() <= std::size_t( std::numeric_limits< uint16_t >::max() ) ) ;
				uint16_t const id_size = uint16_t( identifier.size() ) ;
				write_length_followed_by_data( aStream, id_size, identifier ) ;
			}
			return block_size ;
		}
		
		
		namespace impl {
			void check_for_two_alleles( uint16_t numberOfAlleles ) {
				if( numberOfAlleles != 2 ) {
					//std::cerr << "genfile::bgen::impl::check_for_two_alleles: only biallelic variants are currently supported.\n" ;
					//assert(0) ;
				}
			}
		
			struct TwoAlleleSetter {
				TwoAlleleSetter( std::string* allele1, std::string* allele2 ):
					m_allele1( allele1 ),
					m_allele2( allele2 )
				{
					assert( allele1 != 0 ) ;
					assert( allele2 != 0 ) ;
				}

				void operator()( uint16_t i, std::string const& value ) {
					if( i == 0 ) {
						*m_allele1 = value ;
					} else if( i == 1 ) {
						*m_allele2 = value ;
					} else {
						assert(0) ;
					}
				}
				std::string* m_allele1 ;
				std::string* m_allele2 ;
			} ;
		}
		
		namespace v10 {
			bool read_snp_identifying_data(
				std::istream& aStream,
				Context const& context,
				std::string* SNPID,
				std::string* RSID,
				std::string* chromosome,
				uint32_t* SNP_position,
				std::string* first_allele,
				std::string* second_allele
			) {
				// v1.0-style layout, deprecated
				uint32_t number_of_samples = 0 ;
				unsigned char max_id_size = 0 ;
				unsigned char SNPID_size = 0 ;
				unsigned char RSID_size = 0 ;
				if( aStream ) {
					read_little_endian_integer( aStream, &number_of_samples ) ;
					if( !aStream ) {
						return false ;
					}
					if( number_of_samples != context.number_of_samples ) {
						ThrowBGenError() ;
					}
				}
				if( aStream ) {
					read_little_endian_integer( aStream, &max_id_size ) ;
				}
				if( aStream ) {
					read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
					assert( SNPID_size <= max_id_size ) ;
					aStream.ignore( max_id_size - SNPID_size ) ;
				}
				if( aStream ) {
					read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
					assert( RSID_size <= max_id_size ) ;
					aStream.ignore( max_id_size - RSID_size ) ;
				}
				if( aStream ) {
					unsigned char chromosome_char = 0 ;
					read_little_endian_integer( aStream, &chromosome_char ) ;
					read_little_endian_integer( aStream, SNP_position ) ;

					switch( chromosome_char ) {
						case 1: *chromosome = "01" ; break ;
						case 2: *chromosome = "02" ; break ;
						case 3: *chromosome = "03" ; break ;
						case 4: *chromosome = "04" ; break ;
						case 5: *chromosome = "05" ; break ;
						case 6: *chromosome = "06" ; break ;
						case 7: *chromosome = "07" ; break ;
						case 8: *chromosome = "08" ; break ;
						case 9: *chromosome = "09" ; break ;
						case 10: *chromosome = "10" ; break ;
						case 11: *chromosome = "11" ; break ;
						case 12: *chromosome = "12" ; break ;
						case 13: *chromosome = "13" ; break ;
						case 14: *chromosome = "14" ; break ;
						case 15: *chromosome = "15" ; break ;
						case 16: *chromosome = "16" ; break ;
						case 17: *chromosome = "17" ; break ;
						case 18: *chromosome = "18" ; break ;
						case 19: *chromosome = "19" ; break ;
						case 20: *chromosome = "20" ; break ;
						case 21: *chromosome = "21" ; break ;
						case 22: *chromosome = "22" ; break ;
						case 23: *chromosome = "0X" ; break ;
						case 24: *chromosome = "0Y" ; break ;
						case 253: *chromosome = "XY" ; break ;
						case 254: *chromosome = "MT" ; break ;
						default: *chromosome = "NA" ;
					}
					unsigned char allele1_size = 0, allele2_size = 0 ;
					read_length_followed_by_data( aStream, &allele1_size, first_allele ) ;
					read_length_followed_by_data( aStream, &allele2_size, second_allele ) ;
				}
				return true ;
			}
		}
		
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) {
#if DEBUG_BGEN_FORMAT
			std::cerr << "genfile::bgen::impl::read_snp_identifying_data(): flags = 0x" << std::hex << context.flags << ".\n" ;
#endif
			uint32_t const layout = context.flags & e_Layout ;
			if( layout == e_Layout1 || layout == e_Layout2 ) {
				// forward to v12 version which handles multiple alleles.
				impl::TwoAlleleSetter allele_setter( first_allele, second_allele ) ;
				return bgen::read_snp_identifying_data(
					aStream, context,
					SNPID, RSID, chromosome, SNP_position,
					&impl::check_for_two_alleles,
					allele_setter
				) ;
			} else if( layout == e_Layout0 ) {
				return v10::read_snp_identifying_data( aStream, context, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;
			} else {
				assert(0) ;
			}
			return true ;
		}
		
		namespace v11 {
			namespace impl {
				double get_probability_conversion_factor( uint32_t flags ) {
					uint32_t layout = flags & e_Layout ;
					if( layout == e_Layout0 ) {
						// v1.0-style blocks, deprecated
						return 10000.0 ;
					} else if( layout == e_Layout1 ) {
						// v1.1-style blocks
						return 32768.0 ;
					} else {
						// v1.2 style (or other) blocks, these are treated differently and this function does not apply.
						assert(0) ;
					}
					return -1 ;
				}
			}
		}

		void ignore_genotype_data_block(
			std::istream& aStream,
			Context const& context
		) {
			if( (context.flags & bgen::e_CompressedSNPBlocks) != e_NoCompression ) {
				uint32_t compressed_data_size = 0 ;
				read_little_endian_integer( aStream, &compressed_data_size ) ;
				if( compressed_data_size > 0 ) {
					// gcc std::istream::ignore() has a bug / feature in which
					// it peeks at the next char and sets eof() if you ignore all the bytes in the file.
					// This breaks our expected invariant and means subsequent calls to peekg() fail with -1,
					// which stops us getting file size.
					// We deal with this by simply ignoring one less character here.
					aStream.ignore( compressed_data_size - 1 ) ;
					aStream.get() ;
				}
			}
			else {
				aStream.ignore( 6 * context.number_of_samples ) ;
			}
		}

		void read_genotype_data_block(
			std::istream& aStream,
			Context const& context,
			std::vector< byte_t >* buffer
		) {
			uint32_t payload_size = 0 ;
			if( (context.flags & e_Layout) == e_Layout2 || ((context.flags & e_CompressedSNPBlocks) != e_NoCompression ) ) {
				read_little_endian_integer( aStream, &payload_size ) ;
			} else {
				payload_size = 6 * context.number_of_samples ;
			}
			buffer->resize( payload_size ) ;
			aStream.read( reinterpret_cast< char* >( &(*buffer)[0] ), payload_size ) ;
		}

		void uncompress_probability_data(
			Context const& context,
			std::vector< byte_t > const& compressed_data,
			std::vector< byte_t >* buffer
		) {
			// compressed_data contains the (compressed or uncompressed) probability data.
			uint32_t const compressionType = (context.flags & bgen::e_CompressedSNPBlocks) ;
			if( compressionType != e_NoCompression ) {
				byte_t const* begin = &compressed_data[0] ;
				byte_t const* const end = &compressed_data[0] + compressed_data.size() ;
				uint32_t uncompressed_data_size = 0 ;
				if( (context.flags & e_Layout) == e_Layout1 ) {
					uncompressed_data_size = 6 * context.number_of_samples ;
				} else {
					begin = read_little_endian_integer( begin, end, &uncompressed_data_size ) ;
				}
				buffer->resize( uncompressed_data_size ) ;
				assert( compressionType == e_ZstdCompression );
				zstd_uncompress( begin, end, buffer ) ;
				assert( buffer->size() == uncompressed_data_size ) ;
			}
			else {
				// copy the data between buffers.
				buffer->assign( compressed_data.begin(), compressed_data.end() ) ;
			}
		}

		namespace v12 {
			namespace impl {
				namespace {
					// std::floor seems to sometimes eat cycles
					// so roll our own.
					// All values should be in the range 0...2^31-1 so should fit in a uint32.
					double floor( double v ) {
						return double( long( v )) ;
					}

					double fractional_part( double v ) {
						return( v - floor(v)) ;
					}

					double round( double v ) {
						return floor( v + 0.5 ) ;
					}
					
					struct CompareFractionalPart{
						CompareFractionalPart( double* v, std::size_t n ):
							m_v( v ), m_n( n )
						{}
						CompareFractionalPart( CompareFractionalPart const& other ):
							m_v( other.m_v ), m_n( other.m_n )
						{}
						CompareFractionalPart& operator=( CompareFractionalPart const& other ) {
							m_v = other.m_v ;
							m_n =  other.m_n ;
							return *this ;
						}
						
						bool operator()( std::size_t a, std::size_t b ) const {
							return( fractional_part( m_v[a] ) > fractional_part( m_v[b] )) ;
						}
					private:
						double* m_v ;
						std::size_t m_n ;
					} ;
				}

				void compute_approximate_probabilities( double* p, std::size_t* index, std::size_t const n, int const number_of_bits ) {
					double const scale = ( 0xFFFFFFFFFFFFFFFF >> ( 64 - number_of_bits ) ) ;
					double total_fractional_part = 0.0 ;
					double sum = 0.0 ;
					for( std::size_t i = 0; i < n; ++i ) {
						p[i] *= scale ;
						sum += p[i] ;
						index[i] = i ;
						total_fractional_part += fractional_part( p[i] ) ;
					}
					// Suppose the n input numbers sum to 1
					// Each has rounding error of at most machine epsilon
					// the sum thus has maximum error of n * epsilon * scale.
					assert( sum < (scale * (1 + n * std::numeric_limits< double >::epsilon()) ) ) ;

					// We have n numbers which sum to scale ± delta,
					// where delta is the above rounding error < n*epsilon.
					// sum_i floor(p_i) is an integer by definition.
					// Total fractional part is therefore of the form r ± delta where r is an integer.
					// Since scale = sum_i floor(p_i) + r, rounding up r of the p_i's yields a
					// set of integers summing to scale.
					std::size_t const r = round( total_fractional_part ) ;
					std::sort( index, index + n, CompareFractionalPart( p, n ) ) ;

					for( std::size_t i = 0; i < r; ++i ) {
						p[ index[i] ] = std::ceil( p[ index[i] ] ) ;
					}
					for( std::size_t i = r; i < n; ++i ) {
						p[ index[i] ] = floor( p[ index[i] ] ) ;
					}
				}
				
				byte_t* write_scaled_probs(
					uint64_t* data,
					std::size_t* offset,
					double const* probs,
					std::size_t const n,
					int const number_of_bits,
					byte_t* destination,
					byte_t* const end
				) {
					for( std::size_t i = 0; i < (n-1); ++i ) {
						uint64_t const storedValue = uint64_t( probs[i] ) ;
						*data |= storedValue << (*offset) ;
						(*offset) += number_of_bits ;
						if( (*offset) >= 32 ) {
							assert( (destination+4) <= end ) ;
							destination = std::copy(
								reinterpret_cast< byte_t const* >( data ),
								reinterpret_cast< byte_t const* >( data ) + 4,
								destination
							) ;
							(*offset) -= 32 ;
							(*data) >>= 32 ;
						}
					}
					return destination ;
				}
			}
		}
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include "db/sqlite3.hpp"
#include "genfile/bgen/IndexQuery.hpp"

// #define DEBUG 1

namespace genfile {
	namespace bgen {
		IndexQuery::UniquePtr IndexQuery::create( std::string const& filename, std::string const& table_name ) {
			return IndexQuery::UniquePtr( new SqliteIndexQuery( filename, table_name )) ;
		}

		SqliteIndexQuery::SqliteIndexQuery( std::string const& filename, std::string const& table_name ):
			m_connection( open_connection( filename ) ),
			m_metadata( load_metadata( *m_connection ) ),
			m_index_table_name( table_name ),
			m_initialised( false )
		{
		}
	
		boost::optional< SqliteIndexQuery::FileMetadata > const&
		SqliteIndexQuery::file_metadata() const {
			return m_metadata ;
		}

		void SqliteIndexQuery::initialise( ProgressCallback callback ) {
			db::Connection::StatementPtr stmt = build_query() ;
			if( callback ) {
				callback( 0, boost::optional< std::size_t >() ) ;
			}
			m_positions.reserve( 1000000 ) ;
			std::size_t batch_i = 0 ;
			for( stmt->step() ; !stmt->empty(); stmt->step(), ++batch_i ) {
				int64_t const pos = stmt->get< int64_t >( 0 ) ;
				int64_t const size = stmt->get< int64_t >( 1 ) ;
				assert( pos >= 0 ) ;
				assert( size >= 0 ) ;
				m_positions.push_back( std::make_pair( int64_t( pos ), int64_t( size ))) ;
				if( callback ) {
					callback( m_positions.size(), boost::optional< std::size_t >() ) ;
				}
			}
	#if DEBUG
			std::cerr << "SqliteIndexQuery::initialise(): read positions for " << m_positions.size() << " variants.\n" ;
	#endif
	
			m_initialised = true ;
		}

		std::size_t SqliteIndexQuery::number_of_variants() const {
			assert( m_initialised ) ;
			return m_positions.size() ;
		}

		SqliteIndexQuery::FileRange SqliteIndexQuery::locate_variant( std::size_t index ) const {
	#if DEBUG
			std::cerr << "SqliteIndexQuery::locate_variant(" << index << ")...\n" ;
	#endif
			assert( m_initialised ) ;
			assert( index < m_positions.size() ) ;
			return m_positions[index] ;
		}

		SqliteIndexQuery& SqliteIndexQuery::from_row( int row ) {
			m_query_parts.limit =
				(boost::format( "LIMIT %d OFFSET %d" ) % 1e13 % row).str() ;
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::include_range( GenomicRange const& range ) {
			m_query_parts.inclusion += ((m_query_parts.inclusion.size() > 0) ? " OR " : "" ) + (
				boost::format( "( chromosome == '%s' AND position BETWEEN %d AND %d )" ) % range.chromosome() % range.start() % range.end()
			).str() ;
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::exclude_range( GenomicRange const& range ) {
			m_query_parts.exclusion += ( m_query_parts.exclusion.size() > 0 ? " AND" : "" ) + (
				boost::format( " NOT ( chromosome == '%s' AND position BETWEEN %d AND %d )" )
					% range.chromosome() % range.start() % range.end()
			).str() ;
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::include_rsids( std::vector< std::string > const& ids ) {
			m_connection->run_statement( "CREATE TEMP TABLE IF NOT EXISTS tmpIncludedId( identifier TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = m_connection->get_statement( "INSERT INTO tmpIncludedId( identifier ) VALUES( ? )" ) ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				insert_stmt->bind( 1, ids[i] ).step() ;
				insert_stmt->reset() ;
			}
			if( m_query_parts.join.find( "tmpIncludedId" ) == std::string::npos ) {
				m_query_parts.join += " LEFT OUTER JOIN tmpIncludedId TI ON TI.identifier == V.rsid" ;
				m_query_parts.inclusion += ( m_query_parts.inclusion.size() > 0 ? " OR" : "" ) + std::string( " TI.identifier IS NOT NULL" ) ;
			}
			m_initialised = false ;
			return *this ;
		}

		SqliteIndexQuery& SqliteIndexQuery::exclude_rsids( std::vector< std::string > const& ids ) {
			m_connection->run_statement( "CREATE TEMP TABLE IF NOT EXISTS tmpExcludedId( identifier TEXT NOT NULL PRIMARY KEY ) WITHOUT ROWID" ) ;
			db::Connection::StatementPtr insert_stmt = m_connection->get_statement( "INSERT INTO tmpExcludedId( identifier ) VALUES( ? )" ) ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				insert_stmt->bind( 1, ids[i] ).step() ;
				insert_stmt->reset() ;
			}
			if( m_query_parts.join.find( "tmpExcludedId" ) == std::string::npos ) {
				m_query_parts.join += " LEFT OUTER JOIN tmpExcludedId TE ON TE.identifier == V.rsid" ;
				m_query_parts.exclusion += ( m_query_parts.exclusion.size() > 0 ? " AND" : "" ) + std::string( " TE.identifier IS NULL" ) ;
			}
			m_initialised = false ;
			return *this ;
		}

		db::Connection::UniquePtr SqliteIndexQuery::open_connection( std::string const& filename ) const {
			db::Connection::UniquePtr result ;
			try {
				result = db::Connection::create( "file:" + filename + "?nolock=1", "r" ) ;
				//result = db::Connection::create( filename, "r" ) ;
			} catch( db::ConnectionError const& e ) {
				throw std::invalid_argument( "Could not open the index file \"" + filename + "\"" ) ;
			}
			return result ;
		}

		SqliteIndexQuery::OptionalFileMetadata
		SqliteIndexQuery::load_metadata( db::Connection& connection ) const {
			OptionalFileMetadata result ;
			db::Connection::StatementPtr stmt = connection.get_statement( "SELECT * FROM sqlite_master WHERE name == 'Metadata' AND type == 'table'" ) ;
			stmt->step() ;
			if( !stmt->empty() ) {
				db::Connection::StatementPtr mdStmt = connection.get_statement( "SELECT filename, file_size, first_1000_bytes FROM Metadata" ) ;
				mdStmt->step() ;

				if( mdStmt->empty() ) {
					throw std::invalid_argument( "Index file appears malformed (empty \"Metadata\" table)" ) ;
				}
				FileMetadata metadata ;
				// Get metadata fields for comparison
				metadata.filename = mdStmt->get< std::string >( 0 ) ;
				metadata.size = mdStmt->get< int64_t >( 1 ) ;
				metadata.first_bytes = mdStmt->get< std::vector< uint8_t > >( 3 ) ;
				
				result = metadata ;
			}
			return result ;
		}

		db::Connection::StatementPtr SqliteIndexQuery::build_query() const {
			std::string const select = "SELECT file_start_position, size_in_bytes FROM `"
				+ m_index_table_name + "` V" ;
			std::string const inclusion = ( m_query_parts.inclusion.size() > 0 ) ? ("(" + m_query_parts.inclusion + ")") : "" ;
			std::string const exclusion = ((m_query_parts.inclusion.size() > 0 && m_query_parts.exclusion.size() > 0 ) ? "AND " : "" )
				+ (( m_query_parts.exclusion.size() > 0 ) ? ("(" + m_query_parts.exclusion + ")") : "" ) ;
			std::string const where = (inclusion.size() > 0 || exclusion.size() > 0) ? ("WHERE " + inclusion + exclusion) : "" ;
			std::string const orderBy = "ORDER BY chromosome, position, rsid, allele1, allele2" ;
			std::string const select_sql = select + " " + m_query_parts.join + " " + where + " " + orderBy + " " + m_query_parts.limit ;
	#if DEBUG
			std::cerr << "BgenIndex::build_query(): SQL is: \"" << select_sql << "\"...\n" ;
	#endif
	
			return m_connection->get_statement( select_sql ) ;
		}
	}
}


//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include "genfile/MissingValue.hpp"

namespace genfile {
	std::ostream& operator<<( std::ostream& o, MissingValue const& ) {
		return o << "NA" ;
	}
	
	bool MissingValue::operator==( MissingValue const& ) const {
		return true ;
	}

	bool MissingValue::operator<( MissingValue const& ) const {
		return false ;
	}

	bool MissingValue::operator<=( MissingValue const& ) const {
		return true ;
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/IndexQuery.hpp"
#include "genfile/bgen/View.hpp"

namespace genfile {
	namespace bgen {
		View::UniquePtr View::create( std::string const& filename ) {
			return View::UniquePtr( new View( filename )) ;
		}

		/* View implementation */
		View::View( std::string const& filename ):
			m_filename( filename ),
			m_variant_i(0),
			m_have_sample_ids( false ),
			m_state( e_NotOpen )
		{
			setup( m_filename ) ;
			m_file_position = m_stream->tellg() ;
		}

		std::size_t View::number_of_samples() const {
			return m_context.number_of_samples ;
		}

		void View::set_query( IndexQuery::UniquePtr &query ) {
			m_index_query = std::move(query) ;
			if( m_index_query->number_of_variants() > 0 ) {
				m_stream->seekg( m_index_query->locate_variant(0).first ) ;
			}
			m_file_position = m_stream->tellg() ;
		}

		View::FileMetadata const& View::file_metadata() const {
			return m_file_metadata ;
		}

		genfile::bgen::Context const& View::context() const {
			return m_context ;
		}

		std::streampos View::current_file_position() const {
			return m_file_position ;
		}

		uint32_t View::number_of_variants() const {
			if( m_index_query.get() ) {
				return m_index_query->number_of_variants() ;
			} else {
				return m_context.number_of_variants ;
			}
		}

		std::ostream& View::summarise( std::ostream& o ) const {
			o << "View: bgen file ("
				<< ( m_context.flags & genfile::bgen::e_Layout2 ? "layout = 2" : "layout = 1" )
				<< ", " ;
			std::string compression = "no" ;
			if( (m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZlibCompression ) {
				compression = "zlib" ;
			} else if( (m_context.flags & genfile::bgen::e_CompressedSNPBlocks) == genfile::bgen::e_ZstdCompression ) {
				compression = "zstd" ;
			}
			o << compression << " compression)" ;
			o << " with " 
				<< m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
				<< m_context.number_of_variants << " variants.\n" ;
			if( m_index_query.get() ) {
				o << "IndexQuery: query will return " << m_index_query->number_of_variants() << " variants.\n" ;
			}
			return o ;
		}

		namespace {
			// callbacks to avoid using C++11 lambda functions
			void resize_vector( std::vector< std::string >* target, std::size_t n ) {
				target->resize(n) ;
			}
			
			void set_vector_element( std::vector< std::string >* target, std::size_t i, std::string const& value ) {
				target->at(i) = value ;
			}

			void push_back_vector( std::vector< std::string >* target, std::string const& value ) {
				target->push_back( value ) ;
			}
		}

		bool View::read_variant(
			std::string* SNPID,
			std::string* rsid,
			std::string* chromosome,
			uint32_t* position,
			std::vector< std::string >* alleles
		) {
			assert( m_state == e_ReadyForVariant ) ;

			if( m_index_query.get() ) {
				if( m_variant_i == m_index_query->number_of_variants() ) {
					return false ;
				}
				IndexQuery::FileRange const range = m_index_query->locate_variant( m_variant_i ) ;
				m_stream->seekg( range.first ) ;
			}

			if(
				genfile::bgen::read_snp_identifying_data(
					*m_stream, m_context,
					SNPID, rsid, chromosome, position,
//					[&alleles]( std::size_t n ) { alleles->resize( n ) ; },
					boost::bind( &resize_vector, alleles, _1 ),
//					[&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
					boost::bind( &set_vector_element, alleles, _1, _2 )
				)
			) {
				m_state = e_ReadyForProbs ;
				return true ;
			} else {
				return false ;
			}
		}

		// Read and uncompress genotype probability data, and unpack
		// it into constituent parts, but don't do a full parse.
		// This can lead to more efficient code paths than a full parse for some operations.
		//
		// Currently this function works for 'layout=2' files, e.g. v1.2 and above only.
		// Data is returned in the fields of the supplied 'pack' object, which
		// is defined in bgen.hpp.
		void View::read_and_unpack_v12_genotype_data_block(
			genfile::bgen::v12::GenotypeDataBlock* pack
		) {
			assert( (m_context.flags & genfile::bgen::e_Layout) == genfile::bgen::e_Layout2 ) ;
			std::vector< byte_t > const& buffer = read_and_uncompress_genotype_data_block() ;
			pack->initialise( m_context, &buffer[0], &buffer[0] + buffer.size() ) ;
			++m_variant_i ;
		}

		// Ignore genotype probability data for the SNP just read using read_variant()
		// After calling this method it should be safe to call read_variant()
		// to fetch the next variant from the file.
		void View::ignore_genotype_data_block() {
			assert( m_state == e_ReadyForProbs ) ;
			genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
			m_file_position = m_stream->tellg() ;
			m_state = e_ReadyForVariant ;
			++m_variant_i ;
		}

		// Open the bgen file, read header data and gather metadata.
		void View::setup( std::string const& filename ) {
			m_file_metadata.filename = filename ;

			// Open the stream
			m_stream.reset(
				new std::ifstream( filename.c_str(), std::ifstream::binary )
			) ;
			if( !*m_stream ) {
				throw std::invalid_argument( filename ) ;
			}

			// get file size
			{
				auto origin = m_stream->tellg() ;
				m_stream->seekg( 0, std::ios::end ) ;
				m_file_metadata.size = m_stream->tellg() - origin ;
				m_stream->seekg( 0, std::ios::beg ) ;
			}
			// read first (up to) 1000 bytes.
			{
				m_file_metadata.first_bytes.resize( 1000, 0 ) ;
				m_stream->read( reinterpret_cast< char* >( &m_file_metadata.first_bytes[0] ), 1000 ) ;
				m_file_metadata.first_bytes.resize( m_stream->gcount() ) ;
				m_stream->clear() ;
			}

			m_state = e_Open ;

			// Read the offset, header, and sample IDs if present.
			m_stream->seekg( 0, std::ios::beg ) ;
			genfile::bgen::read_offset( *m_stream, &m_offset ) ;
			genfile::bgen::read_header_block( *m_stream, &m_context ) ;

			if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
				genfile::bgen::read_sample_identifier_block(
					*m_stream, m_context,
					boost::bind( &push_back_vector, &m_sample_ids, _1 )
					//[this]( std::string id ) { m_sample_ids.push_back( id ) ; }
				) ;
				m_have_sample_ids = true ;
			}
	
			// read data up to first data block.
			m_postheader_data.resize( m_offset+4 - m_stream->tellg() ) ;
			m_stream->read( reinterpret_cast< char* >( &m_postheader_data[0] ), m_postheader_data.size() ) ;
			if( m_stream->gcount() != m_postheader_data.size() ) {
				throw std::invalid_argument(
					(
						boost::format(
							"BGEN file (\"%s\") appears malformed - offset specifies more bytes (%d) than are in the file."
						) % filename % m_offset
					).str()
				) ;
			}

			// Jump to the first variant data block.
			// m_stream->seekg( m_offset + 4 ) ;

			// We keep track of state (though it's not really needed for this implementation.)
			m_state = e_ReadyForVariant ;
		}

		// Utility function to read and uncompress variant genotype probability data
		// without further processing.
		std::vector< byte_t > const& View::read_and_uncompress_genotype_data_block() {
			assert( m_state == e_ReadyForProbs ) ;
			genfile::bgen::read_genotype_data_block( *m_stream, m_context, &m_buffer1 ) ;
			m_file_position = m_stream->tellg() ;
			m_state = e_ReadyForVariant ;
			genfile::bgen::uncompress_probability_data( m_context, m_buffer1, &m_buffer2 ) ;
			return m_buffer2 ;
		}
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <zlib.h>
#include "genfile/zlib.hpp"

namespace genfile {
	void zstd_compress(
		uint8_t const* buffer,
		uint8_t const* const end,
		std::vector< uint8_t >* dest,
		std::size_t const offset,
		int const compressionLevel
	) {
		assert( dest != 0 ) ;
		assert( compressionLevel >= 1 && compressionLevel <= 22 ) ;
		std::size_t const source_size = ( end - buffer ) ;
		std::size_t compressed_size = ZSTD_compressBound( source_size ) ;
		dest->resize( compressed_size + offset ) ;
		compressed_size = ZSTD_compress(
			reinterpret_cast< void* >( &( dest->operator[](0) ) + offset ),
			compressed_size,
			reinterpret_cast< void const* >( buffer ),
			source_size,
			compressionLevel
		) ;
		assert( !ZSTD_isError( compressed_size )) ;
		dest->resize( compressed_size + offset ) ;
	}
}
