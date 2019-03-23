
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "db/Connection.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"

namespace db {
	Connection::UniquePtr Connection::create( std::string const& filename, std::string const& mode ) {
		return Connection::UniquePtr( new SQLite3Connection( filename, true, mode )) ;
	}

	void Connection::run_statement( std::string const& SQL ) {
		get_statement( SQL )->step() ;
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "db/SQLite3Statement.hpp"
#include "db/Error.hpp"

namespace db {
	Error::Error( std::string const& caller, std::string const& db_spec, int error, std::string const& sql ):
		m_spec( db_spec ),
		m_error( error ),
		m_sql( sql )
	{
		assert( m_error != SQLite3Statement::Error::OK ) ;
	}
	
	Error::~Error() throw() {}
	
	std::string Error::description() const {
		typedef db::SQLite3Statement::Error Error ;
		std::string result ;
		switch( m_error ) {
			case Error::ERROR: 		result = "SQL error or missing database"; break ;
			case Error::INTERNAL: 	result = "Internal logic error in SQLite"; break ;
			case Error::PERM: 		result = "Access permission denied"; break ;
			case Error::ABORT: 		result = "Callback routine requested an abort"; break ;
			case Error::BUSY: 		result = "The database file is locked"; break ;
			case Error::LOCKED: 	result = "A table in the database is locked"; break ;
			case Error::NOMEM: 		result = "A malloc() failed"; break ;
			case Error::READONLY: 	result = "Attempt to write a readonly database"; break ;
			case Error::INTERRUPT: 	result = "Operation terminated by sqlite3_interrupt()"; break ;
			case Error::IOERR: 		result = "Some kind of disk I/O error occurred"; break ;
			case Error::CORRUPT: 	result = "The database disk image is malformed"; break ;
			case Error::NOTFOUND: 	result = "NOT USED. Table or record not found"; break ;
			case Error::FULL: 		result = "Insertion failed because database is full"; break ;
			case Error::CANTOPEN: 	result = "Unable to open the database file"; break ;
			case Error::PROTOCOL: 	result = "Database lock protocol error"; break ;
			case Error::EMPTY: 		result = "Database is empty"; break ;
			case Error::SCHEMA: 	result = "The database schema changed"; break ;
			case Error::TOOBIG: 	result = "String or BLOB exceeds size limit"; break ;
			case Error::CONSTRAINT: result = "Abort due to constraint violation"; break ;
			case Error::MISMATCH: 	result = "Data type mismatch"; break ;
			case Error::MISUSE: 	result = "Library used incorrectly"; break ;
			case Error::NOLFS: 		result = "Uses OS features not supported on host"; break ;
			case Error::AUTH: 		result = "Authorization denied"; break ;
			case Error::FORMAT: 	result = "Auxiliary database format error"; break ;
			case Error::RANGE: 		result = "2nd parameter to sqlite3_bind out of range"; break ;
			case Error::NOTADB: 	result = "File opened that is not a database file"; break ;
			case Error::ROW: 		result = "sqlite3_step() has another row ready"; break ;
			case Error::DONE: 		result = "sqlite3_step() has finished executing"; break ;
			default: 				assert(0) ;
		}
		result += ", in statement \"" + m_sql + "\"" ;
		return result ;
	}
}


//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <stdint.h>
#include <thread>
#include <chrono>
#include "db/SQLite3Connection.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/Error.hpp"

extern "C" {
	void sqlite3_trace_callback( void* udp, const char* sql ) {
#if DEBUG
		std::cerr << "SQLite3 trace: SQL = \"" << sql << "\".\n" ;
#endif
	}
	
	int sqlite3_busy_callback( void*, int number_of_tries ) {
		if( number_of_tries > 10 ) {
			return 0 ;
		}
		std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) ) ;
		return 1 ;
	}
}

namespace db {
	SQLite3Connection::SQLite3Connection( std::string const& filename, bool overwrite, std::string const& mode ):
		m_filename( filename ),
		m_db_connection(0),
		m_managed( true )
	{
		open_db_connection( filename, overwrite, mode ) ;
	}

	SQLite3Connection::~SQLite3Connection() {
		close_db_connection_if_necessary() ;
	}

	
	SQLite3Connection::StatementPtr SQLite3Connection::get_statement( std::string const& SQL ) {
		return StatementPtr( new SQLite3Statement( this, SQL ) ) ;
	}	

	Connection::RowId SQLite3Connection::get_last_insert_row_id() const {
		uint64_t result = sqlite3_last_insert_rowid( m_db_connection ) ;
		return result ;
	}

	sqlite3_stmt* SQLite3Connection::prepare_sql( std::string const& SQL ) const {
		assert( m_db_connection != 0 ) ;
		sqlite3_stmt* statement ;
		int code = sqlite3_prepare_v2(
			m_db_connection,
			SQL.c_str(),
			static_cast< int >( SQL.size()+1 ), // +1 accounts for null terminating byte.
			&statement,
			0 // ignore pzTail
		) ;
		if( code != SQLITE_OK ) {
			throw StatementPreparationError( "SQLite3Connection::prepare_sql()", get_spec(), code, SQL ) ;
		}
		// SQLite might return 0 if the SQL consisted of comments and whitespace only
		// but I want to treat this as a programmer error.
		assert( statement != 0 ) ;
		return statement ;
	}

	int SQLite3Connection::finalise_statement( sqlite3_stmt* statement ) {
		assert( statement != 0 ) ;
		return sqlite3_finalize( statement ) ;
	}

	int SQLite3Connection::step_statement( sqlite3_stmt* statement ) {
		int code = sqlite3_step( statement ) ;
		if( code != SQLITE_ROW && code != SQLITE_DONE ) {
			throw StatementStepError( "SQLite3Connection::step_statement()", get_spec(), code, std::string( sqlite3_sql( statement ) ) ) ;
		}
		return (code == SQLITE_ROW);
	}
	
	void SQLite3Connection::open_db_connection( std::string const& filename, bool overwrite, std::string const& mode ) {
		int flags = 0 ;
		if( mode == "r" ) {
			flags |= SQLITE_OPEN_READONLY ;
		} else if( mode == "rw" ) {
			flags |= SQLITE_OPEN_READWRITE ;
			if( overwrite ) {
				flags |= SQLITE_OPEN_CREATE ;
			}
		} else {
			assert( 0 ) ;
		}
		int code = sqlite3_open_v2( filename.c_str(), &m_db_connection, flags, NULL ) ;
		if( code != SQLITE_OK ) {
			throw ConnectionError( "SQLite3Connection::open_db_connection()", get_spec(), code ) ;
		}
		// We add a busy handler.  This makes the database more robust by retrying failed transactions.
		sqlite3_busy_handler( m_db_connection, &sqlite3_busy_callback, NULL ) ;
		// Uncomment the next line to trace SQL statements executed.
		//sqlite3_trace( m_db_connection, &sqlite3_trace_callback, NULL ) ;
	}

	void SQLite3Connection::close_db_connection_if_necessary() {
		if( m_managed && m_db_connection != 0 ) {
			// According to SQLite docs, we must finalise any prepared statements
			// before we can close the db connection.
			finalise_prepared_statements() ;
			sqlite3_close( m_db_connection ) ;
			m_db_connection = 0 ;
		}
	}

	void SQLite3Connection::finalise_prepared_statements() {
		assert( m_db_connection != 0 ) ;
#if SQLITE_VERSION_NUMBER > 3006000
		sqlite3_stmt* statement ;
		while(( statement = sqlite3_next_stmt( m_db_connection, 0)) != 0 ) {
			finalise_statement( statement ) ;
		}
#endif
	}
	
	SQLite3Connection::ScopedTransactionPtr SQLite3Connection::open_transaction( double max_seconds_to_wait ) {		
		ScopedTransactionPtr transaction ;
		// we wait 10 milliseconds between attempts.
		
		for( std::size_t count = 0 ; true; ++count ) {
			try {
				transaction.reset( new SQLite3Connection::Transaction( *this ) ) ;
				break ;
			}
			catch( db::StatementStepError const& e ) {
				// Because of the busy handler (see top of this file)
				// each attempt takes ~0.1s anyway
				// We wait an additional 0.1s so that each attempt takes 0.2s in total.
				if( ( count * 0.2 ) > max_seconds_to_wait ) {
#if DEBUG
					std::cerr << "Open transaction: failure count=" << count << " (~" << count*0.2 << "s).  Bailing out.\n" ;
#endif
					std::this_thread::sleep_for( std::chrono::milliseconds( 100 ) ) ;
					throw TransactionError( "SQLite3Connection::open_transaction()", get_spec(), e.error_code(), e.sql() ) ;
				}
			}
		}
		return transaction ;
	}
	
	SQLite3Connection::Transaction::Transaction( SQLite3Connection& connection ):
		m_connection( connection )
	{
		m_connection.run_statement( "BEGIN IMMEDIATE TRANSACTION" ) ;
	}
	
	SQLite3Connection::Transaction::~Transaction() {
		m_connection.run_statement( "COMMIT" ) ;
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <cassert>
#include <string>
#include <exception>
#include <boost/lexical_cast.hpp>
#include "sqlite3/sqlite3.h"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/SQLite3Statement.hpp"
#include "db/SQLite3Error.hpp"

namespace db {
	SQLite3Statement::SQLite3Statement( SQLite3Connection* connection, std::string const& SQL ):
		m_connection( connection ),
		m_have_results( false )
	{
		assert( m_connection ) ;
		m_statement = m_connection->prepare_sql( SQL ) ;
		//std::cerr << "SQLite3Statement::SQLite3Statement(): statement is \"" + SQL + "\".\n" ;
		assert( m_statement != 0 ) ;
	}
	
	SQLite3Statement::~SQLite3Statement() {
		m_connection->finalise_statement( m_statement ) ;
	}

	bool SQLite3Statement::step() {
		m_have_results = m_connection->step_statement( m_statement ) ;
		return m_have_results ;
	}

	bool SQLite3Statement::empty() const {
		return !m_have_results ;
	}

	std::size_t SQLite3Statement::get_number_of_columns() const {
		return std::size_t( get_column_count() ) ;
	}

	std::string SQLite3Statement::get_name_of_column( std::size_t i ) const {
		assert( m_statement != 0 ) ;
		assert( i < get_number_of_columns() ) ;
		return std::string( sqlite3_column_origin_name( m_statement, i )) ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, int32_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint32_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, int64_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int64( m_statement, i, sqlite3_int64( value ) ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint64_t value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_int64( m_statement, i, sqlite3_int64( value ) ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, double value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_double( m_statement, i, value ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, std::string const& value ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_text( m_statement, i, value.c_str(), value.size(), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, char const* buffer, char const* const end ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_blob( m_statement, i, reinterpret_cast< void const* >( buffer ), int( end - buffer ), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind( std::size_t i, uint8_t const* buffer, uint8_t const* const end ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_blob( m_statement, i, reinterpret_cast< void const* >( buffer ), int( end - buffer ), SQLITE_TRANSIENT ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::bind_NULL( std::size_t i ) {
		assert( m_statement != 0 ) ;
		int error = sqlite3_bind_null( m_statement, i ) ;
		if( error != SQLITE_OK ) {
			throw ValueBindError( "SQLite3Statement::bind()", m_connection->get_spec(), error, boost::lexical_cast<std::string>( i ) ) ;
		}
		return *this ;
	}

	SQLite3Statement& SQLite3Statement::reset() {
		assert( m_statement != 0 ) ;
		int error = sqlite3_reset( m_statement ) ;
		if( error != SQLITE_OK ) {
			throw db::Error( "SQLite3Statement::reset()", m_connection->get_spec(), error ) ;
		}
		return *this ;
	}

	std::string SQLite3Statement::get_sql() const {
		assert( m_statement != 0 ) ;
		return std::string( sqlite3_sql( m_statement )) ;
	}

	bool SQLite3Statement::is_null( int column_id ) const {
		return ( sqlite3_column_type( m_statement, column_id ) == SQLITE_NULL ) ;
	}

	int SQLite3Statement::get_column_int( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_int( m_statement, column_id ) ;
	}

	int64_t SQLite3Statement::get_column_int64( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_int64( m_statement, column_id ) ;
	}

	double SQLite3Statement::get_column_double( int column_id ) const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_double( m_statement, column_id ) ;
	}
	
	std::string SQLite3Statement::get_column_string( int column_id ) const {
		assert( m_statement != 0 ) ;
		return reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id ) ) ;
	}

	char SQLite3Statement::get_column_char( int column_id ) const {
		assert( m_statement != 0 ) ;
		char const* p = reinterpret_cast< char const * >( sqlite3_column_text( m_statement, column_id )) ;
		int bytes = sqlite3_column_bytes( m_statement, column_id ) ;
		assert( bytes == 1 ) ;
		return *p ;
	}
	
	std::vector< uint8_t > SQLite3Statement::get_column_blob( int column_id ) const {
		assert( m_statement != 0 ) ;
		uint8_t const* p = reinterpret_cast< uint8_t const* >( sqlite3_column_blob( m_statement, column_id )) ;
		int nBytes = sqlite3_column_bytes( m_statement, column_id ) ;
		return std::vector< uint8_t >( p, p+nBytes ) ;
	}
	
	int SQLite3Statement::get_column_count() const {
		assert( m_statement != 0 ) ;
		return sqlite3_column_count( m_statement ) ;
	}
}

//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <string>
#include <stdint.h>
#include "sqlite3/sqlite3.h"
#include "db/SQLStatement.hpp"

namespace db {
	SQLStatement::~SQLStatement() {}
	
	template<>
	int SQLStatement::get_column< int >( int column_id ) const {
		return this->get_column_int( column_id ) ;
	}

	template<>
	int64_t SQLStatement::get_column< int64_t >( int column_id ) const {
		return this->get_column_int64( column_id ) ;
	}

	template<>
	double SQLStatement::get_column< double >( int column_id ) const {
		return this->get_column_double( column_id ) ;
	}

	template<>
	std::string SQLStatement::get_column< std::string >( int column_id ) const {
		return this->get_column_string( column_id ) ;
	}

	template<>
	char SQLStatement::get_column< char >( int column_id ) const {
		return this->get_column_char( column_id ) ;
	}

	template<>
	std::vector< uint8_t > SQLStatement::get_column< std::vector< uint8_t > >( int column_id ) const {
		return this->get_column_blob( column_id ) ;
	}
	
	SQLStatement& SQLStatement::bind( std::size_t i, char const* value ) {
		bind( i, std::string( value )) ;
		return *this ;
	}
}
