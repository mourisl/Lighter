#ifndef _MOURISL_CLASS_STORE
#define _MOURISL_CLASS_STORE
/**
  The wrapper for bloom filters to store kmers.
*/
#include <stdio.h>
#include <stdint.h>
#include <map>

#include "bloom_filter.hpp"
#include "KmerCode.h"

class Store
{
private:
	uint64_t size ;
	bloom_parameters bfpara ;
	bloom_filter bf ;
	std::map<uint64_t, int> hash ;
	int method ; //0-bloom filter. 1-std:map

	int Put( uint64_t val, int kmerLength, bool testFirst )
	{
		//return 0 ;
		//printf( "%d\n", method ) ;
		//printf( "%llu\n", val ) ;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		//printf( "%llu\n", val ) ;
		//exit(1) ;
		if ( method == 1 )
		{
			hash[ val ] = 1 ;
			return 1 ;
		}
		
		if ( numOfThreads > 1 && testFirst && bf.contains( val ) )
			return 0 ;
		bf.insert( val ) ;
		return 0 ;
	}

	int IsIn( uint64_t val, int kmerLength ) 
	{
		//printf( "1. %llu\n", val ) ;
		val = GetCanonicalKmerCode( val, kmerLength ) ;
		//printf( "2. %llu\n", val ) ;
		if ( method == 1 )
		{
			return ( hash.find( val ) != hash.end() ) ;
		}

		return bf.contains( val ) ;
	}

	int numOfThreads ;
public:
	Store( double fprate = 0.01 ): size( 10000003 ), bfpara( 10000003, fprate ), bf( bfpara )
	{
		numOfThreads = 1 ;
		method = 0 ;
	}

	Store( uint64_t s, double fprate = 0.01 ): size( s ), bfpara( s, fprate ), bf( bfpara )  
	{
		numOfThreads = 1 ;
		method = 0 ;
	}

	~Store() 
	{
	}
	
	double Occupancy()
	{
		return bf.occupancy() ;
	}
	double GetFP()
	{
		if ( method == 1 )
			return 0 ;
		return bf.GetActualFP() ;
	}
	int Put( KmerCode &code, bool testFirst = false ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		Put( code.GetCode(), code.GetKmerLength(), testFirst ) ;
		return 0 ;
	}

	int IsIn( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
		return IsIn( code.GetCode(), code.GetKmerLength() ) ;
	}


	void SetNumOfThreads( int in ) 
	{ 
		numOfThreads = in ;
		bf.SetNumOfThreads( in ) ;
	}

	uint64_t GetCanonicalKmerCode( uint64_t code, int k ) 
	{
		int i ;
		uint64_t crCode = 0ull ; // complementary code
		for ( i = 0 ; i < k ; ++i )
		{
			//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
			//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;

			uint64_t tmp = ( code >> ( 2ull * i ) ) & 3ull ;
			crCode = ( crCode << 2ull ) | ( 3ull - tmp ) ;
		}
		return crCode < code ? crCode : code ;
	}

	void TemporaryOutput( char *file)
	{
		FILE *fp = fopen( file, "w" ) ;
		bf.Output( fp ) ;
		fclose( fp ) ;
	}

	void TemporaryInput( char *file ) 
	{
		FILE *fp = fopen( file, "r" ) ;
		bf.Input( fp ) ;
		fclose( fp ) ;
	}

	int Clear() 
	{
		return 0 ;
	}
} ;
#endif
