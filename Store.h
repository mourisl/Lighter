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

#if MAX_KMER_LENGTH <= 32 
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
#else
	int Put( uint64_t code[], int kmerLength, bool testFirst )
	{
		//return 0 ;
		//printf( "%d\n", method ) ;
		//printf( "%llu\n", val ) ;
		GetCanonicalKmerCode( code, kmerLength ) ;
		//printf( "%llu\n", val ) ;
		//exit(1) ;
		if ( method == 1 )
		{
			//hash[ val ] = 1 ;
			return 1 ;
		}
		
		//printf( "1: %d\n", bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )  ;
		if ( numOfThreads > 1 && testFirst && bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )
			return 0 ;
		//printf( "2: %lld %d\n", code[0], sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) )  ;
		bf.insert( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) ) ;
		//printf( "3: %d\n", bf.contains( (char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32  + 1 ) ) )  ;
		return 0 ;
	}

	int IsIn( uint64_t code[], int kmerLength ) 
	{
		//printf( "1. %llu\n", val ) ;
		GetCanonicalKmerCode( code, kmerLength ) ;
		//printf( "2. %llu\n", val ) ;
		if ( method == 1 )
		{
			//return ( hash.find( val ) != hash.end() ) ;
		}
		return bf.contains( ( char *)code, sizeof( uint64_t ) * ( ( kmerLength - 1 ) / 32 + 1 ) ) ;
	}
#endif 
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
#if MAX_KMER_LENGTH <= 32 
		Put( code.GetCode(), code.GetKmerLength(), testFirst ) ;
#else
		uint64_t c[K_BLOCK_NUM] ;
		code.GetCode( c ) ;
		Put( c, code.GetKmerLength(), testFirst ) ;
#endif
		return 0 ;
	}

	int IsIn( KmerCode &code ) 
	{
		if ( !code.IsValid() )
			return 0 ;
#if MAX_KMER_LENGTH <= 32 
		return IsIn( code.GetCode(), code.GetKmerLength() ) ;
#else
		uint64_t c[K_BLOCK_NUM] ;
		code.GetCode( c ) ;
		return IsIn( c, code.GetKmerLength() ) ;
#endif
	}


	void SetNumOfThreads( int in ) 
	{ 
		numOfThreads = in ;
		bf.SetNumOfThreads( in ) ;
	}

#if MAX_KMER_LENGTH <= 32 
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
#else
	void GetCanonicalKmerCode( uint64_t code[], int k ) 
	{
		int i, j ;
		uint64_t crCode[ K_BLOCK_NUM ] ; // complementary code
		const int largestBlock = ( k - 1 ) / ( sizeof( uint64_t ) * 4 ) ;
		for ( i = 0 ; i < K_BLOCK_NUM ; ++i )
			crCode[i] = 0 ;

		for ( i = 0, j = k - 1 ; i < k ; ++i, --j )
		{
			//uint64_t tmp = ( code >> ( 2ull * (k - i - 1) ) ) & 3ull   ; 
			//crCode = ( crCode << 2ull ) | ( 3 - tmp ) ;
			
			/*int tagI = i >> 5 ;
			int tagJ = j >> 5 ;
			int residualI = i & 31 ;
			int residualJ = j & 31 ;*/
			uint64_t tmp = ( code[ i >> 5 ] >> ( 2ull * ( i & 31 ) ) ) & 3ull ;
			crCode[j >> 5] = crCode[ j >> 5 ] | ( ( 3ull - tmp ) << ( 2ull * ( j & 31 ) ) ) ;
		}

		bool crFlag = false ;
		for ( i = largestBlock ; i >= 0 ; --i )
		{
			if ( crCode[i] < code[i] )
			{
				crFlag = true ;
				break ;
			}
			else if ( crCode[i] > code[i] )
			{
				crFlag = false ;
				break ;
			}
		}
		//printf( "%llu,%llu %llu,%llu: %d\n", code[1], code[0], crCode[1], crCode[0], crFlag ) ;
		if ( crFlag )
		{
			for ( i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = crCode[i] ;
		}
	}
#endif

	void BloomOutput( char *file)
	{
		FILE *fp = fopen( file, "w" ) ;
		if ( fp == NULL )
		{
			printf( "Could not find file %s.\n", file ) ;
			exit( EXIT_FAILURE ) ;
		}
		bf.Output( fp ) ;
		fclose( fp ) ;
	}

	void BloomInput( char *file ) 
	{
		FILE *fp = fopen( file, "r" ) ;
		if ( fp == NULL )
		{
			printf( "Could not find file %s.\n", file ) ;
			exit( EXIT_FAILURE ) ;
		}
		bf.Input( fp ) ;
		fclose( fp ) ;
	}

	int Clear() 
	{
		return 0 ;
	}
} ;
#endif
