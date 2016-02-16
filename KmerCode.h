#ifndef _MOURISL_KMERCODE_HEADER
#define _MOURISL_KMERCODE_HEADER

#include <stdio.h>
#include "utils.h"

extern char nucToNum[26] ;
extern char numToNuc[26] ;

const int K_BLOCK_NUM = ( MAX_KMER_LENGTH - 1 ) / ( sizeof( uint64_t ) * 4 ) + 1 ;


#if MAX_KMER_LENGTH <= 32
class KmerCode
{
	private:
		int kmerLength ;
		int invalidPos ; // The position contains characters other than A,C,G,T in the code
		uint64_t code ;
		uint64_t mask ;
	public: 
		KmerCode() 
		{
		}
		KmerCode( int kl ) 
		{
			int i ;
			kmerLength = kl ;
			code = 0 ;
			invalidPos = -1 ;

			mask = 0 ;
			for ( i = 0 ; i < kmerLength ; ++i )
			{
				mask = mask << 2 ;
				mask = mask | 3 ;
			}
		}

		KmerCode( const KmerCode& in )
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			mask = in.mask ;
			code = in.code ;
		}

		void Restart() { code = 0ull ; invalidPos = -1 ; } 
		uint64_t GetCode() { return code ; } 
		int GetKmerLength() { return kmerLength ; }

		bool IsValid() 
		{
			/*if ( invalidPos != -1 )
			{
				printf( "hi\n") ;
			}*/
			return ( invalidPos == -1 ) ;
		}

		void Append( char c ) 
		{
			if ( invalidPos != -1 )
				++invalidPos ;

			code = ( ( code << 2ull ) & mask ) | 
				( (uint64_t)( nucToNum[ c - 'A' ] & 3 ) ) ;

			if ( nucToNum[c - 'A'] == -1 )
			{
				invalidPos = 0 ;	
			}
			if ( invalidPos >= kmerLength )
				invalidPos = -1 ;
		}

		void Prepend( char c )
		{
			ShiftRight() ;
			if ( nucToNum[c-'A'] == -1 )
			{
				invalidPos = kmerLength - 1 ;
			}
			code = ( code | ( (uint64_t)( nucToNum[c - 'A'] & 3 ) << ( 2ull * ( kmerLength - 1 ) ) ) ) & mask ;
		}

		inline void ShiftRight() 
		{
			if ( invalidPos != -1 )
				--invalidPos ;

			code = ( code >> 2ull ) & ( mask >> 2ull ) ;	

			if ( invalidPos < 0 )
				invalidPos = -1 ;
		}

		KmerCode& operator=( const KmerCode& in ) 
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			mask = in.mask ;
			code = in.code ;

			return *this ;
		}
} ;
#else
class KmerCode
{
	private:
		int kmerLength ;
		int invalidPos ; // The position contains characters other than A,C,G,T in the code
		int largestBlock ;
		uint64_t code[K_BLOCK_NUM] ;
		uint64_t mask ;
	public: 
		KmerCode() 
		{
		}
		KmerCode( int kl ) 
		{
			int i ;
			kmerLength = kl ;
			for ( i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = 0 ;
			invalidPos = -1 ;

			mask = 0 ;
			const int residual = ( kmerLength - 1 ) & ( sizeof( uint64_t ) * 4 - 1 ) ;
			for ( i = 0 ; i <= residual ; ++i )
			{
				mask = mask << 2 ;
				mask = mask | 3 ;
			}
			largestBlock = ( kmerLength - 1 ) / ( sizeof( uint64_t ) * 4 ) ;
		}

		KmerCode( const KmerCode& in )
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			largestBlock = in.largestBlock ;
			mask = in.mask ;
			for ( int i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = in.code[i] ;
		}

		void Restart() 
		{ 
			for ( int i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = 0 ;
			invalidPos = -1 ; 
		} 
		int GetCode( uint64_t outCode[] ) 
		{ 
			for ( int i = 0 ; i < K_BLOCK_NUM ; ++i )
				outCode[i] = code[i] ;
			return K_BLOCK_NUM ;
		} 

		uint64_t GetCodeI( int index )
		{
			return code[index] ;
		}
		int GetKmerLength() { return kmerLength ; }

		bool IsValid() 
		{
			/*if ( invalidPos != -1 )
			{
				printf( "hi\n") ;
			}*/
			return ( invalidPos == -1 ) ;
		}
		// the bits are organsized:
		// 4321 8765 ...
		void Append( char c ) 
		{
			ShiftLeft() ;
			code[0] |= ( (uint64_t)( nucToNum[ c - 'A' ] & 3 ) ) ;
			if ( nucToNum[c - 'A'] == -1 )
			{
				invalidPos = 0 ;
			}
		}

		void Prepend( char c )
		{
			ShiftRight() ;
			if ( nucToNum[c-'A'] == -1 )
			{
				invalidPos = kmerLength - 1 ;
			}
			const uint64_t residual = ( kmerLength - 1 ) % ( sizeof( uint64_t ) * 4 ) ;
			code[largestBlock] = ( code[largestBlock] | ( (uint64_t)( nucToNum[c - 'A'] & 3 ) << ( 2ull * ( residual ) ) ) ) & mask ;
		}

		void ShiftLeft()
		{
			if ( invalidPos != -1 )
				++invalidPos ;
			int i ;
			for ( i = largestBlock ; i >= 0 ; --i )
			{
				code[i] = ( code[i] << (2ull) ) ;
				if ( i > 0 )
				{
					uint64_t head = ( code[i - 1] >> (62ull ) ) & 3ull;
					code[i] |= head ;
				}
			}
			code[ largestBlock ] &= mask ;
			if ( invalidPos >= kmerLength )
				invalidPos = -1 ;
		}

		void ShiftRight()
		{
			if ( invalidPos != -1 )
				--invalidPos ;
			int i ;
			for ( i = 0 ; i <= largestBlock ; ++i )
			{
				if ( i > 0 )
				{
					uint64_t tail = code[i] & 3ull ;
					code[i - 1] |= ( tail << ( 2ull * 31 ) ) ;
				}
				code[i] = ( code[i] >> ( 2ull ) ) ;	
			}
			code[ largestBlock ] &= ( mask >> 2ull ) ;
			if ( invalidPos < 0 )
				invalidPos = -1 ;
		}

		KmerCode& operator=( const KmerCode& in ) 
		{
			kmerLength = in.kmerLength ;
			invalidPos = in.invalidPos ;
			largestBlock = in.largestBlock ;
			mask = in.mask ;
			for ( int i = 0 ; i < K_BLOCK_NUM ; ++i )
				code[i] = in.code[i] ;

			return *this ;
		}
} ;
#endif 

#endif 
