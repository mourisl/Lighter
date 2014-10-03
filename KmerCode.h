#ifndef _MOURISL_KMERCODE_HEADER
#define _MOURISL_KMERCODE_HEADER

#include <stdio.h>
#include "utils.h"

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

		void Append( char c ) ;
		void Prepend( char c ) ;
		void ShiftRight( int k ) ; 

		KmerCode& operator=( const KmerCode& in ) ;	
} ;

#endif 
