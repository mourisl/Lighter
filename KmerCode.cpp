#include "KmerCode.h"
#include "utils.h"

extern char nucToNum[26] ;
extern char numToNuc[26] ;

void KmerCode::Append( char c )
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

void KmerCode::Prepend( char c )
{
	ShiftRight( 1 ) ;
	if ( nucToNum[c-'A'] == -1 )
	{
		invalidPos = kmerLength - 1 ;
	}
	code = ( code | ( (uint64_t)( nucToNum[c - 'A'] & 3 ) << ( 2ull * ( kmerLength - 1 ) ) ) ) & mask ;
}

void KmerCode::ShiftRight( int k )
{
	if ( invalidPos != -1 )
		invalidPos -= k ;

	code = ( code >> ( 2ull * k ) ) & ( mask >> ( 2ull * k ) ) ;	

	if ( invalidPos < 0 )
		invalidPos = -1 ;
}

KmerCode& KmerCode::operator=( const KmerCode& in )
{
	kmerLength = in.kmerLength ;
	invalidPos = in.invalidPos ;
	mask = in.mask ;
	code = in.code ;

	return *this ;
}
