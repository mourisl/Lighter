// A verification program to test the correction result from mason simulator.
// Usage: ./a.out *.fa/fq [OPTIONS]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char buffer[2048] ;	
char id[2048] ;
char oriSeq[2048] ;
char seq[2048] ;
char qual[2048] ;
char origSeq[2048] ;
char cigar[2048] ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

void ReverseComplement( char *s )
{
	char buffer[2048] ;
	int len = strlen( s ), i ;
	for ( i = 0 ; i < len ; ++i )
		buffer[i] = numToNuc[ 3 - nucToNum[ s[len - i - 1] - 'A' ] ] ;
	buffer[i] = '\0' ;
	strcpy( s, buffer ) ;
}

bool StrCompWithTrim( char *ref, char *s )
{
	int i ;
	for ( i = 0 ; s[i] && ref[i] ; ++i )
		if ( s[i] != ref[i] )
			break ;
	return s[i] != '\0' ;
}

char *FindIdColumn( char *id, const char *tag ) 
{
	char *p = strstr( id, tag ) ;
	if ( p == NULL )
		return p ;
	while ( *p && *p != '=' )
		++p ;
	++p ;
	//return ( p + strlen( tag ) + 1 ) ;
	return  p ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;
	int len ;
	char *p ;
	FILE *fp ;
	int FILE_TYPE ; // 0-fasta, 1-fastq
	int correctCount = 0, errorCount = 0 ;
	int sameCount = 0 ;
	int trimCount = 0, trimSum = 0 ;
	bool verbose = false ;
	int baseTP = 0, baseFP = 0, baseFN = 0 ;
	int readTP = 0, readFP = 0, readFN = 0 ;

	bool perBase = false ;
	int perBaseCorrect[2000] ;
	int perBaseTotal[2000] ;
	int largestReadLength = 0 ;

	for ( i = 2 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i],"-v" ) )
			verbose = true ;
		else if ( !strcmp( argv[i], "-per" ) )
			perBase = true ;
	}

	// Decide whether it is FASTQ or FASTA.
	fp = fopen( argv[1], "r" ) ;
	fscanf( fp, "%s", buffer ) ;
	if ( buffer[0] == '>' )
		FILE_TYPE = 0 ;
	else
		FILE_TYPE = 1 ;
	fclose( fp ) ;
	
	memset( perBaseCorrect, 0, sizeof( perBaseCorrect ) ) ;
	memset( perBaseTotal, 0, sizeof( perBaseTotal ) ) ;
	
	fp = fopen( argv[1], "r" ) ;
	while ( fgets( id, sizeof( id ), fp ) != NULL )
	{
		if ( FILE_TYPE == 0 )
		{
			fgets( seq, sizeof( seq ), fp ) ;
		}
		else if ( FILE_TYPE == 1 )
		{
			fgets( seq, sizeof( seq ), fp ) ;
			fgets( buffer, sizeof( buffer ), fp ) ;
			fgets( qual, sizeof( qual ), fp ) ;
		}
		//printf( "%s%s%s", id, seq,qual ) ;
		// Clean the return symbol	
		len = strlen( id ) ;		
		if ( id[len - 1] == '\n')
			id[len - 1] = '\0' ;
		len = strlen( seq ) ;
		if ( seq[len - 1] == '\n' )
			seq[len - 1] = '\0' ;
		if ( qual[len - 1] == '\n' )
			qual[len - 1] = '\0' ;
		
		// Parse the id field
		p = FindIdColumn( id, "haplotype_infix" ) ;
		sscanf( p, "%s", origSeq ) ;
	
		p = FindIdColumn( id, "edit_string" ) ;
		sscanf( p, "%s", cigar ) ;

		p = FindIdColumn( id, "strand=reverse" ) ;
		if ( p != NULL )
		{
			ReverseComplement( origSeq ) ;
		}
		if ( verbose )
			printf( "%s\n", id ) ;
		//printf( "%s %s\n", seq, origSeq ) ;
		if ( StrCompWithTrim( origSeq, seq ) )
		{
			/*if ( verbose )
			{
				printf( "Diff:\n%s\n%s\n", id, seq ) ;
			}*/
			++errorCount ;
			
			for ( i = 0 ; cigar[i] ; ++i )
			{
				if ( cigar[i] != 'M' )
					break ;
			}
			if ( cigar[i] )
			{
				if ( verbose )
					printf( "FN\n" ) ;
				++readFN ;
			}
			else
			{
				if ( verbose )
					printf( "FP\n" ) ;
				++readFP ;
			}
		}
		else
		{
			/*if ( verbose )
			{
				printf( "Same:\n%s\n%s\n", id, seq ) ;
			}*/
			//printf( "S\n" ) ;
			++correctCount ;
			
			for ( i = 0 ; cigar[i] ; ++i )
			{
				if ( cigar[i] != 'M' )
					break ;
			}
			if ( cigar[i] )
			{
				if ( verbose )
					printf( "TP\n" ) ;
				++readTP ;
			}
		}

		for ( i = 0 ; cigar[i] ; ++i )
		{
			if ( cigar[i] != 'M' )
				break ;
		}
		if ( !cigar[i] )
			++sameCount ;


		p = FindIdColumn( id, "trim" ) ;
		if ( p != NULL )
		{
			int tmp = atoi( p ) ;
			//printf( "%d %s\n", tmp, p ) ;
			++trimCount ;
			trimSum += tmp ;
		}

		// Collect information of TP, FP, FN for base level
		for ( i = 0 ; seq[i] ; ++i )
		{
			if ( cigar[i] == 'M' )
			{
				if ( seq[i] != origSeq[i] )
					++baseFP ;
			}
			else if ( cigar[i] == 'E' )
			{
				if ( seq[i] == origSeq[i] )
					++baseTP ;
				else
					++baseFN ;
			}

			if ( seq[i] == origSeq[i] )
				++perBaseCorrect[i] ;
			++perBaseTotal[i] ;
			if ( i > largestReadLength )
				largestReadLength = i ;
		}
 	}
	printf( "correct #: %d\n"
		"error #: %d\n", correctCount, errorCount ) ;
	printf( "Original Correct Reads Count: %d\n", sameCount ) ;
	printf( "Trimmed Reads Count: %d. Average trim length: %lf\n", trimCount, (double)trimSum / trimCount ) ;
	printf( "\nBase level:\n" ) ;
	printf( "TP: %d\nFP: %d\nFN: %d\n", baseTP, baseFP, baseFN ) ;
	double recall = ( double )baseTP/(baseTP+baseFN) ;
	double precision = (double)baseTP/(baseTP+baseFP) ;
	printf( "Recall: %lf\n"
		"Precision: %lf\n"
		"F-score: %lf\n"
		"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
				(double)(baseTP-baseFP)/(baseTP+baseFN) ) ;
	
	printf( "\nRead level:\n" ) ;
	printf( "TP: %d\nFP: %d\nFN: %d\n", readTP, readFP, readFN ) ;
	recall = ( double )readTP/(readTP+readFN) ;
	precision = (double)readTP/(readTP+readFP) ;
	printf( "Recall: %lf\n"
			"Precision: %lf\n"
			"F-score: %lf\n"
			"Gain: %lf\n", recall, precision, 2*recall*precision / ( recall + precision ),
			(double)(readTP-readFP)/(readTP+readFN) ) ;
	printf( "\n") ;

	if ( perBase )
	{
		printf( "\nPer-base evaluation:\n" ) ;
		for ( i = 0 ; i < largestReadLength ; ++i )
			printf( "%d %lf\n", i, (double)perBaseCorrect[i]/perBaseTotal[i] ) ;
	}
	return 0 ;
}
