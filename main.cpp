#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "utils.h"
#include "Store.h"
#include "ErrorCorrection.h"
#include "Reads.h"
#include "KmerCode.h"
#include "GetKmers.h"
#include "pthread.h"

#define READ_BUFFER_PER_THREAD 1024

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

int MAX_CORRECTION ;
bool ALLOW_TRIMMING ;

struct _summary
{
	uint64_t corrCnt ;
	uint64_t trimReadsCnt ;
	uint64_t trimBaseCnt ;
	uint64_t discardReadsCnt ;
	uint64_t totalReads ;
	uint64_t errorFreeReadsCnt ;
} ;

void PrintHelp()
{
	printf( "Usage: ./lighter [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
		"\t-k kmer_length genome_size alpha\n"
		"Other parameters:\n"
		"\t-od: output_file_directory (default: ./)\n"
		"\t-t: number of threads to use (default: 1)\n"
		"\t-trim: allow trimming (default: false)\n"
		"\t-discard: discard unfixable reads. Will LOSE paired-end matching when discarding (default: false)\n"
		"\t-noQual: ignore the quality socre (default: false)\n"
		"\t-stable: sequentialize the sampling stage, output the same result with different runs (default: false)\n"
		"\t-maxcor: the maximum number of correction for within a kmer_length window (default: 4)\n" ) ;
}

uint64_t StringToUint64( char *s )   
{
	int i ;
	uint64_t ret = 0 ;
	for ( i = 0 ; s[i] ; ++i )
	{
		ret = ret * 10 + s[i] - '0' ;
	}
	return ret ;
}

inline void ExtractKmer( char *s, int offset, int kmerLength, char *buffer )
{
	int i ;
	for ( i = offset ; i < offset+ kmerLength ; ++i )
		buffer[i - offset] = s[i] ;
	buffer[i - offset] = '\0' ;
}

void GetCumulativeBinomialDistribution( double F[], int l, double p )
{
	// p is the probability of getting 1.
	int i ;
	double coef = 1 ;
	double exp = pow( 1 - p, l ) ;
	F[0] = pow( 1 - p, l ) ;

	for ( i = 1 ; i <= l ; ++i )
	{
		coef = coef / i * ( l - i + 1 ) ;
		exp = exp / ( 1 - p ) * p ;
		F[i] = F[i - 1] + coef * exp ;		
	}
}

char GetGoodQuality( Reads &reads )
{
	int i ;
	int qualHisto[300] ;
	int totalCnt, cnt ;

	//Reads reads( readFile ) ;
	if ( !reads.HasQuality() )
		return 127 ;

	memset( qualHisto, 0, sizeof( qualHisto ) ) ;
	for ( i = 0 ; i < 1000000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		++qualHisto[ (int)reads.qual[0] ] ;
	}

	totalCnt = i ;
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += qualHisto[i] ;
		if ( cnt > totalCnt * 0.25 )
			break ;
	}
	return (char)i ;
}

char GetBadQuality( Reads &reads )
{
	int i ;
	int qualHisto[300], firstQualHisto[300] ;
	int totalCnt, cnt ;
	int t1, t2 ;
	//Reads reads( readFile ) ;
	if ( !reads.HasQuality() )
		return 0 ;

	memset( qualHisto, 0, sizeof( qualHisto ) ) ;
	memset( firstQualHisto, 0, sizeof( firstQualHisto )) ;
	for ( i = 0 ; i < 1000000 ; ++i )
	{
		if ( !reads.Next() )
			break ;
		++qualHisto[ (int)reads.qual[ strlen( reads.seq ) - 1 ] ] ;
		++firstQualHisto[ (int)reads.qual[0] ] ;
	}

	totalCnt = i ;
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += firstQualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t1 = i - 1 ;
	
	cnt = 0 ;
	for ( i = 0 ; i < 300 ; ++i )
	{
		cnt += qualHisto[i] ;
		if ( cnt > totalCnt * 0.05 )
			break ;
	}
	t2 = i ;
	//printf( "%d %d\n", t1, t2 ) ;
	return (char)( t2 < t1 ? t2 : t1 ) ;
}

void PrintLog( const char *log )
{
	time_t rawtime ;
	struct tm *timeInfo ;
	char buffer[128] ;
	//FILE *fp = fopen( "lighter.log", "a" ) ;
	
	time( &rawtime ) ;
	timeInfo = localtime( &rawtime ) ;
	strftime( buffer, sizeof( buffer ), "%F %H:%M:%S", timeInfo ) ;

	printf( "[%s] %s\n", buffer, log ) ;

	//fclose( fp ) ;
}

void UpdateSummary( char *seq, int correction, int badSuffix, bool paraDiscard, struct _summary &summary ) 
{
	if ( correction == 0 )
		++summary.errorFreeReadsCnt ;			
	else if ( correction > 0 )
		++summary.corrCnt ;	
	else if ( paraDiscard ) // tmp < 0
		++summary.discardReadsCnt ;

	if ( ALLOW_TRIMMING && badSuffix > 0 )
	{
		++summary.trimReadsCnt ;
		summary.trimBaseCnt += badSuffix ;
	}
	++summary.totalReads ;	
}

void PrintSummary( const struct _summary &summary )
{
	printf( "Processed %" PRIu64 " reads:\n"
		"\t%" PRIu64 " are error-free\n"
		"\tCorrected %" PRIu64 " bases(%lf corrections for reads with errors)\n"
		"\tTrimmed %" PRIu64 " reads with average trimmed bases %lf\n"
		"\tDiscard %" PRIu64 " reads\n",
		summary.totalReads, summary.errorFreeReadsCnt,
		summary.corrCnt, 
		summary.totalReads == summary.errorFreeReadsCnt ? 0.0 : 
					(double)summary.corrCnt / ( summary.totalReads - summary.errorFreeReadsCnt ),
		summary.trimReadsCnt, 
		summary.trimReadsCnt == 0 ? 0.0 : (double)summary.trimBaseCnt / summary.trimReadsCnt, 
		summary.discardReadsCnt ) ;	
}

int main( int argc, char *argv[] )
{
	int kmerLength ;
	double alpha = 0.0 ;
	char *readId/**, read, *qual*/ ;
	char buffer[1023] ;
	double untrustF[100][100] ;
	//double trustF[100][100] ;
	int threshold[100] ;
	char goodQuality = '\0', badQuality = '\0' ;
	int badPrefix, badSuffix ;
	bool paraDiscard ;
	bool ignoreQuality, stable ;
	//double bloomFilterFP = 0.0005 ;
	int i, j ;
	//uint64_t kmerCode ;
	//uint64_t mask ;
	uint64_t genomeSize = 0;
	struct _summary summary ;

	struct _SamplePattern *samplePatterns = NULL ;

	// variables for threads
	int numOfThreads ;
	pthread_attr_t pthreadAttr ;
	pthread_t *threads = NULL;
	pthread_mutex_t mutexSampleKmers, mutexStoreKmers ;

	if ( argc == 1 || !strcmp( "-h", argv[1] ) )
	{
		PrintHelp() ;
		exit( 0 ) ;
	}

	Reads reads ;
	/*reads.AddReadFile( argv[1] ) ;
	kmerLength = atoi( argv[2] ) ;
	genomeSize = StringToUint64( argv[3] ) ;
	alpha = (double)atof( argv[4] ) ;*/

	paraDiscard = false ; 
	MAX_CORRECTION = 4 ;
	ALLOW_TRIMMING = false ;
	kmerLength = -1 ;
	numOfThreads = 1 ;
	ignoreQuality = false ;
	stable = false ;
	memset( &summary, 0, sizeof( summary ) ) ;
	
	// Parse the arguments
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-discard", argv[i] ) )
		{
			paraDiscard = true ;
		}
		else if ( !strcmp( "-maxcor", argv[i] ) )
		{
			MAX_CORRECTION = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-r", argv[i] ) )
		{
			//reads.AddReadFile( argv[i + 1 ] ) ;
			++i;
			continue ; // wait to be processed after next round
		}
		else if ( !strcmp( "-k", argv[i] ) )
		{
			if(i + 1 >= argc) {
				printf("Must specify k-mer length, genome size, and alpha after -k\n");
				exit(1);
			}
			kmerLength = atoi( argv[i + 1] ) ;
			if(i + 2 >= argc) {
				printf("Must specify k-mer length, genome size, and alpha after -k\n");
				exit(1);
			}
			genomeSize = StringToUint64( argv[i + 2] ) ;
			if(i + 3 >= argc) {
				printf("Must specify k-mer length, genome size, and alpha after -k\n");
				exit(1);
			}
			alpha = (double)atof( argv[i + 3] ) ;
			i += 3 ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			reads.SetOutputDirectory( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-trim", argv[i] ) )
		{
			ALLOW_TRIMMING = true ;
		}
		else if ( !strcmp( "-t", argv[i] ) )
		{
			numOfThreads = atoi( argv[i + 1] ) ;
			++i ;
		}
		else if ( !strcmp( "-noQual", argv[i] ) )
		{
			ignoreQuality = true ;
		}
		else if ( !strcmp( "-stable", argv[i] ) )
		{
			stable = true ;
		}
		else
		{
			printf( "Unknown argument %s\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	// Go the second round to get the reads files
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-r", argv[i] ) )
		{
			reads.AddReadFile( argv[i + 1 ] ) ;
			++i;
		}
	}
	if ( kmerLength == -1 )
	{
		printf( "Require -k parameter!\n" ) ;
		exit( 1 ) ;
	}

	PrintLog( "=============Start====================" ) ;
	KmerCode kmerCode( kmerLength ) ;
	reads.SetDiscard( paraDiscard ) ;	

	// Prepare data structures and other data.
	//Store kmers(1000000000ull) ;
	//Store trustedKmers(1000000000ull) ;
	Store kmers((uint64_t)genomeSize * 1.5, 0.01 ) ;
	Store trustedKmers((uint64_t)genomeSize * 1.5, 0.0005 ) ;


	if ( numOfThreads > 1 )
	{
		// Initialized pthread variables
		pthread_attr_init( &pthreadAttr ) ;
		pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;
		threads = ( pthread_t * )malloc( sizeof( pthread_t ) * numOfThreads ) ;
		pthread_mutex_init( &mutexSampleKmers, NULL ) ;
		pthread_mutex_init( &mutexStoreKmers, NULL ) ;

		//kmers.SetNumOfThreads( numOfThreads ) ;
		trustedKmers.SetNumOfThreads( numOfThreads ) ;
	}
	
	//goodQuality = GetGoodQuality( reads ) ;
	//reads.Rewind() ;
	if ( ignoreQuality == false )
		badQuality = GetBadQuality( reads ) ;
	if ( badQuality != '\0' )
	{
		sprintf( buffer, "Bad quality threshold is %c", badQuality ) ;
		PrintLog( buffer ) ;
	}
	else
	{
		PrintLog( "No quality score used." ) ;
	}
	reads.Rewind() ;

	//printf( "%c\n", badQuality ) ;
	//exit( 1 ) ;
		/*for ( i = 1 ; i <= kmerLength ; ++i )
	{
		for ( j = 0 ; j <= i ; ++j )
		{
			printf( "%.10lf\t", trustF[i][j] ) ;
		}
		printf( "\n" ) ;
	}
	exit( 1 ) ;*/
	//Store kmers((uint64_t)50000000 * 4, 0.001 ) ;
	//Store trustedKmers((uint64_t)50000000 * 2, 0.001 ) ;
		
	
	// Step 1: Sample the kmers 
	//printf( "Begin step1. \n" ) ; fflush( stdout ) ;
	srand( 17 ) ;
	// Build the patterns for sampling
	if ( numOfThreads > 1 && stable == false )
	{
		samplePatterns = ( struct _SamplePattern *)malloc( sizeof( *samplePatterns ) * SAMPLE_PATTERN_COUNT ) ;

		for ( i = 0 ; i < SAMPLE_PATTERN_COUNT ; ++i )
		{
			int k ;
			for ( k = 0 ; k < MAX_READ_LENGTH / 8 ; ++k )
				samplePatterns[i].tag[k] = 0 ;
			for ( k = 0 ; k < MAX_READ_LENGTH ; ++k )
			{
				double p = rand() / (double)RAND_MAX;
				if ( p < alpha )
				{
					samplePatterns[i].tag[ k / 8 ] |= ( 1 << ( k % 8 ) ) ; // Notice within the small block, the order is reversed
				}
			}
		}
	}
	// It seems serialization is faster than parallel
	if ( numOfThreads == 1 || stable == true )
	{
		while ( reads.Next() != 0 )
		{
			SampleKmersInRead( reads.seq, reads.qual, kmerLength, alpha, kmerCode, &kmers ) ;
		}
	}
	else
	{
		struct _SampleKmersThreadArg arg ;
		void *pthreadStatus ;
		kmers.SetNumOfThreads( numOfThreads ) ;

		arg.kmerLength = kmerLength ;
		arg.alpha = alpha ;
		arg.kmers = &kmers ;
		arg.samplePatterns = samplePatterns ;
		// Since there is no output, so we can just directly read in the
		// sequence without considering the order.
		arg.reads = &reads ;
		arg.lock = &mutexSampleKmers ;
		//arg.lockPut = &mutexSampleKmersPut ;

		//numOfThreads = 1 ;
		for ( i = 0 ; i < numOfThreads / 2 ; ++i )
		{
			pthread_create( &threads[i], &pthreadAttr, SampleKmers_Thread, (void *)&arg ) ;	
		}

		for ( i = 0 ; i < numOfThreads / 2 ; ++i )
		{
			pthread_join( threads[i], &pthreadStatus ) ;
		}
	}
	if ( numOfThreads > 1 && stable == false )
		free( samplePatterns ) ;

	// Update the bloom filter's false positive rate.
	// Compute the distribution of the # of sampled kmers from untrusted and trusted position
	double tableAFP = kmers.GetFP() ;
	for ( i = 1 ; i <= kmerLength ; ++i )
	{
		int d = (int)( 0.1 / alpha * 2 );
		double p ;
		if ( d < 2 )
			d = 2 ;
		p = 1 - pow( ( 1 - alpha ), d ) ;
		//else
		//	p = 1 - pow( 1 - 0.05, 2 ) ;
		//double p = 1 - pow( 1 - 0.05, 2 ) ;
		//p = 0 ;
		GetCumulativeBinomialDistribution( untrustF[i], i, p + tableAFP - p * tableAFP ) ;
		//GetCumulativeBinomialDistribution( untrustF[i], i, tableAFP ) ;
		//GetCumulativeBinomialDistribution( untrustF[i], i, alpha ) ;
		//GetCumulativeBinomialDistribution( trustF[i], i, 1 - pow( ( 1 - alpha ), 20 ) ) ;
	}

	/*for ( i = 1 ; i <= kmerLength ; ++i )
	{
		for ( j = 0 ; j <= i ; ++j )
		{
			printf( "%.20lf\t", untrustF[i][j] ) ;
		}
		printf( "\n" ) ;
	}
	printf( "===============\n" ) ;
	exit( 1 ) ;*/


	for ( i = 1 ; i <= kmerLength ; ++i )
	{
		for ( j = 0 ; j <= i ; ++j )
		{
			if ( untrustF[i][j] >= 1 - 0.5 * 1e-2 )
			{
				threshold[i] = j ;
				//if ( threshold[i] <= i / 2 )
				//	threshold[i] = i / 2 ;
				break ;
			}
		}
	}
	/*for ( i = 1 ; i <= kmerLength ; ++i )
	{
		printf( "%d %d\n", threshold[i] + 1, i ) ;
	}
	exit( 1 ) ;*/
	PrintLog( "Finish sampling kmers" ) ;
	
	sprintf( buffer, "Bloom filter A's error rate: %lf", tableAFP ) ;
	PrintLog( buffer ) ;
	// Step 2: Store the trusted kmers
	//printf( "Begin step2.\n") ; fflush( stdout ) ;
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() )
		{
			StoreTrustedKmers( reads.seq, reads.qual, kmerLength, badQuality, threshold,
					kmerCode, &kmers, &trustedKmers ) ;
		}
	}
	else
	{
		struct _StoreKmersThreadArg arg ;
		void *pthreadStatus ;

		arg.kmerLength = kmerLength ;
		arg.threshold = threshold ;
		arg.kmers = &kmers ;
		arg.trustedKmers = &trustedKmers ;
		arg.reads = &reads ;
		arg.goodQuality = goodQuality ;
		arg.badQuality = badQuality ;
		arg.lock = &mutexStoreKmers ;
		
		for ( i = 0 ; i < numOfThreads ; ++i )
		{
			pthread_create( &threads[i], &pthreadAttr, StoreKmers_Thread, (void *)&arg ) ;
		}

		for ( i = 0 ; i < numOfThreads ; ++i )
		{
			pthread_join( threads[i], &pthreadStatus ) ;
		}
	}
	PrintLog( "Finish storing trusted kmers" ) ;
	// Step 3: error correction
	//printf( "%lf %lf\n", kmers.GetFP(), trustedKmers.GetFP() ) ;
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() )
		{
			readId = reads.id ;
			//read = reads.seq ;
			/*kmerCode = 0 ;
			  for ( i = 0 ; i < kmerLength ; ++i )
			  {
			  kmerCode = kmerCode << (uint64_t)2 ;
			  kmerCode = kmerCode | (uint64_t)nucToNum[ read[i] - 'A' ] ;
			  }
			  if ( !trustedKmers.IsIn( kmerCode ) )
			  {
			//printf( "- %d %lld\n", i, kmerCode ) ;
			printf( "%s\n%s\n", readId, read ) ;
			continue ;
			}
			for ( ; read[i] ; ++i )
			{
			kmerCode = ( kmerCode << (uint64_t)2 ) & mask ;
			kmerCode = kmerCode | (uint64_t)nucToNum[ read[i] - 'A' ] ;

			if ( !trustedKmers.IsIn( kmerCode ) )
			{
			printf( "%s\n%s\n", readId, read ) ;
			break ;
			}
			}
			continue ;*/
			int tmp = ErrorCorrection_Wrapper( reads.seq, reads.qual, kmerCode, badQuality, &trustedKmers, badPrefix, badSuffix ) ;

			//if ( reads.HasQuality() )
			//	
			//else
			//	readId[0] = '>' ;
			UpdateSummary( reads.seq, tmp, badSuffix, paraDiscard, summary ) ;			

			reads.Output( tmp, badPrefix, badSuffix, ALLOW_TRIMMING ) ;
		}
	}
	else
	{
		int maxBatchSize = READ_BUFFER_PER_THREAD * numOfThreads ;
		int batchSize ;
		struct _ErrorCorrectionThreadArg arg ;
		pthread_mutex_t errorCorrectionLock ;
		void *pthreadStatus ;

		struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		pthread_mutex_init( &errorCorrectionLock, NULL ) ;

		arg.kmerLength = kmerLength ;
		arg.trustedKmers = &trustedKmers ;
		arg.readBatch = readBatch ;
		arg.lock = &errorCorrectionLock ;
		arg.badQuality = badQuality ;

		while ( 1 )
		{
			batchSize = reads.GetBatch( readBatch, maxBatchSize, true, true ) ;
			if ( batchSize == 0 )
				break ; 
			//printf( "batchSize=%d\n", batchSize ) ;
			arg.batchSize = batchSize ;
			arg.batchUsed = 0 ;
			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_create( &threads[i], &pthreadAttr, ErrorCorrection_Thread, (void *)&arg ) ;	

			for ( i = 0 ; i < numOfThreads ; ++i )
				pthread_join( threads[i], &pthreadStatus ) ;
			
			for ( i = 0 ; i < batchSize ; ++i )
				UpdateSummary( readBatch[i].seq, readBatch[i].correction, readBatch[i].badSuffix, paraDiscard, summary ) ;			
			reads.OutputBatch( readBatch, batchSize, ALLOW_TRIMMING ) ;
		}

		free( readBatch ) ;	
	}
	

	PrintLog( "Finish error correction" ) ;
	PrintSummary( summary ) ;
	return 0 ;
}
