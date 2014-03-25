#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

#include "utils.h"
#include "Store.h"
#include "ErrorCorrection.h"
#include "Reads.h"
#include "KmerCode.h"
#include "GetKmers.h"
#include "pthread.h"

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

int MAX_CORRECTION ;
bool ALLOW_TRIMMING ;

void PrintHelp()
{
	printf( "Usage: ./lighter [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
		"\t-k kmer_length genome_size alpha\n"
		"Other parameters:\n"
		"\t-od output_file_directory (default: ./)\n"
		"\t-t number of threads to use (default: 1)\n"
		"\t-trim allow trimming (default: false)\n"
		"\t-all: output all the reads including those unfixable (default: false)\n"
		"\t-maxcor: the maximum number of correction for a read (default: 10)\n" ) ;
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

void PrintLog( const char *log )
{
	time_t rawtime ;
	struct tm *timeInfo ;
	char buffer[128] ;
	FILE *fp = fopen( "EC.log", "a" ) ;
	
	time( &rawtime ) ;
	timeInfo = localtime( &rawtime ) ;
	strftime( buffer, sizeof( buffer ), "%F %H:%M:%S", timeInfo ) ;

	fprintf( fp, "[%s] %s\n", buffer, log ) ;

	fclose( fp ) ;
}

int main( int argc, char *argv[] )
{
	int kmerLength ;
	double alpha = 0.0 ;
	char *readId/**, read, *qual*/ ;
	//char buffer[1023] ;
	double untrustF[100][100] ;
	//double trustF[100][100] ;
	int threshold[100] ;
	char goodQuality ;
	int badPrefix, badSuffix ;
	bool paraOutputAllReads ;
	//double bloomFilterFP = 0.0005 ;
	int i, j ;
	//uint64_t kmerCode ;
	//uint64_t mask ;
	uint64_t genomeSize = 0;

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

	paraOutputAllReads = false ; 
	MAX_CORRECTION = 10 ;
	ALLOW_TRIMMING = false ;
	kmerLength = -1 ;
	numOfThreads = 1 ;
	// Parse the arguments
	for ( i = 1 ; i < argc ; ++i )
	{
		if ( !strcmp( "-all", argv[i] ) )
		{
			paraOutputAllReads = true ;
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

	KmerCode kmerCode( kmerLength ) ;
	
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
	
	goodQuality = GetGoodQuality( reads ) ;
	reads.Rewind() ;

	//printf( "%c\n", goodQuality ) ;
	/*for ( i = 1 ; i <= kmerLength ; ++i )
	{
		for ( j = 0 ; j <= i ; ++j )
		{
			printf( "%.10lf\t", untrustF[i][j] ) ;
		}
		printf( "\n" ) ;
	}
	printf( "===============\n" ) ;
	for ( i = 1 ; i <= kmerLength ; ++i )
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
		
	PrintLog( "=============Start====================" ) ;

	// Step 1: Sample the kmers 
	//printf( "Begin step1. \n" ) ; fflush( stdout ) ;
	srand( 17 ) ;
	// It seems serialization is faster than parallel
	if ( 1 ) //numOfThreads == 1 )
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

		arg.kmerLength = kmerLength ;
		arg.alpha = alpha ;
		arg.kmers = &kmers ;
		// Since there is no output, so we can just directly read in the
		// sequence without considering the order.
		arg.reads = &reads ;
		arg.lock = &mutexSampleKmers ;
		//numOfThreads = 1 ;
		for ( i = 0 ; i < numOfThreads ; ++i )
		{
			pthread_create( &threads[i], &pthreadAttr, SampleKmers_Thread, (void *)&arg ) ;	
		}

		for ( i = 0 ; i < numOfThreads ; ++i )
		{
			pthread_join( threads[i], &pthreadStatus ) ;
		}
	}

	// Update the bloom filter's false positive rate.
	// Compute the distribution of the # of sampled kmers from untrusted and trusted position
	double tableAFP = kmers.GetFP() ;
	for ( i = 1 ; i <= kmerLength ; ++i )
	{
		int d = (int)( 0.05 / alpha * 2 );
		if ( d < 2 )
			d = 2 ;
		double p = 1 - pow( ( 1 - alpha ), d ) ;
		GetCumulativeBinomialDistribution( untrustF[i], i, p + tableAFP - p * tableAFP ) ;
		//GetCumulativeBinomialDistribution( untrustF[i], i, alpha ) ;
		//GetCumulativeBinomialDistribution( trustF[i], i, 1 - pow( ( 1 - alpha ), 20 ) ) ;
	}

	for ( i = 1 ; i <= kmerLength ; ++i )
	{
		for ( j = 0 ; j <= i ; ++j )
		{
			if ( untrustF[i][j] >= 1 - 0.5 * 1e-2 )
			{
				threshold[i] = j ;
				break ;
			}
		}
	}

	PrintLog( "Finish sampling kmers" ) ;
	// Step 2: Store the trusted kmers
	//printf( "Begin step2.\n") ; fflush( stdout ) ;
	reads.Rewind() ;
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() )
		{
			StoreTrustedKmers( reads.seq, reads.qual, kmerLength, goodQuality, threshold,
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
			int tmp = ErrorCorrection_Wrapper( reads.seq, kmerCode, &trustedKmers, badPrefix, badSuffix ) ;

			//if ( reads.HasQuality() )
			//	
			//else
			//	readId[0] = '>' ;

			reads.Output( tmp, badPrefix, badSuffix, ALLOW_TRIMMING ) ;
		}
	}
	else
	{
		int maxBatchSize = 128 * numOfThreads ;
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

			reads.OutputBatch( readBatch, batchSize, ALLOW_TRIMMING ) ;
		}

		free( readBatch ) ;	
	}
	PrintLog( "Finish error correction" ) ;

	return 0 ;
}
