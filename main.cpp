#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sys/stat.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "utils.h"
#include "Store.h"
#include "ErrorCorrection.h"
#include "Reads.h"
#include "KmerCode.h"
#include "GetKmers.h"
#include "pthread.h"


char LIGHTER_VERSION[] = "Lighter v1.1.2" ;

char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;

int MAX_CORRECTION ;
bool ALLOW_TRIMMING ;
int SET_NEW_QUAL ;
bool zlibVersionChecked = false ; 

struct _summary
{
	uint64_t corrCnt ;
	uint64_t trimReadsCnt ;
	uint64_t trimBaseCnt ;
	uint64_t discardReadsCnt ;
	uint64_t totalReads ;
	uint64_t errorFreeReadsCnt ;
} ;

struct _OutputThreadArg
{
	struct _summary *summary ;
	struct _Read *readBatch ;
	int batchSize ;
	bool paraDiscard ;

	Reads *reads ;
	int fileInd ;
} ;

void PrintHelp()
{
	printf( "Usage: ./lighter [OPTIONS]\n"
		"OPTIONS:\n"
		"Required parameters:\n"
		"\t-r seq_file: seq_file is the path to the sequence file. Can use multiple -r to specifiy multiple sequence files\n"
	        "\t             The file can be fasta and fastq, and can be gzip\'ed with extension *.gz.\n"      
		"\t             When the input file is *.gz, the corresponding output file will also be gzip\'ed.\n"
		"\t-k kmer_length genome_size alpha: (see README for information on setting alpha)\n"
		"\t\t\t\t\tor\n"
		"\t-K kmer_length genom_size: in this case, the genome size should be relative accurate.\n"
		"Other parameters:\n"
		"\t-od output_file_directory: (default: ./)\n"
		"\t-t num_of_threads: number of threads to use (default: 1)\n"
		"\t-maxcor INT: the maximum number of corrections within a 20bp window (default: 4)\n"
		"\t-trim: allow trimming (default: false)\n"
		"\t-discard: discard unfixable reads. Will LOSE paired-end matching when discarding (default: false)\n"
		"\t-noQual: ignore the quality socre (default: false)\n"
		"\t-newQual ascii_quality_score: set the quality for the bases corrected to the specified score (default: not used)\n"
		//"\t-stable: sequentialize the sampling stage, output the same result with different runs (default: false)\n"
		"\t-saveTrustedKmers file: save the trusted kmers to specified file then stop (default: not used)\n"
		"\t-loadTrustedKmers file: directly get solid kmers from specified file (default: not used)\n"
		"\t-zlib compress_level: set the compression level(0-9) of gzip (default: 1)\n"
		"\t-h: print the help message and quit\n"
		"\t-v: print the version information and quit\n") ;
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
	return (char)( t2 < t1 ? t2 : t1 ) ;
}

double InferAlpha( Reads &reads, uint64_t genomeSize ) 
{
	uint64_t totalLen = 0 ;

	while ( reads.Next() )
		totalLen += strlen( reads.seq ) ;		
	
	return 7.0 / ( (double)totalLen / genomeSize ) ;
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

	fprintf( stderr, "[%s] %s\n", buffer, log ) ;

	//fclose( fp ) ;
}

void UpdateSummary( char *seq, int correction, int badSuffix, bool paraDiscard, struct _summary &summary ) 
{
	if ( correction == 0 )
		++summary.errorFreeReadsCnt ;			
	else if ( correction > 0 )
		summary.corrCnt += correction ;	
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
	fprintf( stderr, "Processed %" PRIu64 " reads:\n"
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


void *Output_Thread( void *arg )
{
	int i ;
	struct _OutputThreadArg *myArg = ( struct _OutputThreadArg *)arg ;
	struct _Read *readBatch = myArg->readBatch ;
	int batchSize = myArg->batchSize ;
	
	//fprintf( stderr, "hi\n"  ) ;
	for ( i = 0 ; i < batchSize ; ++i )
		UpdateSummary( readBatch[i].seq, readBatch[i].correction, readBatch[i].badSuffix, myArg->paraDiscard, *( myArg->summary ) ) ;			
	myArg->reads->OutputBatch( readBatch, batchSize, ALLOW_TRIMMING, myArg->fileInd ) ;
	pthread_exit( NULL ) ;
	return NULL ;
}

int main( int argc, char *argv[] )
{
	int kmerLength ;
	double alpha = -1 ;
	//char *readId/**, read, *qual*/ ;
	char buffer[1023] ;
	double untrustF[100][100] ;
	//double trustF[100][100] ;
	int threshold[MAX_KMER_LENGTH+1] ;
	char goodQuality = '\0', badQuality = '\0' ;
	int badPrefix, badSuffix ;
	bool paraDiscard ;
	bool ignoreQuality, inferAlpha ; //stable ;
	int zlibLevel ;
	char *saveTrustedKmers, *loadTrustedKmers ;

	//double bloomFilterFP = 0.0005 ;
	int i, j ;
	//uint64_t kmerCode ;
	//uint64_t mask ;
	uint64_t genomeSize = 0;
	struct _summary summary ;

	struct _SamplePattern *samplePatterns = NULL ;
	bool setMaxCor ;

	// variables for threads
	int numOfThreads ;
	pthread_attr_t pthreadAttr ;
	pthread_t *threads = NULL;
	pthread_mutex_t mutexSampleKmers, mutexStoreKmers ;

	if ( argc == 1 )
	{
		PrintHelp() ;
		exit( EXIT_FAILURE ) ;
	}

	Reads reads ;
	/*reads.AddReadFile( argv[1] ) ;
	kmerLength = atoi( argv[2] ) ;
	genomeSize = StringToUint64( argv[3] ) ;
	alpha = (double)atof( argv[4] ) ;*/

	paraDiscard = false ; 
	MAX_CORRECTION = 4 ;
	setMaxCor = false ;
	ALLOW_TRIMMING = false ;
	SET_NEW_QUAL = -1 ;
	kmerLength = -1 ;
	numOfThreads = 1 ;
	ignoreQuality = false ;
	//stable = false ;
	inferAlpha = false ;
	loadTrustedKmers = NULL ;
	saveTrustedKmers = NULL ;
	zlibLevel = 1 ;
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
			setMaxCor = true ;
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
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			kmerLength = atoi( argv[i + 1] ) ;
			if(i + 2 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			genomeSize = StringToUint64( argv[i + 2] ) ;
			if(i + 3 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size, and alpha after -k\n");
				exit( EXIT_FAILURE );
			}
			alpha = (double)atof( argv[i + 3] ) ;
			i += 3 ;
		}
		else if ( !strcmp( "-K", argv[i] ) ) 
		{
			if(i + 1 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size after -K\n");
				exit( EXIT_FAILURE );
			}
			kmerLength = atoi( argv[i + 1] ) ;
			if(i + 2 >= argc) {
				fprintf( stderr, "Must specify k-mer length, genome size after -K\n");
				exit( EXIT_FAILURE );
			}
			genomeSize = StringToUint64( argv[i + 2] ) ;

			inferAlpha = true ;
			i += 2 ;
		}
		else if ( !strcmp( "-od", argv[i] ) )
		{
			mkdir( argv[i + 1], 0700 ) ;
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
		else if ( !strcmp( "-newQual", argv[i] ) )
		{
			SET_NEW_QUAL = (int)argv[i + 1][0] ;
			++i ;
		}
		else if ( !strcmp( "-stable", argv[i] ) )
		{
			//stable = true ;
		}
		else if ( !strcmp( "-saveTrustedKmers", argv[i] ) )
		{
			saveTrustedKmers = argv[i + 1] ;
			++i ;
		}
		else if ( !strcmp( "-loadTrustedKmers", argv[i] ) )
		{
			loadTrustedKmers = argv[i + 1] ;
			++i ;
		}
		else if ( !strcmp( "-zlib", argv[i] ) )
		{
			zlibLevel = atoi( argv[i+1] ) ;
			reads.SetCompressLevel( zlibLevel ) ;
			++i ;
		}
		else if ( !strcmp( "-h", argv[i] ) )
		{
			PrintHelp() ;
			exit( 0 ) ;
		}
		else if ( !strcmp( "-v", argv[i] ) )
		{
			printf( "%s\n", LIGHTER_VERSION ) ;
			exit( 0 ) ;
		}
		else
		{
			fprintf( stderr, "Unknown argument %s\n", argv[i] ) ;
			exit( EXIT_FAILURE ) ;
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
		fprintf( stderr, "Require -k or -K parameter!\n" ) ;
		exit( EXIT_FAILURE ) ;
	}
	if ( kmerLength > MAX_KMER_LENGTH )
	{
		fprintf( stderr, "K-mer length must be no larger than %d. You can adjust the MAX_KMER_LENGTH constraints in utils.h.\n", MAX_KMER_LENGTH ) ;
		exit( EXIT_FAILURE ) ;
	}
	
	if ( alpha != -1 && inferAlpha == true )
	{
		fprintf( stderr, "Can not use both -k and -K.\n" ) ;
		exit( EXIT_FAILURE ) ;
	}

	if ( loadTrustedKmers != NULL && saveTrustedKmers != NULL )
	{
		fprintf( stderr, "Can't use both -saveTrustedKmers and -loadTrustedKmers at the same time.\n" ) ;
		exit( EXIT_FAILURE ) ;
	}

	PrintLog( "=============Start====================" ) ;
	KmerCode kmerCode( kmerLength ) ;
	reads.SetDiscard( paraDiscard ) ;	

	if ( inferAlpha && loadTrustedKmers == NULL )
	{
		PrintLog( "Scanning the input files to infer alpha(sampling rate)" ) ;
		alpha = InferAlpha( reads, genomeSize ) ;
		
		sprintf( buffer, "Average coverage is %.3lf and alpha is %.3lf", 7.0 / alpha, alpha ) ;
		PrintLog( buffer ) ;
		
		reads.Rewind() ;
	}

	// Prepare data structures and other data.
	//Store kmers(1000000000ull) ;
	//Store trustedKmers(1000000000ull) ;
	Store kmers((uint64_t)( genomeSize * 1.5 ), 0.01 ) ;
	Store trustedKmers((uint64_t)( genomeSize * 1.5 ), 0.0005 ) ;
	

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
		sprintf( buffer, "Bad quality threshold is \"%c\"", badQuality ) ;
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
		
	if ( loadTrustedKmers != NULL )
	{
		trustedKmers.BloomInput( loadTrustedKmers ) ;
		sprintf( buffer, "Finish loading trusted kmers from file %s.", loadTrustedKmers ) ;
		PrintLog( buffer ) ;
	}
if ( loadTrustedKmers == NULL ) // a very long if state-ment, I avoid the indent here to regard this a macro.
{
	// Step 1: Sample the kmers 
	//printf( "Begin step1. \n" ) ; fflush( stdout ) ;
	srand( 17 ) ;
	// Build the patterns for sampling
	if ( numOfThreads > 1 )//&& stable == false )
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
	// It seems serialization is faster than parallel. NOT true now!
	if ( numOfThreads == 1 ) //|| stable == true )
	{
		while ( reads.Next() != 0 )
		{
			SampleKmersInRead( reads.seq, reads.qual, kmerLength, alpha, kmerCode, &kmers ) ;
		}
	}
	else //if ( 0 ) 
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
	if ( numOfThreads > 1 ) //&& stable == false )
		free( samplePatterns ) ;
	//kmers.BloomInput( "sample_bf.out" ) ;
	//kmers.BloomOutput( "sample_bf.out" ) ;

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
	
	sprintf( buffer, "Bloom filter A's false positive rate: %lf", tableAFP ) ;
	PrintLog( buffer ) ;

	if ( setMaxCor == false && tableAFP > 0.1 )
	{
		++MAX_CORRECTION ;
		if ( badQuality != '\0' )
		{
			++badQuality ;
			sprintf( buffer, "The error rate is high. Lighter adjusts -maxcor to %d and bad quality threshold to \"%c\".", MAX_CORRECTION, badQuality ) ;
		}
		else
			sprintf( buffer, "The error rate is high. Lighter adjusts -maxcor to %d.", MAX_CORRECTION ) ;
		PrintLog( buffer ) ;
	}
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
	else //if ( 0 )
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
	if ( saveTrustedKmers != NULL )
	{
		trustedKmers.BloomOutput( saveTrustedKmers ) ;
		sprintf( buffer, "The trusted kmers are saved in file %s.", saveTrustedKmers ) ;
		PrintLog( buffer ) ;
		return 0 ;
	}
}
	//trustedKmers.BloomInput( "bf.out ") ;
	//trustedKmers.BloomOutput( "bf.out ") ;

	// Step 3: error correction
	//printf( "%lf %lf\n", kmers.GetFP(), trustedKmers.GetFP() ) ;
	reads.Rewind() ;
	// Different ways of parallel depending on the number of threads.
	if ( numOfThreads == 1 )
	{
		while ( reads.Next() )
		{
			//readId = reads.id ;
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
			int info ;
			int tmp = ErrorCorrection_Wrapper( reads.seq, reads.qual, kmerCode, badQuality, &trustedKmers, badPrefix, badSuffix, info ) ;

			//if ( reads.HasQuality() )
			//	
			//else
			//	readId[0] = '>' ;
			UpdateSummary( reads.seq, tmp, badSuffix, paraDiscard, summary ) ;			

			reads.Output( tmp, badPrefix, badSuffix, info, ALLOW_TRIMMING ) ;
		}
	}
	else if ( numOfThreads == 2 )
	{
		int maxBatchSize = READ_BUFFER_PER_THREAD * numOfThreads ;
		int batchSize ;
		int fileInd ;
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
			batchSize = reads.GetBatch( readBatch, maxBatchSize, fileInd, true, true ) ;
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
			reads.OutputBatch( readBatch, batchSize, ALLOW_TRIMMING, fileInd ) ;
		}

		free( readBatch ) ;	
	}
	else 
	{
		int maxBatchSize = READ_BUFFER_PER_THREAD * ( numOfThreads - 1 ) ;
		int batchSize[3] ;
		bool init = true, canJoinOutputThread = false ;
		int tag = 2, prevTag ;
		int fileInd[3] ;

		int useOutputThread = 0 ;
		pthread_t outputThread ;
		struct _OutputThreadArg outputArg ;

		struct _ErrorCorrectionThreadArg arg ;
		pthread_mutex_t errorCorrectionLock ;
		void *pthreadStatus ;
		
		struct _Read *readBatch[3] ;
		readBatch[0] = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		readBatch[1] = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		readBatch[2] = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;
		
		pthread_mutex_init( &errorCorrectionLock, NULL ) ;

		arg.kmerLength = kmerLength ;
		arg.trustedKmers = &trustedKmers ;
		//arg.readBatch = readBatch ;
		arg.lock = &errorCorrectionLock ;
		arg.badQuality = badQuality ;
		arg.batchSize = 0 ;
		arg.batchFinished = 0 ;
		
		if ( numOfThreads >= 6 )
			useOutputThread = 1 ;

		while ( 1 )
		{
			prevTag = tag ;
			tag = ( tag + 1 > 2 ) ? 0 : ( tag + 1 ) ;

			batchSize[tag] = reads.GetBatch( readBatch[tag], maxBatchSize, fileInd[tag], true, true ) ;

			// Wait for the previous batch finish
			if ( !init )
			{
				if ( canJoinOutputThread  ) // wait for the finish of the previous previous batch's output
					pthread_join( outputThread, &pthreadStatus ) ;
				for ( i = 0 ; i < numOfThreads - 1 - useOutputThread ; ++i )
					pthread_join( threads[i], &pthreadStatus ) ;
			}

			// Start current batch
			if ( batchSize[tag] != 0 )
			{

				//printf( "batchSize=%d\n", batchSize ) ;
				arg.batchSize = batchSize[tag] ;
				arg.readBatch = readBatch[tag] ;
				arg.batchUsed = 0 ;
				arg.batchFinished = 0 ;
				for ( i = 0 ; i < numOfThreads - 1 - useOutputThread ; ++i )
					pthread_create( &threads[i], &pthreadAttr, ErrorCorrection_Thread, (void *)&arg ) ;	

				//for ( i = 0 ; i < numOfThreads - 1 ; ++i )
				//	pthread_join( threads[i], &pthreadStatus ) ;
			}

			// Output previous batch 
			if ( !init )
			{
				// Create another thread to output previous batch
				if ( !useOutputThread )
				{
					for (  i = 0 ; i < batchSize[prevTag] ; ++i )
						UpdateSummary( readBatch[prevTag][i].seq, readBatch[prevTag][i].correction, readBatch[prevTag][i].badSuffix, paraDiscard, summary ) ;			
					reads.OutputBatch( readBatch[prevTag], batchSize[prevTag], ALLOW_TRIMMING, fileInd[prevTag] ) ;
				}
				else
				{	
					outputArg.readBatch = readBatch[ prevTag ] ;
					outputArg.batchSize = batchSize[prevTag] ;
					outputArg.summary = &summary ;
					outputArg.paraDiscard = paraDiscard ;
					outputArg.reads = &reads ;
					outputArg.fileInd = fileInd[ prevTag] ;

					pthread_create( &outputThread, &pthreadAttr, Output_Thread, (void *)&outputArg ) ;	
					canJoinOutputThread = true ;
				}
			}

			if ( batchSize[tag] == 0 )
				break ;

			init = false ;
		}

		if ( canJoinOutputThread )
			pthread_join( outputThread, &pthreadStatus ) ;
		//fprintf( stderr, "jump out\n" ) ;
		free( readBatch[2] ) ;
		free( readBatch[1] ) ;	
		free( readBatch[0] ) ;
	}

	PrintLog( "Finish error correction" ) ;
	PrintSummary( summary ) ;
	return 0 ;
}
