#include "GetKmers.h"

struct _SampleKmersPutThreadArg
{
	Store *kmers ;
	KmerCode *kmerCodes ;
	int kmerCodesCnt ;

	//pthread_mutex_t *lockPut ;
} ;

void *SampleKmers_PutThread( void *arg )
{
	struct _SampleKmersPutThreadArg *myArg = ( struct _SampleKmersPutThreadArg * )arg ;
	Store *kmers = myArg->kmers ;
	int i ;
	//pthread_mutex_lock( myArg->lockPut ) ;
	for ( i = 0 ; i < myArg->kmerCodesCnt ; ++i )
		kmers->Put( myArg->kmerCodes[i], true ) ;
	//pthread_mutex_unlock( myArg->lockPut ) ;

	pthread_exit( NULL ) ;
	return NULL ;
}

void *SampleKmers_Thread( void *arg )
{
	struct _SampleKmersThreadArg *myArg = ( struct _SampleKmersThreadArg *)arg ; 	
	//char read[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], id[MAX_ID_LENGTH] ;
	int i, tmp ;
	int fileInd ;
	struct _Read *readBatch = ( struct _Read *)malloc( sizeof( *readBatch ) * 128 ) ;
	struct _SamplePattern *samplePatterns = myArg->samplePatterns ;
	int kmerLength = myArg->kmerLength ;
	KmerCode kmerCode( kmerLength ) ;

	//double p ;
	Store *kmers = myArg->kmers ;
	//double factor = 1.0 ;
	const int bufferSizeFactor = 9 ;
	KmerCode *kmerCodeBuffer[2] ; //[ ( bufferSizeFactor + 1 )* MAX_READ_LENGTH] ;
	int bufferTag ;
	int kmerCodeBufferUsed = 0 ;
	
	for ( bufferTag = 0 ; bufferTag < 2 ; ++bufferTag )
	{
		kmerCodeBuffer[ bufferTag ] = ( KmerCode * )malloc( sizeof( KmerCode ) * ( bufferSizeFactor + 1 ) * MAX_READ_LENGTH ) ;
		for ( i = 0 ; i < ( bufferSizeFactor + 1 ) * MAX_READ_LENGTH ; ++i )
			kmerCodeBuffer[bufferTag][i] = kmerCode ;
	}

	void *pthreadStatus ;
	pthread_t putThread ;
	pthread_attr_t pthreadAttr ;
	struct _SampleKmersPutThreadArg putArg ;
	bool threadInit = false ;
	pthread_attr_init( &pthreadAttr ) ;
	pthread_attr_setdetachstate( &pthreadAttr, PTHREAD_CREATE_JOINABLE ) ;
	bufferTag = 0 ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		//tmp = myArg->reads->NextWithBuffer( id, read, qual, false ) ;
		tmp = myArg->reads->GetBatch( readBatch, 128, fileInd, false, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;
		int k ;
		if ( tmp != 0 )
		{
			for ( k = 0 ; k < tmp ; ++k )
			{
				int len ; 
				char *id = readBatch[k].id ;
				char *read = readBatch[k].seq ;
				char *qual = readBatch[k].qual ;
				/*len = strlen( id ) ;		
				  if ( id[len - 1] == '\n')
				  id[len - 1] = '\0' ;*/

				//int tag = rand() % SAMPLE_PATTERN_COUNT ; // Decide to use which pattern
				int tag = fileInd ;

				len = (int)strlen( read ) ;
				if ( read[len - 1] == '\n' )
					read[len - 1] = '\0' ;

				if ( qual[0] != '\0' )
				{
					if ( qual[len - 1] == '\n' )
						qual[len - 1] = '\0' ;
				}

				for ( i = 0 ; id[i] ; ++i )	
					tag = tag * 17  + ( id[i] - 'A' ) ;

				for ( i = 0 ; read[i] ; ++i )
					tag = tag * 7 + ( read[i] - 'A' )  ;
				
				for ( i = 0 ; qual[i] ; ++i )
					tag = tag * 17+ ( qual[i] - 'A' )  ;
					
				tag %= SAMPLE_PATTERN_COUNT ;
				if ( tag < 0 )
					tag += SAMPLE_PATTERN_COUNT ;
				//tag = rand() % SAMPLE_PATTERN_COUNT ; // Decide to use which pattern
				//printf( "%d\n", tag ) ;
				//printf( "%s\n", read ) ;
				//SampleKmersInRead( read, qual, myArg->kmerLength, myArg->alpha, 
				//		kmerCode, myArg->kmers ) ;

				if ( len - 1 < kmerLength )
					continue ;

				for ( i = 0 ; i < kmerLength ; ++i )
				{
					kmerCode.Append( read[i] ) ;
				}
				//printf( "%d %d %d\n", tag, (int)samplePatterns[tag].tag[2], (int)samplePatterns[tag].tag[3] ) ;	
				if ( samplePatterns[tag].tag[ ( i - 1 ) / 8 ] & ( 1 << ( ( i - 1 ) % 8 ) ) )
				{
					//kmers->Put( kmerCode ) ;
					kmerCodeBuffer[ bufferTag ][ kmerCodeBufferUsed ] = kmerCode ;
					++kmerCodeBufferUsed ;
				}

				for ( ; read[i] ; ++i )
				{
					kmerCode.Append( read[i] ) ;

					if ( samplePatterns[tag].tag[ i / 8 ] & ( 1 << ( i % 8 ) ) )
					{
						//kmers->Put( kmerCode ) ;
						kmerCodeBuffer[bufferTag][ kmerCodeBufferUsed ] = kmerCode ;
						++kmerCodeBufferUsed ;
					}
				}
				
				if ( kmerCodeBufferUsed >= bufferSizeFactor * MAX_READ_LENGTH )
				{
					if ( threadInit )
						pthread_join( putThread, &pthreadStatus ) ;
					
					putArg.kmers = kmers ;
					putArg.kmerCodes = ( KmerCode *)kmerCodeBuffer[ bufferTag ] ;
					putArg.kmerCodesCnt = kmerCodeBufferUsed ;
					//putArg.lockPut = myArg->lockPut ;
					//printf( "%d %d\n", bufferTag, kmerCodeBufferUsed ) ;
					pthread_create( &putThread, &pthreadAttr, SampleKmers_PutThread, ( void *)&putArg ) ;

					kmerCodeBufferUsed = 0 ;
					bufferTag = 1 - bufferTag ;

					threadInit = true ;
				}
			}
		}
		else
			break ;
	}
	
	//put the remaining
	if ( threadInit )
		pthread_join( putThread, &pthreadStatus ) ;

	putArg.kmers = kmers ;
	putArg.kmerCodes = ( KmerCode *)kmerCodeBuffer[ bufferTag ] ;
	putArg.kmerCodesCnt = kmerCodeBufferUsed ;
	//putArg.lockPut = myArg->lockPut ;

	pthread_create( &putThread, &pthreadAttr, SampleKmers_PutThread, ( void *)&putArg ) ;
	kmerCodeBufferUsed = 0 ;
	bufferTag = 1 - bufferTag ;
	threadInit = true ;

	pthread_join( putThread, &pthreadStatus ) ;

	free( readBatch ) ;
	for ( i = 0 ; i < 2 ; ++i )
		free( kmerCodeBuffer[i] ) ;
	pthread_attr_destroy( &pthreadAttr ) ;
	pthread_exit( NULL ) ;
	return NULL ;
}

void SampleKmersInRead( char *read, char *qual, int kmerLength, double alpha, KmerCode &kmerCode, Store *kmers )
{
	int i ;
	double p ;
	double factor = 1 ;
	/*int badQualPartialCount[MAX_READ_LENGTH]*/ ;

	// Get the partial counting sum of not good qualies
	// NOTICE: we shift one position here
	/*badQualPartialCount[0] = 0 ; 
	for ( i = 0 ; qual[i] ; ++i )
	{
		if ( qual[0] != '\0' )
			badQualPartialCount[i + 1] = badQualPartialCount[i] + 
				( qual[i] < goodQuality ? 1 : 0 ) ;
		else
			badQualPartialCount[i + 1] = i + 1 ;
	}*/
	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength && read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}

	if ( i < kmerLength )
		return ;

	p = rand() / (double)RAND_MAX ;
	/*if ( badQualPartialCount[kmerLength] - badQualPartialCount[0] == 0 )
		factor = 1 ;
	else
		factor = 1 ;*/
	//factor = 1 ;
	//printf( "%lf %lf\n", p, alpha ) ;
	if ( p < alpha * factor )
	{
		//printf( "%lf %lf\n", p, alpha ) ;
		kmers->Put( kmerCode ) ;
	}

	for ( ; read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;

		/*if ( badQualPartialCount[i + 1] - badQualPartialCount[i - kmerLength + 1] == 0 )
			factor = 1 ;
		else
			factor = 1 ;*/
		//factor = 1 ;
		p = rand() / (double)RAND_MAX ;
		if ( p < alpha * factor )
		{
			/*if ( !strcmp( read, "CCAGTACCTGAAATGCTTACGTTGCCGTTGCTGGCTCATCCTGCCCAGAG" ) )
			  {
			  ExtractKmer( read, i - kmerLength + 1, kmerLength, buffer ) ;
			  printf( "Put %s\n", buffer ) ;
			  }*/
			//printf( "%lf %lf\n", p, alpha ) ;
			kmers->Put( kmerCode ) ;
		}
	}
}
void *StoreKmers_Thread( void *arg )
{
	struct _StoreKmersThreadArg *myArg = ( struct _StoreKmersThreadArg *)arg ; 	
	int i ;
	int tmp, fileInd ;
	int maxBatchSize = READ_BUFFER_PER_THREAD ; 
	struct _Read *readBatch = ( struct _Read *)malloc( sizeof( struct _Read ) * maxBatchSize ) ;

	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		//tmp = myArg->reads->( id, read, qual, false ) ;
		tmp = myArg->reads->GetBatch( readBatch, maxBatchSize, fileInd, false, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;

		if ( tmp == 0 )
			break ;
		
		for ( i = 0 ; i < tmp ; ++i )
		{
			/*id = readBatch[i].id ;
			int len = (int)strlen( id ) ;
			if ( id[len - 1] == '\n')
				id[len - 1] = '\0' ;*/
			char *read = readBatch[i].seq ;
			char *qual = readBatch[i].qual ;

			int len = (int)strlen( read ) ;
			if ( read[len - 1] == '\n' )
				read[len - 1] = '\0' ;

			if ( qual[0] != '\0' )
			{
				if ( qual[len - 1] == '\n' )
					qual[len - 1] = '\0' ;
			}

			StoreTrustedKmers( read, qual, myArg->kmerLength, myArg->badQuality, myArg->threshold, 
				kmerCode, myArg->kmers, myArg->trustedKmers ) ;
		}
	}

	free( readBatch ) ;

	pthread_exit( NULL ) ;
	return NULL ;
}

void StoreTrustedKmers( char *read, char *qual, int kmerLength, char badQuality, int *threshold,
	KmerCode &kmerCode, Store *kmers, Store *trustedKmers )
{
	bool occur[MAX_READ_LENGTH] ;
	bool trustedPosition[MAX_READ_LENGTH] ;
	int i ;

	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	if ( kmers->IsIn( kmerCode) )
		occur[i - kmerLength] = true ;
	else
		occur[i - kmerLength] = false ;

	for ( ; read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;

		if ( kmers->IsIn( kmerCode ) )
			occur[i - kmerLength + 1] = true ;
		else
		{
			occur[i - kmerLength + 1] = false ;
		}
	}

	int readLength = i ;
	int occurCnt = readLength - kmerLength + 1 ;
	int zeroCnt = 0, oneCnt = 0 ;

	/*printf( "%s\n", read ) ;
	for ( i = 0 ; i < occurCnt ; ++i )
	  printf( "%d", occur[i] ) ;
	printf( "\n" ) ;*/

	// Set the trusted positions
	for ( i = 0 ; i < readLength ; ++i )
	{
		if ( i >= kmerLength )
		{
			if ( occur[i - kmerLength] )
				--oneCnt ;
			else
				--zeroCnt ;
		}

		if ( i < occurCnt )
		{
			if ( occur[i] )
				++oneCnt ;
			else
				++zeroCnt ;
		}
		//printf( "%d %d\n",i, oneCnt ) ;
		int sum = oneCnt + zeroCnt ;	
		int adjust = 0 ;
		/*if ( i - kmerLength + 1 >= kmerLength - 1 && i + kmerLength - 1 < readLength &&
			qual[i] >= goodQuality )
		{
			adjust -= 1 ;
		}*/
		/*if ( i + kmerLength / 2 - 1 >= readLength || i < kmerLength / 2 )
		{
			adjust += 1 ;
		}*/
	
		//TODO: play with the parameters here
		//if ( oneCnt > alpha * sum && 
		//	( oneCnt - alpha * sum ) * ( oneCnt - alpha * sum ) > 
		//		8 * alpha * (1 - alpha ) * sum + 1 )  

		if ( qual[0] != '\0' && qual[i] <= badQuality )
		{
			trustedPosition[i] = false ;
			continue ;
		}

		if ( oneCnt > threshold[sum] + adjust )
		{
			//if ( !strcmp( read, "CCAGTACCTGAAATGCTTACGTTGCCGTTGCTGGCTCATCCTGCCCAGAG" ) )
			//	printf( "trustedPosition: %d %d\n", i, oneCnt ) ;
			trustedPosition[i] = true ;
		}
		else
		{
			trustedPosition[i] = false ;
		}
		//trustedPosition[i] = true ;
	}

	oneCnt = 0 ;
	kmerCode.Restart() ;

	//printf( "!! %s\n", readId ) ;
	//printf( "!! %s\n!! ", read ) ;
	/*for ( i = 0 ; i < readLength ; ++i )
	  printf( "%d", trustedPosition[i] ) ;
	 printf( "\n" ) ; */

	for ( i = 0 ; i < readLength ; ++i )
	{
		if ( trustedPosition[i] )
			++oneCnt ;
		if ( i >= kmerLength )
		{
			if ( trustedPosition[i - kmerLength] ) 
				--oneCnt ;
		}

		kmerCode.Append( read[i] ) ;
		if ( oneCnt == kmerLength 
		   )//|| ( i - kmerLength + 1 >= kmerLength - 1 && i + kmerLength - 1 < readLength && oneCnt >= kmerLength - 1 ) )
		   {
			   //ExtractKmer( read, i - kmerLength + 1, kmerLength, buffer ) ;
			   //if ( !strcmp( read,"CCAGTACCTGAAATGCTTACGTTGCCGTTGCTGGCTCATCCTGCCCAGAG" ) ) 
			   //if ( !strcmp( read, "ATCATCAGAGGGTCTCGTGTAGTGCTCCAGTACCTGAAATGCTTACGTTG" ) )
			   //	printf( "Into table B: %d %d\n", i, oneCnt ) ;
			   //printf( "%d %lld\n", i, kmerCode ) ;
			   trustedKmers->Put( kmerCode, true ) ;
		   }
	}
}
