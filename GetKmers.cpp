#include "GetKmers.h"

void *SampleKmers_Thread( void *arg )
{
	struct _SampleKmersThreadArg *myArg = ( struct _SampleKmersThreadArg *)arg ; 	
	//char read[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], id[MAX_ID_LENGTH] ;
	int i, tmp ;
	struct _Read readBatch[128] ;
	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		//tmp = myArg->reads->NextWithBuffer( id, read, qual, false ) ;
		tmp = myArg->reads->GetBatch( readBatch, 128, false, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;
		
		if ( tmp != 0 )
		{
			for ( i = 0 ; i < tmp ; ++i )
			{
				int len ;
				char *read = readBatch[i].seq ;
				char *qual = readBatch[i].qual ;
				/*len = strlen( id ) ;		
				  if ( id[len - 1] == '\n')
				  id[len - 1] = '\0' ;*/
				
				len = (int)strlen( read ) ;
				if ( read[len - 1] == '\n' )
					read[len - 1] = '\0' ;

				if ( qual[0] != '\0' )
				{
					if ( qual[len - 1] == '\n' )
						qual[len - 1] = '\0' ;
				}
				//printf( "%s\n", read ) ;
				SampleKmersInRead( read, qual, myArg->kmerLength, myArg->alpha, 
						kmerCode, myArg->kmers ) ;
			}
		}
		else
			break ;
	}
	pthread_exit( NULL ) ;	
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
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}

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
	char read[MAX_READ_LENGTH], qual[MAX_READ_LENGTH], id[MAX_ID_LENGTH] ;
	int tmp ;
	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		tmp = myArg->reads->NextWithBuffer( id, read, qual, false ) ;
		pthread_mutex_unlock( myArg->lock ) ;
		if ( tmp != 0 )
		{
			int len = (int)strlen( id ) ;
			if ( id[len - 1] == '\n')
				id[len - 1] = '\0' ;
			len = (int)strlen( read ) ;
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
		else
			break ;
	}
	pthread_exit( NULL ) ;	
}

void StoreTrustedKmers( char *read, char *qual, int kmerLength, char badQuality, int *threshold,
	KmerCode &kmerCode, Store *kmers, Store *trustedKmers )
{
	bool occur[MAX_READ_LENGTH] ;
	bool trustedPosition[MAX_READ_LENGTH] ;
	int i, k ;

	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	k = 0 ;
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
	  printf( "\n" ) ;*/

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
