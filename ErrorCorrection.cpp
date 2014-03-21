#include "ErrorCorrection.h"
#include <string.h>

extern char nucToNum[26] ; 
extern char numToNuc[26] ;
extern int MAX_CORRECTION ;
extern bool ALLOW_TRIMMING ;

void *ErrorCorrection_Thread( void *arg )
{
	int ind ;
	int correction, badPrefix, badSuffix ;
	struct _ErrorCorrectionThreadArg *myArg = ( struct _ErrorCorrectionThreadArg *)arg ; 	
	
	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		ind = myArg->batchUsed ;
		++myArg->batchUsed ;
		pthread_mutex_unlock( myArg->lock ) ;
		//printf( "%d %d\n", ind, myArg->batchSize ) ;
		if ( ind >= myArg->batchSize )
			break ;
		correction = ErrorCorrection_Wrapper( myArg->readBatch[ind].seq, kmerCode, myArg->trustedKmers, 
				badPrefix, badSuffix ) ;
		myArg->readBatch[ind].correction = correction ;
		myArg->readBatch[ind].badPrefix = badPrefix ;
		myArg->readBatch[ind].badSuffix = badSuffix ;
	}

	pthread_exit( NULL ) ;
}


int ErrorCorrection( char *read, KmerCode& kmerCode, int maxCorrection, Store *kmers, int &badPrefix, int &badSuffix )
{
	int i, j, k ;	
	bool storedKmer[1000] ;
	int kmerCnt = 0 ;
	int readLength = 0 ;
	int fix[1000] ; // The fixed character for each untrusted position.
	int tag ;
	int from, to ;
	int trimStart = -1 ;
	
//	KmerCode kmerCode( inCode ) ;
	KmerCode tmpKmerCode( 0 ) ;
	badPrefix = badSuffix = 0 ;

	int kmerLength = kmerCode.GetKmerLength() ;
	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	if ( !kmers->IsIn( kmerCode ) )
		storedKmer[0] = false ;
	else
		storedKmer[0] = true ;
	kmerCnt = 1 ;
	for (  ; read[i] ; ++i, ++kmerCnt )
	{
		kmerCode.Append( read[i] ) ;

		if ( !kmers->IsIn( kmerCode ) )
			storedKmer[ kmerCnt ] = false ;
		else
			storedKmer[ kmerCnt ] = true ;
	}

	readLength = i ;
	trimStart = readLength ;
	for ( i = 0 ; i < readLength ; ++i )
		fix[i] = -1 ;
	//printf( "%d %d %d\n", kmerLength, kmerCnt, i ) ;	
	
	/*printf( "%s\n", read ) ;
	for ( i = 0 ; i < kmerCnt ; ++i )
		printf( "%d", storedKmer[i] ) ;
	printf( "\n" ) ;*/
	//exit( 1 ) ;
	// All the kmers are reliable. 
	for ( i = 0 ; i < kmerCnt ; ++i )
	{
		if ( storedKmer[i] == false )
			break ;
	}
	/*if ( !strcmp( read, "CCAGTACCTGAAATGCTTACGTTGCCGTTGCTGGCTCATCCTGCCCAGAG" ) )
	//if ( !strcmp( read, "ATCATCAGAGGGTCTCGTGTAGTGCTCCAGTACCTGAAATGCTTACGTTG" ) )
	{
		printf( "## %d\n", i ) ;
	}*/
	if ( i >= kmerCnt )
		return 0 ;

	// Find the xxfirst trusted kmer
	int longestStoredKmerCnt = 0, storedKmerCnt = 0 ;
	tag = -1 ;
	for ( i = 0 ; i < kmerCnt ; ++i )
	{
		if ( storedKmer[i] )
		{
			++storedKmerCnt ;
		}
		else
		{
			if ( storedKmerCnt > longestStoredKmerCnt )
			{
				longestStoredKmerCnt = storedKmerCnt ;
				tag = i - storedKmerCnt ; 
			}
			storedKmerCnt = 0 ;
		}
	}
	if ( storedKmerCnt > longestStoredKmerCnt )
	{
		longestStoredKmerCnt = storedKmerCnt ;
		tag = i - storedKmerCnt ; 
	}

	if ( longestStoredKmerCnt >= kmerCnt )
	{
		//printf( "0%s\n", read ) ;
		return 0 ;
	}

	//if ( tag > 0 )
	//	return 0 ;
	for ( k = tag + 1 ; k < kmerCnt ; ++k )
		if ( !storedKmer[k] )
			break ;

	// Scan towards right
	kmerCode.Restart() ;
	if ( k >= kmerCnt )
		// The kmers are all correct after tag, so force skip the correction.
		i = readLength + 1 ;
	else
	{
		for ( i = k - 1 + 1 ; i < k - 1 + kmerLength - 1 + 1 ; ++i  )
		{
			kmerCode.Append( read[i] ) ;
		}
	}
	/*for ( i = k ; i < k + kmerLength - 1 ; ++i  )
	{
		kmerCode = ( kmerCode << (uint64_t)2 ) & mask ;
		kmerCode = kmerCode | (uint64_t)nucToNum[ read[i] - 'A' ] ;
	}*/

	for ( ; i < readLength ; )
	{
		KmerCode fixedKmerCode( kmerCode ) ;
		int maxTo = -1 ;
		int maxChange = -1 ;
		int maxCnt = 0 ;
		from = i + 1 ;
		to = ( i + kmerLength - 1 < readLength ) ? i + kmerLength - 1 : readLength - 1 ; 
		for ( j = 0 ; j < 4 ; ++j )
		{
			//kmerCode = ( fixedKmerCode << (uint64_t)2 ) & mask ;
			//kmerCode = kmerCode | (uint64_t)j ;
			kmerCode = fixedKmerCode ;
			kmerCode.Append( numToNuc[j] ) ;
			//printf( "?code=%llu\n", kmerCode.GetCode() ) ;
			if ( !kmers->IsIn( kmerCode ) )
				continue ;

			if ( maxTo == -1 )
				maxTo = i ;
			// How many kmers this change can fix
			for ( k = from ; k <= to ; ++k )
			{
				kmerCode.Append( read[k] ) ;
				if ( !kmers->IsIn( kmerCode ) )
					break ;
			}
			
			// Try to extend 1 position
			if ( k > to && to == readLength - 1 )
			{
				int l, m ;
				for ( m = 0 ; m < kmerLength - 1 - ( to - from + 1 ) ; ++m )
				{
					for ( l = 0 ; l < 4 ; ++l )
					{
						KmerCode tmpKmerCode( kmerCode ) ;
						tmpKmerCode.Append( numToNuc[l] ) ;
						if ( kmers->IsIn( tmpKmerCode ) )
							break ;
					}

					if ( l < 4 )
					{
						//printf( "hi\n" ) ;
						kmerCode.Append( numToNuc[l] ) ;
						++k ;
					}
				}
			}

			//printf( "hi %d %d\n", k, maxTo ) ;

			if ( k > maxTo )
			{
				/*if ( maxTo - from + 1 >= kmerLength / 2 )
				{
					++maxCnt ;
				}
				else*/
				maxCnt = 1 ;
				maxTo = k ;
				maxChange = j ;
				tmpKmerCode = kmerCode ;
			}
			/*else if ( k - from + 1 >= kmerLength / 2 )
			{
				++maxCnt ;
			}*/
			else if ( k == maxTo )
			{
				++maxCnt ;
				if ( k == i + 1 && j == nucToNum[ read[i] - 'A' ] )
				{
					maxCnt = 1 ;
					maxChange = j ;
					tmpKmerCode = kmerCode ;
				}
				else if ( k == i + 1 && maxChange == nucToNum[ read[i] - 'A' ] )
				{
					maxCnt = 1 ;
				}
			}
		}
		if ( maxTo == -1 || maxCnt > 1  )
		{
			//printf( "+%s\n", read ) ;

			// trim
			//if ( maxCnt > 1 )
			//{
			trimStart = i ;
			break ;
			//}
			//else
			//return 0 ;
		}
		//printf( "hi %d: %d %d=>%d, (%d): code=%llu\n", i, maxTo, to, maxChange, maxCnt, tmpKmerCode.GetCode() ) ;	
		fix[i] = maxChange ;

		if ( maxTo >= readLength )
			break ;
		
		if ( maxTo <= to )
		{
			// There are multiple errors in the region

			if ( 0 ) //maxTo > i + 1 )
			{
				// The next start search can start one place earlier
				//uint64_t tmp ;
				/*kmerCode.ShiftRi  ;
				if ( fix[maxTo - kmerLength + 1] != -1 )
					kmerCode |= ( (uint64_t)fix[maxTo - kmerLength + 1] << ( 2ull * ( kmerLength - 2 ) )  ) ;
				else
					kmerCode |= ( (uint64_t)nucToNum[ read[maxTo - kmerLength + 1] - 'A' ] << ( 2ull * ( kmerLength -2 ) )  ) ;
				i = maxTo - 1 ;*/
			}
			else
			{
				kmerCode = tmpKmerCode ;
				kmerCode.ShiftRight( 1 ) ;
				i = maxTo ;
			}
			continue ;
		}
		else
		{
			// Search for next error.
			for ( k = to - kmerLength + 2 ; k < kmerCnt ; ++k )
				if ( !storedKmer[k] )
					break ;

			if ( k >= kmerCnt )
				break ;

			kmerCode.Restart() ;
			for ( i = k - 1 + 1 ; i < k - 1 + kmerLength - 1 + 1 ; ++i  )
			{
				if ( fix[i] == -1 )
					kmerCode.Append( read[i] ) ;
				else
					kmerCode.Append( numToNuc[ fix[i] ] ) ;
			}
			continue ;
		}
	}

	// Scan towards left
	kmerCode.Restart() ;
	for ( i = tag + 1 - 1 ; i < tag + kmerLength - 1 ; ++i )
	{
		kmerCode.Append( read[i] ) ;	
	}
	kmerCode.Append( 'A' ) ;
	//printf( "%d\n", tag ) ;
	if ( tag == 0 )
	{
		// Force skip
		tag = -1 ;
	}
	//kmerCode = ( kmerCode << (uint64_t)2 ) & mask ;
	for ( i = tag - 1 ; i >= 0 ;  )
	{
		KmerCode fixedKmerCode( kmerCode ) ;
		int minTo = readLength + 1 ;
		int minChange = -1 ;
		int minCnt = 0 ;
		from = i - 1 ;
		to = ( i - kmerLength + 1 < 0 ) ? 0 : ( i - kmerLength + 1 ) ; 
		for ( j = 0 ; j < 4 ; ++j )
		{
			//kmerCode = ( fixedKmerCode >> (uint64_t)2 ) & ( mask >> 2ull ) ;
			//printf( "1.%llu\n", kmerCode.GetCode() ) ;
			//kmerCode = kmerCode | ( ((uint64_t)j) << (uint64_t)(2 * kmerLength - 2 ) ) ;

			kmerCode = fixedKmerCode ;
			kmerCode.Prepend( numToNuc[j] ) ;
			
			//printf( "2.%llu\n", kmerCode.GetCode() ) ;
			//printf( "%d %lld %lld\n", j, (long long int)kmerCode, (long long int)fixedKmerCode ) ;
			if ( !kmers->IsIn( kmerCode ) )
				continue ;
			if ( minTo == tag + 1 )
				minTo = i ;
			// How many kmers this change can fix
			for ( k = from ; k >= to ; --k )
			{
				kmerCode.Prepend( read[k] ) ;	
				if ( !kmers->IsIn( kmerCode ) )
					break ;
			}

			// try extension
			if ( k < to && to == 0 )
			{
				int l, m ;
				for ( m = 0 ; m < kmerLength - 1 - ( from - to + 1 ) ; ++m )
				{
					for ( l = 0 ; l < 4 ; ++l )
					{
						KmerCode tmpKmerCode( kmerCode ) ;
						tmpKmerCode.Prepend( numToNuc[l] ) ;
						if ( kmers->IsIn( tmpKmerCode ) )
							break ;
					}

					if ( l < 4 )
					{
						//printf( "hi\n" ) ;
						kmerCode.Prepend( numToNuc[l] ) ;
						--k ;
					}
				}
			}

			if ( k < minTo )
			{
				minCnt = 1 ;
				minTo = k ;
				minChange = j ;
				tmpKmerCode = kmerCode ;
			}
			else if ( k == minTo )
			{
				++minCnt ;
				if ( k == i - 1 && j == nucToNum[ read[i] - 'A' ] )
				{
					minCnt = 1 ;
					minChange = j ;
					tmpKmerCode = kmerCode ;
				}
				else if ( k == i -1 && minChange == nucToNum[ read[i] - 'A' ] )
				{
					minCnt = 1 ;
				}
			}
		}
		//printf( "-hi %d: %d %d=>%d, (%d)\n", i, minTo, to, minChange, minCnt ) ;	
		if ( minTo == readLength + 1 || minCnt > 1  )
		{
			//printf( "-%s\n", read ) ;
			//return -1 ;
			badPrefix = i + 1 ;
			break ;
		}
		//printf( "---%d %d\n", minChange,  nucToNum[read[0] - 'A'] ) ;
		fix[i] = minChange ;

		if ( minTo < 0 )
			break ;
		if ( minTo >= to )
		{
			
			/*if ( maxTo > i + 1 )
			{
				// The next start search can start one place earlier
				kmerCode = tmpKmerCode >> 2ull ;
				i = maxTo - 1 ;
			}
			else*/
			//kmerCode = tmpKmerCode << 2ull ;
			kmerCode = tmpKmerCode ;
			kmerCode.Append( 'A' ) ;
			//printf( "i=>minTo: %d=>%d\n", i, minTo ) ;
			i = minTo ;
			continue ;
		}
		else
		{
			// Search for next error.
			for ( k = to - 1 ; k >= 0 ; --k )
				if ( !storedKmer[k] )
					break ;

			if ( k < 0 )
				break ;

			kmerCode.Restart() ;
			for ( i = k + 2  - 1 ; i < k + 2 + kmerLength - 1 - 1 ; ++i  )
			{
				if ( fix[i] == -1 )
					kmerCode.Append( read[i] ) ;
				else
					kmerCode.Append( numToNuc[ fix[i] ] ) ;
			}
			i = k + 1 - 1 ;
			kmerCode.Append( 'A' ) ;
			continue ;
		}
		/*{
			// This happens when the distance between two errors is just of kmerLength
			//printf( "%llu\n", tmpKmerCode.GetCode() ) ;
			kmerCode = tmpKmerCode ;
			i = minTo ;
		}*/
	}
	
	/*printf( "fix: ") ;
	for ( i = 0 ; i < readLength; ++i )
		printf( "%d", fix[i] ) ;
	printf( "\n" ) ;*/

	int ret = 0 ;
	int correctCnt = 0 ;
	//printf( "%s\n%d\n", read, ret ) ;
	for ( i = badPrefix ; i < trimStart ; ++i )
	{
		if ( fix[i] == -1 || read[i] == numToNuc[ fix[i] ] || read[i] == 'N' )
			continue ;
		++correctCnt ;

		if ( correctCnt > maxCorrection )
		{
			// There are too many corrections, adjust the trimStart
			for ( j = i - 1 ; j >= badPrefix ; --j )
			{
				if ( fix[j] == -1 || read[j] == numToNuc[ fix[j] ] )
					break ;
			}
			trimStart = j + 1 ;
			if ( trimStart == badPrefix )
				return -1 ;
			
			correctCnt = 0 ;
			for ( i = badPrefix ; i < trimStart ; ++i )
			{
				if ( fix[i] == -1 || read[i] == numToNuc[ fix[i] ] )
					continue ;
				++correctCnt ;
			}
			if ( correctCnt > maxCorrection - 1 )
				return -1 ;

			break ;
		}
	}

	/*if ( correctCnt > MAX_CORRECTION )
	{
		// Find the point where the correct count become to increase rapidly.
		//return 0 ;
		int correctPartialSum[1000] ; 
		correctPartialSum[0] = 0 ;
		for ( i = 0 ; i < trimStart ; ++i )
		{
			correctPartialSum[i + 1] = correctPartialSum[i - 1] + 
				( ( fix[i] == -1 || read[i] == numToNuc[ fix[i] ] ) ? 0 : 1 ) ;
		}

		for ( i = trimStart - 1 ; i >= 0 ; --i )
		{
		}
		
		//int cnt = 0 ;
		//for ( i = 0 ; i < trimStart ; ++i )
		
	}*/

	for ( i = badPrefix ; i < trimStart ; ++i )
	{
		if ( fix[i] == -1 )
			continue ;
		if ( read[i] != numToNuc[ fix[i] ] )
		{
			read[i] = numToNuc[ fix[i] ] ;
			//if ( read[i] != 'N' )
			++ret ;
		}
	}
	
	//printf( "%d %d\n", ret, trimStart ) ;

	//badPrefix = 0 ;
	badSuffix = readLength - trimStart ;
		
	//if ( ALLOW_TRIMMING ) // else we do partial correction
	//	read[ trimStart ] = '\0' ; 
	//ret += trimmed ;
	return ret ;
}

int ErrorCorrection_Wrapper( char *read, KmerCode& kmerCode, Store *kmers, int &badPrefix, int &badSuffix )
{
	int correction ;
	int tmpBadPrefix, tmpBadSuffix ;
	int len = (int)strlen( read ) ;
	int interim = 0 ;

	correction = ErrorCorrection( read, kmerCode, MAX_CORRECTION, kmers,
			tmpBadPrefix, tmpBadSuffix ) ;

	if ( correction == -1 )
	{
		badPrefix = badSuffix = 0 ;
		return -1 ;
	}

	badPrefix = tmpBadPrefix ;
	badSuffix = tmpBadSuffix ;
	if ( badPrefix > len / 2 || badPrefix > kmerCode.GetKmerLength() * 2 )
	{
		int tmp ;
		char c = read[badPrefix] ;
		read[ badPrefix ] = '\0' ;
		tmp = ErrorCorrection( read, kmerCode, MAX_CORRECTION - correction, kmers, 
				tmpBadPrefix, tmpBadSuffix ) ;

		read[ badPrefix ] = c ;
		if ( tmp != -1 )
		{
			badPrefix = tmpBadPrefix ;
			correction += tmp ;
			interim += tmpBadSuffix ;
		}
	}

	if ( badSuffix > len / 2 || badSuffix > kmerCode.GetKmerLength() * 2 )
	{
		int tmp ;
		tmp = ErrorCorrection( read + len - badSuffix, kmerCode, MAX_CORRECTION - correction, kmers,
				tmpBadPrefix, tmpBadSuffix ) ;	
		if ( tmp != -1 )
		{
			badSuffix = tmpBadSuffix ;
			correction += tmp ;
			interim += tmpBadPrefix ;
		}
	}

	if ( correction == 0 && interim > 0 )
	{
		correction = -1 ;
		badSuffix = badPrefix = 0 ;
	}
	else if ( ALLOW_TRIMMING )
	{
		read[len - badSuffix] = '\0' ;
	}

	return correction ;
}
