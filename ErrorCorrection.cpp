#include "ErrorCorrection.h"
#include <string.h>

//#define DEBUG 

extern char nucToNum[26] ; 
extern char numToNuc[26] ;
extern int MAX_CORRECTION ;
extern bool ALLOW_TRIMMING ;
extern int SET_NEW_QUAL ;

void *ErrorCorrection_Thread( void *arg )
{
	int ind ;
	int correction, badPrefix, badSuffix, info ;
	struct _ErrorCorrectionThreadArg *myArg = ( struct _ErrorCorrectionThreadArg *)arg ; 	
	bool init = true ;
	KmerCode kmerCode( myArg->kmerLength ) ;
	while ( 1 )
	{
		pthread_mutex_lock( myArg->lock ) ;
		if ( !init )
			++myArg->batchFinished ;
		ind = myArg->batchUsed ;
		++myArg->batchUsed ;
		pthread_mutex_unlock( myArg->lock ) ;
		//printf( "%d %d\n", ind, myArg->batchSize ) ;
		if ( ind >= myArg->batchSize )
			break ;
		correction = ErrorCorrection_Wrapper( myArg->readBatch[ind].seq, myArg->readBatch[ind].qual, kmerCode, 
					myArg->badQuality, myArg->trustedKmers, badPrefix, badSuffix, info ) ;
		myArg->readBatch[ind].correction = correction ;
		myArg->readBatch[ind].badPrefix = badPrefix ;
		myArg->readBatch[ind].badSuffix = badSuffix ;
		myArg->readBatch[ind].info = info ;
		init = false ;
	}

	pthread_exit( NULL ) ;
	return NULL ;
}

// If we could not find an anchor, it is mostly likely there are errors in every kmer.
// We will try to find 1 position that can create a stretch of stored kmers.
// return-which position is changed. -1: failed
int CreateAnchor( char *read, char *qual, int *fix, bool *storedKmer, KmerCode &kmerCode, Store *kmers )
{
	int readLength = strlen( read ) ;	
	int kmerLength = kmerCode.GetKmerLength() ; 
	int i, j, k ;
	//if ( readLength < 2 * kmerLength )
  	//	return -1 ;
	int maxLen = 0 ;
	int maxLenStats[2] = {0, 0}; // 0-position, 1-which nucleutide changed to

	// This stored the partial result of substitution, we use two arrays
	// to keep the result for the optimal plan.
	int tag = 0 ;
	int stored[2][MAX_READ_LENGTH] ; 
	int scnt = 0 ;

	int start = readLength / 2 - kmerLength ;
	int end = readLength / 2 + kmerLength + 1 ;
	/*if ( start < 0 )
		start = 0 ;
	if ( end > readLength )
		end = readLength ;*/
	start = 0 ;
	end = readLength ;
	for ( i = start ; i < end ; ++i )
	{
		int from = i - kmerLength + 1 ;
		if ( from < 0 )
			from = 0 ;
		for ( j = 0 ; j < 4 ; ++j )
		{
			if ( numToNuc[j] == read[i] )
				continue ;
			char c = read[i] ;
			read[i] = numToNuc[j] ;
			
			// For efficiency, use one kmer to test whether we want to choose this candidate
			// But for bigger kmer, we need more effort

			// kind of test 1/32 of the kmers. Since the sequencing error are less likely on left-side,
			// we shift more on left
			int loop = ( kmerLength - 1 ) / 32 + 1 ;
			int t ;
			for ( t = 0 ; t < loop ; ++t )
			{
				int l = i - ( 2 * loop - t - 1 ) * kmerLength / ( 2 * loop ) + 1 ;
				if ( l < 0 )
					l = 0 ;
				if ( l + kmerLength - 1 >= readLength )
					l = readLength - 1 - kmerLength + 1 ;
				kmerCode.Restart() ;
				for ( k = l ; k <= l + kmerLength - 1 ; ++k )
					kmerCode.Append( read[k] ) ;
				if ( kmers->IsIn( kmerCode ) )
				{
					break ;
				}
			}
			if ( t >= loop )
			{
				read[i] = c ;
				continue ;
			}
			
			scnt = 0 ;
			for ( k = from ; k < from + kmerLength - 1 ; ++k )
				kmerCode.Append( read[k] ) ;
			scnt = 0 ;
			int missCnt = 0 ;
			int hitCnt = 0 ;
			for ( ; k < readLength && k < i + kmerLength ; ++k, ++scnt )
			{
				kmerCode.Append( read[k] ) ;
				if ( kmers->IsIn( kmerCode ) )
				{
					stored[tag][scnt] = 1 ;
					++hitCnt ;
				}
				else
				{
					stored[tag][scnt] = 0 ;
					++missCnt ;
					
					if ( kmerLength - missCnt < maxLen )
						break ;
				}
			}

			int sum = 0 ;
			int max = 0 ;
			for ( k = 0 ; k < scnt ; ++k )
			{
				if ( stored[tag][k] == 0 )
				{
					if ( sum > max )
						max = sum ;
					sum = 0 ;
				}
				else
					++sum ;
			}
			if ( sum > max )
				max = sum ;

			if ( max > maxLen )
			{
				maxLen = max ;
				maxLenStats[0] = i ;
				maxLenStats[1] = j ;
				/*for ( int l = 0 ; l < scnt ; ++l )
					printf( "%d ", stored[tag][l] ) ;
				printf( "\n" ) ;*/
				tag = 1 - tag ; 
				
				//read[i] = c ;
				//break ;
			}
			else if ( max > 0 && max == maxLen && qual[0] != '\0' && qual[i] < qual[ maxLenStats[0] ] )
			{
				maxLenStats[0] = i ;
				maxLenStats[1] = j ;
				tag = 1 - tag ;
			}
			
			read[i] = c ;
		}
	}

	if ( maxLen == 0 )
	{
		//printf( "%s\n", read ) ;
		return -1 ;
	}

	// Pass the effect of substitution
	fix[ maxLenStats[0] ] = maxLenStats[1] ;
	i = maxLenStats[0] ;
	//printf( "%d: %d %d\n", maxLen, i, maxLenStats[1] ) ;
	int from = i - kmerLength + 1 ;
	if ( from < 0 )
		from = 0 ;
	for ( j = 0, k = from ; k + kmerLength - 1 < readLength && k <= i ; ++k, ++j )
	{
		//printf( "%d %d %d\n", j, k, stored[1 - tag][j] ) ;
		storedKmer[k] = ( stored[1 - tag][j] == 1 ) ? true : false ;
	}

	return maxLenStats[0] ;
}

int ErrorCorrection( char *read, char *qual, KmerCode& kmerCode, int maxCorrection, char badQuality, Store *kmers, int &badPrefix, int &badSuffix, int &info )
{
	int i, j, k ;	
	bool storedKmer[MAX_READ_LENGTH] ;
	int kmerCnt = 0 ;
	int readLength = 0 ;
	int fix[MAX_READ_LENGTH] ; // The fixed character for each untrusted position.
	bool trusted[MAX_READ_LENGTH] ; // Do not correct these makred positions.
	int tag ;
	int from, to ;
	int trimStart = -1 ;
	int ambiguousCnt = 0 ;
	int alternativeCnt = 0 ;
	bool hasAnchor = false ;
	int createAnchorPos = -1 ;
	
	bool noCandidateKmer = false ;

//	KmerCode kmerCode( inCode ) ;
	KmerCode tmpKmerCode( 0 ) ;
	badPrefix = badSuffix = 0 ;
	info = 0 ;

	int kmerLength = kmerCode.GetKmerLength() ;


	kmerCode.Restart() ;
	for ( i = 0 ; i < kmerLength && read[i] ; ++i )
	{
		kmerCode.Append( read[i] ) ;
	}
	if ( !kmers->IsIn( kmerCode ) )
		storedKmer[0] = false ;
	else
	{
		storedKmer[0] = true ;
		hasAnchor = true ;
	}
	kmerCnt = 1 ;
	for (  ; read[i] ; ++i, ++kmerCnt )
	{
		kmerCode.Append( read[i] ) ;

		if ( !kmers->IsIn( kmerCode ) )
			storedKmer[ kmerCnt ] = false ;
		else
		{
			storedKmer[ kmerCnt ] = true ;
			hasAnchor = true ;
		}
	}

	readLength = i ;
	if ( readLength < kmerLength )
		return 0 ;
	
	trimStart = readLength ;
	for ( i = 0 ; i < readLength ; ++i )
		fix[i] = -1 ;

	
	for ( i = 0 ; i < readLength ; ++i )
		trusted[i] = false ;
	tag = -1 ;
	for ( i = 0 ; i < kmerCnt ; ++i )
	{
		if ( storedKmer[i] )
		{
			if ( tag == -1 )
				tag = i ;
		}
		else
		{
			if ( tag != -1 )// && i - tag >= 3 )
			{
				for ( j = tag ; j < i + kmerLength - 1 ; ++j )
					trusted[j] = true ;
			}
			tag = -1 ;
		}
	}

	if ( !hasAnchor )
	{
		createAnchorPos = CreateAnchor( read, qual, fix, storedKmer, kmerCode, kmers ) ;
		if ( createAnchorPos == -1 )
		{
			info = 3 ;
			return -1 ;
		}

		++alternativeCnt ;
	}

	//printf( "%d %d %d\n", kmerLength, kmerCnt, i ) ;	

#ifdef DEBUG
	printf( "%s\n", read ) ;
	for ( i = 0 ; i < kmerCnt ; ++i )
		printf( "%d", storedKmer[i] ) ;
	printf( "\n" ) ;
#endif
	//exit( 1 ) ;
	// All the kmers are reliable. 
	for ( i = 0 ; i < kmerCnt ; ++i )
	{
		if ( storedKmer[i] == false )
			break ;
	}
	/*if ( !strcmp( read, "TTGCGTAAGATGGGGGTACCCACGTGGTGTCAGAGTGTCTCTCATTTCGGTTTGATCTACGCAGATCTACAAAAAATGCGGGAGAATAGACGCAGAGTTCT" ) )
	{
		printf( "## %d\n", i ) ;
		exit( 1 ) ;
	}*/
	if ( i >= kmerCnt )
		return 0 ;

	// Find the first trusted kmer
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

	if ( longestStoredKmerCnt == 0 )
	{
		info = 3 ;
		return -1 ;
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
		char backupC = '\0' ;
		if ( createAnchorPos != -1 )
		{
			backupC = read[ createAnchorPos ] ;
			read[ createAnchorPos ] = numToNuc[ fix[ createAnchorPos ] ] ;
		}

		if ( longestStoredKmerCnt < kmerLength )
		{
			for ( i = k - 1 + 1 ; i < k - 1 + kmerLength - 1 + 1 ; ++i  )
			{
				kmerCode.Append( read[i] ) ;
			}
		}
		else
		{
			// Adjust the anchor if necessary
			// Test whether the current anchor is good
			for ( i = k ; i < k + kmerLength ; ++i )
				kmerCode.Append( read[i] ) ;
			int c ;
			for ( c = 0 ; c < 4 ; ++c )
			{
				if ( numToNuc[c] == read[i - 1] )
					continue ;

				if ( numToNuc[c] == read[i - 1] )
					continue ;
				tmpKmerCode = kmerCode ;
				tmpKmerCode.ShiftRight() ;
				tmpKmerCode.Append( numToNuc[c] ) ;
				if ( kmers->IsIn( tmpKmerCode ) ) 
				{	
					// Test whether this branch makes sense
					int t = 0 ;
					for ( t = 0 ; t < kmerLength && read[i + t] ; ++t ) // and it is should be a very good fix 
					{
						tmpKmerCode.Append( read[i + t] ) ;
						if ( !kmers->IsIn( tmpKmerCode ) )
							break ;
					}
					if ( !read[i + t] || t >= kmerLength )
						break ;
				}
			}

			if ( c < 4 )
			{
				//kmerCode.ShiftRight( 1 ) ; Seems wrong
				for ( i = k - 1 ; i < k + kmerLength - 1 ; ++i )
					kmerCode.Append( read[i] ) ;
			}
			else
			{

				// Adjust the right side
				for ( j = kmerLength / 2 - 1 ; j >= 0 ; --j )
				{
					int c ;
					KmerCode tmpKmerCode( kmerLength ) ;

					for ( i = k - j - 1 ; i < k - j + kmerLength - 1 ; ++i )
						kmerCode.Append( read[i] ) ;
					//printf( "%d %d %c\n", j, i, read[i - 1] ) ;
					for ( c = 0 ; c < 4 ; ++c )
					{
						if ( numToNuc[c] == read[i - 1] )
							continue ;
						tmpKmerCode = kmerCode ;
						tmpKmerCode.ShiftRight() ;
						tmpKmerCode.Append( numToNuc[c] ) ;
						if ( kmers->IsIn( tmpKmerCode ) ) 
						{	
							++alternativeCnt ;
							// Test whether this branch makes sense
							int t = 0 ;
							for ( t = 0 ; t <= kmerLength / 2 && read[i + t] ; ++t ) // and it is should be a very good fix 
							{
								tmpKmerCode.Append( read[i + t] ) ;
								if ( !kmers->IsIn( tmpKmerCode ) )
									break ;
							}
							if ( t > kmerLength / 2 )
								break ;
						}
					}
					if ( c < 4 )
					{
						// adjust the anchor
						--i ;
						kmerCode.ShiftRight() ;
						break ;
					}
				}
			}
		}

		if ( createAnchorPos != -1 )
		{
			read[ createAnchorPos ] = backupC ;
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
		int testCnt = 0 ;
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
			++testCnt ; 
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
			//printf( "%d\n", j ) ;
		}
#ifdef DEBUG
		printf( "+hi %d: %d %d=>%d, (%d)\n", i, maxTo, to, maxChange, maxCnt ) ;
#endif
		if ( testCnt > 1 )
			alternativeCnt += ( testCnt - 1 ) ;
		// TODO: if maxTo is far from i, then we may in a repeat. Try keep this base unfixed
		//       see whether the next fixing makes sense.	
		if ( maxTo == -1 || ( maxCnt > 1 && ( maxTo <= to || to - i + 1 < kmerLength ) ) )
		{
			//printf( "+%s\n", read ) ;

			// trim
			//if ( maxCnt > 1 )
			//{
			if ( maxTo == -1 )
				noCandidateKmer = true ;
			trimStart = i ;
			info = 2 ;
			break ;
			//}
			//else
			//return 0 ;
		}
		fix[i] = maxChange ;
		if ( maxCnt > 1 )
		{
			// Remove the effect of the ambiguous fixing
			fix[i] = -2 ;
			++ambiguousCnt ;
		}

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
				kmerCode.ShiftRight() ;
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
				if ( fix[i] < 0 )
					kmerCode.Append( read[i] ) ;
				else
					kmerCode.Append( numToNuc[ fix[i] ] ) ;
			}
			continue ;
		}
	}

	// Scan towards left
	//printf( "%d\n", tag ) ;
	if ( tag == 0 )
	{
		// Force skip
		tag = -1 ;
	}
	else
	{
		char backupC = '\0' ;
		if ( createAnchorPos != -1 )
		{
			backupC = read[ createAnchorPos ] ;
			read[ createAnchorPos ] = numToNuc[ fix[ createAnchorPos ] ] ;
		}

		if ( longestStoredKmerCnt < kmerLength )
		{
			kmerCode.Restart() ;
			for ( i = tag + 1 - 1 ; i < tag + kmerLength - 1 ; ++i )
			{
				kmerCode.Append( read[i] ) ;	
			}
			kmerCode.Append( 'A' ) ;
		}
		else
		{
			// Test whether the current anchor is good
			int c ;
			j = -1 ;
			for ( i = tag -1 ; i < tag - 1 + kmerLength ; ++i )
				kmerCode.Append( read[i] ) ;

			for ( c = 0 ; c < 4 ; ++c )
			{
				if ( numToNuc[c] == read[tag + j] )
					continue ;
				tmpKmerCode = kmerCode ;
				tmpKmerCode.Append( 'A' ) ;
				tmpKmerCode.Prepend( numToNuc[c] ) ;
				if ( kmers->IsIn( tmpKmerCode ) ) 
				{	
					// Test whether this branch makes sense
					int t = 0 ;
					for ( t = 0 ; t <= kmerLength - 1 && tag + j - t - 1 >= 0 ; ++t ) // and it is should be a very good fix 
					{
						tmpKmerCode.Prepend( read[tag + j - t - 1] ) ;
						if ( !kmers->IsIn( tmpKmerCode ) )
							break ;
					}
					if ( t > kmerLength - 1 || tag + j -t -1 < 0 )
						break ;
				}
			}

			if ( c < 4 )
			{
				j = 0 ;
				kmerCode.Restart() ;
				for ( i = tag + j ; i < tag + j + kmerLength ; ++i )
					kmerCode.Append( read[i] ) ;
			}
			else
			{
				// Adjust the left side of the anchor
				for ( j = kmerLength / 2 - 1 ; j >= 0 ; --j )
				{
					int c ;
					KmerCode tmpKmerCode( kmerLength ) ;

					kmerCode.Restart() ;
					for ( i = tag + j ; i < tag + j + kmerLength ; ++i )
						kmerCode.Append( read[i] ) ;
					//printf( "%d %d %c\n", j, i, read[i - 1] ) ;
					for ( c = 0 ; c < 4 ; ++c )
					{
						if ( numToNuc[c] == read[tag + j] )
							continue ;
						tmpKmerCode = kmerCode ;
						tmpKmerCode.Append( 'A' ) ;
						tmpKmerCode.Prepend( numToNuc[c] ) ;
						if ( kmers->IsIn( tmpKmerCode ) ) 
						{
							++alternativeCnt ;
							// Test whether this branch makes sense
							int t = 0 ;
							for ( t = 0 ; t <= kmerLength / 2 && tag + j - t - 1 >= 0 ; ++t ) // and it is should be a very good fix 
							{
								tmpKmerCode.Prepend( read[tag + j - t - 1] ) ;
								if ( !kmers->IsIn( tmpKmerCode ) )
									break ;
							}
							if ( t > kmerLength / 2 )
								break ;
						}
					}
					if ( c < 4 )
					{
						// adjust the anchor
						tag = tag + j + 1 ;
						kmerCode.Append( 'A' ) ;
						break ;
					}
				}
			}
		}

		if ( createAnchorPos != -1 )
		{
			read[ createAnchorPos ] = backupC ;
		}
	}
	//kmerCode = ( kmerCode << (uint64_t)2 ) & mask ;
	for ( i = tag - 1 ; i >= 0 ;  )
	{
		KmerCode fixedKmerCode( kmerCode ) ;
		int minTo = readLength + 1 ;
		int minChange = -1 ;
		int minCnt = 0 ;
		int testCnt = 0 ;
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
			++testCnt ;

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
#ifdef DEBUG
		printf( "-hi %d: %d %d=>%d, (%d)\n", i, minTo, to, minChange, minCnt ) ;	
#endif
		if ( testCnt > 1 )
			alternativeCnt += ( testCnt - 1 ) ;
		
		if ( minTo == readLength + 1 || ( minCnt > 1 && ( minTo >= to || i - to + 1 < kmerLength ) ) )
		{
			//printf( "-%s\n", read ) ;
			//return -1 ;
			if ( minTo == readLength + 1 )
				noCandidateKmer = true ;
			badPrefix = i + 1 ;
			info = 2 ;
			break ;
		}
		//printf( "---%d %d\n", minChange,  nucToNum[read[0] - 'A'] ) ;
		fix[i] = minChange ;
		if ( minCnt > 1 )
		{
			fix[i] = -2 ;
			++ambiguousCnt ;
		}

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
				if ( fix[i] < 0 )
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

#ifdef DEBUG
	printf( "fix: ") ;
	for ( i = 0 ; i < readLength; ++i )
	{
		if ( fix[i] == -1  )
			printf( "5" ) ;
		else if ( fix[i] == -2  )
			printf( "6" ) ;
		else
			printf( "%d", fix[i] ) ;
	}
	printf( "\n" ) ;
#endif

	int ret = 0 ;
	double correctCnt = 0 ;
	//printf( "%s\n%d\n", read, ret ) ;
	/*for ( i = badPrefix ; i < trimStart ; ++i )
	{
		if ( trusted[i] && fix[i] != -1 )
			return -1 ;
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
	}*/

	bool overCorrected = false ;
	int adjustMaxCor = 0 ;
	for ( i = 0 ; i < readLength ; ++i )	
		if ( trusted[i] && fix[i] != -1 )
			break ;
	if ( alternativeCnt == 0 && i >= readLength )
		adjustMaxCor = 1 ;

	for ( i = 0 ; i < readLength ; ++i )
	{
		int overCorrectWindow = 20 ;
		if ( i >= overCorrectWindow && ( fix[i - overCorrectWindow] >= 0 && read[i - overCorrectWindow] != 'N' ) )
		{
			if ( qual[0] != '\0' && qual[i - overCorrectWindow] <= badQuality )
				correctCnt -= 0.5 ;
			else
				--correctCnt ;
		}
		if ( fix[i] >= 0 && read[i] != 'N' )
		{
			if ( qual[0] != '\0' && qual[i] <= badQuality )
				correctCnt += 0.5 ;
			else
				++correctCnt ;
		}
		int tmp = maxCorrection ;
		if ( i >= overCorrectWindow  && i + overCorrectWindow - 1 < readLength )
			tmp += adjustMaxCor ;
		if ( correctCnt > tmp ) //maxCorrection )
		{
			if ( fix[i] >= 0 )
			{
				fix[i] = 4 ;
				overCorrected = true ;
			}
		}
	}

	if ( overCorrected )
	{
		for ( i = 0 ; i < readLength ; ++i )	
		{
			if ( fix[i] == 4 )
			{
				fix[i] = -1 ; 
				int tag = i ;
				for ( j = i - 1 ; j >= 0 && j >= tag - kmerLength + 1 ; --j )
					if ( fix[j] >= 0 )
					{
						fix[j] = -1 ;
						tag = j ;
					}

				tag = i ;
				for ( j = i + 1 ; j < readLength && j <= tag + kmerLength - 1 ; ++j )
					if ( fix[j] >= 0 )
					{
						fix[j] = -1 ;
						tag = j ;
					}
				i = j - 1 ; // the minus 1 here is to compensate for the ++i in the outer-loop
			}
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

#ifdef DEBUG
	printf( "prefix=%d suffix=%d\n", badPrefix, trimStart ) ;
	printf( "fix: ") ;
	for ( i = 0 ; i < readLength; ++i )
	{
		if ( fix[i] == -1  )
			printf( "5" ) ;
		else if ( fix[i] == -2  )
			printf( "6" ) ;
		else
			printf( "%d", fix[i] ) ;
	}
	printf( "\n" ) ;
#endif

	for ( i = badPrefix ; i < trimStart ; ++i )
	{
		if ( fix[i] < 0 )
			continue ;
		if ( read[i] != numToNuc[ fix[i] ] )
		{
			read[i] = numToNuc[ fix[i] ] ;
			//if ( read[i] != 'N' )
			if ( qual[0] != '\0' && SET_NEW_QUAL >= 0 )
				qual[i] = (char)SET_NEW_QUAL ;
			++ret ;
		}
	}
	
	//printf( "%d %d\n", ret, trimStart ) ;

	//badPrefix = 0 ;
	badSuffix = readLength - trimStart ;
		
	//if ( ALLOW_TRIMMING ) // else we do partial correction
	//	read[ trimStart ] = '\0' ; 
	//ret += trimmed ;
	
	/*if ( !strcmp( read, "GTAAACGCCTTATCCGGCCTACGGAGGGTGCGGGAATTTGTAGGCCTGATAAGACGCGCAAGCGTCGCATCAGGCAGTCGGCACGGTTGCCGGATGCAGCG" ) )
	{
		printf( "## %d\n", i ) ;
		exit( 1 ) ;
	}*/
	if ( ret == 0 && badPrefix == 0 && badSuffix == 0 && ambiguousCnt > 0 )
	{
		info = 2 ;		
		ret = -1 ;
	}
	else if ( ret == 0 && badPrefix == 0 && badSuffix == 0 && overCorrected )
	{
		info = 1 ;
		ret = -1 ;
	}
	if ( noCandidateKmer )
		info = 3 ;
	/*if ( ret == 0 && badSuffix > 0 )	
	{
		printf( "%d\n", info ) ;
	}*/
	return ret ;
}

int ErrorCorrection_Wrapper( char *read, char *qual, KmerCode& kmerCode, char badQuality, Store *kmers, int &badPrefix, int &badSuffix, int &info )
{
	int correction ;
	int tmpBadPrefix, tmpBadSuffix, tmpInfo ;
	int len = (int)strlen( read ) ;
	int interim = 0 ;

	info = 0 ;
	correction = ErrorCorrection( read, qual, kmerCode, MAX_CORRECTION, badQuality, kmers, 
			tmpBadPrefix, tmpBadSuffix, tmpInfo ) ;
	
	if ( correction == -1 )
	{
		badPrefix = badSuffix = 0 ;
		info = tmpInfo ;
		return -1 ;
	}

	badPrefix = tmpBadPrefix ;
	badSuffix = tmpBadSuffix ;
	info = tmpInfo ;
	if ( badPrefix > len / 2 || badPrefix > kmerCode.GetKmerLength() * 2 )
	{
		int tmp ;
		char c = read[badPrefix] ;
		read[ badPrefix ] = '\0' ;
		tmp = ErrorCorrection( read, qual, kmerCode, MAX_CORRECTION, badQuality, kmers,   
				tmpBadPrefix, tmpBadSuffix, tmpInfo ) ;

		read[ badPrefix ] = c ;
		if ( tmp != -1 )
		{
			badPrefix = tmpBadPrefix ;
			correction += tmp ;
			interim += tmpBadSuffix ;
		}
		if ( tmpInfo > info )
			info = tmpInfo ;
	}

	if ( badSuffix > len / 2 || badSuffix > kmerCode.GetKmerLength() * 2 )
	{
		int tmp ;

		tmp = ErrorCorrection( read + len - badSuffix, qual + len - badSuffix, kmerCode, MAX_CORRECTION, badQuality, kmers, 
				tmpBadPrefix, tmpBadSuffix, tmpInfo ) ;	
		if ( tmp != -1 )
		{
			badSuffix = tmpBadSuffix ;
			correction += tmp ;
			interim += tmpBadPrefix ;
		}

		if ( tmpInfo > info )
			info = tmpInfo ;
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
