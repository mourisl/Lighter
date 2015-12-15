#ifndef _MOURI_ERROR_CORRECTION
#define _MOURI_ERROR_CORRECTION

#include <stdio.h>

#include "Store.h"
#include "KmerCode.h"
#include "Reads.h"

struct _ErrorCorrectionThreadArg
{
	int kmerLength ;
	Store *trustedKmers ;
	struct _Read *readBatch ;
	int batchSize ;
	int batchUsed ;
	int batchFinished ;
	char badQuality ;

	pthread_mutex_t *lock ;
} ;


void *ErrorCorrection_Thread( void *arg ) ;

//@ return:0: this read is correct. -1-this read is unfixable. Otherwise, return the number 
// of corrected positions.
int ErrorCorrection_Wrapper( char *read, char *qual, KmerCode& kmerCode, char badQuality, Store *kmers, int &badPrefix, int &badSuffix, int &info ) ;

#endif
