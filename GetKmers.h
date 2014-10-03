// Functions used to smaple kmers and get trusted kmers
#ifndef _MOURISL_GETKMERS
#define _MOURISL_GETKMERS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "Reads.h"
#include "Store.h"
#include "KmerCode.h"
#include "utils.h"

// The pattern of sampling for each read. 
// Marks the bits to indicate whether to sample the 
// kmer end at that position.
#define SAMPLE_PATTERN_COUNT 65536 

struct _SamplePattern
{
	char tag[ MAX_READ_LENGTH / 8 ] ;	
} ;

// Strcture for handling multi-threads
struct _SampleKmersThreadArg
{
	int kmerLength ;
	double alpha ;
	//int batchSize ;
	Store *kmers ; 
	Reads *reads ;

	struct _SamplePattern *samplePatterns ;
	
	pthread_mutex_t *lock ;
	pthread_mutex_t *lockPut ;
} ;

struct _StoreKmersThreadArg
{
	int kmerLength ;
	//int batchSize ;
	int *threshold ;
	Store *kmers ;
	Store *trustedKmers ;
	Reads *reads ;
	char goodQuality ;
	char badQuality ;

	pthread_mutex_t *lock ;
} ;

void *SampleKmers_Thread( void *arg ) ;
void SampleKmersInRead( char *read, char *qual, int kmerLength, double alpha, KmerCode &kmerCode, Store *kmers ) ;

void *StoreKmers_Thread( void *arg ) ;
void StoreTrustedKmers( char *read, char *qual, int kmerLength, char badQuality, int *threshold,  
	KmerCode &kmerCode, Store *kmers, Store *trustedKmers ) ;
#endif
