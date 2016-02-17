#ifndef _MOURISL_CORRECTOR_UTIL_HEADER
#define _MOURISL_CORRECTOR_UTIL_HEADER

#include <stdint.h>


#define MAX_READ_LENGTH 1024
#define MAX_ID_LENGTH 512

// By changing this definition, the user can use longer kmer length.
// Using longer kmer does NOT mean better accuracy!!
#define MAX_KMER_LENGTH 32

#define READ_BUFFER_PER_THREAD 1024
/*char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2, 
	-1, -1, -1, -1, -1, -1, -1,
	-1, -1, -1, -1, -1, 3,
	-1, -1, -1, -1, -1, -1 } ;

char numToNuc[26] = {'A', 'C', 'G', 'T'} ;*/


#endif
