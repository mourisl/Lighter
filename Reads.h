// The class for reading reads from a file
#ifndef _MOURISL_READS
#define _MOURISL_READS

#include <stdio.h>
#include <string.h>
#include "utils.h"

#define MAX_READ_FILE 100

struct _Read
{
	char id[MAX_ID_LENGTH] ;
	char seq[MAX_READ_LENGTH] ;
	char qual[MAX_READ_LENGTH] ;

	int correction ;
	int badPrefix, badSuffix ;
} ;

class Reads
{
	private:
		FILE *fp[MAX_READ_FILE] ;
		FILE *outputFp[MAX_READ_FILE] ;
		int FILE_TYPE[MAX_READ_FILE] ; // 0-FASTA, 1-FASTQ
		int fpUsed ;
		int currentFpInd ;
		char outputDirectory[256] ;

		void GetFileName( char *in, char *out ) 
		{
			int i, j ;
			int len = strlen( in ) ;
			for ( i = len ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
				;
			for ( j = len ; j >= 0 && in[j] != '/' ; --j )
				;
			if ( i >= 0 && in[i] == '.' )
			{
				in[i] = '\0' ;
				strcpy( out, in + j + 1 ) ;
				in[i] = '.' ;
			}
			else
			{
				strcpy( out, in + j + 1 ) ;
			}
		}

	public:
		char id[MAX_ID_LENGTH] ;
		char seq[MAX_READ_LENGTH] ;
		char qual[MAX_READ_LENGTH] ;

		Reads(): fpUsed(0), currentFpInd(0)
		{
			strcpy( outputDirectory, "./" ) ;
		}

		~Reads()
		{
			int i ;
			for ( i = 0 ; i < fpUsed ; ++i )
			{
				fclose( fp[i] ) ;
				fclose( outputFp[i] ) ;
			}
		}

		void AddReadFile( char *file )
		{
			if ( fpUsed >= MAX_READ_FILE )
			{
				printf( "The number of read files exceeds the limit %d.\n", MAX_READ_FILE ) ;
				exit( 1 ) ;
			}
			char buffer[1024], fileName[1024] ;
			fp[ fpUsed ] = fopen( file, "r" ) ;
			fscanf( fp[ fpUsed ], "%s", buffer ) ;
			if ( buffer[0] == '>' )
			{
				FILE_TYPE[ fpUsed ] = 0 ;
				qual[0] = '\0' ;
			}
			else
			{
				FILE_TYPE[ fpUsed ] = 1 ;
			}
			fclose( fp[fpUsed] ) ;
			fp[fpUsed] = fopen( file, "r" ) ;

			GetFileName( file, fileName ) ;
			if ( FILE_TYPE[ fpUsed ] == 1 )
				sprintf( buffer, "%s/%s.cor.fq", outputDirectory, fileName ) ;
			else
				sprintf( buffer, "%s/%s.cor.fa", outputDirectory, fileName ) ;
			outputFp[ fpUsed ] = fopen( buffer, "w" ) ;

			++fpUsed ;
		}

		void SetOutputDirectory( char *d )
		{
			strcpy( outputDirectory, d ) ;
		}

		bool HasQuality()
		{
			return ( FILE_TYPE[ currentFpInd ] != 0  ) ;
		}

		void Rewind() 
		{
			int i ;
			for ( i = 0 ; i < fpUsed ; ++i )
			{
				rewind( fp[i] ) ;
				rewind( outputFp[i] ) ;
			}
			currentFpInd = 0 ;
		}

		int Next() 
		{
			int len ;
			char buffer[2048] ;
			FILE *lfp ;
			while ( currentFpInd < fpUsed && ( fgets( id, sizeof( id ), fp[ currentFpInd ] ) == NULL ) )
			{
				++currentFpInd ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;
			lfp = fp[ currentFpInd ] ;
			if ( FILE_TYPE[ currentFpInd ] == 0 )
			{
				fgets( seq, sizeof( seq ), lfp ) ;
			}
			else if ( FILE_TYPE[ currentFpInd ] == 1 )
			{
				fgets( seq, sizeof( seq ), lfp ) ;
				fgets( buffer, sizeof( buffer ), lfp ) ;
				fgets( qual, sizeof( qual ), lfp ) ;
			}
			// Clean the return symbol	
			len = strlen( id ) ;		
			if ( id[len - 1] == '\n')
				id[len - 1] = '\0' ;
			len = strlen( seq ) ;
			if ( seq[len - 1] == '\n' )
				seq[len - 1] = '\0' ;

			if ( FILE_TYPE[ currentFpInd ] == 1 )
			{
				if ( qual[len - 1] == '\n' )
					qual[len - 1] = '\0' ;
			}
			return 1 ;
		}

		int NextWithBuffer( char *id, char *seq, char *qual, bool removeReturn = true, bool stopWhenFileEnds = false ) 
		{
			int len ;
			char buffer[2048] ;
			FILE *lfp ;
			while ( currentFpInd < fpUsed && ( fgets( id, sizeof( char ) * MAX_ID_LENGTH, 
					fp[ currentFpInd ] ) == NULL ) )
			{
				++currentFpInd ;
				if ( stopWhenFileEnds )
					return -1 ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;
			lfp = fp[ currentFpInd ] ;
			if ( FILE_TYPE[ currentFpInd ] == 0 )
			{
				fgets( seq, sizeof( char ) * MAX_READ_LENGTH, lfp ) ;
				qual[0] = '\0' ;
			}
			else if ( FILE_TYPE[ currentFpInd ] == 1 )
			{
				fgets( seq, sizeof( char ) * MAX_READ_LENGTH, lfp ) ;
				fgets( buffer, sizeof( buffer ), lfp ) ;
				fgets( qual, sizeof( char ) * MAX_READ_LENGTH, lfp ) ;
			}
			// Clean the return symbol	
			if ( removeReturn )
			{
				len = strlen( id ) ;		
				if ( id[len - 1] == '\n')
					id[len - 1] = '\0' ;
				len = strlen( seq ) ;
				if ( seq[len - 1] == '\n' )
					seq[len - 1] = '\0' ;

				if ( FILE_TYPE[ currentFpInd ] == 1 )
				{
					if ( qual[len - 1] == '\n' )
						qual[len - 1] = '\0' ;
				}
			}
			return 1 ;
		}

		void Output( int correction, int badPrefix, int badSuffix, bool allowTrimming )
		{
			char buffer[1024] ;
			if ( correction == 0 && badPrefix == 0 && badSuffix == 0 )
			{
				fprintf( outputFp[ currentFpInd ], "%s\n%s\n", id, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
					fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
			}
			else if ( correction == -1 )
			{
				/*if ( !strcmp( readId + 1, read ) ) 
				  {
				  printf( "%s\n%s\n", readId, read ) ;
				  }*/
				fprintf( outputFp[ currentFpInd], "%s unfixable_error\n%s\n", id, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
					fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
				//printf( "%s\n%s\n", readId, read ) ;
			}
			else
			{
				char buffer1[20] = "" ;
				char buffer2[20] = "" ;
				char buffer3[20] = "" ;
				/*if ( allowTrimming )
					strcpy( buffer2, "trimmed" ) ;
				else
					strcpy( buffer2, "bad_suffix" ) ;

				if ( badSuffix == 0 )
					strcpy( buffer, "cor" ) ;
				else if ( trimmed > 0 && correction == 0 )
					sprintf( buffer, "%s=%d", buffer2, trimmed ) ;
				else 
					sprintf( buffer, "cor %s=%d", buffer2, trimmed ) ;*/

				if ( correction > 0 )
					strcpy( buffer1, " cor" ) ;
				if ( badPrefix > 0 )
					sprintf( buffer2, " bad_prefix=%d", badPrefix ) ;
				if ( badSuffix > 0 )
				{
					if ( allowTrimming )
						sprintf( buffer3, " trimmed=%d", badSuffix ) ;
					else
						sprintf( buffer3, " bad_suffix=%d", badSuffix ) ;
				}

				fprintf( outputFp[ currentFpInd ], "%s%s%s%s\n%s\n", id, buffer1, buffer2, buffer3, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
				{
					if ( allowTrimming )
						qual[ strlen( qual ) - badSuffix ] = '\0' ;
					fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
				}
			}
		}

		// Get a batch of reads, it terminates until the buffer is full or 
		// the file ends.
		int GetBatch( struct _Read *readBatch, int maxBatchSize, bool trimReturn, bool stopWhenFileEnds  )
		{
			int batchSize = 0 ;
			while ( batchSize < maxBatchSize ) 
			{
				int tmp = NextWithBuffer( readBatch[ batchSize].id, readBatch[batchSize].seq,
							readBatch[batchSize].qual, trimReturn, stopWhenFileEnds ) ;
				if ( tmp == -1 && batchSize > 0 )
				{
					--currentFpInd ;
					return batchSize ; // Finished read current file.
				}
				else if ( tmp == -1 && batchSize == 0 )
					continue ; // The current read file is empty
				else if ( tmp == 0 && batchSize == 0 )
					return 0 ; // Finished reading	

				++batchSize ;
			}
			return batchSize ;
		}

		void OutputBatch( struct _Read *readBatch, int batchSize, bool allowTrimming )
		{
			int i ;
			for ( i = 0 ; i < batchSize ; ++i )
			{
				char *id = readBatch[i].id ;
				char *seq = readBatch[i].seq ;
				char *qual = readBatch[i].qual ;
				int correction = readBatch[i].correction ;
				int badPrefix = readBatch[i].badPrefix ;
				int badSuffix = readBatch[i].badSuffix ;
				
				char buffer[1024] ;
				if ( correction == 0 && badPrefix == 0 && badSuffix == 0 )
				{
					fprintf( outputFp[ currentFpInd ], "%s\n%s\n", id, seq ) ;
					if ( FILE_TYPE[ currentFpInd ] != 0 )
						fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
				}
				else if ( correction == -1 )
				{
					/*if ( !strcmp( readId + 1, read ) ) 
					  {
					  printf( "%s\n%s\n", readId, read ) ;
					  }*/
					fprintf( outputFp[ currentFpInd], "%s unfixable_error\n%s\n", id, seq ) ;
					if ( FILE_TYPE[ currentFpInd ] != 0 )
						fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
					//printf( "%s\n%s\n", readId, read ) ;
				}
				else
				{
					char buffer1[20] = "" ;
					char buffer2[20] = "" ;
					char buffer3[20] = "" ;

					if ( correction > 0 )
						strcpy( buffer1, " cor" ) ;
					if ( badPrefix > 0 )
						sprintf( buffer2, " bad_prefix=%d", badPrefix ) ;
					if ( badSuffix > 0 )
					{
						if ( allowTrimming )
							sprintf( buffer3, " trimmed=%d", badSuffix ) ;
						else
							sprintf( buffer3, " bad_suffix=%d", badSuffix ) ;
					}

					fprintf( outputFp[ currentFpInd ], "%s%s%s%s\n%s\n", id, buffer1, buffer2, buffer3, seq ) ;
					if ( FILE_TYPE[ currentFpInd ] != 0 )
					{
						if ( allowTrimming )
							qual[ strlen( qual ) - badSuffix ] = '\0' ;
						fprintf( outputFp[ currentFpInd ], "+\n%s\n", qual ) ;
					}
				}
			}
		}
} ;

// The class handling read in a batch of reads
/*class ReadBatch
{
} ;*/

#endif
