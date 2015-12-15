// The class for reading reads from a file
#ifndef _MOURISL_READS
#define _MOURISL_READS

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "utils.h"
#include "File.h"

#define MAX_READ_FILE 100

struct _Read
{
	char id[MAX_ID_LENGTH] ;
	char seq[MAX_READ_LENGTH] ;
	char qual[MAX_READ_LENGTH] ;

	int correction ;
	int badPrefix, badSuffix ;
	int info ;
} ;

class Reads
{
	private:
		File fp[MAX_READ_FILE] ;
		File outputFp[MAX_READ_FILE] ;
		int FILE_TYPE[MAX_READ_FILE] ; // 0-FASTA, 1-FASTQ
		int fpUsed ;
		int currentFpInd ;
		char outputDirectory[256] ;
		bool discard ;
		int compressLevel ;

		void GetFileName( char *in, char *out ) 
		{
			int i, j ;
			int len = (int)strlen( in ) ;
			for ( i = len ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
				;
			if ( i >= 0 && !strcmp( &in[i], ".gz" ) )
			{
				int tmp = i ;
				for ( i = i - 1 ; i >= 0 && in[i] != '.' && in[i] != '/' ; --i )
					;
				in[tmp] = '\0' ;
				if ( i >= 0 && ( !strcmp( &in[i], ".fastq" ) || !strcmp( &in[i], ".fasta" ) ||
					!strcmp( &in[i], ".fq" ) || !strcmp( &in[i], ".fa" ) ) )
				{
					;
				}
				else
				{
					i = tmp ;
				}
				in[tmp] = '.' ;
			}

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
			compressLevel = 1 ;
			strcpy( outputDirectory, "./" ) ;
		}

		~Reads()
		{
			int i ;
			for ( i = 0 ; i < fpUsed ; ++i )
			{
				fp[i].Close() ;
				outputFp[i].Close() ;
			}
		}

		void SetDiscard( bool d )
		{
			discard = d ;
		} 

		void SetCompressLevel( int cl )
		{
			compressLevel = cl ;
		}

		void AddReadFile( char *file )
		{
			if ( fpUsed >= MAX_READ_FILE )
			{
				fprintf( stderr, "The number of read files exceeds the limit %d.\n", MAX_READ_FILE ) ;
				exit( 1 ) ;
			}
			char buffer[1024], fileName[1024] ;
			fp[ fpUsed ].Open( file, "r" ) ;
			
			fp[ fpUsed ].Gets( buffer, sizeof( buffer ) ) ;
			if ( buffer[0] == '>' )
			{
				FILE_TYPE[ fpUsed ] = 0 ;
				qual[0] = '\0' ;
			}
			else if ( buffer[0] == '@' )
			{
				FILE_TYPE[ fpUsed ] = 1 ;
			}
			else
			{
				fprintf( stderr, "\"%s\"'s format is wrong.\n", file ) ;
				exit( 1 ) ;
			}

			fp[fpUsed].Close() ;

			fp[fpUsed].Open( file, "r" ) ;

			GetFileName( file, fileName ) ;

			int len = strlen( file ) ;
			if ( file[ len - 2] == 'g' && file[ len - 1 ] == 'z' && compressLevel > 0 )
			{
				if ( FILE_TYPE[ fpUsed ] == 1 )
					sprintf( buffer, "%s/%s.cor.fq.gz", outputDirectory, fileName ) ;
				else
					sprintf( buffer, "%s/%s.cor.fa.gz", outputDirectory, fileName ) ;
			}
			else
			{
				if ( FILE_TYPE[ fpUsed ] == 1 )
					sprintf( buffer, "%s/%s.cor.fq", outputDirectory, fileName ) ;
				else
					sprintf( buffer, "%s/%s.cor.fa", outputDirectory, fileName ) ;
			}
			
			outputFp[ fpUsed ].SetCompressLevel( compressLevel ) ;
			outputFp[ fpUsed ].Open( buffer, "w" ) ;
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
				fp[i].Rewind() ;
				outputFp[i].Rewind() ;
			}
			currentFpInd = 0 ;
		}

		int Next() 
		{
			int len ;
			char buffer[2048] ;
			while ( currentFpInd < fpUsed && ( fp[currentFpInd].Gets( id, sizeof( id ) ) == NULL ) )
			{
				++currentFpInd ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;

			File &lfp = fp[currentFpInd] ;
			if ( FILE_TYPE[ currentFpInd ] == 0 )
			{
				lfp.Gets( seq, sizeof( seq ) ) ;
			}
			else if ( FILE_TYPE[ currentFpInd ] == 1 )
			{
				lfp.Gets( seq, sizeof( seq ) ) ;
				lfp.Gets( buffer, sizeof( buffer ) ) ;
				lfp.Gets( qual, sizeof( qual ) ) ;
			}
			// Clean the return symbol	
			len = (int)strlen( id ) ;
			if ( id[len - 1] == '\n')
				id[len - 1] = '\0' ;
			len = (int)strlen( seq ) ;
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
			while ( currentFpInd < fpUsed && 
				( fp[currentFpInd].Gets( id, sizeof( char ) * MAX_ID_LENGTH ) == NULL ) )
			{
				++currentFpInd ;
				if ( stopWhenFileEnds )
					return -1 ;
			}
			if ( currentFpInd >= fpUsed )
				return 0 ;
			File &lfp = fp[ currentFpInd ] ;
			if ( FILE_TYPE[ currentFpInd ] == 0 )
			{
				lfp.Gets( seq, sizeof( char ) * MAX_READ_LENGTH ) ;
				qual[0] = '\0' ;
			}
			else if ( FILE_TYPE[ currentFpInd ] == 1 )
			{
				lfp.Gets( seq, sizeof( char ) * MAX_READ_LENGTH ) ;
				lfp.Gets( buffer, sizeof( buffer ) ) ;
				lfp.Gets( qual, sizeof( char ) * MAX_READ_LENGTH ) ;
			}
			// Clean the return symbol	
			if ( removeReturn )
			{
				len = (int)strlen( id ) ;
				if ( id[len - 1] == '\n')
					id[len - 1] = '\0' ;
				len = (int)strlen( seq ) ;
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

		void Output( int correction, int badPrefix, int badSuffix, int info, bool allowTrimming )
		{
			/* char buffer[1024] ; */
			char failReason[4] = "" ;
			if ( info == 1 )
				strcpy( failReason, " oc" ) ;
			else if ( info == 2 )
				strcpy( failReason, " ak" ) ;
			else if ( info == 3 )
				strcpy( failReason, " lc" ) ;

			if ( correction == 0 && badPrefix == 0 && badSuffix == 0 )
			{
				outputFp[currentFpInd].Printf( "%s\n%s\n", id, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
					outputFp[currentFpInd].Printf( "+\n%s\n", qual ) ;
			}
			else if ( correction == -1 )
			{
				/*if ( !strcmp( readId + 1, read ) ) 
				  {
				  printf( "%s\n%s\n", readId, read ) ;
				  }*/

				if ( discard )
					return ;
									
				outputFp[ currentFpInd].Printf( "%s unfixable_error%s\n%s\n", id, failReason, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
					outputFp[ currentFpInd ].Printf( "+\n%s\n", qual ) ;
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

				outputFp[ currentFpInd ].Printf( "%s%s%s%s%s\n%s\n", id, buffer1, buffer2, buffer3, failReason, seq ) ;
				if ( FILE_TYPE[ currentFpInd ] != 0 )
				{
					if ( allowTrimming )
						qual[ strlen( qual ) - badSuffix ] = '\0' ;
					outputFp[ currentFpInd ].Printf( "+\n%s\n", qual ) ;
				}
			}
		}

		// Get a batch of reads, it terminates until the buffer is full or 
		// the file ends.
		int GetBatch( struct _Read *readBatch, int maxBatchSize, int &fileInd, bool trimReturn, bool stopWhenFileEnds  )
		{
			int batchSize = 0 ;
			while ( batchSize < maxBatchSize ) 
			{
				int tmp = NextWithBuffer( readBatch[ batchSize].id, readBatch[batchSize].seq,
							readBatch[batchSize].qual, trimReturn, stopWhenFileEnds ) ;
				if ( tmp <= 0 && batchSize > 0 )
				{
					--currentFpInd ;
					fileInd = currentFpInd ;
					return batchSize ; // Finished read current file.
				}
				else if ( tmp == -1 && batchSize == 0 )
					continue ; // The current read file is empty
				else if ( tmp == 0 && batchSize == 0 )
				{
					fileInd = currentFpInd ;
					return 0 ; // Finished reading	
				}

				++batchSize ;
			}

			fileInd = currentFpInd ;
			return batchSize ;
		}

		void OutputBatch( struct _Read *readBatch, int batchSize, bool allowTrimming, int fileInd = -1 )
		{
			int i ;
			if ( fileInd == -1 )
				fileInd = currentFpInd ;
			for ( i = 0 ; i < batchSize ; ++i )
			{
				char *id = readBatch[i].id ;
				char *seq = readBatch[i].seq ;
				char *qual = readBatch[i].qual ;
				int correction = readBatch[i].correction ;
				int badPrefix = readBatch[i].badPrefix ;
				int badSuffix = readBatch[i].badSuffix ;
				int info = readBatch[i].info ;
				
				char failReason[4] = "" ;
				if ( info == 1 )
					strcpy( failReason, " oc" ) ;
				else if ( info == 2 )
					strcpy( failReason, " ak" ) ;
				else if ( info == 3 )
					strcpy( failReason, " lc" ) ;

				if ( correction == 0 && badPrefix == 0 && badSuffix == 0 )
				{
					outputFp[ fileInd ].Printf( "%s\n%s\n", id, seq ) ;
					if ( FILE_TYPE[ fileInd ] != 0 )
						outputFp[ fileInd ].Printf( "+\n%s\n", qual ) ;
				}
				else if ( correction == -1 )
				{
					/*if ( !strcmp( readId + 1, read ) ) 
					  {
					  printf( "%s\n%s\n", readId, read ) ;
					  }*/
					if ( discard )
						continue ;
					outputFp[ fileInd ].Printf( "%s unfixable_error%s\n%s\n", id, failReason, seq ) ;
					if ( FILE_TYPE[ fileInd ] != 0 )
						outputFp[ fileInd ].Printf( "+\n%s\n", qual ) ;
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

					outputFp[ fileInd ].Printf( "%s%s%s%s%s\n%s\n", id, buffer1, buffer2, buffer3, failReason, seq ) ;
					if ( FILE_TYPE[ fileInd ] != 0 )
					{
						if ( allowTrimming )
							qual[ strlen( qual ) - badSuffix ] = '\0' ;
						outputFp[ fileInd ].Printf( "+\n%s\n", qual ) ;
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
