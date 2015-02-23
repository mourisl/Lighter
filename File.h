// The wrapper to handle reading and writingregular files and compressed files

#ifndef _MOURISL_FILE
#define _MOURISL_FILE

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

#define UNCOMPRESSED_FILE 0
#define COMPRESSED_FILE 1

extern bool zlibVersionChecked ; 

class File
{
private:
	bool type ; 
	FILE *fp ;
	gzFile gzFp ;
	int gzCompressLevel ;

	bool opened ;
public:
	File() { opened = false ; } 
	~File() 
	{
		if ( !opened )
			return ;
		if ( type == COMPRESSED_FILE )
			gzclose( gzFp ) ;
		else if ( type == UNCOMPRESSED_FILE )
			fclose( fp ) ;
		opened = false ;
	}

	void Open( char *fileName, const char *mode )
	{
		opened = true ;
		// Test it is gz or normal file by looking at the last to bit
		int len = strlen( fileName ) ;
		if ( fileName[len - 2] == 'g' && fileName[len - 1] == 'z' )
			type = COMPRESSED_FILE ;
		else
			type = UNCOMPRESSED_FILE ;

		if ( type == COMPRESSED_FILE )
		{
			char modeBuffer[5] ;
			strcpy( modeBuffer, mode ) ;

			if ( modeBuffer[0] == 'w' )
			{
				modeBuffer[1] = gzCompressLevel + '0' ;
				modeBuffer[2] = '\0' ;
			}

			gzFp = gzopen( fileName, modeBuffer ) ;
			if ( gzFp == Z_NULL )
			{
				fprintf( stderr, "ERROR: Could not access file %s\n", fileName ) ;
				exit( 1 ) ;
			}

			if ( zlibVersionChecked == false )
			{
				zlibVersionChecked = true ;
#ifdef ZLIB_VERNUM	
				if ( ZLIB_VERNUM < 0x1240 )
					fprintf( stderr, "WARNING: zlib version on your system is %s(< 1.2.4). "
							"Newer veresion (>=1.2.4) is much faster.\n", ZLIB_VERSION ) ;
#else
					fprintf( stderr, "WARNING: Unknown zlib version. Newer veresion (>=1.2.4) is much faster.\n" ) ;
#endif
			}
		}
		else if ( type == UNCOMPRESSED_FILE )
		{
			fp = NULL ;
			fp = fopen( fileName, mode ) ;
			if ( fp == NULL )
			{
				fprintf( stderr, "ERROR: Could not access file %s\n", fileName ) ;
				exit( 1 ) ;
			}
		}
	}
	
	void Close()
	{
		if ( !opened )
			return ;
		if ( type == COMPRESSED_FILE )
			gzclose( gzFp ) ;
		else if ( type == UNCOMPRESSED_FILE )
			fclose( fp ) ;
		opened = false ;
	}

	char *Gets( char *buf, int len )
	{
		if ( type == COMPRESSED_FILE )
		{
			return gzgets( gzFp, buf, len ) ;
		}
		else if ( type == UNCOMPRESSED_FILE )
		{
			return fgets( buf, len, fp ) ;	
		}
		return NULL ;
	}
	
	int Puts( char *buf )
	{
		if ( type == COMPRESSED_FILE )
		{
			return gzwrite( gzFp, buf, strlen( buf ) ) ;
		}
		else if ( type == UNCOMPRESSED_FILE )
		{
			return fputs( buf, fp ) ;	
		}
		return 0 ;
	}

	int Printf( const char *fmt, ... )
	{
		char buffer[1024] ;
		va_list args ;
		va_start( args, fmt ) ;
		vsprintf( buffer, fmt, args ) ;
		if ( type == COMPRESSED_FILE )
		{
			return gzwrite( gzFp, buffer, strlen( buffer ) ) ;
		}
		else if ( type == UNCOMPRESSED_FILE )
		{
			return fputs( buffer, fp ) ;
		}
		return 0 ;
	}

	void Rewind()
	{
		if ( type == COMPRESSED_FILE )
		{
			gzrewind( gzFp ) ;
		}
		else if ( type == UNCOMPRESSED_FILE )
		{
			rewind( fp ) ;
		}
	}

	void SetCompressLevel( int cl )
	{
		if ( cl < 0 || cl > 9 )
		{
			fprintf( stderr, "Compress level must be 0-9.\n" ) ;
			exit( 1 ) ;
		}
		gzCompressLevel = cl ;
	}
} ;

#endif
