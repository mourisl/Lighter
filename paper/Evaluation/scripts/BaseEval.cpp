#include <stdio.h>
#include <memory.h>

// Get the matched ratio for each position
// usage: a.out < xxx.sam

char line[10000] ;

//ERR022075.20	81	gi|49175990|ref|NC_000913.2|	4088567	42	100M	=	4088139	-528	ACAGCGGTTGTTGCTTTTGCTTTTCCGTTAACGCATGAGGCGCGACGGAACCAATCACACCAGGGATTTTCACTCCCTTGTGTGTGCGTATGGTTACCCN	###################################################################################################!	AS:i:-1	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:99G0	YS:i:-2	YT:Z:DP

int main()
{
	char id[1000], chr[100], cigar[1000], md[1000] ;
	int tag ;
	int start ;
	int qual ;
	char *p ;
	int cnt[1000] ;
	int totalCnt ;
	int i, j ;
	int num, pos ;

	memset( cnt, 0, sizeof( cnt ) ) ;
	while ( fgets( line, sizeof( line ), stdin ) )
	{
		sscanf( line, "%s %d %s %d %d %s\n", id, &tag, chr, &start, &qual, cigar ) ;
		if ( strcmp( cigar, "100M" ) )
			continue ;
		p = strstr( line, "MD:Z:" ) ;
		if ( p == NULL || ( tag & 0x40 ) == 0  )
			continue ;
		p += 5 ;
		sscanf( p, "%s", md ) ;
		
		pos = 0, num = 0 ;
		for ( i = 0 ; md[i] ; ++i )
		{
			if ( md[i] >= '0' && md[i] <= '9' )
			{
				num = num * 10 + md[i] - '0' ;
			}
			else
			{
				if ( num != 0 )
				{
					for ( j = pos ; j < pos + num ; ++j )
						++cnt[j] ;
					pos += num ;
					num = 0 ;
				}
				++pos ;
			}
		}
		if ( num != 0 )
		{
			for ( j = pos ; j < pos + num ; ++j )
				++cnt[j] ;
			pos += num ;
		}
		/*if ( pos != 100 )
		{
			printf( "hi\n" ) ;
		}*/
		++totalCnt ;
	}

	for ( i = 0 ; i < 100 ; ++i )
	{
		printf( "%lf\n", (double)cnt[i] / totalCnt * 100.0 ) ;
	}
	return 0 ;
}
