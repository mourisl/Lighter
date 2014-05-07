#!/bin/perl

my $ref ;
my $i ;
my $line ;
open FP2, ">mutate.pos" ;
print ">ref2\n" ;
my $pos = 0 ;
my @alphabet ;
$line = <> ;
push @alphabet, "A" ;
push @alphabet, "C" ;
push @alphabet, "G" ;
push @alphabet, "T" ;
srand( 17 ) ;
while ( <> )
{
	chomp ;
	$line = $_ ;
	my $c, $nc ;
	for ( $i = 0 ; $i < length( $line ) ; ++$i )
	{
		if ( rand() < 0.001 )
		{
			$c = substr( $line, $i, 1 ) ;
			do
			{
				$nc = $alphabet[ int(rand() * 4) ] ;
				#print $nc, "\n" ;
			} while ( $nc eq $c ) ;
			substr( $line, $i, 1, $nc ) ;

			print FP2 $pos + $i, "\n" ;
		}
	}
	print $line, "\n" ;
	$pos += length( $line ) ;
}

