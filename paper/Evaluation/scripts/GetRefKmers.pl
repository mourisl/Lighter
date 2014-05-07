#!/bin/perl

use strict ;

my $k = $ARGV[0] ;
my $ref ;
my $i ;
my %kmers ;
my $key ;

$ref = <STDIN> ;
$ref="" ;
while ( <STDIN> )
{
	chomp ;
	$ref = $ref.$_ ;
}

for ( $i = $k - 1 ; $i < length( $ref ) ; ++$i )
{
	$key = substr( $ref, $i - $k + 1, $k ) ;
	$kmers{ $key } = 1 ;
}

foreach $key (keys %kmers )
{
	print $key, "\n" ;
}
