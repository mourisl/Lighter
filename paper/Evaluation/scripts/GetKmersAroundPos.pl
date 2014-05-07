#!/bin/perl

#usage: a.pl kmer_length xxx.ref1 yyy.ref2 pos_file

use strict ;

my $k = $ARGV[0] ;
my $i ;
my $j ;
open FP1, $ARGV[1] ;
open FP2, $ARGV[2] ;
open FP3, $ARGV[3] ;

my %kmers ;
my @mutates ;
my $ref ;
my $key ;

while ( <FP3> )
{
	chomp ;
	push @mutates, $_ ;
}

$ref = <FP1> ;
$ref = "" ;
while ( <FP1> )
{
	chomp ;
	$ref = $ref.$_ ;
}

for ( $i = 0 ; $i < scalar( @mutates ) ; ++$i )
{
	my $pos = $mutates[$i] ;
	for ( $j = $pos - $k + 1 ; $j <= $pos ; ++$j )
	{
		next if ( $j < 0 || $j + $k - 1>= length( $ref ) ) ;
		$key = substr( $ref, $j, $k ) ;
		$kmers{ $key } = 1 ;
	}
}

$ref = <FP2> ;
$ref = "" ;
while ( <FP2> )
{
	chomp ;
	$ref = $ref.$_ ;
}

for ( $i = 0 ; $i < scalar( @mutates ) ; ++$i )
{
	my $pos = $mutates[$i] ;
	for ( $j = $pos - $k + 1 ; $j <= $pos ; ++$j )
	{
		next if ( $j < 0 || $j + $k - 1>= length( $ref ) ) ;
		$key = substr( $ref, $j, $k ) ;
		$kmers{ $key } = 1 ;
	}
}

foreach $key ( keys %kmers )
{
	print $key, "\n" ;
}
