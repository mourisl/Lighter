#!/bin/perl

my $line ;

while (<>)
{
	$line = $_ ;
	print $line ;
	
	# trim the last 2 base
	$line = <> ;
	chomp $line ;
	substr( $line, -2, 2, "\n" ) ;
	print $line ;
}
