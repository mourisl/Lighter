#!/bin/perl

#ARGV[0]=0.077 when building 75x coverage

my $line ;

srand(17) ;

while ( <STDIN> )
{
	if ( rand() < $ARGV[0] )
	{
		print $_ ;
		$line = <STDIN> ;
		print $line ;
		$line = <STDIN> ;
		print $line ;
		$line = <STDIN> ;
		print $line ;
	}
	else
	{
		$line = <STDIN> ;
		$line = <STDIN> ;
		$line = <STDIN> ;
	}
}
