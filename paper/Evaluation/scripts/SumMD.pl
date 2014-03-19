#!/bin/perl
use strict ;

my $sumAS = 0 ;
my $sum = 0 ;
my @cols = 0 ;
my $i ;

while(<>)
{
	if ( /MD:Z:(.+?)\s/ )
	{
		#print $1, "\n" ;
		my $sum1 = 0 ;
		@cols = split /[A-Z\^]+/, $1 ;
		#print @cols, "\n" ;
		
		for ( $i = 0 ; $i < scalar(@cols) ; ++$i )
		{
			#print $cols[$i], " " ;
			$sum1 += $cols[$i] ;
		}
		$sumAS += $sum1 ;
		
		#print "\n" ;
		# then the penalty
		#$sum1 = 0 ;
		#for ( $i = 0 ; $i < length( $1 ) ; ++$i )
		#{
		#	++$sum1 if ( substr( $1, $i, 1 ) =~ /[A-Z]/ ) ;
		#}
		#print "$sum1\n" ;
		#$sumAS -= $sum1 ;
	}
	++$sum ;
}
print $sum, " ",$sumAS, "\n"
