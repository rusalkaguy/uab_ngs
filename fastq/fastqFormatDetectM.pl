#!/usr/bin/perl
# 
#
# Original Author: Martin Dahlo
# http://www.uppmax.uu.se/userscript/check-fastq-quality-score-format
#
# 2013-04-23 Modified by Curtis Hendrickson @ UAB
#   * handle a list of files
#   * -q only prints 32/64 to make more script friendly
#   * handle .gz files
#   * change lower threshold to > 74, so that Illumina 1.8+ come out as Sanger
#
# Usage:  perl scriptname.pl -q <infile> [<infiles]
# ex.
#      perl scriptname.pl reads.fq

use warnings;
use strict;


=pod

Used to detect the format of a fastq file. In its current state,
it can only differentiate between sanger and solexa/illumina.
If need arises, checking for different versions of illumina formats
could easily be implemented. ( Please upload an update if you implement this )

Can easily be copy/pasted into any other script and altered to do other
things than die when it has determined the format.

Pseudo code

* Open the fastq file
* Look at each quality ASCII char and convert it to a number
* Depending on if that number is above or below certain thresholds,
  determine the format.


=cut


# get variables
my $usage = "Usage:  perl scriptname.pl [-q] <infile> [<infile>] \n\t-q prints eitehr 64(Solexa/Illumina1.3) or 33(Sanger/Illumina1.8)\n";
my $quiet = 0;
if( $ARGV[0] eq "-q" ) {
    $quiet = 1;
    shift;
}

if( scalar(@ARGV) < 1 ) { die $usage; }
foreach my $fq (@ARGV) {

    # open the files
    # open FQ, "<", $fq or die $!;
    my $rc = 0;
    if( $fq =~ m/\.gz$/ ) {
	$rc = open( FQ, "zcat $fq |");
    } else {
	$rc = open( FQ, "<", $fq );
    }
    if( ! $rc ){ 
	print "$fq\t$!\n"; 
	next ; 
    }


    # initiate
    my @line;
    my $l;
    my $number;


FILE: 
    # go thorugh the file
    while(<FQ>){

	# if it is the line before the quality line
	if($_ =~ /^\+/){

	    $l = <FQ>;		# get the quality line
	    @line = split(//,$l); # divide in chars
	    for(my $i = 0; $i <= $#line; $i++){	# for each char
		$number = ord($line[$i]); # get the number represented by the ascii char

		# check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
		if($number > 74){ # if solexa/illumina
		    #die "This file is solexa/illumina format\n"; # print result to terminal and die
		    if( $quiet ) {
			print "64\n";
		    } else {
			print "$fq\t64\tIllumina 1.3+/Solexa (Phred/Solexa+64)\n";
		    }
		    last FILE;
		}elsif($number < 59){ # if sanger
                    #die "This file is sanger format\n";	# print result to terminal and die
		    if( $quiet ) {
			print "32\n";
		    } else {
			print "$fq\t32\tSanger/Illumina 1.8+ (Phred+33)\n";
		    }
		    last FILE;
		}
	    } # next QUAL char 
	} # if pre-QUAL line
    } # next line
    close(FQ);
} # next ARG (file)
