#!/usr/bin/perl
#
# custom multi-column filtering that ANNOVAR does not support
#
use strict;
use Data::Dumper;

# 
# argument parsing (basic)
#
my $debug=0;
if( $ARGV[0] eq "-v" ) {
    $debug = 1;
    shift @ARGV;
}
my $col_name  =shift @ARGV;
my $col_index = "NA";
my $col_func  =shift @ARGV;
my $col_value =shift @ARGV;

my %CHROM;
my $parsed_header = 0;

# stats
my %counts = (
    linenum => 0,
    header_lines => 0,
    data_lines => 0,
    data_dropped => 0,
    data_kept => 0,
);

#
# data on STDIN
#
while(<>) {
    $counts{linenum}++;
    if( !$parsed_header ) {
	#
	# parse the header line to get column names
	#
	$counts{header_lines}++;
	print "[$counts{linenum}] HEADER: $_" if($debug);
	if( $_ =~ m/^#CHROM/ || $_ =~ m/^Chr\tStart\tEnd\tRef\tAlt/ ) {
	    print "[$counts{linenum}] #CHROM header\n" if($debug);
	    my $index = 0;
	    my @header_keys = split(m/[\t\n\r]/,$_);
	    foreach my $key ( @header_keys ) {
		$CHROM{$key} = $index; 
		print "\t#CHROM[$index]: $key -> $index\n" if($debug);
		$index++;
	    }
	    $parsed_header = 1;
	    print join("\t", @header_keys), "\n";

	    # check validity of request
	    $col_index = $CHROM{$col_name};
	    if( ! defined($col_index) ) {
		print STDERR "ERROR: filter column '$col_name' not found in ",join(",",@header_keys),"\n";
		exit(1);
	    }
	    print "col_index=",$col_index+1,"\n" if( $debug );
	    next;
	}
	print $_; # print non-column header comment lines
    } else {
	$counts{data_lines}++;
	print "[$counts{linenum}]DATA $_\n" if($debug);

	#
	# parse and FILTER data lines
	# 
	my @fields = split(m/[\t\n\r]/,$_);
	my $keep = 0;

	# string functions
	$keep=1 if( "$col_func" eq "eq" && $fields[$col_index] eq $col_value );
	$keep=1 if( "$col_func" eq "ne" && $fields[$col_index] ne $col_value );

	# numeric functions
	if( "NA" ne $fields[$col_index] ) {
	    # have a number
	    print "NUMTEST $fields[$col_index] $col_func $col_value\n" if( $debug );
	    if( "$col_func" eq "==" ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if( $fields[$col_index] == $col_value );
	    }
	    if( "$col_func" eq ">" ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if( $fields[$col_index] > $col_value );
	    }
	    if( "$col_func" eq ">="  ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if( $fields[$col_index] >= $col_value );
	    }
	    if( "$col_func" eq "<"   ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if(  $fields[$col_index] < $col_value );
	    }
	    if( "$col_func" eq "<="  ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if( $fields[$col_index] <= $col_value );
	    }
	    if( "$col_func" eq "!="  ) {
		print "TEST $fields[$col_index] $col_func $col_value\n" if( $debug );
		$keep=1 if( $fields[$col_index] != $col_value );
	    }
	}	
	print "[$counts{linenum}] keep=$keep ($fields[$col_index] $col_func $col_value)\n" if( $debug );

	# actually filter
	if( !$keep ) {
	    $counts{data_dropped}++;
	    next;
	} else {
	    $counts{data_kept}++; 
	    print join("\t",@fields),"\n";
	}
    }
}

# summary at end

print STDERR "# stats ${col_name}[",($col_index+1),"] $col_func $col_value\n";
foreach my $k ( sort keys %counts ) {
    print STDERR "#         $k=$counts{$k}\n";
}
