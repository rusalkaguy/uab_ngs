#!/usr/bin/perl
#
# Split compound columns into multiple columns
#
# example: CADD annotation with -otherinfo puts 2 CSV scores in the column: raw,scaled
# this splits that into two TSV columns, with appropriate headers.
#
use strict;
use Data::Dumper;

my %CHROM;
my $parsed_header = 0;
my $debug=0;

my %split_keys = (
    "caddgt10" => { sep => ",", out_keys => ["caddgt10_raw", "caddgt10"], present=>0 },
    "caddgt20" => { sep => ",", out_keys => ["caddgt20_raw", "caddgt20"], present=>0 },
    "cadd" => { sep => ",", out_keys => ["cadd_raw", "cadd"], present=>0 },
);
#
# data on STDIN
#
my $linenum=0;
while(<>) {
    $linenum++;
    if( !$parsed_header ) {
	#
	# parse the header line to get column names
	#
	print "[$linenum] HEADER: $_" if($debug);
	if( $_ =~ m/^#CHROM/ || $_ =~ m/^Chr\tStart\tEnd\tRef\tAlt/ ) {
	    my @keys = split(m/[\t\n\r]/,$_);
	    my @out_keys;
	    print "[$linenum] #CHROM header: ",scalar(@keys)," columns\n" if($debug);
	    my $index = 0;
	    foreach my $key ( @keys ) {
		$CHROM{$key} = $index; 
		print "\t#CHROM[$index]: $key -> $index\n" if($debug);
		$index++;
		# add column header for columns we will split
		if( $split_keys{$key} ) {
		    push @out_keys, @{$split_keys{$key}->{out_keys}};
		    $split_keys{$key}->{present} = 1;
		} else {
		    push @out_keys, $key;
		}
	    }
	    $parsed_header = 1;
	    print join("\t",@out_keys)."\n";

	} else {
	    # pass through other header lines
	    print $_;
	}
    } else {
	#
	# parse and fix data lines
	# 
	my @fields = split(m/[\t\n\r]/,$_);

	# strip first value in CSV in caddgt10 (-otherinfo): raw_score,phred_scaled_score
	for my $split_key ( keys %split_keys ) {
	    # if observered in header
	    if ( $split_keys{$split_key}->{present} ) {

		# if value is illegal - error out entirely
		if ( $fields[$CHROM{$split_key}] ne "NA" && $fields[$CHROM{$split_key}] !~ m/^[0-9.-]+,[0-9.]+$/ ) {
		    print "ERROR: line $linenum : $split_key (col=$CHROM{$split_key}]) isn't a 2-value CSV: '$fields[$CHROM{$split_key}]' \n";
		    exit(1);
		} else {
		    # split column 
		    if( $debug) { print "split [$split_key,$CHROM{$split_key}]: '$fields[$CHROM{$split_key}]' => "; }
		    if( $fields[$CHROM{$split_key}] ne "NA" ) {
			# split CSV
			$fields[$CHROM{$split_key}] = join("\t",split($split_keys{$split_key}->{sep},$fields[$CHROM{$split_key}]));
		    } else {
			# replicate NA value
			$fields[$CHROM{$split_key}] = "NA\tNA";
		    }
		    if( $debug) { print "'$fields[$CHROM{$split_key}]'\n"; }
		}
	    }
	}
	print "[$linenum] DATA[",scalar(@fields),"]: "  if($debug);
	print join("\t",@fields),"\n";
    }
    
}

