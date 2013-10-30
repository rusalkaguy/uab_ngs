#!/usr/bin/env perl
use strict;
use File::Basename; 
use Data::Dumper;

my $dest_genome = $ARGV[0];
my $in_fname = $ARGV[1];
my $debug = 0;

# build simple chr name map
my %genome_map = (
    "hg19" => { src_assembly => "b37", "MT" => "chrM" },
    "b37" => { src_assembly => "hg19", "chrM" => "MT" },
);
foreach my $chr ( 1..23, "X", "Y" ) {
    $genome_map{hg19}->{$chr}= "chr".$chr;
    $genome_map{b37}->{"chr".$chr}= $chr;
}

# construct out filenames
$in_fname = "in_fname.vcf" if( ! $in_fname );
my ($in_name,$in_path, $in_ext) = fileparse($in_fname, (".vcf" ));
my $out_fname = join("", $in_path, $in_name, ".", $dest_genome,"", $in_ext);
my $error_fname = join("", $in_path, $in_name, ".", $dest_genome, "_unmapped", $in_ext);

# check conditions
if ( ! $genome_map{$dest_genome} ) {
    print "SYNTAX: ".basename($0)." dest_genome in_fname\n";
    print "dest_genome = [".join(",", keys %genome_map )."]\n";
    print "out_fname = $out_fname\n";
    print "error_fname = $error_fname\n";
    exit(1);
}

open( SRC, "<", $in_fname ) || die "$in_fname: $!\n";
open( DEST, ">", $out_fname )|| die "$out_fname: $!\n";
open( UNMAPPED, ">", $error_fname )|| die "$error_fname: $!\n";

my $error_header = "";
my $line_num = 0;
my $first_contig_header = 1;
my%stats = ( 
	     contig_mapped => "", contig_unmapped => "", 
	     variant_mapped => "", variant_unmapped => "",
	     variant_unparsable => "",
	     );

while( my $line = <SRC> ) {
    $line_num++;
    # check first line for fileformat=VCFv4.1
    if( $line_num == 1 ) {
	if( $line !~ m/^..fileformat=VCFv4.1/ ) {
	    die "$in_fname: not a VFC v4.1\n";
	}
    }
    if( $line =~ m/^#/ ) {

	#
	# parse headers
	#
	if( $line =~ m/(^..contig=<ID=)([^,]+)(,.*)assembly=([a-z0-9._-]+)(.*)/i ) {
	    #
            # contig header
	    # 
	    my ($contig_left,$src_contig,$contig_len,$src_assembly,$contig_right)=($1,$2,$3,$4,$5);

	    #  output our commandline 
	    if( $first_contig_header ) {
		# put outselves on the header
		print DEST "##MapFastAssemblyCmd=\"".basename($0)." $ARGV[0] $ARGV[1]\"\n";
		$first_contig_header = 0;
	    }

	    # check source assmebly
	    if( $src_assembly != $dest_genome ) {
		die "${in_fname}:$line_num: assembly=$src_assembly, not src_assembly=$genome_map{$dest_genome}->{src_assembly}\n";
	    }
	    # look for mapped contig name
	    my $dest_contig = $genome_map{$dest_genome}->{$src_contig};
	    if( $dest_contig ) {
		# map the contig
		print DEST join("", $contig_left,$dest_contig,$contig_len,"assembly=",$dest_genome, $contig_right), "\n";
		print "[$line_num] mapped contig header $src_contig -> $dest_contig\n" if( $debug );
		$stats{contig_mapped}++;
	    } else {
		# reject the contig
		$error_header .= $line;
		print "[$line_num] UNmapped contig header $src_contig\n" if( $debug );
		$stats{contig_unmapped}++;
	    }
	} else {
	    #
	    # other header 
	    #
	    print DEST $line;
	    $error_header .= $line;
	    print "[$line_num] header\n" if( $debug );
	}
    } else {
	# parse variant line
	if( $line =~ m/(^[a-z0-9_.-]+)(\t.*)/i ) {
	    # valid line
	    my( $src_contig, $variant_right ) = ($1, $2);
	    my $dest_contig = $genome_map{$dest_genome}->{$src_contig};
	    if( $dest_contig ) {
		# mappable variant
		print DEST $dest_contig,$variant_right, "\n";
		print "[$line_num] mapped variant $src_contig -> $dest_contig\n" if( $debug );
		$stats{variant_mapped}++;
	    } else {
		# un-mappable
		if( $error_header ) {
		    print UNMAPPED $error_header;
		    $error_header = "";
		}
		print UNMAPPED $line;
		$stats{variant_unmapped}++;
		print "[$line_num] UNmapped variant $src_contig\n" if( $debug );
	    }
	} else {
	    # unparsable
	    if( $error_header ) {
		print UNMAPPED $error_header;
		$error_header = "";
	    }
	    print UNMAPPED $line;
	    $stats{variant_unparsable}++;
	    print STDERR "${in_fname}:$line_num PARSE ERROR $line";
	    print "[$line_num] unparsable\n" if( $debug );
	}
    }
}
close( SRC );
close( DEST );
close( UNMAPPED );

# if didn't use unmapped output, delete the file
if( $error_header ) {
    unlink( $error_fname );
}

foreach my $key ( sort keys %stats ) {
    print "$stats{$key}\t$key\n";
}
