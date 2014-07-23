#!/usr/bin/perl
#
# extract names columns from ANNOVAR multi-anno output
# compute statistics on overlap between variant impact predictors
# add is_novel and max_aaf columns to right, based on existing column values. 
#
# (evolved out of annovar_multianno_filter.pl, with the filtering removed)
use strict;
use Data::Dumper;


my %CHROM;
my $parsed_header = 0;
my $debug=0;

# stats
my %counts = (
    variant => 0,
    snp => 0,
    indel => 0,
    snp_pass => 0,
    cadd_hit => 0,
    pphen_hit => 0,
    sift_hit => 0,
    sift_and_pphen_hit => 0,
    sift_or_pphen_hit => 0,
    cadd_and_sift_and_pphen_hit => 0,
);
my %and_counts;
my %or_counts;

my $cadd_obs_min = 50; 
my $cadd_obs_max= -50;
my $cadd_pass_min = 50; 
my $cadd_pass_max= -50;


#
# thresholds 
#
# 
# http://www.openbioinformatics.org/annovar/annovar_filter.html#ljb23
my $pphenHdiv2_min = 0.452;
my $sift2_maxeq = 0.05;
# http://cadd.gs.washington.edu/info
# we would suggest to put a cutoff somewhere between 10 and 20. 
# Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants
# 10 => 10% highest
# 15 => mean for splice sites/nonsynon 
# 20 => 1% highest
my $cadd_min = 10;
my $req_novel = 0;

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
	    print "[$linenum] #CHROM header\n" if($debug);
	    my $index = 0;
	    my @header_keys = split(m/[\t\n\r]/,$_);
	    push @header_keys, "max_aaf";
	    push @header_keys, "is_novel";
            push @header_keys, "is_snp";
	    foreach my $key ( @header_keys ) {
		$CHROM{$key} = $index; 
		print "\t#CHROM[$index]: $key -> $index\n" if($debug);
		$index++;
	    }
	    $parsed_header = 1;
	    print join("\t", @header_keys), "\n";
	    next;
	}
	print $_; # print non-column header comment lines
    } else {
	$counts{variant}++;
	#
	# parse and FILTER data lines
	# 
	my @fields = split(m/[\t\n\r]/,$_);

	# how do we tell if it's a INDEL?
	my $is_snp = 1;
	my $ref_seq = $fields[$CHROM{Ref}]; $ref_seq =~ s/-//;
	my $alt_seq = $fields[$CHROM{Alt}]; $alt_seq =~ s/-//;
	if( length($ref_seq) != 1 && length($alt_seq) != 1 ) {
	    $is_snp = 0;
	    $counts{indel}++;
	} else {
	    $counts{snp}++;
	}

	# compute novel - check snp137 & all reference populations
	my $is_novel = 1;
	my $max_maf = "0";
	# dbSNP
	my $snp137 = $fields[$CHROM{"snp137"}];
	my $snp138 = $fields[$CHROM{"snp138"}];
	$is_novel = 0 if( "NA" ne $snp137 || "NA" ne $snp138);
	# compute max observer allele freq
	foreach my $pop ( "1000g2012apr_all", "1000g2012apr_eur", "esp6500si_all", "esp6500si_ea", "cg46" ) {
	    my $maf= $fields[$
