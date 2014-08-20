#!/usr/bin/env perl
#
# 
use strict;

my $cohort_filename;
my %sample_groups;

#
# parse SAMPLE/GROUP map, if available
#
if( $ARGV[0] eq "--groups" ) {
    $cohort_filename = @ARGV[1];
    pop @ARGV;
    pop @ARGV;

    open(GRP, "<", $cohort_filename) || die "ERROR: ${cohort_filename}: $!\n";
    my $linenum;
    while(<GRP>) {
	$linenum++;
	next if( $_ =~ m/^\#/ ); # skip comments; 
	my ($sample_id, $sample_group) = split(/[\t\n]/,$_);
	$sample_groups{$sample_id} = $sample_group;
    }
    close(GRP);
}

my $vcf_filename="STDIN";
my $first_format; # format for first observed variant
my %formats;  # hash of format abbreviations
my %infos; # hash of info abbreviations
my @chrom; # hash of #CHROM line (sample list)


 
# output header, part1
print join("\t","pos","score","qc","sample","group");

# parse input VCF
my $linenum=0;
while (my $line=<> ) {
    $linenum++;
    if( $linenum == 1 && $line !~ m/^##fileformat=VCFv4/ ) {
	print STDERR "ERROR: ${vcf_filename}: not a VCF file.\n";
	exit(1);
    }
    # capture format lines
    if( $line =~/^##FORMAT=<ID=([^,]+),.*Description="([^"]*)"/ ) {
	$formats{_index}++;
	$formats{$1}=(id=>$1,description=>$2,index=>$formats{_index}+0);
	next;
    }
    # capture info lines
    if( $line =~/^##INFO=<ID=([^,]+),.*Description="([^"]*)"/ ) {
	$infos{_index}++;
	$infos{$1}=(id=>$1,description=>$2,index=>$infos{_index}+0);
	next;
    }
    # capture CHROM line
    if( $line =~ /^#CHROM/ ) {
	@chrom = split(m/[\t\n]/,$line);
	next;
    }
    # skip other headers
    if( $line =~ /^#/ ) {
	next;
    }

    # PARSE variants
    my @variant = split(m/[\t\n]/, $line);
    my $v_chr =$variant[0];
    my $v_start =$variant[1];
    my $pos="${v_chr}:${v_start}";
    my $ref="$variant[3]";
    my @alleles=($ref, split(",",$variant[4]));
    my $score =$variant[5];
    my $qc =$variant[6];
    my $info=$variant[7]; 
    my $format=$variant[8];
    if( $first_format && $first_format ne $format ) {
	print "ERROR: ${vcf_filename}:${linenum}: format string not consistent: $format != $first_format\n";
	exit(1);
    } elsif(!$first_format) {
	$first_format=$format;
	# header, part2
	print "\t",join("\t", split(":",$format)), "\n";
    }
	
    # EXPORT samples per variant
    for(my $i=9; $i < scalar(@chrom); $i++) {
	my @v_info = split(/:/, $variant[$i]);

	# lookup sample grou
	my $group = $sample_groups{$chrom[$i]};
	# convert index genotype to real genotype
	my $genotype;
	if( $v_info[0] =~ m/([.0-9]+)(.)([.0-9]+)(.*)/ ) {
	    $genotype = join("", ($1 eq ".")?".":$alleles[$1], $2,  ($3 eq ".")?".":$alleles[$3], $4);
	    if( $4 ) {
		print STDERR "WARNING: ${vcf_filename}:${linenum}:${i}: genotype info has more than 2 alleles per sample: $chrom[$i] => $v_info[0]\n";
	    }		
	} else {
	    print STDERR "ERROR: ${vcf_filename}:${linenum}:${i}: can't parse genotype info: $v_info[0]\n";
	    exit 1;
	}
	print join("\t", 
		   #"[${linenum}:$i]",
		   # include these for TABIX indexing
		   $v_chr, $v_start, 
		   # these are for Excel Pivot Table
		   $pos,
		   $score,
		   $qc,
		   #@alleles, # debugging
		   # sample info
		   $chrom[$i],
		   $group,
		   # INFO stuff here
		   # variant data
		   $genotype,
		   @v_info,
		   ),
	"\n";
    }
}


