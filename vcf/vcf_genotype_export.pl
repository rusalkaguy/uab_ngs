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
print join("\t","#CHROM", "start", "pos","id","score","qc","sample","group");

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
    my $v_id =$variant[2];
    my $v_ref=$variant[3];
    my @alleles=($v_ref, split(",",$variant[4]));
    my $v_score =$variant[5];
    my $v_qc =$variant[6];
    my $v_info=$variant[7]; 
    my $v_format=$variant[8];
    if( $first_format && $first_format ne $v_format ) {
	print "ERROR: ${vcf_filename}:${linenum}: format string not consistent: $v_format != $first_format\n";
	exit(1);
    } elsif(!$first_format) {
	$first_format=$v_format;
	# header, part2
	print "\t",join("\t","GTclass", "GTsubclass", split(":",$v_format)), "\n";
    }
	
    # EXPORT samples per variant
    for(my $i=9; $i < scalar(@chrom); $i++) {
	my @v_info = split(/:/, $variant[$i]);

	# lookup sample group
	my $group = $sample_groups{$chrom[$i]};
	# classify genotype
	my $genoclass = "other";
	my $genosubclass = "other";

	# convert index genotype to real genotype
	my $genotype;
	if( $v_info[0] =~ m/([.0-9]+)(.)([.0-9]+)(.*)/ ) {
	    my($gt1, $phased, $gt2, $gtx)=($1,$2,$3,$4);
	    
	    $genotype = join("", ($gt1 eq ".")?".":$alleles[$gt1], $phased,  ($gt2 eq ".")?".":$alleles[$gt2], $gtx);
	    if( $gtx ) {
		print STDERR "WARNING: ${vcf_filename}:${linenum}:${i}: genotype info has more than 2 alleles per sample: $chrom[$i] => $v_info[0]\n";
	    }		
	    # genotype classification
	    if( $gt1 eq "0" && $gt2 eq "0" ) { $genoclass="HomRef"; $genosubclass=$genoclass; }
	    if( $gt1 eq "." || $gt2 eq "." ) { $genoclass="NoData"; $genosubclass=$genoclass; }
	    $genosubclass="HomAlt" if( $gt1 eq $gt2 && $gt1 ne "." && $gt1 ne "0" );
	    $genosubclass="HetRef" if( $gt1 ne $gt2 && ($gt1 eq "0" || $gt2 eq "0") && ($gt1 ne "." && $gt2 ne "."));
	    $genosubclass="HetAlt" if( $gt1 ne $gt2 && ($gt1 ne "0" && $gt2 ne "0") && ($gt1 ne "." && $gt2 ne "."));
	    #print STDERR "[${linenum}:${i}] $gt1$phased$gt2 = $genoclass, $genosubclass\n";
	} else {
	    print STDERR "ERROR: ${vcf_filename}:${linenum}:${i}: can't parse genotype info: $v_info[0]\n";
	    exit 1;
	}

	# replace index based GT
	$v_info[0] = $genotype;


	print join("\t", 
		   #"[${linenum}:$i]",
		   # include these for TABIX indexing
		   $v_chr, $v_start, 
		   # these are for Excel Pivot Table
		   $pos,
		   $v_id,
		   $v_score,
		   $v_qc,
		   #@alleles, # debugging
		   # sample info
		   $chrom[$i],
		   $group,
		   # INFO stuff here
		   # variant data
		   #$genotype,
		   $genoclass,
		   $genosubclass,
		   @v_info,
		   ),
	"\n";
    }
}


