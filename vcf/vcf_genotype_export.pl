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

#
# parse ANNOVAR output for Gene.refGene (WARNING: stores array in memory!)
# please use SnpEff instead!
#
my @annovar_variant_gene;
if( $ARGV[0] eq "--annovar_gene" ) {
    my $annovar_filename = @ARGV[1];
    pop @ARGV;
    pop @ARGV;

    open(GNAME, "<", $annovar_filename) || die "ERROR: ${annovar_filename}: $!\n";
    #print STDERR "reading $annovar_filename...\n";
    my $annovar_linenum;
    my $annovar_variant_num; 
    my $col_index; # for Gene.refGene
    while(<GNAME>) {
	$annovar_linenum++;
	next if( $_ =~ m/^\#\#/ ); # skip comments; 
	if( $_ =~ m/^#CHROM/ ) {
	    my @col_names = split(/[\t\n]/,$_);
	    # get index for our column name: http://www.perlmonks.org/?node_id=75660
	    # print STDERR "\t",$_,"=>",$col_names[$_],"\n";
	    my @gr = grep { $col_names[$_] =~ m/Gene.refGene/i } 0..$#col_names;
	    $col_index = $gr[0];
	    #print STDERR "\tcol_index=$col_index\n";
	    next;
	}
	my @var_cols = split(/[\t\n]/,$_);
	$annovar_variant_num++;
	my $gene_name = $var_cols[$col_index];
	# clean up 
	$gene_name =~ s/[\(;].*//;  
	$annovar_variant_gene[$annovar_variant_num] = $gene_name;
	#if( $annovar_variant_num < 10 ) {print STDERR "\t annovar[$annovar_variant_num] $var_cols[0] $var_cols[1] $var_cols[$col_index] (idx=$col_index)\n";}
    
    }
    close(GNAME);
}

my $vcf_filename="STDIN";
my $first_format; # format for first observed variant
my %formats;  # hash of format abbreviations
my %infos; # hash of info abbreviations
my @chrom; # hash of #CHROM line (sample list)


 
# output header, part1
print join("\t","#CHROM", "start", "pos","id","score","qc","gene","sample","group");

# parse input VCF
my $linenum=0;
my $varnum=0;
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
    $varnum++;
    my @variant = split(m/[\t\n]/, $line);
    my $v_chr =$variant[0];
    my $v_start =$variant[1];
    my $pos="${v_chr}:${v_start}";
    my $v_id =$variant[2];
    my $v_ref=$variant[3];
    my @alleles=($v_ref, split(",",$variant[4]));
    my $v_score =$variant[5];
    my $v_qc =$variant[6];
    my $v_gene = $annovar_variant_gene[$varnum];
    my $v_info=$variant[7]; 
    my $v_format=$variant[8];


    # handle info fields

    # check formatting of variants
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
	my @s_info = split(/:/, $variant[$i]);

	# lookup sample group
	my $group = $sample_groups{$chrom[$i]};
	# classify genotype
	my $genoclass = "other";
	my $genosubclass = "other";

	# convert index genotype to real genotype
	my $genotype;
	if( $s_info[0] =~ m/([.0-9]+)(.)([.0-9]+)(.*)/ ) {
	    my($gt1, $phased, $gt2, $gtx)=($1,$2,$3,$4);
	    
	    $genotype = join("", ($gt1 eq ".")?".":$alleles[$gt1], $phased,  ($gt2 eq ".")?".":$alleles[$gt2], $gtx);
	    if( $gtx ) {
		print STDERR "WARNING: ${vcf_filename}:${linenum}:${i}: genotype info has more than 2 alleles per sample: $chrom[$i] => $s_info[0]\n";
	    }		
	    # genotype classification
	    if( $gt1 eq "0" && $gt2 eq "0" ) { $genoclass="HomRef"; $genosubclass=$genoclass; }
	    if( $gt1 eq "." || $gt2 eq "." ) { $genoclass="NoData"; $genosubclass=$genoclass; }
	    $genosubclass="HomAlt" if( $gt1 eq $gt2 && $gt1 ne "." && $gt1 ne "0" );
	    $genosubclass="HetRef" if( $gt1 ne $gt2 && ($gt1 eq "0" || $gt2 eq "0") && ($gt1 ne "." && $gt2 ne "."));
	    $genosubclass="HetAlt" if( $gt1 ne $gt2 && ($gt1 ne "0" && $gt2 ne "0") && ($gt1 ne "." && $gt2 ne "."));
	    #print STDERR "[${linenum}:${i}] $gt1$phased$gt2 = $genoclass, $genosubclass\n";
	} else {
	    print STDERR "ERROR: ${vcf_filename}:${linenum}:${i}: can't parse genotype info: $s_info[0]\n";
	    exit 1;
	}

	# replace index based GT
	$s_info[0] = $genotype;


	print join("\t", 
		   #"[${linenum}:$i]",
		   # include these for TABIX indexing
		   $v_chr, $v_start, 
		   # these are for Excel Pivot Table
		   $pos,
		   $v_id,
		   $v_score,
		   $v_qc,
		   $v_gene,
		   #@alleles, # debugging
		   # sample info
		   $chrom[$i],
		   $group,
		   # INFO stuff here
		   # variant data
		   #$genotype,
		   $genoclass,
		   $genosubclass,
		   @s_info,
		   ),
	"\n";
    }
}


