#!/usr/bin/perl
#
# NOT YET TESTED
#
# custom multi-column filtering that ANNOVAR does not support
#
# (evolved out of annovar_multianno_filter.pl, with stats removed, just filtering now)
#
use strict;
use Data::Dumper;

my %counts = (
	      variants => 0,
	      );

my %CHROM;
my $parsed_header = 0;
my $debug=0;

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
	my $is_indel = 0;
	my $ref_seq = $fields[$CHROM{Ref}]; $ref_seq =~ s/-//;
	my $alt_seq = $fields[$CHROM{Alt}]; $alt_seq =~ s/-//;
	if( length($ref_seq) != 1 && length($alt_seq) != 1 ) {
	    $is_indel = 1;
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
	    my $maf= $fields[$CHROM{$pop}];
	    if( "NA" ne $maf ) { 
		$is_novel=0; 
		if("NA" eq $max_maf || $max_maf < $maf ) { 
		    $max_maf=$maf; 
		}
	    }
	}

	# actually filter novel
	if( !$is_novel)	{ 
	    $counts{variant_known}++; 
	    $counts{indel_known}++ if( $is_indel); 
	    if( $req_novel ) {
		print "[$linenum] DROP DATA[#col=",scalar(@fields),",is_indel=?, is_novel=$is_novel : $snp137,$max_maf\n"  if($debug);
		next;
	    }
	} else {
	    $counts{variant_novel}++ ;
	    $counts{indel_novel}++ if( $is_indel); 
	};

	#
	# SNP scoring
	# see http://www.openbioinformatics.org/annovar/annovar_filter.html#ljb23
	#

	# CADD score (scaled)
	my $cadd_phred=$fields[$CHROM{cadd}];
	if( "NA" eq $fields[$CHROM{cadd}] ) {
	    $cadd_phred = 0;  # replace missing data with sentinal value 0
	}

	# PolyPhen 
	my $pphenHdiv2=$fields[$CHROM{LJB2_PolyPhen2_HDIV}];
	if( "NA" eq $fields[$CHROM{LJB2_PolyPhen2_HDIV}] ) {
	    $pphenHdiv2 = 0;  # replace missing data with sentinal value 0
	}
	
	# SIFT
	my $sift2=$fields[$CHROM{LJB2_SIFT}];
	if( "NA" eq $fields[$CHROM{LJB2_PolyPhen2_HDIV}] ) {
	    $pphenHdiv2 = 0;  # replace missing data with sentinal value 0
	}
	    

	# eleminate if it is a SNP and fails all three impact scores. 
	if( !$is_indel ) {
	    my $pphen_hit = ($pphenHdiv2 > $pphenHdiv2_min)?1:0;
	    my $sift_hit = ($sift2 <= $sift2_maxeq)?1:0;
	    my $cadd_hit = ($cadd_phred > $cadd_min)?1:0;

	    #print "DROP KEY $hit_key $and_counts{$hit_key} $or_counts{$hit_key}\n";
	    if( ! ($pphen_hit || $sift_hit || $cadd_hit) ) {
		# reject snp as un-interesting
		print "[$linenum] DROP DATA[#col=",scalar(@fields),",is_indel=$is_indel]\n"  if($debug);
		next;
	    } else {
		# pass as possibly interesting 
		$counts{snp_pass}++;
	    }
	}
	print "[$linenum] PASS DATA[#col=",scalar(@fields),",is_indel=$is_indel]\n"  if($debug);
	print join("\t",@fields),"\n";
    }
    
}

# summary at end

print STDERR "# stats\n";
print STDERR "# variant_count=$counts{variant}\n";
print STDERR "#         sift2_maxeq=$sift2_maxeq\n";
print STDERR "#         pphenHdiv2_min=$pphenHdiv2_min\n";
print STDERR "#         cadd_min=$cadd_min\n";
print STDERR "#     indel_count=$counts{indel}\n";
print STDERR "#     snp_count=$counts{snp}\n";
print STDERR "#     snp_pass=$counts{snp_pass}\n";
