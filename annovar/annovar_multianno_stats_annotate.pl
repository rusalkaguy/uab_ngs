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
	my $id = $fields[$CHROM{"ID"}];
	my $snp137 = $fields[$CHROM{"snp137"}];
	my $snp138 = $fields[$CHROM{"snp138"}];
	$is_novel = 0 if( "NA" ne $snp137 || "NA" ne $snp138 || ("." ne $id && "" ne $id) );
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
	# add max_aaf as column
	$fields[$CHROM{"max_aaf"}] = $max_maf;
	$fields[$CHROM{"is_novel"}] = $is_novel;
	$fields[$CHROM{"is_snp"}] = $is_snp;

	# count at levels
	$counts{variant_maf05}++ if($max_maf <= 0.05);
	$counts{indel_maf05}++ if($max_maf <= 0.05 && !$is_snp);
	$counts{variant_maf03}++ if($max_maf <= 0.03);
	$counts{indel_maf03}++ if($max_maf <= 0.03 && !$is_snp);
	$counts{variant_maf01}++ if($max_maf <= 0.01);
	$counts{indel_maf01}++ if($max_maf <= 0.01 && !$is_snp);
	# actually filter novel
	if( !$is_novel)	{ 
	    $counts{variant_known}++; 
	    $counts{indel_known}++ if( ! $is_snp ); 
	    if( $req_novel ) {
		print "[$linenum] DROP DATA[#col=",scalar(@fields),",is_snp=?, is_novel=$is_novel : $snp137,$max_maf\n"  if($debug);
		next;
	    }
	} else {
	    $counts{variant_novel}++ ;
	    $counts{indel_novel}++ if( !$is_snp ); 
	};

	#
	# SNP scoring
	# see http://www.openbioinformatics.org/annovar/annovar_filter.html#ljb23
	#

	# CADD score (scaled)
	my $cadd_phred=$fields[$CHROM{cadd}];
	if( "NA" ne $fields[$CHROM{cadd}] ) {
	    $cadd_obs_min = $cadd_phred if( $cadd_phred < $cadd_obs_min);
	    $cadd_obs_max = $cadd_phred if( $cadd_phred > $cadd_obs_max);
	} else {
	    $cadd_phred = 0;  # replace missing data with sentinal value 0
	}

	# PolyPhen 
	my $pphenHdiv2=$fields[$CHROM{LJB2_PolyPhen2_HDIV}];
	
	# SIFT
	my $sift2=$fields[$CHROM{LJB2_SIFT}];
	    

	# eleminate if it is a SNP and fails all three impact scores. 
	if( $is_snp ) {
	    my $pphen_hit = ($pphenHdiv2 > $pphenHdiv2_min)?1:0;
	    my $sift_hit = ($sift2 <= $sift2_maxeq)?1:0;
	    my $cadd_hit = ($cadd_phred > $cadd_min)?1:0;
	    my $hit_key = "$cadd_hit$pphen_hit$sift_hit";

	    $counts{cadd_hit}++ if( $cadd_hit);
	    $counts{sift_hit}++ if( $sift_hit);
	    $counts{pphen_hit}++ if( $pphen_hit);
	    $counts{sift_and_pphen_hit}++ if( $pphen_hit && $sift_hit);
	    $counts{sift_or_pphen_hit}++ if( $pphen_hit || $sift_hit);
	    $counts{cadd_sift_and_pphen_hit}++ if( $cadd_hit && $pphen_hit && $sift_hit);
	    $counts{cadd_sift_or_pphen_hit}++ if( $cadd_hit || $pphen_hit || $sift_hit);

	    $and_counts{$hit_key}++;
	    $or_counts{$hit_key}++ if( $cadd_hit || $pphen_hit || $sift_hit);
	    #print "DROP KEY $hit_key $and_counts{$hit_key} $or_counts{$hit_key}\n";
#	    if( ! ($pphen_hit || $sift_hit || $cadd_hit) ) {
#		# reject snp as un-interesting
#		print "[$linenum] DROP DATA[#col=",scalar(@fields),",is_snp=$is_snp]\n"  if($debug);
#		next;
#	    } else {
		# pass as possibly interesting 
		$counts{snp_pass}++;
#	    }
	}
	# track cadd scores for pphen/sift hits
	if( "NA" ne $fields[$CHROM{cadd}] ) {
	    $cadd_pass_min = $cadd_phred if( $cadd_phred < $cadd_pass_min);
	    $cadd_pass_max = $cadd_phred if( $cadd_phred > $cadd_pass_max);
	}

	print "[$linenum] PASS DATA[#col=",scalar(@fields),",is_snp=$is_snp]\n"  if($debug);
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
foreach my $k ( sort keys %counts ) {
    print STDERR "#         $k=$counts{$k}\n";
}
#print STDERR "AND=",Dumper(\%and_counts),"\n";
#print STDERR "OR=",Dumper(\%or_counts),"\n";
print STDERR join("\t", "#", "AND",  "cadd", "pphen", "sift", "pp+sift" ), "\n";
print STDERR join("\t", "#", "tot",   $counts{cadd_hit},$counts{pphen_hit},$counts{sift_hit},$counts{sift_and_pphen_hit}), "\n";
print STDERR join("\t", "#", "cadd",  $and_counts{"100"},$and_counts{"110"},$and_counts{"101"},$and_counts{"111"}), "\n";
print STDERR join("\t", "#", "pphen", $and_counts{"110"},$and_counts{"010"},$and_counts{"011"},$and_counts{"011"}), "\n";
print STDERR join("\t", "#", "sift",  $and_counts{"101"},$and_counts{"011"},$and_counts{"001"},$and_counts{"011"}), "\n";
print STDERR join("\t", "#", "OR", "cadd|pphen", "cadd|sift", "pphen|sift", "cadd|pphen|sift" ), "\n";
print STDERR join("\t", "#", "tot",
		  $or_counts{"110"}+$or_counts{"010"}+$or_counts{"100"},
		  $or_counts{"101"}+$or_counts{"001"}+$or_counts{"100"},
		  $or_counts{"011"}+$or_counts{"001"}+$or_counts{"010"},
		  $or_counts{"111"}+$or_counts{"001"}+$or_counts{"010"}+$or_counts{"100"}+$or_counts{"110"}+$or_counts{"101"}+$or_counts{"011"}
    ), "\n";
#print STDERR "# counts=", Dumper(\%counts),"\n";

print STDERR "#\n" ;
    print STDERR "# extra stats\n";
print STDERR "# cadd_obs_min=$cadd_obs_min\n";
print STDERR "# cadd_obs_max=$cadd_obs_max\n";
print STDERR "# cadd_pass_min=$cadd_pass_min\n";
print STDERR "# cadd_pass_max=$cadd_pass_max\n";
