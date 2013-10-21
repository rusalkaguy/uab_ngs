#!/usr/bin/perl
# 
# input must be sorted read name, then bitscore (desc)
#
# sort -k1,1 -k12,12nr
#
my $min_len = 1;
my $ln = 50;
my $cur_read;
my $cur_score;
my $read_stats = {}; 
my $genotype = {};

while(my $line=<>) {
    # count lines
    $ln++;
    # parse blast 2-column out
    my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split("\t", $line);
    # min length constraint
    if( $length < $min_len ) {
	# skip this read
	next;
    }
    # parse genotype
    my ($gene,$subtype) = split(m/\./, $sseqid);
    #print "[$ln] $line";
    #print "[$ln] $qseqid, $sseqid, $gene, $subtype\n";
    # already in this read
    if( $cur_read eq $qseqid ) {
	# worse score, ignore
	if($bitscore < $cur_score) {
	    # skip this read
	    next;
	}
	# accounting for this read
	$read_stats->{$gene}->{$subtype}++;
    } else {
	# process prev read info
	if( $cur_read ) {
	    # for each gene
	    foreach my $g ( keys %{$read_stats} ) {
		# for each genotype
		my $composite_subtype = join(":", sort( keys %{$read_stats->{$g}} ) );
		$genotype->{$g}->{$composite_subtype}->{readcount}++;
		$genotype->{$g}->{$composite_subtype}->{bitscore}+=$cur_score;
	    }
	}	    
	# move to new read
	$cur_read = $qseqid;
	$cur_score = $bitscore;
	$read_stats = {};
    }
}
# output
foreach my $g ( sort keys %{$genotype} ) {
    foreach my $cs ( sort keys %{$genotype->{$g}} ) {
	print join("\t", 
		   "$g.$cs", # gene name & composit subtype
		   $genotype->{$g}->{$cs}->{readcount}, 
		   $genotype->{$g}->{$cs}->{bitscore},
		   "\n" );
    }
}
		   

    
