#!/usr/bin/env perl 
#
# Multi-file join based on gene symbol
# 
use strict;
use Data::Dumper;
Data::Dumper::Terse(1);

my %genes; #=>  { name, chr, min, max, values => { filename=>{line,value} }
my @header = ("#CHROM","start","end","gene");
my $missing = 0;

# zero-based column indices
my $col_chr  =0;
my $col_start=1;
my $col_end  =$col_start;
my $col_gene =2;
my $col_value=3;

# 
# argument parsing (basic)
#
my $debug=0;
if( $ARGV[0] eq "-v" ) {
    $debug = 1;
    shift @ARGV;
}
my @filenames = @ARGV;
print "# FILES: ",join(",",@filenames),"\n" if($debug);

foreach my $filename ( @filenames ) {
    open(INFILE, "<", $filename ) || die "$filename: $!\n";
    print "# OPENED($filename)\n" if ($debug);
    my $linenum = 0;
    while( <INFILE> ) {
	# parse and count
	$linenum++;
	my @cols = split(/[ \t\n]+/,$_);

	# handle header lines
	if( $_ =~ m/^#CHROM/ ) {
	    # already have a header, add, the new sample name
      	    push @header, $cols[$col_value];
	    print "# [${filename}:$linenum] SAMPLE_NAME $cols[$col_value]\n" if($debug);
	    next;
	} elsif( $_ =~ m/^#/ ) {
	    # some other header - ignore
	    print "# [${filename}:$linenum] SKIP HEADERE $_\n" if($debug);
	    next; 
	} else {
	    print "# [${filename}:$linenum] COLS: ", join("|",@cols),"\n" if($debug);
	    # look up gene record
	    my $gene = $genes{$cols[$col_gene]};
	    if( ! $gene ) {
		# create if does not exist
		$gene = { 
		    name  => $cols[$col_gene],
		    chr   => $cols[$col_chr],
		    start => $cols[$col_start],
		    end   => $cols[$col_end],
		};
		$genes{$cols[$col_gene]} = $gene;
		print "# NEW GENE = ", DumperOne($gene),"\n" if($debug);
	    }
	    # update gene record
	    $gene->{start} = $cols[$col_start] if( $cols[$col_start] < $gene->{start});
	    $gene->{end}   = $cols[$col_end]   if( $cols[$col_end]   > $gene->{end});
	    # add or create value record
	    if( $gene->{values}->{$filename} ) {
		# add values (hope they're numbers
		$gene->{values}->{$filename}->{value} += $cols[$col_value];
	    } else {
		# create value record
		$gene->{values}->{$filename} = {
		    linenum => $linenum,
		    value => $cols[$col_value],
		};
	    }
	    print "# UPDATED GENE = ", DumperOne($gene),"\n" if($debug);
	}
    }
    close(INFILE);
    print "# CLOSED($filename)\n" if ($debug);
}

#
# output the mess
#

# header
print join("\t", @header ), "\n";

# unsorted gene list
foreach my $gene ( values %genes ) {

    # gene info
    print join("\t", $gene->{chr}, $gene->{start}, $gene->{end}, $gene->{name} );

    # sample info
    foreach my $filename ( @filenames ) {
	print "\t",
	defined($gene->{values}->{$filename}->{value})?
	    $gene->{values}->{$filename}->{value}:$missing
	    ;
	    
    }
    # end of line/record
    print "\n";
}


##
# DumperOne
##
sub DumperOne
{
    my $str = Dumper($_[0]);
    $str =~ s/[\t\n]/ /g;
    $str =~ s/ [ ]+/ /g;
    return($str);
}

    
