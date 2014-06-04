#!/usr/bin/perl
#
# Merge galaxy generated taxonomy reports across samples
#
use strict;
my $debug = 0;
my $verbose = 1;
my $min_perc = 0.02; # filters tax based on total reads

my $target = "order";
my %taxname = (
	       $target => "target_level", # level to extract
	       "root" => "root",  # number that had blast hits
	       "number_of_reads_to_sample" => "total", # total reads
	       "total" => "number_of_reads_to_sample", # total reads
);

print "ARGV ", join(",", @ARGV), "\n" if($debug);
my @sample_files = @ARGV;
print "SAMPLE_FILES ", join(",", @ARGV), "\n" if($debug);
if( ! scalar(@ARGV) ) {
    print STDERR "Syntax: script tax_summary1.txt [tax_summary2.txt....]\n";
    exit(1);
}

my @sample_names;
my %taxon; # key is taxon_name, then sample_name
foreach my $sample_file (@sample_files) {
    print "sample_file=$sample_file\n" if($debug);
    open(SIN, "<", $sample_file) || die "$sample_file: $!\n";
    print "OPENED $sample_file\n" if($debug);
    my $sample_name;
    my $count;
    while(<SIN>){
	my @col = split(m/[\n\t]/,$_);
	print "    SCAN ",join("|",@col),"\n" if($debug);
	# fist line should have sample name
	if ("sample_name" eq $col[0] ) { 
	    $sample_name = $col[2]; 
	    # hack to pull sample name from filename
	    $sample_name = $sample_file;
	    $sample_name =~ s/_(.*)__.*/$1/;
	    $sample_name =~ s/(.*)_v.*/$1/;
	    push @sample_names, $sample_name;
	    next; 
	}
	# handle numeric data
	if( $taxname{$col[0]} ) {
	    $count++;
	    my $taxon_name=$col[1];
	    my $read_count=$col[2];
	    print "$col[0] $taxon_name $read_count\n" if($debug);
	    $taxon{$taxon_name}->{read_count} += $read_count;
	    $taxon{$taxon_name}->{$sample_name}->{read_count} += $read_count;
	    if( $taxon{$taxname{total}}->{$sample_name}->{read_count} == 0  ){
		$taxon{$taxname{total}}->{$sample_name}->{read_count} = 5000;
	    }
	    $taxon{$taxon_name}->{$sample_name}->{read_perc} = $read_count/$taxon{$taxname{total}}->{$sample_name}->{read_count};
	    #print STDERR"  PERC [$sample_name|$taxon_name] = $read_count / $taxon{$taxname{total}}->{$sample_name}->{read_count} = $taxon{$taxon_name}->{$sample_name}->{read_perc}\n";
	}
    }
    close(SIN);
    print STDERR "$sample_file : processed $count ${target}s\n" if($verbose);
}

#
# output summary
#

# column headers
print join("\t", "total_read_count",$target);
foreach my $sample_name ( @sample_names ) {
    print "\t", $sample_name;
}
foreach my $sample_name ( @sample_names ) {
    print "\t", $sample_name ."_raw";
}
print "\n";

# data: rows=taxa, cols=samples
foreach my $taxon_name ( sort keys %taxon ) {
    #next if( $taxname{$taxon_name} ); # skip totals, admin stuff
    # filter out rare taxa
    if( $min_perc ) {
	next if( $taxon{$taxon_name}->{read_count}/5000 < $min_perc );
    }
	    
    # print taxon header colums
    print join("\t", $taxon{$taxon_name}->{read_count}, $taxon_name);

    # print per-sample PERCENT values
    foreach my $sample_name ( @sample_names ) {
	my $sample_value = 0;
	if( $taxon{$taxon_name}->{$sample_name} ) {
	    $sample_value = $taxon{$taxon_name}->{$sample_name}->{read_perc};
	}
	print "\t", $sample_value;
    }
    # print per-sample RAW values
    foreach my $sample_name ( @sample_names ) {
	my $sample_value = 0;
	if( $taxon{$taxon_name}->{$sample_name} ) {
	    $sample_value = $taxon{$taxon_name}->{$sample_name}->{read_count};
	}
	print "\t", $sample_value;
    }
    print "\n";
}
#use Data::Dumper;
#print Dumper(\%taxon);
#print Dumper(\%taxname);

exit(0);
