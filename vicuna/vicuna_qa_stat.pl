#!/usr/bin/perl
#
# Parse VICUNA QA_StatsDetailed.txt
#
use JSON -support_by_pp;
#use LWP::Simple;
use Getopt::Long;
use Data::Dumper;
#
# parse params 
#
my %opt= {
    stats_filename => undef,
    sample_name => "sample",
};
# Declare our arguments, and parse the actual inputs
GetOptions(
	   'in=s'		=>\$opt{stats_filename}, 	
	   'sample_name=s'	=>\$opt{sample_name},
	   'no_header'  	=>\$opt{no_header},
	   'verbose'		=>\$opt{verbose},	# flag, no args
	   'help|?'		=> sub { usage(); },	# flag, no args, call subroutine
	   ) || &usage("missing value");
&usage("-in REQUIRED") if( !$opt{stats_filename} );
&usage("-sample_name REQUIRED") if( !$opt{sample_name} );

sub usage {
    print STDERR "SYNTAX $0 [-v] -sample sam1 -in QA_StatsDeftailed.txt\n";
    exit 1;
}

my @sections = (
		# vfat stats
		"Assembly QC Coverage Stats",
		"Assembly QC General Stats",
		"Reference data",
		"Files used:",
		# flagstats 
		"Mapping to Reference",
		"Mapping to Assembly",
		"Variants vs Assembly",
		);

# parse VFAT files
print "SAMPLE: $opt{sample_name}\n" if( $opt{verbose} );
print "PARSE: $opt{stats_filename}\n" if( $opt{verbose} );
$stats = parse_vfat_stats_file($opt{stats_filename});

# parse samtools flagstat files
# REF
$opt{sample_vs_ref_stats}=$opt{stats_filename};
$opt{sample_vs_ref_stats}=~ s|(.*)_QA_StatsDetailed.txt$|$1_alignVsRef_flagstat.txt|;
$stats->{$sections[4]} = parse_flagstats_file($opt{sample_vs_ref_stats});
# ASSEMBLY
$opt{sample_vs_asm_stats}=$opt{stats_filename};
$opt{sample_vs_asm_stats}=~ s|(.*)_QA_StatsDetailed.txt$|$1_alignVsAssembly_flagstat.txt|;
$stats->{$sections[5]} = parse_flagstats_file($opt{sample_vs_asm_stats});

# parse VCFs
# ASSMEBLY
$opt{variants_vs_asm_vcf}=$opt{stats_filename};
$opt{variants_vs_asm_vcf}=~ s|(.*)_QA_StatsDetailed.txt$|$1_alignVsAssembly_sort.bam.vcf|;
$stats->{$sections[6]} = parse_vcf_file($opt{variants_vs_asm_vcf});


#
# output 
#


# CSV
my $csv_filename=$opt{stats_filename};
$csv_filename =~ s/\.txt$/.csv/g;
open($csv, ">", $csv_filename ) || die "${csv_filename}:$!\n";
{
    # sample name
    my @headers = ("sample_name");
    my @values = ($opt{sample_name});
    # sections
    foreach my $sec_name (@sections) {
	# keys sorted by line number
	foreach $key_name (
			   # values sorted by line number
			   sort { $stats->{$sec_name}->{$a}->{line_num} <=> $stats->{$sec_name}->{$b}->{line_num} }
			   keys(%{$stats->{$sec_name}})
			   ) {
	    push @headers, $key_name;
	    push @values, $stats->{$sec_name}->{$key_name}->{value};
	}
    }
    # headers
    if( ! $opt{no_header} ) {
	print $csv join(",", @headers), "\n";
	print join(",", @headers), "\n"  if( $opt{verbose} );
    }
    # values
    print $csv join(",", @values), "\n";
    print join(",", @values), "\n"  if( $opt{verbose} );
}
close($csv);
print "Wrote $csv_filename\n" if( $opt{verbose} );

# Perl
my $pl_filename=$opt{stats_filename};
$pl_filename =~ s/\.txt$/.pl/g;
open($pl, ">", $pl_filename ) || die "${pl_filename}:$!\n";
$stats->{sample_name}=$opt{sample_name};
print $pl Dumper($stats);
close($pl);
print "Wrote $pl_filename\n" if( $opt{verbose} );

# JSON
my $json_filename=$opt{stats_filename};
$json_filename =~ s/\.txt$/.json/g;
open($json, ">", $json_filename ) || die "${json_filename}:$!\n";
my $json_text = to_json($stats, {utf8 => 1, pretty => 1});
print $json $json_text;
close($json);
print "Wrote $json_filename\n" if( $opt{verbose} );

sub parse_vcf_file
{
    my $filename=$_[0];

    if( ! -e $filename ) {
	print STDERR "WARNING: $filename doesn't exist\n";
	return {};
    }
    # open input
    open(my $src, "<", $filename)|| die "$filename : $!";
    #print "OPENED: $filename\n";
    #
    # parse the file
    # 
    my $line_num=0;
    my $section={ 
	variant_count => {line_num=>1,value=>0}, 
	multi_alt_count=> {line_num=>2,value=>0},
	indel_count=>     {line_num=>3,value=>0},
	align_to_n=>     {line_num=>4,value=>0},
	#comment=>     {line_num=>5,value=>0},
	#vcf_line=>     {line_num=>6,value=>0},
    };
    while(my $line=<$src>) {
	$line_num++;
	#$section->{vcf_line}->{value}++;
	chomp($line);
	print "[vcf:$line_num] $line\n";

	# skip comments
	if( $line =~ m/^#/ ) {
	    #$section->{comment}->{value}++;
	    print "SKIP\n"
	    next;
	}

	# parse columns
	my @cols = split(m/\t/,$line);

	# skip lines with N in the ref
	if( $cols[3] eq "N" ) {
	    $section->{align_to_n}->{value}++;
	    print "N_REF\n"
	    next;
	}

	$section->{variant_count}->{value}++;


	# more than one observed SNP
	if( $cols[4] =~ m/,/ ) { 
	    $section->{multi_alt_count}->{value}++;
	    print "MULTI\n"
	}

	# check indels
	if( $cols[7] =~ m/INDEL/ ) {
	    $section->{indel_count}->{value}++;
	    print "INDEL\n"
	}
    }
    close($src);

    #print Dumper($section);
    return $section;
}
sub parse_flagstats_file
{
    my $filename=$_[0];

    if( ! -e $filename ) {
	print STDERR "WARNING: $filename doesn't exist\n";
	return {};
    }
    # open input
    open(my $src, "<", $filename)|| die "$filename : $!";

    #
    # parse the file
    # 
    my $line_num=0;
    my $section={};
    while(my $line=<$src>) {
	$line_num++;
	chomp($line);

	# parse num + num title extra
	if( $line =~ m/^([0-9]+)\s+\+\s+[0-9]+\s+([a-z0-9 ]+)(.*)/ # newer flagstat output - double numbers
	    || $line =~ m/^([0-9]+)\s+([a-z0-9 ]+)(.*)/ # older flagstat output
	    ) {
	    my $key_value=$1;
	    my $key_value_perc='';
	    my $key_name=$2;
	    my $extra = $3;
	    
	    if( $extra =~ m/(.mapQ>=5.)/ ) {
		$key_name .= $1;
	    } elsif ( $extra =~ m/\(([0-9.]+)/ ) {
		$key_value_perc=$1;
	    }
	    #print "KV: flagstat: $key_name = $key_value [$key_value_perc] [$extra]\n";
	    # clean lead/trailing whitespace
	    $key_name=~s/^\s+//;
	    $key_name=~s/\s+$//;

	    # store
	    $section->{$key_name} = { name=>$key_name, value=>$key_value, perc_value=>$key_value_perc, line_num=>$line_num };
	    if( $key_value_perc ) {
		$section->{${key_name}."_perc"} = { name=>$key_name."_perc", value=>$key_value_perc, line_num=>($line_num+0.5) };
	    }
	}
    }
    close($src);

    #print Dumper($section);
    return $section;
}

sub parse_vfat_stats_file
{
    my $stats_filename=$_[0];
    # open input
    open(my $src, "<", $stats_filename)|| die "$stats_filename : $!";

    my %sections = (
		    "Files used:" => {name=>"files"},
		    "Reference data" => {name=>"ref_data"},
		    "Assembly QC General Stats" => {name=>"assembly_stats"},
		    "Assembly QC Coverage Stats" => {name=>"assembly_coverage"}
		    );

    #print "SECTIONS: ", join(",", keys(%sections)), "\n";

    #
    # parse the file
    # 
    my $line_num=0;
    my $section_name=undef;
    my $section=undef;
    while(my $line=<$src>) {
	$line_num++;
	chomp($line);

	# recognize section headers
	if( $sections{$line} ) {
	    $section=$sections{$line};
	    $section_name = $section->{name};
	    #print "[$line_num] ENTERING: $section_name\n";
	    next;
	}

	# ignore blanks and visual separateors
	if( ! $line || $line =~ m/^-+$/ ) {
	    next;
	}

	# parse name-value-units triple
	if( $line =~ m/^([^:(]+)\s*(\([^)]\))*:\s*([^%]+)(.*)/ ) {
	    my $key_name=$1;
	    my $key_comment=$2;
	    my $key_value=$3;
	    my $key_units=$4;
	    #print "KV: $section_name : $key_name = $key_value [$key_units]\n";
	    # clean lead/trailing whitespace
	    $key_name=~s/^\s+//;
	    $key_name=~s/\s+$//;
	    # store
	    $section->{$key_name} = { name=>$key_name, comment=>$key_comment, value=>$key_value, units=>$key_units, value_str=>"$key_value$key_units", line_num=>$line_num };
	}
    }
    close( $src );
    return \%sections;
}

