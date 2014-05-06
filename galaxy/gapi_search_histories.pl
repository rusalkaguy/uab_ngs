#!/usr/bin/env perl
use strict;
use LWP::Simple qw(get);
use Data::Dumper;
use Getopt::Std;
use lib '/home/curtish/lib/perl5/site_perl/5.8.8/';
use JSON        qw(from_json);

# defaults
my $galaxy_url = "http://galaxy.uabgrid.uab.edu";
my $key_env = "GALAXY_API_KEY";
my $dataset_status = "ok";

my %opts;  # option are in %opts
if( !getopts('vHqubg:k:h:d:s:t:l:c:f:', \%opts) ) {
    print <<EOF
$0:  
# --- connection info ---
# k = API key (can alsop be set in $key_env env var [$ENV{$key_env}])
# g = galaxy URL [$galaxy_url]
# --- output ---
# v = verbose [default: no]
# H = history metadata
# q = DONT print header line
# u = DONT print obj URLs
# b = Bracket fields
# f = field separateor [default: tab]
#  --- filters ---
# h = history name REGEX [default: .*]
# d = dataset name REGEX [default: none]
# s = dataset state REGEX [default: $dataset_status]
# t = dataset data_type REGEX [default: .*]
# --- actions ---
# l = symlink datasets into DIR [default: no action]
# c = COPY datasets into DIR [default: no action]
EOF
;
    exit(1);
}
# defaults
$opts{g}= $galaxy_url if( !$opts{g} );
$opts{k}= $ENV{$key_env} if( !$opts{k} );
if( !$opts{k} ) {
    die("The api key is mandatory (-k or $key_env envvar)");
}
my $key_url = "?key=".$opts{k};
$opts{s} = $dataset_status if( !$opts{s} );
my ($field_start,$field_sep,$field_end)= ("","\t","");
if( $opts{f} ) {
    $field_sep=$opts{f};
}
if( $opts{b} ) { 
    # bracket the fields
    ($field_start,$field_end)= ("[","]");
}
my $field_join=$field_end.$field_sep.$field_start;

# data_type -> extension mapping
my $ext_map = {
    "fastqsanger" => "fastq",
    "tabular"     => "txt",
    };

# header line
if( !$opts{q} ) {
    print $field_start, join($field_join,
	       "object.type",
	       "history.name",
	       "object.url",
	       "history.nice_size",
	       "history.deleted",
	       "history.published",
	       "dataset.name",
	       "dataset.type",
	       "dataset.state",
	       "dataset.data_type",
	       "dataset.file_size",
	       "dataset.file_name" 
	       ), $field_end, "\n";
}

	       
	       
# get history list
my $histories_url = $opts{g}."/api/histories".$key_url;
if( $opts{v} ) {print "# Fetching $histories_url\n"; }
my $histories = get($histories_url);
if( $opts{v} ) {print "# got ", length($histories), " bytes\n"; }
my $j = from_json($histories);
if( $opts{v} ) {print "# got ", scalar(@{$j}), " histories\n"; }


# filter NAMES of histories
my @matches = @{$j};
if( $opts{h} ) {
    @matches = grep { $_->{name} =~ m/$opts{h}/i } @{$j};
    if( $opts{v}  ) {print "# matches[$opts{h}] ", scalar(@matches), " histories\n"; }
}
# list matching histories
foreach my $e ( @matches ) {
    # check if we need details on the history
    my $h_info;
    if( $opts{H} ) {
	if( $opts{v} ) {print "# get ", $opts{g}.$e->{url},"\n"; }
	$h_info = from_json(get($opts{g}.$e->{url}.$key_url));    
	#print Dumper($h_info);
    }

    # format output
    $e->{tsv} = join( $field_join, "history", $e->{name},
		      # optional URL
		      $opts{u}?"":($opts{g}.$e->{url}),
		      # optional history info size
		      $h_info->{nice_size},
		      $h_info->{deleted},
		      $h_info->{published},
		      );
    if( !$opts{d} ) {
	print $field_start, $e->{tsv}, $field_end, "\n";
    }
}

# process CONTENTS of histories
if( $opts{d} || $opts{t} ) {
    foreach my $h (@matches) {
	# get history CONTENTS
	if( $opts{v} ) {print "#",$h->{tsv},"\n"; }
	my $hist_url=$opts{g}.$h->{url}."/contents";
	if( $opts{v} ) {print "# get ",$hist_url,"\n"; }
	my $contents = get($hist_url.$key_url);
	if( $opts{v} ) {print "# got ", length($contents), " bytes\n"; }
	my $jcontents = from_json($contents);
	if( $opts{v} ) {print "# got ", scalar(@{$jcontents}), " datasets\n"; }


	if( $opts{d} ) {
	    # filter datasets by name (optional)

	    foreach my $d ( grep {$_->{name} =~ m/$opts{d}/i } @{$jcontents} ) {
		# get details on each data set that matched
		my $jinfo = from_json(get($opts{g}.$d->{url}.$key_url));
		# filter dataset by state (optional)
		if( !$opts{s} || $jinfo->{state} =~ m/$opts{s}/i ) {
		    # filter dataset by data_type (optional)
		    if( !$opts{t} || $jinfo->{data_type} =~ m/$opts{t}/i ) {
			print $field_start, join($field_join, "dataset",
						 $h->{name},
						 $opts{u}?"":($opts{g}.$d->{url}),
						 $d->{name},
						 $d->{type}, 
						 $jinfo->{state}, 
						 $jinfo->{data_type}, 
						 $jinfo->{file_size},
						 $jinfo->{file_name}
						 ), $field_end, "\n";
			# check if should symlink
			if( $opts{c} || $opts{l} ) {
			    # map dataset data_type to file extension
			    my $file_ext = $jinfo->{data_type};
			    if( $ext_map->{$file_ext} ) { 
				if($opts{v}) { print "# mapped $file_ext to $ext_map->{$file_ext}\n"; }
				$file_ext=$ext_map->{$file_ext};
			    }
			    # clean up dataset name to be a file name
			    my $clean_name = $jinfo->{name};
			    $clean_name =~ s/[^A-Za-z0-9\-\.]/_/g;
			    $clean_name .= ".".$file_ext;
			    # copy files out
			    if($opts{c}) {
				my @args = ("cp",$jinfo->{file_name}, $opts{c}."/".$clean_name);
				if($opts{v}) {print "# ",join(" ",@args),"\n";}
				my $status = system(@args);
				if( ($status>>=8) != 0 ) {
				    print "# ERROR: ",join(" ", @args), ": ", $?, "\n";
				}
			    }
			    # link files out
			    if($opts{l}) {
				my @args = ("ln","-sf",$jinfo->{file_name}, $opts{l}."/".$clean_name);
				if($opts{v}) {print "# ",join(" ",@args),"\n";}
				my $status = system(@args);
				if( ($status>>=8) != 0 ) {
				    print "# ERROR: ",join(" ", @args), ": ", $?, "\n";
				}
			    }
			}		    
		    }
		}
	    }
	}
    }
}
