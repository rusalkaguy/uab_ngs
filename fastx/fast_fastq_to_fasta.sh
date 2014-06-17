#!/bin/bash
#
# dumb fastq to fasta conversion - just loose the quality!
#
perl -e 'my $ln=0; while(my $line=<>) { if($ln%4==0) { $line =~ s/^@//; $line =~ s/ /_/g; print ">$line"; } elsif($ln%4==1){print $line;}else {;} $ln++;}' 
exit $? 2>&1 > /dev/null
ZZ
