#!/bin/bash
perl -e 'while(<>) {if( $ln % 4 == 0 ) { s/^@(.*)/>$1/; print $_; } elsif( $ln % 4 == 1 ) { print $_; } $ln++;} print stderr $ln/4, "\n";' < $1
exit $?

