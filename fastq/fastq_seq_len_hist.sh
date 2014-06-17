#!/bin/bash
awk '(2==NR%4){print length($0);}' $* \
    | sort -n \
    | uniq -c \
    | awk 'BEGIN{OFS="\t";print "length","count"; } {print $2,$1;}'
