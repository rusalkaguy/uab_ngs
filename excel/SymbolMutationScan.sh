#! /bin/csh
# written by b zeeberg 11/18/02

foreach arg ($*)
echo $arg
#The following line is supposed to be split in two for some arcane reason in sed
tr "\r" "\n" < $arg | sed '1,$s/\^M/\\
/g' | gawk '(/[0-9]\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)/ || /[0-9]\.[0-9][0-9]E\+[[0-9][0-9]/) {print NR,$0}'
end