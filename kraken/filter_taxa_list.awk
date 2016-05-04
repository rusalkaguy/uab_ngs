BEGIN { 
    FS="\t"; OFS="\t";
    #debug=0
}
# load taxa-of-interest table: c1 = taxid
(tax_filename == "") {
  if(1==debug){print "Loading " FILENAME}
  tax_filename=FILENAME
}
(FILENAME==tax_filename && ! /^#/ ) {
    taxarr[$1]=1;
    if(1==debug){print "\tloaded["$1"] => " $0}
    next;
}
(tax_filename != FILENAME && tax_filename != "" && filter_start == "") {
  if(1==debug){print "Filtering " FILENAME}
  filter_start = FILENAME
}
# get seqIds from kreken.txt for matching taxa
(FILENAME!=tax_filename && "C"==$1 && 1==taxarr[$3]) {
  if(1==debug){
    print $2, $3, taxarr[$3] # seqID, taxid, flag
  } else {
    print $2  # print fasta seq ID
  }
}
END {
  if(1==debug){print "done"}
}
