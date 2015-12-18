# sum_fpkm.awk 
# 
# ARGUMENTS
#     GENE_TXN_MAP := name of mapping file [required]
#     DEBUG := boolean [default 0]
#     GENE_NAMES := boolean [default 0] - prints gene names/ids
#     SAMPLE_NAME := name of sample to print in output [default ""]
#     COLUMN := transpose output: one gene per column 
# EXAMPLE
#
# # one gene per row
# awk -v GENE_NAMES=1 -v SAMPLE_NAME=SL62334  -v GENE_TXN_MAP=gene_table_format2.txt  -f sum_txns2.awk SL62334.isoforms.fpkm_tracking
#
# # one gene per column 
# awk -v GENE_NAMES=1 -v COLUMN=1 -v SAMPLE_NAME=SL62334  -v GENE_TXN_MAP=gene_table_format2.txt  -f sum_txns2.awk SL62334.isoforms.fpkm_tracking
#
# # just sums in a row with sample_id
# awk -v COLUMN=1 -v SAMPLE_NAME=SL62334  -v GENE_TXN_MAP=gene_table_format2.txt  -f sum_txns2.awk SL62334.isoforms.fpkm_tracking
#
# # just sums in a column
# awk -v GENE_TXN_MAP=gene_table_format2.txt  -f sum_txns2.awk SL62334.isoforms.fpkm_tracking
#

BEGIN {
  OFS="\t";
  # ----------------------------------------------------------------------
  # 
  # load tab-separated mapping file (tripples): gene_id gene_name txn_id
  #
  # ----------------------------------------------------------------------
  # creates associative arrays:
  #     gene[1..gene_count] = gene_name
  #     gene_list[gene_name] = gene_id
  #     id2gene[gene_id|txn_id] = gene_name
  if(!GENE_TXN_MAP) {
    print "ERROR: GENE_TXN_MAP must be set to the filename of the mapping file";
    print "mapping file should contain tab-separated tripples: gene_id gene_name txn_id"
    exit(1);
  }
  if(DEBUG){print "##### loading mapping file" GENE_TXN_MAP};
  while (getline < GENE_TXN_MAP) {
    if(DEBUG){print "line ", $0};
    # parse mapping file line
    ct=split($0,a,"\t");
    gene_id=a[1];
    gene_name=a[2];
    txn_id=a[3];
    if( $0 ~ /^#/ ) { continue; } # skip comments
    if(gene_list[gene_name] == "") {
      # new gene
      gene_count++;
      gene[gene_count]=gene_name;
      gene_list[gene_name]=gene_id;
    }
    # associate gene_id with gene_name
    id2gene[gene_id]=gene_name; 
    if(DEBUG){print gene_id "=>" gene_name}
    # associate txn_id with gene_name
    id2gene[txn_id]=gene_name;
    if(DEBUG){print txn_id "=>" gene_name}
  }
  if(DEBUG){print "##### mapping file loaded"}

  #
  # column mappings for FPKM file
  #
  cTRACKING_ID  = 1
  cCLASS_CODE  = 2
  cNEAREST_REF_ID  = 3
  cGENE_ID  = 4
  cGENE_SHORT_NAME  = 5
  cTSS_ID  = 6
  cLOCUS  = 7
  cLENGTH  = 8
  cCOVERAGE  = 9
  cFPKM  = 10
  cFPKM_CONF_LO  = 11
  cFPKM_CONF_HI  = 12
  cFPKM_STATUS  = 13
}
# ----------------------------------------------------------------------
# 
# parse fpkm file (gene or iso)
# 
# ----------------------------------------------------------------------
{
  # strip suffix off accession
  split($cTRACKING_ID,a,".");tracking_acc=a[1]
}
(id2gene[tracking_acc] != "") {
  # sum fpkm value for gene symbol associated with this accession
  gene_name=id2gene[tracking_acc]
  summed[gene_name]=summed[gene_name]+$cFPKM
  if(DEBUG){print $1 "=>" tracking_acc "=>" gene_name " = "   summed[gene_name]}
}
# ----------------------------------------------------------------------
# 
# output final per-gene list
#
# ----------------------------------------------------------------------
END {
  if(DEBUG){print "##### END - dump results"}
  if( ! COLUMN ) {
    # ROW OUTPUT

    # one gene per row
    for(g=1; g<= gene_count; g++) {
      gene_name=gene[g];
      gene_id=gene_list[gene_name];
      if( SAMPLE_NAME ) {
	printf SAMPLE_NAME OFS;
      }
      if( GENE_NAMES ) {
	printf gene_name OFS gene_id OFS;
      }
      print summed[gene_name];
    }
  } else {  
    # COLUMN output
    # one gene per col
    if( GENE_NAMES ) {
      printf "#gene_name"
      # gene_names row
      for(g=1; g<= gene_count; g++) {
	printf OFS gene[g]
      }
      print "";
      # gene_ids row
      printf "#gene_id"
      for(g=1; g<= gene_count; g++) {
	printf OFS gene_list[gene[g]]
      }
      print "";
    }
    # sample & data  row
    if( SAMPLE_NAME ) {
      printf SAMPLE_NAME;
    }else {
      printf "summed_fpkm"
    }
    for(g=1; g<= gene_count; g++) {
      printf OFS summed[gene[g]]
    }
    print "";
  }
}
  