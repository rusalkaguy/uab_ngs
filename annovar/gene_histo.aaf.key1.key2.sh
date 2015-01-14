# 1:#CHROM
# 2:start
# 3:end
# 4:gene
# 5:accession
# 6:ESRD_HOM
# 7:ESRD_HET
# 8:ESRD_TOT       ESRD
# 9:SLE_HOM
# 10:SLE_HET
# 11:SLE_TOT       SLE
# 12:HEALTHY_HOM
# 13:HEALTHY_HET
# 14:HEALTHY_TOT   NORM
# 15:ESRD-SLE
# 16:SLE-NORM
# 17:ESRD_TOT:SLE_TOT:HEALTHY_TOT


AAF=00;  if [ "$1" != "" ]; then AAF=$1;  fi; echo AAF=$AAF
KEY1=15; if [ "$2" != "" ]; then KEY1=$2; fi; echo KEY1=$KEY1
KEY2=16; if [ "$3" != "" ]; then KEY2=$3; fi; echo KEY2=$KEY2
INFILE=COLS/bwamem_pad200_hg19.fs500.adp.dp13.gq30.miss.pass.left.dbsnp.caseControl.annovar.data.cols.filt.aaf$AAF.sum.ALL.txt
echo INFILE=$INFILE
OUTFILE=COLS/HIST/bwamem_pad200_hg19.fs500.adp.dp13.gq30.miss.pass.left.dbsnp.caseControl.annovar.data.cols.filt.aaf$AAF.sum.ALL.k${KEY1}k${KEY2}.txt
echo OUTFILE=$OUTFILE
mkdir -p `dirname $OUTFILE`

wc $INFILE
cat $INFILE \
| awk 'BEGIN{OFS="\t";ESRD=8;SLE=11;NORM=14}\
{if(1==NR){ESRD_SLE="ESRD-SLE";SLE_NORM="SLE-NORM"}else{ESRD_SLE=$ESRD-$SLE;SLE_NORM=$SLE-$NORM};\
print $0,ESRD_SLE,SLE_NORM,$ESRD":"$SLE":"$NORM}' \
| awk -v key1=$KEY1 -v key2=$KEY2 -v label1=4 -v label2=17 -v sep=, \
'BEGIN{OFS="\t"} \
{
 key=$key1 ":" $key2; count_arr[key]++; \
 if(1==NR){header_key=key;header_key1=$key1;header_key2=$key2}; \
 if(label1_arr[key]==""){label1_arr[key]=$label1;label2_arr[key]=$label2} \
 else{label1_arr[key]=label1_arr[key] sep $label1;label2_arr[key]=label2_arr[key] sep $label2};
};\
END{
  count_arr[header_key]="count_key"; \
  for(key in label1_arr){ \
    kv=key;sub(/:/,"\t",kv); \
    print key, kv, count_arr[key], label1_arr[key], label2_arr[key] \
  } \
}\
' \
| sort -k4,4n \
> $OUTFILE

wc $OUTFILE

echo "cat $OUTFILE | column -t | less"

exit 0
