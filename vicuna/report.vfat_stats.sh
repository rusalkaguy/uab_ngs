#!/bin/bash
#
# parse $SAMPLE/vfat/{REF}-QA

# parse stats files
for vfat_stats in `\ls */vfat/*_QA_StatsDetailed.txt`; do
    D1=`dirname $vfat_stats`
    SAMPLE=`dirname $D1`
    echo "Parsing $SAMPLE"
    ~/uab_ngs/vicuna/vicuna_qa_stat.pl \
	-sample_name $SAMPLE \
	-in $vfat_stats
done

# concat results
MERGE_FILE=vfat_QA_StatsDetailed.csv
TMP_FILE=`mktemp`
for vfat_stats in `\ls */vfat/*_QA_StatsDetailed.txt`; do
    CSV=`dirname $vfat_stats`/`basename $vfat_stats .txt`.csv
    cat $CSV >> $TMP_FILE
done
head -n 1 $TMP_FILE > $MERGE_FILE
grep -v "^sample_name," $TMP_FILE >> $MERGE_FILE

# tar up individual files
echo "Creating stats.tar.gz ...."
tar czf stats.tar.gz $MERGE_FILE */vfat/*_QA_StatsDetailed.*
echo "done"
