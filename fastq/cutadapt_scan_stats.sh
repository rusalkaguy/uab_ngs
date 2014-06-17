echo -n "sample,"
echo -n "R1reads,"
echo -n "R1bases,"
echo -n "R2trimmed,R2trimmed,"
echo -n "R1PercAdpt,R2PercAdpt,"
echo -n "R1PercLowQ,R2PercLowQ,"
echo -n "R1goodBases,R2goodBases"
echo ","
for x in SR*; do 
    R1="$x/${x}_R1.cutadapt.fastq.gz";
    R2=`echo $R1|perl -pe 's/_R1/_R2/;'`;
    S1="$x/${x}_R1.cutadapt.stats";
    S2=`echo $S1|perl -pe 's/_R1/_R2/;'`;
    echo -n "$x,";
    echo -n `grep Processed.reads $S1 | perl -pe 's/.* ([0-9.]+).*/$1/;s/\n/,/;'`;
    echo -n `grep Processed.bases $S1 | perl -pe 's/.* ([0-9.]+) bp.*/$1/;s/\n/,/;'`;
    echo -n `grep Trimmed.reads $S1 $S2| perl -pe 's/.*\(([0-9.]+%).*/$1/;s/\n/,/;'`;
    echo -n `grep Trimmed.bases $S1 $S2 | perl -pe 's/.*\(([0-9.]+%).*/$1/;s/\n/,/;'`;
    echo -n `grep Quality-trimmed $S1 $S2| perl -pe 's/.*\(([0-9.]+%).*/$1/;s/\n/,/;'`;
    echo -n `~/uab_ngs/fastq/fastq_bases.sh -q $R1 $R2 | perl -pe 's/\n/,/;'`;
    echo ""
done

