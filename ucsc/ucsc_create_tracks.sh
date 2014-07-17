#
# create a UCSC tracks.txt file
# that loads all VCFs in this directory
#
MAP=name_map.txt

PREFIX_URL=$1; shift
if [[ (-z "$PREFIX_URL" || "$PREFIX_URL" == "MAP") && -e "$MAP" ]]; then
    PREFIX_URL=`grep "^#URL_BASE" $MAP | cut -f 2`
fi
if [ -z "$PREFIX_URL" ]; then
    echo "ERROR: no PREFIX_URL specified" >2 
    echo "SYNTAX: $0 PREFIX_URL" >2 
    exit 1
fi
FILES=$*
if [ -z "$FILES" ]; then
    FILES=*.vcf.gz
fi
if [[ -e "$MAP" && "$FILES"=="MAP" ]]; then
    FILES=`grep -v "^#" $MAP | cut -f 1`
fi

# location header
(echo "browser position chr10:49,872,280-50,199,447"; echo "")

# add each VCF.gz
for X in $FILES; do 
    if [ ! -e $X.tbi ]; then echo "ERROR: missing: $X.tbi" >2; exit 1; fi
    
    FNAME=`basename $X .vcf.gz`
    NAME=`grep -w "^$FNAME" name_map.txt | cut -f 2`
    if [ -z "$NAME" ]; then NAME=`basename $X .vcf.gz|cut -d . -f 4-`; fi
    COLOR=`grep -w "^$FNAME" name_map.txt | cut -f 3`
    if [ -z "$COLOR" ]; then COLOR="0,0,0"; fi

    echo -n "track "
    echo -n "type=vcfTabix "
    echo -n "name=$NAME description=\"$NAME\" "
    echo -n "color=$COLOR " 
    echo -n "visibility=dense "
    echo "bigDataUrl=$PREFIX_URL/$X"; 
    echo "" 
done

