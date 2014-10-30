#!/bin/bash
#
# Install R libraries GATK depends on 
# 
# technique from 
#     http://www.stat.osu.edu/computer-support/mathstatistics-packages/installing-r-libraries-locally-your-home-directory
#
module load ngs-ccts/GATK/3.3-0

GATK_TAG=`echo $GATK_VER  | cut -d - -f 1`

# AnalyzeCovariates requires graphing libraries
# http://biosupport.se/questions/706/gatk-analyzecovariates-error
AC_SCRIPT_URL=https://raw.githubusercontent.com/broadgsa/gatk/${GATK_TAG}/public/gatk-tools-public/src/main/resources/org/broadinstitute/gatk/utils/recalibration/BQSR.R
AC_LIBS=` \
    wget -q --no-check-certificate -O - $AC_SCRIPT_URL \
    | grep "^library" \
    | perl -pe 's/^library[("]+([^)"]+).*/$1/;' \
`
echo `basename $AC_SCRIPT_URL`" requires $AC_LIBS"
LIB_LIST="$LIB_LIST $AC_LIBS"

#
# setup GATK-specific R library location
#
export R_LIBS=$PWD/R_libs

chmod u+w .
mkdir -p $R_LIBS
chmod -R u+w $R_LIBS

#
# install list of libraries
#
for RLIB in $LIB_LIST ; do
    Rscript -e " install.packages('$RLIB', repos=\"http://cran.r-project.org\")"
done

#
# re-secure directory
#
chmod -R -w . 
