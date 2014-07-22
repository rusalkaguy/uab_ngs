#!/bin/bash
#
# compress and tabix a VCF
#
#------- qsub---------
#$ -N bgzip_tabix
#$ -S /bin/bash
#$ -cwd # remember what dir I'm launched in 
#$ -V   # need this for parameter passing from MASTER to SLAVE!
#$ -m beas  #email at:  Begining, End, Abort, Suspend
# *** output logs ***
#$ -j y # merge stderr into stdout
#$ -l h_rt=5:00:00 -l s_rt=4:55:00  # 5 hr run time
#$ -l vf=1.8G -l h_vmem=2G          # 1 CPUs, 2G
#*** END QSUB ****
module load ngs-ccts/tabix/0.2.6 

# arg parsing
while [[ "$1" == -* ]]; do 
    FLAG=$1; shift 1
    if [[ "$FLAG" == -qsub ]]; then DO_QSUB=yes; continue; fi
    FLAGS="$FLAGS $FLAG"
done
if [[ `basename $0` == qsub* ]]; then DO_QSUB=yes; fi

# iterate over targets
while [[ -n "$1" ]]; do
    # get target from ARGS
    if [ -n "$1" ]; then TARGET=$1; shift; fi
    echo "# TARGET=$TARGET"

    # check existance of input
    if [ ! -e "$TARGET" ]; then
	echo "ERROR: no such TARGET: $TARGET"
	exit 1
    fi

    # figure out compress named
    if [[ "$TARGET" == *.gz  ]]; then 
	GZ=$TARGET
    else 
	GZ="${TARGET}.gz"
    fi

    # compress if not compressed
    if [[ -e "$GZ" && ("$GZ" -nt "$TARGET" || "$TARGET" == "$GZ" ) ]]; then 
	echo "SKIP: Found up-to-date TARGET.gz - hope it's bgzipped..."
    else
        # if auto-qsub, do that here, per sample
	if [ -n "$DO_QSUB" ]; then echo -n "qsub job_id="; qsub -terse -N `basename $0 .sh`-`basename $TARGET` -M $USER@uab.edu $0 $TARGET; continue; fi

	# do work in-line
	echo "bgzip -c $TARGET"
	bgzip -c $TARGET > $GZ
	RC=$?
	if [ "$RC" != 0 ]; then 
	    echo "ERROR: RC=$RC; bgzip -c $TARGET > $GZ"
	    exit $RC
	fi
    fi

    # tabix index - setup flags
    if [[ -e "${GZ}.tbi" && "${GZ}.tbi" -nt "$GZ" ]]; then 
	echo "SKIP: Found up-to-date TARGET.gz.tbi ..."
    else
        # if auto-qsub, do that here, per sample
	if [ -n "$DO_QSUB" ]; then echo -n "qsub job_id="; qsub -terse -N `basename $0 .sh`-`basename $TARGET` -M $USER@uab.edu $0 $TARGET; continue; fi

	# do work in-line
	# setup flags from extension
	if [[ "$GZ" == *.vcf.gz && -z "$FLAGS" ]]; then FLAGS="-p vcf"; fi
	if [[ "$GZ" == *.sam.gz && -z "$FLAGS" ]]; then FLAGS="-p sam"; fi
	if [[ "$GZ" == *.gff.gz && -z "$FLAGS" ]]; then FLAGS="-p gff"; fi
	if [[ "$GZ" == *.bed.gz && -z "$FLAGS" ]]; then FLAGS="-p bed"; fi

	CMD="tabix $FLAGS -f $GZ"
	echo $CMD
	$CMD
	RC=$?
	if [ "$RC" != 0 ]; then 
	    echo "ERROR: RC=$RC; $CMD"
	    exit $RC
	fi
    fi
done
