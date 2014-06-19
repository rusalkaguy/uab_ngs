#!/bin/bash
#
# compress and tabix a VCF
#
#------- qsub---------
#$ -N tabix_bgzip
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
	# add self-qsub here!
    FLAGS="FLAGS $1"
    shift 1
done

# iterate over targets
while [[ -n "$1" ]]; do
    # get target from ARGS
    TARGET=$1; shift
    echo "# TARGET=$TARGET"

    # check existance of input
    if [ ! -e "$TARGET" ]; then
	echo "ERROR: no such TARGET: $TARGET"
	exit 1
    fi

    # figure out compress named
    if [[ "$VCF" == *.vcf.gz  ]]; then 
	GZ=$TARGET
    else 
	GZ="${TARGET}.gz"
    fi

    # compress if not compressed
    if [[ -e "${TARGET}.gz" && "$GZ" -nt "$TARGET" ]]; then 
	echo "SKIP: Found up-to-date TARGET.gz - hope it's bgzipped..."
    else
	echo "bgzip -c $TARGET"
	bgzip -c $TARGET > $GZ
	RC=$?
	if [ "$RC" != 0 ]; then 
	    echo "ERROR: RC=$RC; bgzip -c $TARGET > $GZ"
	    exit $RC
	fi
    fi

    # tabix index - setup flags
    if [[ "$GZ" == *.vcf.gz && -z "$FLAGS" ]]; then FLAGS="-p vcf"; fi
    if [[ "$GZ" == *.sam.gz && -z "$FLAGS" ]]; then FLAGS="-p sam"; fi
    if [[ "$GZ" == *.gff.gz && -z "$FLAGS" ]]; then FLAGS="-p gff"; fi
    if [[ "$GZ" == *.bed.gz && -z "$FLAGS" ]]; then FLAGS="-p bed"; fi
    if [[ -e "${GZ}.tbi" || "${GZ}.tbi" -nt "$GZ" ]]; then 
	echo "SKIP: Found up-to-date TARGET.gz.tbi ..."
    else
	CMD="tabix $FLAGS $GZ"
	echo $CMD
	$CMD
	RC=$?
	if [ "$RC" != 0 ]; then 
	    echo "ERROR: RC=$RC; $CMD"
	    exit $RC
	fi
    fi
    TARGET=""
done
