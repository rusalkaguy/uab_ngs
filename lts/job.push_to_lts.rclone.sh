#!/bin/bash
#
# use RCLONE to transfer files to LTS
#
# FILES: *.bam[.bai]
#
#SBATCH --job-name=push_to_lts-kimberly-ics1474-ccs-wgs-bam
#SBATCH --output=%x.%J.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
# special partition for LTS xfers
# like long, has 6d-6h max time:
#  scontrol show partition sciencedmz | egrep -C 100 --color=yes 'MaxTime=[^ ]+'
#SBATCH --partition=sciencedmz
#SBATCH --time=6-06:00:00
##SBATCH --mem-per-cpu=500M # failed on PacBio raw subreads.bam
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=curtish@uab.edu

export DIR="."
export UPDIR=$(echo $DIR | perl -pe 's/[^\/]+/../g')

export RCLONE_SRC="."
export BUCKET="kimberly-ics1474-ccs-wgs-bam"
export RCLONE_DEST="lts:/${BUCKET}/output_minimap_unfiltered/"
#export RCLONE_INCLUDE=" --filter-from <(awk -e 'BEGIN{FS=\"/\"}(\$7!=\"\"){print \"+ \"\$7}END{print \"- *\"}' ../$UPDIR/$SRC_FILE_LIST) "
export RCLONE_INCLUDE=" --include '*.sort.bam'  --include '*.sort.bam.bai' --include 'README.txt'  --include 'README.md'  --include 'push_to_lts.sh'"
# copy only things that haven't changed, ie where the globus transfer has already finished.
#export RCLONE_FLAGS="--min-age 1d " 
#export RCLONE_FLAGS="--min-age 2h " 
echo "RCLONE_SRC=$RCLONE_SRC"
echo RCLONE_DEST=$RCLONE_DEST
echo RCLONE_INCLUDE=$RCLONE_INCLUDE

if [ -z "$(which rclone 2>/dev/null)" ]; then 
    echo module load rclone ...
    module load rclone
fi

#
# copy up 
#
# use "copy" it will avoid re-coyping files that haven't changed, and will delete nothing
# if we use "sync", it will delete things from the bucket we have deleted from the FS (bad)
#
# CPU: normally, rclone defaults to max 4 CPU. If running as a SLURM job, override that with 
# number of job CPUs on node
RCLONE_CPU=" --transfers=${SLURM_CPUS_ON_NODE:-4} --multi-thread-streams ${SLURM_CPUS_ON_NODE:-4}"


echo  rclone $* copy --progress $RCLONE_CPU $RCLONE_SRC $RCLONE_DEST $RCLONE_INCLUDE $RCLONE_FLAGS
eval "rclone $* copy --progress $RCLONE_CPU $RCLONE_SRC $RCLONE_DEST $RCLONE_INCLUDE $RCLONE_FLAGS"

#
# get final md5 list
#

#echo rclone md5sum $RCLONE_DEST \> s3.${BUCKET}.md5.txt
#rclone md5sum $RCLONE_DEST > s3.${BUCKET}.md5.txt


