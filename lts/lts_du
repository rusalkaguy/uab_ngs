#!/usr/bin/env bash
#
# lts_du [lts:/]
#
# concise list of s3 bucket size
#
# uses rclone

# parse arges
export SERVICE="lts"
if [ ! -z "$1" ]; then export SERVICE="$1";fi
echo "# SERVICE=$SERVICE"

# load rclone
module load rclone

# get list of top-level buckets
export BUCKETS=$(rclone lsd $SERVICE:/ | awk '{printf $5" "}')
echo "# BUCKETS=$BUCKETS"

# get sizes for all buckets
(
echo "BUCKET NAME OBJS BYES SIZE SIZE_U"
( for BUCKET in $BUCKETS; do 
    echo $(rclone size $SERVICE:/$BUCKET) |awk -v NAME=$BUCKET 'BEGIN{OFS="\t";COUNT=3;SIZE=6;SIZEU=7;BYTES=8}{B=$BYTES;sub(/[^0-9]/,"",B);print "bucket",NAME, $COUNT, B, $SIZE, $SIZEU}'
done ) | sort -k3,3nr 
) | rev | column -t | rev


