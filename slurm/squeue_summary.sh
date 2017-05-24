#!/bin/bash

(
	echo "COUNT USER NAME ST PARTITION NODES"\
	; 
	squeue $* \
	| awk '($1 != "USER"){print $4, $3, $5, $2, $7}' \
	| sort -k1,1 -k2,2 -k3,3 \
	| uniq -c
) \
| column -t
