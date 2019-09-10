#!/bin/bash

trap ctrl_c SIGINT
function ctrl_c() {
    echo "ABORTED"
	mv "$CURRENT" "$CURRENT.fail"
	exit
}

OUTDIR="./volumes/special_felix"
HOSTNAME=`hostname`

for f in `ls -Sr results/special_felix/*.json.gz.tar`; do
	bf=`basename "$f"`
	of="$OUTDIR/$bf-best.npz"
	if [ -f "$of" ]; then
		continue
	fi
	lf="$of.lock"

	if [ -f "$lf" ]; then
		echo "skipping $bf"
		continue
	fi
	echo "$HOSTNAME" >> "$lf"
	CURRENT="$lf"

	echo "running $bf"
	nice -19 ./postprocess-3dplot.py --outfile "$of" "$f"
	rm "$lf"
done
