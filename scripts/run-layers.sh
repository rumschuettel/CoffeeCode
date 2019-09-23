#!/bin/bash

trap ctrl_c SIGINT
function ctrl_c() {
    echo "ABORTED"
	mv "$CURRENT" "$CURRENT.fail"
	exit
}

OUTDIR="./volumes/concat_codes"
HOSTNAME=`hostname`

sleep $[ ( $RANDOM % 50 )+1 ]s

for f in `ls -Sr results/concat_codes/*-in-*.json.gz.tar`; do
	bf=`basename "$f"`
	of="$OUTDIR/$bf-rate-layers.npz"
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
	nice -19 ./postprocess-layerplot.py --external=True --outfile "$of" "$f"
	rm "$lf"
done
