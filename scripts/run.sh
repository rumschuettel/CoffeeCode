#!/bin/bash

trap ctrl_c SIGINT


INFILE=$1
INFILEBASE=${INFILE##*/}
OUTDIR="/results/jo/all-env-processed/$INFILEBASE"

function ctrl_c() {
        echo "ABORTED"
	echo "manually clean up $OUTDIR"
	exit
}

# skip if directory exists, otherwise create
if [ -d "$OUTDIR" ]
then
	echo "skipping $INFILE"
	exit
fi
mkdir -p "$OUTDIR"
echo "processing $INFILE"

rungraph() {
	local kSys=$1
	local kEnv=$2
	local adjmat=$3
	local outfile=$4

	../build/release_full_arbch/CoffeeCode.$kSys.$kEnv <<< $adjmat > $outfile
}
export -f rungraph

while read -r kSys kEnv graph adjmat; do
	outfile="$OUTDIR/$adjmat"
       	rungraph "$kSys" "$kEnv" "$adjmat" "$outfile"
done < $INFILE 

echo "done $INFILE"

