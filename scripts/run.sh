#!/bin/bash

trap ctrl_c SIGINT


INFILE="$1"
INFILEWOPATH="${INFILE##*/}"
INFILEBASE="${INFILEWOPATH%.*}"
INTERMEDIATEDIR="./cc-runs/$INFILEBASE"
OUTDIR="/results/jo/all-env-processed/$INFILEBASE"

function ctrl_c() {
        echo "ABORTED"
	echo "manually clean up $OUTDIR"
	exit
}

# skip if directory exists, otherwise create
MATFILE=`find ./out-all-env-adjm -name "*$INFILEBASE*" -print`
MATFILEWOPATH="${MATFILE##*/}"
MATFILEBASE="${MATFILEWOPATH%.*}"
OUTMATFILE="$OUTDIR/$MATFILEWOPATH.gz"
if [ -f "$OUTMATFILE" ]
then
	echo "skipping $INFILEBASE"
	exit
fi
# check mat file exists
if [ ! -f "$MATFILE" ]
then
	echo "adjm file missing for $INFILEBASE"
	exit
fi



echo "processing $MATFILE"


# create intermediate directory to store CC runs in
# the subfolders are to circumvent a dir_index bug
mkdir -p "$INTERMEDIATEDIR"
cd "$INTERMEDIATEDIR"
seq 10000 | xargs mkdir
cd -

rungraph() {
	local kSys=$1
        local kEnv=$2
	local adjmat=$3

	local rndfolder="$((1 + RANDOM % 10000))"
	local outfile="$INTERMEDIATEDIR/$rndfolder/$adjmat.json.gz"
	echo "../build/release_full_arbch/CoffeeCode.$kSys.$kEnv <<< $adjmat | gzip --best > $outfile"
}

runall() {
	while read kSys kEnv graph adjmat; do
		rungraph $kSys $kEnv $adjmat
	done < $MATFILE
}

runall | parallel -j 128 --pipe --block 250k bash


# pack up everything and archive on persistent store
mkdir -p "$OUTDIR"

echo "zipping folder $INTERMEDIATEDIR/"
tar -cf - "$INTERMEDIATEDIR/" | pigz --best > "$OUTDIR/$MATFILEBASE.lambdas.tar.gz"
rm -r "$INTERMEDIATEDIR" &

echo "zipping and moving $MATFILE"
pigz -c "$MATFILE" > "$OUTMATFILE"
rm "$MATFILE"

echo "done $MATFILE"

