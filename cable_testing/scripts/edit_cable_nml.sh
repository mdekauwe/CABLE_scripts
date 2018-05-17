#!/bin/bash

# Edits CABLE namelist files in-place.
#
# TODO: Currently only edits options that already exist. If adding options is
# required, perhaps see https://stackoverflow.com/questions/15965073
#
# TODO: Not sure how to deal with case-sensitivity
#
# Author: ned haughton <ned@nedhaughton.com>

set -eux

SCRIPT=$0
BASENAME=$(basename $SCRIPT)

function usage {
    echo "usage: $SCRIPT [-o OUTFILE] INFILE KEY=VALUE [KEY=VALUE...]"
    exit 1
}

OUTFILE=""

while getopts "o:" flag; do
    case $flag in
        o) OUTFILE=$OPTARG ; shift 2;;
    esac
done

[ 2 -gt "$#" ] && { usage; }

INFILE=$1
shift

CHANGES=$@

function sed_pattern {
    PARTS=(${1//=/ })
    NOTE="! [edited by $BASENAME]"
    echo "s#^\s*${PARTS[0]}\s*=\s*[^!]*\s*\(!.*\)*\( $NOTE\)*\$#   ${PARTS[0]} = ${PARTS[1]}  \1 $NOTE#"
}

SED_CMDS=""
for c in $CHANGES ; do
    SED_PAT=$(sed_pattern $c)
    SED_CMDS="$SED_PAT; $SED_CMDS"
done

if [ -z "$OUTFILE" ] ; then
    sed -i "$SED_CMDS" -i $INFILE
    echo "Edited $INFILE in-place"
else
    sed "$SED_CMDS" $INFILE > $OUTFILE
    echo "Edited version of $INFILE saved to $OUTFILE"
fi
