#!/bin/bash

# Parse the args. We use the key to set a variable named that way.
for arg in "$@"; do
    key=${arg%%=*}
    val=${arg#*=}
    eval "$key"='$val'
done

METHOD=$method
if [[ "$METHOD" == "karlsson" ]]; then
    METHODSTR="--karlsson-fst"
else
    METHODSTR=""
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
POPOOL="${SCRIPT_DIR}/popoolation2"

mkdir -p fst
mkdir -p logs
# rm fst/*

echo "Start `date`"
START=$(date +%s.%N)

perl ${POPOOL}/fst-sliding.pl \
    --input ${fileid} \
    ${METHODSTR} \
    --output "fst/${outid}.fst" \
    --window-size ${windowsize} \
    --step-size ${windowsize} \
    --pool-size ${poolsizes} \
    --min-count 2 \
    --min-coverage 4 \
    --max-coverage 100 \
    --min-covered-fraction 0 \
    > logs/fst-${outid}.log 2>&1

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

echo "End `date`"
echo "Internal time: ${DIFF}"
