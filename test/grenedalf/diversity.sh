#!/bin/bash

# Parse the args. We use the key to set a variable named that way.
for arg in "$@"; do
    key=${arg%%=*}
    val=${arg#*=}
    eval "$key"='$val'
done

# Set the args that we need here
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
GRENEDALF="${SCRIPT_DIR}/../../bin/grenedalf"

mkdir -p diversity
mkdir -p logs

echo "Start `date`"
START=$(date +%s.%N)

${GRENEDALF} diversity \
    --pileup-path ${fileid} \
    --window-type interval \
    --window-interval-width ${windowsize} \
    --filter-sample-min-count 2 \
    --filter-sample-min-read-depth 2 \
    --pileup-min-base-qual 0 \
    --pool-sizes ${poolsize} \
    --window-average-policy "window-length" \
    --tajima-d-denominator-policy "provided-min-read-depth" \
    --out-dir "diversity" \
    --file-suffix "-${outid}" \
    --na-entry nan \
    --no-extra-columns \
    --allow-file-overwriting \
    > "logs/diversity-${outid}.log" 2>&1

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

echo "End `date`"
echo "Internal time: ${DIFF}"
