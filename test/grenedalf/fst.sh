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

mkdir -p fst
mkdir -p logs

echo "Start `date`"
START=$(date +%s.%N)

${GRENEDALF} fst \
    --sync-path ${fileid} \
    --window-type interval \
    --window-interval-width ${windowsize} \
    --pool-sizes ${poolsizes} \
    --method ${method} \
    --window-average-policy "window-length" \
    --no-extra-columns \
    --out-dir "fst" \
    --file-suffix "-${outid}" \
    --na-entry nan \
    --allow-file-overwriting \
    > "logs/fst-${outid}.log" 2>&1

END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)

echo "End `date`"
echo "Internal time: ${DIFF}"
