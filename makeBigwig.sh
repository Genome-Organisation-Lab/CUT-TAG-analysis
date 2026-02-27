#!/usr/bin/bash

set -e

# call this script like this:
# $ ./makeBigwig.sh ./path/to/sample/dirs/ SAMPLE_1,SAMPLE_2,SAMPLE_3
# the sample argument should be a comma-separated list with no spaces between names

# split the sample list into an array
IFS=',' read -r -a samples <<< "$1"

# specify effective genome size
genomeSize=135000000

# run the following for each given sample
for sample in ${samples[@]}
do
  bamCoverage -b "${sample}".bam -o "${sample}".bigwig -of bigwig --binSize 10 --normalizeUsing CPM --effectiveGenomeSize $genomeSize
done

