#!/usr/bin/bash

set -e

# call this script like this:
# $ ./makeMergedBigwig.sh ./path/to/replicate/dirs/ REPLICATE_1,REPLICATE_2,REPLICATE_3
# the sample argument should be a comma-separated list with no spaces between names

# the directory for the genotype, including the trailing /
directoryPrefix=$1

# split the sample list into an array
IFS=',' read -r -a samples <<< "$2"

# specify effective genome size
genomeSize=135000000

# run the following for each given sample
for sample in ${samples[@]}
do
  samtools merge -o "${directoryPrefix}"/"${sample}"_mergedBam/"${sample}"_merged.bam "${directoryPrefix}"/"${sample}"_1/"${sample}"_1.bam "${directoryPrefix}"/"${sample}"_2/"${sample}"_2.bam
  samtools index "${directoryPrefix}"/"${sample}"_mergedBam/"${sample}"_merged.bam
  bamCoverage -b "${directoryPrefix}"/"${sample}"_mergedBam/"${sample}"_merged.bam -o "${directoryPrefix}"/"${sample}"_bigwig/"${sample}"_merged.bigwig -of bigwig --binSize 10 --normalizeUsing CPM --effectiveGenomeSize $genomeSize
done

