#!/usr/bin/bash

# Make sure that the shell script stops running if it hits an error
set -e

# Provide path to reference genome index files
PATH_TO_REFERENCE_GENOME="ColCEN/ColCEN_index_files"

# Set parameters
picardCMD="java -jar picard.jar"
minQualityScore=2
binLen=500

# Read sample file path if not already set in a shell parameter
 [ -z "$1" ] && read -p "Provide sample file path --> " sampleFilePath || sampleFilePath=$1

# Read sample name if not already set in a shell parameter
 [ -z "$2" ] && read -p "Provide sample name --> " sampleName || sampleName=$2

echo "Running $sampleName ($sampleFilePath)"

# Read trimming
trim_galore --output_dir "$sampleFilePath" --length 30 --quality 20 --stringency 1 -e 0.1 --fastqc --paired "$sampleFilePath"/"$sampleName"_1.fq "$sampleFilePath"/"$sampleName"_2.fq 

rm "$sampleFilePath"/"$sampleName"_1.fq
rm "$sampleFilePath"/"$sampleName"_2.fq

# Read alignment
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x "$PATH_TO_REFERENCE_GENOME" -1 "$sampleFilePath"/"$sampleName"_1_val_1.fq -2 "$sampleFilePath"/"$sampleName"_2_val_2.fq -S "$sampleFilePath"/"$sampleName".sam &> "$sampleFilePath"/"$sampleName"_bowtie2.txt

rm "$sampleFilePath"/"$sampleName"_1_val_1.fq
rm "$sampleFilePath"/"$sampleName"_2_val_2.fq
rm "$sampleFilePath"/"$sampleName"_1_val_1_fastqc.zip
rm "$sampleFilePath"/"$sampleName"_2_val_2_fastqc.zip

# Sort .sam file and add read group (RG) info
$picardCMD SortSam -I "$sampleFilePath"/"$sampleName".sam -O "$sampleFilePath"/"$sampleName".sorted.sam -SORT_ORDER coordinate
rm "$sampleFilePath"/"$sampleName".sam

$picardCMD AddOrReplaceReadGroups I="$sampleFilePath"/"$sampleName".sorted.sam O="$sampleFilePath"/"$sampleName".newSorted.sam RGLB="$sampleFilePath"/"$sampleName" RGSM="$sampleFilePath"/"$sampleName" RGPL=Illumina RGPU=X204SC24095257-Z01-F001
rm "$sampleFilePath"/"$sampleName".sorted.sam

# Deduplicate
$picardCMD MarkDuplicates -I "$sampleFilePath"/"$sampleName".newSorted.sam -O "$sampleFilePath"/"$sampleName".sorted.dupMarked.sam -METRICS_FILE "$sampleFilePath"/"$sampleName".dupMark.txt
$picardCMD MarkDuplicates -I "$sampleFilePath"/"$sampleName".newSorted.sam -O "$sampleFilePath"/"$sampleName".sorted.rmDup.sam -ASSUME_SORTED false -REMOVE_DUPLICATES true -METRICS_FILE "$sampleFilePath"/"$sampleName".rmDup.txt
rm "$sampleFilePath"/"$sampleName".sorted.dupMarked.sam

# Get fragment length information
samtools view -F 0x04 "$sampleFilePath"/"$sampleName".newSorted.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > "$sampleFilePath"/"$sampleName"_fragmentLen.txt
rm "$sampleFilePath"/"$sampleName".newSorted.sam

samtools view -F 0x04 "$sampleFilePath"/"$sampleName".sorted.rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > "$sampleFilePath"/"$sampleName"_fragmentLen.rmDup.txt

# Filter for high quality reads
samtools view -hbS -F 0x04 -q $minQualityScore "$sampleFilePath"/"$sampleName".sorted.rmDup.sam > "$sampleFilePath"/"$sampleName".bam
samtools index "$sampleFilePath"/"$sampleName".bam

# Convert to .bed file
samtools sort -n -o "$sampleFilePath"/"$sampleName"_sortedByName.bam "$sampleFilePath"/"$sampleName".bam
bedtools bamtobed -i "$sampleFilePath"/"$sampleName"_sortedByName.bam -bedpe > "$sampleFilePath"/"$sampleName".bed
rm "$sampleFilePath"/"$sampleName"_sortedByName.bam

awk '$1==$4 && $6-$2 < 1000 {print $0}' "$sampleFilePath"/"$sampleName".bed > "$sampleFilePath"/"$sampleName".clean.bed

cut -f 1,2,6 "$sampleFilePath"/"$sampleName".clean.bed | sort -k1,1 -k2,2n -k3,3n  > "$sampleFilePath"/"$sampleName".cleanFragments.bed

awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' "$sampleFilePath"/"$sampleName".cleanFragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | sort -k1,1V -k2,2n > "$sampleFilePath"/"$sampleName".cleanFragmentsCount.bin$binLen.bed

echo "Done with $sampleName."
