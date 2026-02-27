# CUT-TAG-analysis

A set of command line tools for CUT&Tag data processing and analysis, adapted from the [tutorial](https://yezhengstat.github.io/CUTTag_tutorial/) by Zheng Y *et al.* (2020) Protocol.io.

## Pre-requisites and dependencies
* [trim-galore](https://github.com/FelixKrueger/TrimGalore)
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [picard](https://github.com/broadinstitute/picard)
* [samtools](https://github.com/samtools/samtools)
* [bedtools](https://github.com/arq5x/bedtools2)
* [deeptools](https://github.com/deeptools/deepTools)
* [DANPOS3](https://github.com/boenc28-cmyk/DANPOS)

## Usage
The following batch script can be used to execute the read processing pipeline. The sample argument should be a comma-separated list with no spaces between names and should not include *.fastq* extension.
```
$ ./readProcessing_batch.sh /path/to/sample/dirs/ SAMPLE_1,SAMPLE_2,SAMPLE_3
```

The following batch scripts can be used to generate bigwig files for each sample or merged across replicates. Adapt the script to specify the effective genome size.
```
$ ./makeBigwig.sh /path/to/sample/dirs/ SAMPLE_1,SAMPLE_2,SAMPLE_3
$ ./makeMergedBigwig.sh /path/to/replicate/dirs/ REPLICATE_1,REPLICATE_2,REPLICATE_3
```

Peakcalling...

## Output
The following output files will be saved to */path/to/sample/dirs/*:
* alignment, deduplication and replicate reproducibility summary files
* coordinate-sorted and deduplicated *.sam* file for each sample
* coordinated-sorted *.bam* and *.bai* files for each sample and merged replicates
* *.bigwig* files for each sample and merged replicates
  
## Credits
This pipeline was assembled by [Jessica Taylor](https://github.com/jessica-a-taylor).
