#!/bin/bash

BAM="$1"
THREADS="${2:-48}"

while read -r CHUNK CHROMS; do
    samtools view -@ $THREADS -b $BAM $CHROMS -o $(basename $BAM .bam).${CHUNK}.bam && samtools index -@ $THREADS $(basename $BAM .bam).${CHUNK}.bam &
done < /projects/DGX_BASECALLING/aligned/chromosome_chunks.txt
wait
