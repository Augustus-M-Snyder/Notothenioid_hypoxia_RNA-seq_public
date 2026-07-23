#!/bin/bash

for bam in ./*_align.s.bam; do
    samtools index "$bam"
done

