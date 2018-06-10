#!/bin/bash

while IFS='' read -r samp || [[-n "$samp"]]; do
	echo "Processing sample ${samp}"
	salmon quant -i transcriptome_index -l A \
			-1 fastq/${samp}_1.fastq.gz \
			-2 fastq/${samp}_2.fastq.gz \
			-p 4 -o quants/${samp}_quant \
			--gcBias --biasSpeedSamp 5
done < sample_list.txt
