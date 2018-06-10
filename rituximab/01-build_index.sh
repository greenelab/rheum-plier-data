#!/bin/bash

# obtain transcriptome from Ensembl
curl ftp://ftp.ensembl.org/pub/release-92/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o hsapiens.fa.gz
curl ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz -o Homo_sapiens.GRCh38.92.gtf.gz
gunzip Homo_sapiens.GRCh38.92.gtf.gz

# build a salmon index -- these are 50bp reads, so we'll set k to 23
salmon index \
	-t hsapiens.fa.gz \
	-i transcriptome_index \
	--type quasi -k 23
