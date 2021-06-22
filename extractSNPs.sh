#!/bin/bash
NAME=${1?Error: Species required...}
samtools faidx maela3K.fasta $NAME | perl -pe 'chomp unless /^>/' | tail -1
