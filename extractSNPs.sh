#!/bin/bash
NAME=${1?Error: Species required...}
# maela3K.fasta - maela alignment in fasta format, extract the whole sequence with label $NAME 

samtools faidx maela3K.fasta $NAME | perl -pe 'chomp unless /^>/' | tail -1
