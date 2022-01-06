#! /bin/bash
# extract the list of sequences from a fasta file
NAME=${1?Error: Filename required...}

sed -n '/>/p' $NAME | cut -c2- > $NAME".seqs"

