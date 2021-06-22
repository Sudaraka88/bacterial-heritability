#! /bin/bash
NAME=${1?Error: Filename required...}

sed -n '/>/p' $NAME | cut -c2- > $NAME".seqs"

