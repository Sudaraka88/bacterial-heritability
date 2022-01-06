#!/bin/bash
# Running this may take a long time, the 4 output tree files are provided

# cd/mic_isolates.fasta - isolates in fasta format
# CD/MICtree.treefile - output from iqtree in newick format


# Build initial trees
iqtree -s cd_isolates.fasta -T AUTO -m GTR+G --prefix CDtree_
iqtree -s mic_isolates.fasta -T AUTO -m GTR+G --prefix MICtree_

# ClonalFrameML 
ClonalFrameML CDtree.treefile cd_isolates.fasta cfrm_cd_
ClonalFrameML MICtree.treefile mic_isolates.fasta cfrm_MIC_
