####################################################################
##                                                                ##
##  Script to run CITE-seq-Count on the cluster for one sample    ##
##  of the project IFNa-BM-Niche                            	  ##
##                                                                ##
####################################################################

#!/usr/bin/bash

module load python/3.6.1 


# File locations
#CITE-seq-Count=/home/gruensch/.local/bin/CITE-seq-Count
READ1=/dir/X_S1_L001_R1_001.fastq.gz
READ2=/dir/X_S1_L001_R2_001.fastq.gz
TAGS=/dir/CITEseqCount_input/tags.csv
WHITELIST=/dir/CITEseqCount_input/barcodes.csv
OUTPATH=/dir/CITE_Seq_Count/$1

# Barcoding parameters
CB_FIRST=1
CB_LAST=16
UMI_FIRST=17
UMI_LAST=28
EXPTECTED_CELLS=10000
N_ERRORS_BC=1
N_ERRORS_UMI=1
MAX_ERROR_HTO=2
N_BASES=0

/.local/bin/CITE-seq-Count \
    -R1 $READ1 \
    -R2 $READ2 \
    -t $TAGS \
    -cbf $CB_FIRST \
    -cbl $CB_LAST \
    -umif $UMI_FIRST \
    -umil $UMI_LAST \
    -cells $EXPTECTED_CELLS \
    --bc_collapsing_dist $N_ERRORS_BC \
    --umi_collapsing_dist $N_ERRORS_UMI \
    -wl $WHITELIST \
    --max-error $MAX_ERROR_HTO \
    -trim $N_BASES \
    -o $OUTPATH
