#!/bin/bash

# Load required modules
module load PYTHON/2.7.5
module load lims/1.2

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project

subproject="POST_12"
/scratch/project/production/DAT/apps/LIMSQ/limsq -sp $subproject | sed 's/;/\t/g' > data/${subproject}_info.tsv
