#!/bin/bash

# Submit one cellranger job per library,
# where cellranger maps the fastq files to a reference
# and calls the feature-barcode matrices.

project_path=$(pwd)
jobs_path="${project_path}/jobs"
cd "$jobs_path"
for dir in ./*; do
  cd "$dir"
  sbatch "${dir}.cmd"
  cd "$jobs_path"
done

