#! /usr/bin/env python3

"""
This script initializes the filesystem of this subproject [cellranger multi 3'GEX+CMO (HTO multiplexed)]:
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cmoreference

It creates a "jobs" folder with as many subdirectories as samples it has
For each sample directory, it creates the following files/folders:
1. fastq_paths.csv -- file containing fastq paths
2. jobs -- dir with libraries.csv and feature_reference.csv files
    2.1. feature_reference.csv -- file indicating correspondence between HTO and sample_id  [if hashed libraries]
    2.2. (sample_id)_config.csv -- configuration file needed to run cellranger multi
    2.3. (sample_id).cmd -- jobscript to compute the features-barcode matrix with cellranger (SC3Pv3 chemistry)
3. fastq -- dir with the symlinks pointing to the fastq files
4. logs -- dir which contains standard error and output of cellranger
"""

# Load packages
import os
import subprocess
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import config as cfg


def main():
    # Read data
    subproject = cfg.subproject
    subproject_path = cfg.subproject_path

    # Metadata -- manually curated metadata file
    metadata_df = pd.read_csv(cfg.metadata_path, sep=",", header=0)
    mask = (metadata_df["subproject"] == subproject)
    libraries = metadata_df.loc[mask, "library_id"]
    libraries = list(libraries)
    # LIMS -- data from Illumina sequencing
    lims = pd.read_csv(cfg.infofile_path, sep="\t", header=0)
    lims = lims.loc[lims.id.isin(libraries)]
    # Feature_reference -- manually curated metadata file
    if cfg.CMO_ref_path != None:
        CMO_df = pd.read_csv(cfg.CMO_path, sep=",", header=0)

    # Assemble fastq paths combining flowcell, lane and index
    fastq_path = cfg.fastq_path
    fastq_path_list_r1 = []
    fastq_path_list_r2 = []
    for idx in lims.index:
        fc = lims.loc[idx, "flowcell"]
        lane = lims.loc[idx, "lane"]
        index = lims.loc[idx, "index"]
        fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(
            fastq_path, fc, lane, fc, lane, index)
        fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(
            fastq_path, fc, lane, fc, lane, index)
        fastq_path_list_r1.append(fastq_path_r1)
        fastq_path_list_r2.append(fastq_path_r2)
    library_id_l = list(lims["id"].append(lims["id"]))
    p_l = "P" * len(fastq_path_list_r1)
    indx_l = list(range(1, len(fastq_path_list_r1) + 1))
    pair_id = [p_l[x] + str(indx_l[x]) for x in range(len(indx_l))]
    fastq_path_list_r1.extend(fastq_path_list_r2)
    pair_id.extend(pair_id)
    fastq_path_l = fastq_path_list_r1
    read_l = (["R1"] * lims.shape[0]) + (["R2"] * lims.shape[0])
    fastq_dict = {"library_id": library_id_l, "fastq_path": fastq_path_l, "read": read_l, "pair_id": pair_id}
    fastq_df = pd.DataFrame(fastq_dict)

    # Write dataframe to csv
    fastq_df.to_csv("{}/data/fastq_paths.csv".format(subproject_path), header=True, index=False)

    # Define paths
    if cfg.reference == "human":
        ref_path = cfg.Hsapiens_path
    elif cfg.reference == "mouse":
        ref_path = cfg.Mmus_path


    # Initialize directory
    jobs_path = "{}/jobs".format(subproject_path)
    if not os.path.exists(jobs_path):
        os.mkdir(jobs_path)

    # For each GEM sample, create directories and jobscript
    gem_id_list = np.unique(metadata_df["gem_id"])
    for gem_id in gem_id_list:
        # Define and create directories
        gemid_dir = "{}/{}".format(jobs_path, gem_id)
        fastq_dir = "{}/fastq".format(gemid_dir)
        logs_dir = "{}/logs".format(gemid_dir)
        for new_dir in [gemid_dir, fastq_dir, logs_dir]:
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)

        # Define variables and subset dataframes
        filt = (metadata_df["gem_id"] == gem_id)
        metadata_df_filt = metadata_df.loc[filt]
        lib_hashing = metadata_df_filt.loc[filt, "hashing"]
        lib_hashing = list(set(lib_hashing))[0]
        library_id = metadata_df_filt.loc[filt, "library_id"]
        fastq_sub_df = fastq_df.loc[fastq_df["library_id"].isin(library_id), :]
         
        for lib in library_id:
            # Create symmlinks to fastq files
            lib_type = metadata_df_filt.loc[metadata_df_filt["library_id"] == lib, "type"]
            lib_type = lib_type.values[0] # LIB_TYPE can be cDNA, HTO, CPS, CMO

            # Creating folders for each library type
            if not os.path.exists("{}/{}".format(fastq_dir, lib_type)):
                os.mkdir("{}/{}".format(fastq_dir, lib_type))

            symlink_path = "{}/{}".format(fastq_dir, lib_type)
            create_fastq_symlink(gem_id, lib, lib_type, fastq_sub_df, symlink_path)

        # List libraries information: file indicating fastq_id, fastq_path and feature_type
        libraries_data = lists_libraries_for_config_file(gem_id, gemid_dir)

        # Create config.csv file
        # [only for hashed samples] Create feature_reference.csv: file indicating correspondence between HTO and sample id
        if lib_hashing == "CellPlex":
            filt_row = ((CMO_df["subproject"] == subproject) & (CMO_df["gem_id"] == gem_id))
            filt_col = np.invert(CMO_df.columns.isin(["subproject", "gem_id"]))
            CMO_filt = CMO_df.loc[filt_row, filt_col]
            CMO_data = lists_CMOinfo_for_config_file(CMO_filt)
            expect_cells = 20000
            CMO_path = cfg.CMO_ref_path
            write_cellranger_config_csv_hashed(gem_id, gemid_dir, ref_path, CMO_path, expect_cells, libraries_data, CMO_data) 

        elif lib_hashing == "not_CellPlex":
            expect_cells = 6000
            write_cellranger_config_csv_not_hashed(gem_id, gemid_dir, ref_path, expect_cells, libraries_data) 
            
        # Create cellranger script
        make_cellranger_clustermode_cmd(gem_id, gemid_dir)





# FUNCTIONS

def create_fastq_symlink(gem_id, library, lib_type, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation for cell-hashed samples.
    Keyword arguments:
        gem_id -- identifier of the Gelbeads-in-Emulsion (GEM)
        library -- library id
        lib_type -- type of library (cDNA or HTO)
        fastq_path_df -- pandas dataframe with the fastq paths for that gem_id
        symlink_path -- string specifying where to create the symlinks (fastq/xxx)
    Returns:
        None
    """
    fastq_path_sub = fastq_path_df.loc[fastq_path_df["library_id"] == library, :]
    pair_ids = np.unique(fastq_path_sub["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            lane = str(i + 1)
            symlink_path_lane = "{}/lane{}".format(symlink_path, lane)
            if not os.path.exists(symlink_path_lane):
                os.mkdir(symlink_path_lane)
            fastq_path = pair_df.loc[j, "fastq_path"]
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            if lib_type == "cDNA":
                gem_id_sp = gem_id
            else:
                gem_id_sp = "{}_{}".format(gem_id, lib_type)

            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path_lane, gem_id_sp, lane, read)])


def lists_libraries_for_config_file(gem_id, gem_id_path):
    """Creates a string (in csv format) with library information required by cellranger multi.
        Each line corresponds with a new library, and contains fastq_id,fastq_path,feature_types
    Keyword arguments:
        gem_id -- identifier of the Gelbeads-in-Emulsion (GEM) well
        gem_id_path -- absolute path to gem_id-specific directory.
    Returns:
        libraries_data -- A string in csv format containing library information
    """
    fastq_dirs = os.listdir("{}/fastq".format(gem_id_path))
    libraries_data = ""
    for d in fastq_dirs:
        if d == "cDNA":
            gem_id_sp = gem_id
            feature_type = "Gene Expression"
        else:
            gem_id_sp = "{}_{}".format(gem_id, d)
            feature_type = "Multiplexing Capture"

        fastq_sub_dirs = "{}/fastq/{}".format(gem_id_path, d)
        libraries_data = libraries_data + "{},{},{}\n".format(gem_id_sp, fastq_sub_dirs, feature_type)
    return libraries_data


def lists_CMOinfo_for_config_file(CMO_filt):
    """Creates a string (in csv format) with CMO information required by cellranger multi.
        Each line corresponds with a multiplexed sample, and contains sample_id,cmo_ids, 
        where CMO is the cell multiplexing oligo IDs used to multiplex this sample (separated by |). 
    Keyword arguments:
        CMO_filt -- Cell multiplexing dataframe information
    Returns:
        libraries_data -- A string in csv format containing library information
    """
    sample_id_list =  CMO_filt["sample_id"].tolist()
    CMO_id_list = CMO_filt["CMO_id"].tolist()
    CMO_data = ""
    for idx, sample_id in enumerate(sample_id_list):
        # This CMO_data will work unless multiple CMO were used to tag the same sample! 
        # (e.g. overloading with unique condition)
        CMO_data = CMO_data + "{},{}\n".format(sample_id, CMO_id_list[idx])
    return CMO_data


def write_cellranger_config_csv_not_hashed(gem_id, gem_id_path, ref_path, expect_cells, libraries_data):
    """Creates the file "config.csv" which is required by cellranger multi for hashed samples.
    Keyword arguments:
        gem_id -- identifier of the Gelbeads-in-Emulsion (GEM) well
        gem_id_path -- absolute path to gem_id-specific directory.
        ref_path -- absolute path to reference genome
        expect_cells -- expected number of recovered cells
        libraries_data -- string (in csv format) containing library information
    Returns:
        None
    """
    config_file = open("{}/{}_config.csv".format(gem_id_path, gem_id), "w")
    config_info = f"""[gene-expression]
reference,{ref_path}
expect-cells,{expect_cells}
no-secondary,true
no-bam,true
[libraries]
fastq_id,fastqs,feature_types
{libraries_data}
"""
    config_file.write(config_info)
    config_file.close()


def write_cellranger_config_csv_hashed(gem_id, gem_id_path, ref_path, CMO_path, expect_cells, libraries_data, CMO_data):
    """Creates the file "config.csv" which is required by cellranger multi for hashed samples.
    Keyword arguments:
        gem_id -- identifier of the Gelbeads-in-Emulsion (GEM) well
        gem_id_path -- absolute path to gem_id-specific directory.
        ref_path -- absolute path to reference genome
        expect_cells -- expected number of recovered cells
        libraries_data -- string (in csv format) containing library information
    Returns:
        None
    """
    config_file = open("{}/{}_config.csv".format(gem_id_path, gem_id), "w")
    config_info = f"""[gene-expression]
reference,{ref_path}
cmo-set,{CMO_path}
expect-cells,{expect_cells}
no-secondary,true
no-bam,true
[libraries]
fastq_id,fastqs,feature_types
{libraries_data}
[samples]
sample_id,cmo_ids
{CMO_data}
"""
    config_file.write(config_info)
    config_file.close()


def make_cellranger_clustermode_cmd(gem_id, gem_id_path):
    """Creates a jobscript with a call to cellranger multi for a GEM well
    to be run in a cluster mode!
    Keyword arguments:
        gem_id -- identifier of the Gelbeads-in-Emulsion (GEM) well
        gem_id_path -- path to save the jobscript
    Returns:
        None
    """
    job_script_file = open("{}/{}.cmd".format(gem_id_path, gem_id), "w")
    job_script = """#!/usr/bin/env bash

#SBATCH --job-name="{}"
#SBATCH --workdir=.

#SBATCH --error=./logs/slurm_%x_%J.err
#SBATCH --output=./logs/slurm_%x_%J.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=normal
#SBATCH --partition=genB,main

#SBATCH --mail-type=END       
#SBATCH --mail-user=laura.jimenez@cnag.crg.eu

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

ulimit -n 16000
export HDF5_USE_FILE_LOCKING="FALSE"

{} multi --id {} \\
    --csv {}/{}_config.csv \\
    --localcores 2 \\
    --jobmode {}

echo [`date "+%Y-%m-%d %T"`] finished job
""".format(
    gem_id, cfg.cellranger_path, gem_id, gem_id_path, gem_id, cfg.slurmtemplate_path)
    job_script_file.write(job_script)
    job_script_file.close()



if __name__ == '__main__':
    exit(main())

