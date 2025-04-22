#! /usr/bin/env python

"""
Author: "Laura Jiménez Gracia"
metrics_dfe: 2021-07-26

This script merges performance metrics of cellranger multi [GEX with CellPlex technology]
for all libraries (GEM ids) in this subproject into 
a single file cellranger_mapping_metrics_CellPlex.csv located in /results folder.
"""

# Load packages
import os
import pandas as pd
import numpy as np
import config as cfg

# Get metrics_dfa paths
subproject = cfg.subproject
subproject_path = cfg.subproject_path
metadata_path = cfg.metadata_path
CMO_path = cfg.CMO_path

# Define list of libraries to merge 'metrics'
metadata_df = pd.read_csv(metadata_path, sep=",", header=0)
mask = (metadata_df["subproject"] == subproject)
libraries = metadata_df.loc[mask, "gem_id"]
libraries = list(set(libraries))

CMO_df = pd.read_csv(CMO_path, sep=",", header=0)

## FOR GEX + CellPlex libraries
summary_gex_dfs = pd.DataFrame()
summary_cmo_dfs = pd.DataFrame()
summary_demux_dfs = pd.DataFrame()
summary_sample_dfs = pd.DataFrame()

for gem_id in libraries:
    library_path=f"{subproject_path}/jobs/{gem_id}/{gem_id}/outs/per_sample_outs"
    samples = os.listdir(library_path)

    for sample_name in samples:
        sample_path=f"{library_path}/{sample_name}/metrics_summary.csv"

        CMO_name = CMO_df[(CMO_df['gem_id'] == gem_id) & (CMO_df['sample_id'] == sample_name)]['CMO_id'].values[0]

        if os.path.exists(sample_path):
            metrics_df = pd.read_csv(sample_path)
            metrics_df.columns = metrics_df.columns.str.replace(' ','_')
            for i in metrics_df.columns:
                metrics_df[i] = metrics_df[i].str.replace(' ','_')

            # Subsetting GEX / CMO LIBRARY
            gex_df = metrics_df[(metrics_df['Library_or_Sample'] == 'Library') &
                (metrics_df['Library_Type'] == 'Gene_Expression')][['Metric_Name','Metric_Value']]
            gex_df = gex_df.set_index('Metric_Name')
            gex_df.index.names = [None]
            gex_df.columns = [gem_id]
            gex_df_t = gex_df.transpose()
            gex_df_t.reset_index(inplace=True)
            gex_df_t = gex_df_t.rename(columns = {'index':'gem_id'})
            summary_gex_dfs = summary_gex_dfs.append(gex_df_t)

            cmo_df = metrics_df[(metrics_df['Library_or_Sample'] == 'Library') &
                (metrics_df['Library_Type'] == 'Multiplexing_Capture') &
                ((metrics_df['Grouped_By'] == 'Fastq_ID') | 
                (metrics_df['Grouped_By'] == np.nan) | 
                (metrics_df['Grouped_By'] == 'Physical_library_ID'))][['Metric_Name','Metric_Value']]
            cmo_df = cmo_df.set_index('Metric_Name')
            cmo_df = cmo_df.set_index(cmo_df.index.astype(str) + '_CMO')
            cmo_df.index.names = [None]
            cmo_df.columns = [gem_id]
            cmo_df_t = cmo_df.transpose()
            cmo_df_t.reset_index(inplace=True)
            cmo_df_t = cmo_df_t.rename(columns = {'index':'gem_id'})
            summary_cmo_dfs = summary_cmo_dfs.append(cmo_df_t)

            # Subsetting DEMULTIPLEXING / SAMPLE
            demux_df = metrics_df[(metrics_df['Library_or_Sample'] == 'Library') &
                (metrics_df['Library_Type'] == 'Multiplexing_Capture') &
                (metrics_df['Group_Name'] == CMO_name)][['Metric_Name','Metric_Value']]
            demux_df = demux_df.set_index('Metric_Name')
            demux_df.index.names = [None]
            demux_df.columns = [CMO_name]
            demux_df_t = demux_df.transpose()
            demux_df_t.reset_index(inplace=True)
            demux_df_t = demux_df_t.rename(columns={'index': 'CMO_name'})
            demux_df_t.insert(0, "sample_name", sample_name)
            demux_df_t.insert(0, "gem_id", gem_id)
            summary_demux_dfs = summary_demux_dfs.append(demux_df_t)        

            # Subsetting SAMPLE
            sample_df = metrics_df[(metrics_df['Library_or_Sample'] == 'Sample') &
                                   (metrics_df['Library_Type'] == 'Gene_Expression')][['Metric_Name', 'Metric_Value']]
            sample_df = sample_df.set_index('Metric_Name')
            sample_df.index.names = [None]
            sample_df.columns = [sample_name]
            sample_df_t = sample_df.transpose()
            sample_df_t.reset_index(inplace=True)
            sample_df_t = sample_df_t.rename(columns={'index': 'sample_name'})
            sample_df_t.insert(0, "gem_id", gem_id)
            summary_sample_dfs = summary_sample_dfs.append(sample_df_t)

# Formatting dataframes
summary_gex_dfs.drop_duplicates(subset = "gem_id", inplace=True)
summary_gex_dfs.drop(
    columns=['Number_of_reads'], 
    inplace=True)
summary_cmo_dfs.drop_duplicates(subset = "gem_id", inplace=True)
summary_cmo_dfs = summary_cmo_dfs.T.drop_duplicates().T

summary_libraries = summary_gex_dfs.merge(summary_cmo_dfs, how = "left", on="gem_id")


summary_demultiplexing = summary_demux_dfs.merge(summary_sample_dfs, how="left", left_on=['gem_id','sample_name'], right_on=['gem_id','sample_name'])
summary_gex_dfs.drop(
    columns=["Cell-associated_barcodes_identified_as_multiplets", 
    "Cell-associated_barcodes_not_assigned_any_CMOs"], 
    inplace=True)
    
# Save combined metrics
# Create results directory
results_directory = "{}/results".format(subproject_path)
if not os.path.exists(results_directory):
    os.mkdir(results_directory)

if not summary_libraries.empty:
    summary_libraries.to_csv(os.path.join( f"{subproject_path}/results/cellranger_mapping_metrics_libraries.csv"), header = True, index = None)

if not summary_demultiplexing.empty:
    summary_demultiplexing.to_csv(os.path.join( f"{subproject_path}/results/cellranger_mapping_metrics_demultiplexing.csv"), header = True, index = None)

