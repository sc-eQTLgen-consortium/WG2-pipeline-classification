#!/usr/bin/env python
import pandas as pd
import os

# Add trailing /.
if not config["inputs"]["repo_dir"].endswith("/"):
    config["inputs"]["repo_dir"] += "/"
if not config["refs"]["ref_dir"].endswith("/"):
    config["refs"]["ref_dir"] += "/"
if not config["outputs"]["output_dir"].endswith("/"):
    config["outputs"]["output_dir"] += "/"

# Check if the singularity image exists.
if not os.path.exists(config["inputs"]["singularity_image"]):
    logger.info("Error, the singularity image does not exist.\n\nExiting.")
    exit("MissingSIFFile")

# Check if the query RDS exists.
if not os.path.exists(config["inputs"]["query_rds"]):
    logger.info("Error, the query RDS does not exist.\n\nExiting.")
    exit("MissingQueryRDSFile")

BATCH = None
BATCHES = ["All"]
if config["inputs"]["batch_info"] is not None:
    logger.info("Loading the batch info")
    BATCH_DF = pd.read_csv(config["inputs"]["batch_info"], sep="\t", dtype=str)
    if BATCH_DF.shape[1] != 1:
        logger.info("\tError, Expect one column in {} file.\n\nExiting.".format(config["inputs"]["batch_info"]))
        exit("InvalidBatchInfoFile")
    if BATCH_DF.shape[0] <= 1:
        logger.info("\tWarning, Expect more than one row in {} file.\n\tIgnoring batch_info.".format(config["inputs"]["batch_info"]))
    else:
        BATCH = BATCH_DF.columns[0]
        BATCHES = BATCH_DF.iloc[:, 0].tolist()

# Check if the references exists.
for reference in ["azimuth", "hierscpred"]:
    if config["settings"][reference] and (not os.path.exists(config["refs"]["ref_dir"] + config["refs_extra"][reference]) or config["refs_extra"][reference] == ""):
        logger.info("Error, could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["refs"]["ref_dir"] + config["refs_extra"]["azimuth"]))
        exit("MissingReferenceFile")

# Check if at least one method is used.
if not config["settings"]["azimuth"] and not config["settings"]["hierscpred"]:
    logger.info("Error, expect at least one classification method to be selected.\n\nExiting.")
    exit("NoClassificationMethod")

def get_input_files(wildcards):
    if config["settings"]["azimuth"] and config["settings"]["hierscpred"]:
        return config["outputs"]["output_dir"] + "compare/comp_contingency_table.tsv"
    elif len(BATCHES) == 1 and config["settings"]["azimuth"]:
        return expand(config["outputs"]["output_dir"] + "map/azimuth_{batch}.RDS", batch=BATCHES)
    elif len(BATCHES) == 1 and config["settings"]["hierscpred"]:
        return expand(config["outputs"]["output_dir"] + "map/hier_scpred_{batch}.RDS", batch=BATCHES)
    else:
        return config["outputs"]["output_dir"] + "reduce/reduced_data.RDS"


rule all:
    input:
        files = get_input_files

# Import individual rules
include: "includes/classification.smk"
