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

# Check if the poolsheet exists.
if not os.path.exists(config["inputs"]["poolsheet_path"]):
    logger.info("Error, the poolsheet file does not exist.\n\nExiting.")
    exit("MissingPoolSheetFile")

# Check if the references exists.
for reference in ["azimuth", "hierscpred"]:
    run_method = config["settings"]["run_" + reference]
    reference_file = config["refs_extra"][reference + "_reference"]
    if config["settings"]["run_" + reference] and (reference_file is None or not os.path.exists(config["refs"]["ref_dir"] + reference_file)):
        logger.info("Error, could not find the {} file. Please check that the file exists.\n\nExiting.".format(config["refs"]["ref_dir"] + reference_file))
        exit("MissingReferenceFile")

# Check if at least one method is used.
if not config["settings"]["run_azimuth"] and not config["settings"]["run_hierscpred"]:
    logger.info("Error, expect at least one classification method to be selected.\n\nExiting.")
    exit("NoClassificationMethod")

logger.info("Loading the input poolsheet")
POOL_DF = pd.read_csv(config["inputs"]["poolsheet_path"], sep="\t", dtype=str)
POOL_DF.fillna("NA", inplace=True)
POOL_DF.index = POOL_DF["Pool"]

required_columns = {
    'Pool': (True, False, False),
    'Counts': (True, False, True),
    'Barcodes': (True, False, True)
}
if config["inputs"]["query_rds"] is not None:
    required_columns = {
        'Pool': (True, False, False)
    }

# Check for missing columns.
missing_columns = [column for column in required_columns.keys() if not column in POOL_DF.columns]
if len(missing_columns) > 0:
    logger.info("\tError, missing columns {} in poolsheet file for the selected methods.".format(", ".join(missing_columns)))
    exit()

# Check if the columns are valid.
poolsheet_is_valid = True
for column, (must_be_unique, must_be_numeric, must_exist) in required_columns.items():
    if must_be_unique and not POOL_DF[column].is_unique:
        logger.info("\tYour {} column contains duplicates, please make sure all values are unique.".format(column))
        poolsheet_is_valid = False

    if must_be_numeric:
        for value in POOL_DF[column]:
            if not str(value).isnumeric():
                logger.info("\tYour {} column is not numeric, please make sure there are only numeric values in this column.".format(column))
                poolsheet_is_valid = False
                break

    if must_exist:
        for fpath in POOL_DF[column]:
            if not os.path.exists(fpath) or not os.path.isfile(fpath):
                logger.info("\tYour {} column contains a file {} that does not exist, please make sure all input files exist.".format(column, os.path.basename(fpath)))
                poolsheet_is_valid = False
                break

if not poolsheet_is_valid:
    logger.info("\n\nExiting.")
    exit("InvalidPoolSheet")

logger.info("\tValid.")
POOLS = POOL_DF["Pool"]

if not config["settings"]["split"]:
    POOLS = ["all"]

input_files = []
if config["settings"]["run_azimuth"]:
    input_files.append(config["outputs"]["output_dir"] + "visualise/visualise.azimuth.done")

if config["settings"]["run_hierscpred"]:
    input_files.append(config["outputs"]["output_dir"] + "visualise/visualise.hier_scpred.done")

if config["settings"]["run_azimuth"] and config["settings"]["run_hierscpred"]:
    input_files.append(config["outputs"]["output_dir"] + "compare/compare.done")

rule all:
    input:
        files = input_files

# Import individual rules
include: "includes/classification.smk"
