# Changelog

All notable changes implemented in this branch compared to the main branch of [WG2-pipeline-classification](https://github.com/sc-eQTLgen-consortium/WG2-pipeline-classification) (2021-11-25, SHA: `9ae2a7565da81fec9c2f791ca1cb4dd291412969`) are documented in this file. 

Note that this branch is in beta and version 2.0.0 is not yet ready for release.

## [2.0.0] - Classification - 2023-11-24

#### Additions
- Build snakemake, similar to other workgroups, around pre-existing scripts
- Added `map_azimuth.R` script for non CITE-seq annotation

#### Fixes
- Fixed issue in split.R where it looks for unknown variable `xaxis` instead of `opt$batch'`
- Fixed issue where a wrong calculation for `future.globals.maxSize` was used

#### Changes
- Aligned `README.md` with the other workgroups
- Moved reference files out of the image and changed the paths to command line arguments
- Moved scripts outside of image for flexibility
- Switched from singularity to Docker
- Updated all software versions
