#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }
which gunzip >/dev/null || { echo "gunzip: command not found"; exit 1; }

mkdir -p data
cd data || exit


# 9bc996d9e8a81016869aa948064599ef  query.RDS
# b2976e2e623db2637b4ed2b8fa78c894  reference.RDS
# 885a5a8fa6ce3cdc1a7b3568c2076efe  cell_classification.sif
echo "Downloading WP2_singularity.zip"
wget https://www.dropbox.com/sh/ekvdocei8r45jxq/AACA17z7PFNbVkeuSavFvjPNa --content-disposition \
  && md5sum -c - <<<"570933608d60d5b26b91a327b71b4c5b  WP2_singularity.zip" \
  && unzip WP2_singularity.zip \
  || rm WP2_singularity.zip
