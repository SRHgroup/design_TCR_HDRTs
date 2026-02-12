### Design TCR HDRTs

Tools for assembling T-cell receptor (TCR) sequences and designing HDR templates using Stitchr.

This repository provides scripts to process TCR sequence data and generate HDR template designs in a reproducible workflow.

### Overview

# The workflow:

Install Stitchr and dependencies

Prepare input TCR data

Run assembly and HDR template design

Obtain output files ready for downstream analysis or synthesis

### Quick Start

# Clone the repository:

git clone https://github.com/SRHgroup/design_TCR_HDRTs.git
cd design_TCR_HDRTs

# Create and activate a conda environment:

conda create -n stitchr python=3.10 -y
conda activate stitchr

# Install dependencies:

pip install -r requirements.txt

# Download gene data for Stitchr:

stitchrdl -s human

# Test the installation:

stitchr -h
thimble -h
Example Usage

# Run the included example:

thimble -in example/example_input.tsv -o example/my_output.tsv -r AB -s human

The file my_output.tsv should match the expected output.

# To generate HDR templates:

python build_TCR_HDRTs.py --use_D112K example/example_input.tsv -o example_results

# Input Data

Input files should be provided as tab-separated tables containing TCR sequence information compatible with Stitchr/Thimble formatting.

An example input file is provided:

example/example_input.tsv


HDR template designs

Intermediate files depending on parameters

Output is written to the directory specified with -o.
