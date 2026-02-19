# Design TCR HDRTs

Tools for assembling T-cell receptor (TCR) sequences and designing HDR templates using Stitchr.

This repository provides scripts to process TCR sequence data and generate HDR template designs in a reproducible workflow.

# Quick Start

### Clone the repository:

Open the terminal and move to the location where you want the program to be.

```
cd /path/to/folder
```

Then clone the repository using the following command: 

```
git clone https://github.com/SRHgroup/design_TCR_HDRTs.git
cd design_TCR_HDRTs
```

### Install Anaconda:

If you do not have Anaconda installed, then install it using the link below. 

https://repo.anaconda.com/archive/

### Create and activate a conda environment:

```
conda create -n stitchr python=3.10 -y
conda activate stitchr
```

### Install dependencies:

```
pip install -r requirements.txt
```

### Download gene data for Stitchr:

```
stitchrdl -s human
```

### Test the installation:

#### Test Stitchr
```
stitchr -h
thimble -h
```
#### Test the build_TCR_HDRTs.py design script
```
python build_TCR_HDRTs.py -h 
```

# Example Usage

### To generate TCR HDR templates:

If the conda environment is not already active, then activate conda environment `conda activate stitchr`

```
python build_TCR_HDRTs.py --use_D112K example/example_input.tsv -o example_results
```

### Input Data

Input files should be provided as tab-separated tables containing TCR sequence information compatible with Stitchr/Thimble formatting.

An example input file is provided in: example/example_input.tsv


