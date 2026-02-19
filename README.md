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

If Anaconda is not installed on your system, then install it using the link below.

https://repo.anaconda.com/archive/

### Create and activate a conda environment:

```
conda create -n stitchr python=3.10 -y
conda activate stitchr
```

### Install STITCHR:

For more information about STITCHR - https://jamieheather.github.io/stitchr

```
pip install -r requirements.txt
```

### Download gene data for STITCHR:

```
stitchrdl -s human
```

### Test the installation:

#### Test STITCHR
```
stitchr -h
thimble -h
```
#### Test build_TCR_HDRTs.py 
```
python build_TCR_HDRTs.py -h 
```

# Example Usage

### To generate TCR HDR templates:

#### Run using D112K TRBC variant
For additional details - https://www.nature.com/articles/s41587-024-02531-6
```
conda activate stitchr
python build_TCR_HDRTs.py --use_D112K example/example_input.tsv -o example_results
```

#### Run using wildtype TRBC variant
```
conda activate stitchr
python build_TCR_HDRTs.py example/example_input.tsv -o example_results
```

### Input Data

Input files should be provided as tab-separated tables containing TCR sequence information compatible with STITCHR formatting (https://jamieheather.github.io/stitchr).

An example input file is provided in: example/example_input.tsv


