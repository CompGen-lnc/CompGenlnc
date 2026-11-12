# CompGenlnc

A Feature Vector Generation Package for lncRNA-miRNA Interaction Probability

<!---Links--->

## Introduction

CompGenlnc is designed to generate feature vectors for lncRNA-miRNA interactions
from the lncRNA and miRNA sequencesand 2D structures, and the subsequences and
2D substructures obtained from the interactions.
The package works as a pipeline, where the steps are as follows:

1. Select base files with the sequences, 2D structures, interactions and
   precursor mature miRNA associations.
3. Calculate missing base files and interaction binding zones.
4. Generate feature vectors for each sequence, 2D structure, subsequence and 2D
   substructure.
5. Index the feature vectors by interactions.

## Installation

### Requirements

Before installing, you must be sure to have installed miRanda from Anaconda. For
this you will need a Unix-like system (Linux or macOS), if you are using
Windows, you first need to have installed WSL and install miRanda there.

For this, you can install miRanda from Anaconda by running the following command
in your terminal:

```bash
conda install bioconda::miranda
```

Be sure to be installed

### Basic Installation

Install the package from PyPI (not yet available, but will be available soon):

```bash
pip install compgenlnc
```

As there are some dependencies out of PyPI, you may need to install them
manually. Instead, use the next command to install all dependencies:

```bash
compgenlnc-install-deps
```

### Development Installation

Clone the repository on your local machine:

```bash
git clone https://github.com/CompGen-lnc/CompGenlnc.git
cd CompGenlnc
```

Before installing the package, it is recommended to create and activate a
virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

Lastly, install the package in edit mode and its dependencies:

```bash
# In the directory of the project
pip install -e .

# For dynamic dependencies
compgenlnc-install-deps
```