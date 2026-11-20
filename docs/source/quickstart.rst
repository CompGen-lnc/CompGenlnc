==========
Quickstart
==========

Install the package
-------------------

To use the package, you need to install it first. You can do this using pip:

.. code:: bash

  pip install compgenlnc

As there are some dependencies that cannot be installed normally using pip, you
have to install them separately. For this, you can run the following:

.. code:: bash

  compgenlnc-install-deps

Additionally, you must have installed `miRanda <https://anaconda.org/channels/bioconda/packages/miranda/overview>`__
from Anaconda to use the full functionality of the package. For this, install
conda and run the following:

.. code:: bash

  conda install bioconda::miranda

If you are using Windows, you can use WSL to run the command above.

To verify the installation, run the following in a Python shell or script:

.. code:: python

  import compgenlnc as cgl
  print(cgl.__version__)

This should print the version of the package. If a warning is displayed, it
means that the dependencies were not installed correctly.

Basics
------

The package provides a set of functions to load the source data and to generate
feature vectors for lncRNA-miRNA interactions. The main functions are
for loading and calculating the source data, in addition to generating and
loading feature vectors.

Base files
^^^^^^^^^^

The package store and load all base files before generating feature vectors.
The ones that user provide are:

- ``lncRNA.fa``: a FASTA file with the sequences of the lncRNAs.
- ``lncRNA_2D.fa``: a FASTA file with the 2D structures of the lncRNAs. This
  can be generated from ``lncRNA.fa``.
- ``miRNA.fa``: a FASTA file with the sequences of the precursor miRNAs.
- ``miRNA_2D.fa``: a FASTA file with the 2D structures of the precursor miRNAs.
  This can be generated from ``miRNA.fa``.
- ``mature_miRNA.fa``: a FASTA file with the sequences of the mature miRNAs.
- ``data_miRNA.csv``: a CSV file with the mature and precursor miRNA
  associations.
  This can be generated from ``miRNA.fa`` and ``mature_miRNA.fa``.
- ``interactions.csv``: a CSV file with the interactions between, including the
  RNA IDs and a mark if they are likely or not. This can be generated from
  ``postive_pairs.csv`` and ``negative_pairs.csv``.
  
  - ``postive_pairs.csv``: a CSV file with the positive interactions.
  - ``negative_pairs.csv``: a CSV file with the negative interactions.

There are also files that the package generate from base files provided by the
user:

- seq+ss files: files with the sequences and 2D structures of lncRNAs and
  precursor miRNAs.
- loops files: files with the loops of lncRNAs and precursor miRNAs.
- miranda files: files with the predicted binding zones of interactions.
- binding zone files: files with the selected binding zones of interactions.
- ``index.npy``: a numpy array with the indices of the interactions.

Using the package
^^^^^^^^^^^^^^^^^

The package provide many functions that can be used to generate the database
from the interactions between lncRNAs and miRNAs. For loading the source_data
you can do:

.. code:: python

  import compgenlnc as cgl

  # Loading files from a same directory
  dir_with_files = "path/to/files" # Files should have the same name
  cgl.load_source_data(dir_with_files)

  # Loading files from others directory
  lnc_fasta = "path/to/lncRNA.fa"
  mirna_fasta = "path/to/miRNA.fa"
  cgl.load_source_data(
    lnc_seq_fasta=lnc_fasta,
    mir_seq_fasta=mirna_fasta,
  )

  # Generating the rest of the base files
  only_indexes = True # If you want to process only RNAs in the indexes
  keep_data = False # If you do want to recalculate the basefiles
  cgl.gen_source_data(filter_indexes=only_indexes, keep_data=keep_data)

  # Generate the database
  cgl.gen_database()

Once the database is generated, you can use it to train machine learning models
to predict the interaction probability or effects in some diseases. The package
also provides functions to get them:

>>> import compgenlnc as cgl
>>> features = cgl.get_interaction_features()
>>> type(features)
<class 'pandas.DataFrame'>
>>> features.shape
(2053, 52)
>>> feat_vectors = cgl.get_feature_vectors()
>>> type(feat_vectors)
<class 'numpy.ndarray'>
>>> feat_vectors.shape
(2053, 4838)
