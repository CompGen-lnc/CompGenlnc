.. CompGenlnc documentation master file, created by
   sphinx-quickstart on Wed Nov 18 16:47:36 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========
CompGenlnc
==========

CompGenlnc processes the lncRNA and miRNA sequences and 2D structures to extract
feature vectors from the interactions between these molecules considering the
original chains and those derived from the binding zone. It uses the miRanda
algorithm to predict the binding zone of interactions and provide a database
to train machine learning models to predict the interaction probability between
a lncRNA and a miRNA or the implicatation that the interaction could have.

.. toctree::
   :caption: Contents
   :maxdepth: 2
   
   Quickstart <quickstart>
   API reference <api/index>

About
-----

This package has been developed by the CompGen-lnc team at the `Costa Rica Institute of Technology <https://www.tec.ac.cr>`__
(TEC). The source code is hosted on `GitHub <https://github.com/CompGen-lnc/CompGenlnc>`__,
you are welcome to contribute to the project and look at the codebase.

