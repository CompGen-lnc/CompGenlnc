========
Database
========

The package generate feature vectors for lncRNA-miRNA interactions from the
sequences, 2D structures, subsequences and 2D substructures. The database
structure is organized as follows.


Sequence Characterization
-------------------------
.. csv-table:: Sequence Features
   :file: data/sequence_table.csv
   :header-rows: 1
   :align: center

.. [#] A→1, C→2, G→3, U→4
.. [#] N = mean RNA length in the set

k-mers
^^^^^^

K-mers are fixed-length subsequences (length k) extracted using a sliding window
of step 1 across RNA sequences. For each :math:`k` value, the frequency of each
possible tuple from the alphabet {A, U, G, C} is normalized by the total number
of available positions. The normalized frequency of a k-mer :math:`w` in
sequence :math:`S` of length :math:`L` is:

.. math::
   f(w, S) = \frac{|{i : S[i:i+k] = w}|}{L - k + 1}

Frequency vectors for each order k are concatenated in ascending order to form
the final descriptor. With a four-nucleotide alphabet, this produces:

- 4\ :sup:`1`  = 4 monomers
- 4\ :sup:`2`  = 16 dinucleotides
- 4\ :sup:`3`  = 64 trinucleotides
- 4\ :sup:`4`  = 256 tetranucleotides

CTD
^^^

The CTD descriptor encodes physicochemical properties of the sequence in a
vector. The procedure groups each position in the sequence into one of two
classes (C :sub:`1` or C :sub:`2`) based on three biochemical properties of the
nucleotide at that position:

.. list-table::
   :header-rows: 1

   * - Property
     - C :sub:`1`
     - C :sub:`2`
   * - Ring Structure	
     - Purines (A, G)
     - Pyrimidines (U, C)
   * - Functional Group
     - Amino (A, C)
     - Keto (G, U)
   * - H-bond Strength
     - Weak bond (A, U; 2 H-bonds)
     - Strong bond (G, C; 3 H-bonds)

From this ternary classification, three components are calculated:

#. Composition (C): Ratio of sequence positions belonging to each type of
   nucleotide.
#. Transition (T): Ratio of adjacent position pairs (i, i+1) where both
   positions belong to different nucleotide types, for each type pair.
#. Distribution (D): Normalized relative position, for each nucleotide type of:

   - First and last position.
   - 25%, 50% and 75% percentiles.

Doc2Vec
^^^^^^^

The Doc2Vec model generates a dense representation for each RNA sequence. The
process consists of two stages: tokenization and model training.

.. rubric:: Tokenization
Each sequence is split into overlapping trinucleotides with step 1, generating
:math:`L-2` tokens for a sequence of length :math:`L`. Example: AUGCGA → AUG,
UGC, GCG, CGA.

This scheme is analogous to a k-mer with k=3 but maintains positional order,
allowing the model to capture positional dependencies.

.. rubric:: Training
One independent model is trained for each RNA type, learning:

- Word vectors (trinucleotides)
- Document vectors (sequences)

Role2Vec
^^^^^^^^

Role2Vec generates embeddings that encode the structural role of each sequence
within a global similarity graph. Unlike previous descriptors operating on
isolated sequences, Role2Vec captures similarity relationships between molecules
in the complete dataset.Four-Stage Procedure:

.. rubric:: Similarity Graph Construction
Build a KDTree structure from k-mer vectors of all sequences. For each sequence,
identify its :math:`k=10` nearest neighbors in Euclidean distance. Add edges
between each sequence and its neighbors, encoding similarity.

.. rubric:: Node Structural Characterization
For each node, calculate local centrality features:

- Degree (number of edges)
- Clustering coefficient
- Histogram of neighbor degrees

These metrics define the “role” of each node in graph topology; two nodes with
the same structural profile share an equivalent role despite being distant in
the graph.

.. rubric:: Biased Random Walks
Generate random walks on the graph that favor visiting nodes with similar
structural profiles to the origin node. Walks are grouped into a corpus of
sequence strings serving as input for the embedding model.

.. rubric:: Embedding Model Training
Process walks using Word2Vec Skip-Gram model, during 16 epochs. Resulting vector
for each node is a 256-dimensional role embedding.

Recoded and Normalized Sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequences are recoded with the following table:

.. list-table::
   :header-rows: 1

   * - Symbol
     - Code
   * - A
     - 1
   * - C
     - 2
   * - G
     - 3
   * - U (or T)
     - 4
   * - Ambiguous (N, etc.)
     - 5
   * - Padding
     - 0

This encoding assigns positive integers to the four canonical nucleotides,
reserves 5 for ambiguous characters, and 0 for padding positions. The resulting
integer vector represents the complete sequence with individual nucleotide
resolution.

.. rubric:: Normalization

Before model usage:

If a sequence exceeds maximum length → truncate at 3’ end up to limit.
If shorter than maximum → pad with zeros (padding) on the right to reach target
length.

Then, recoded sequences are normalized using z-score by column on the training
set:

.. math::
   \hat{x}_{ij} = \frac{x_{ij} - \mu_j}{\sigma_j}\text{, }\sigma_j\leftarrow 1.0 \text{ if }\sigma_j\lt10^{-10}

This adjustment prevents division by zero in columns with near-zero standard
deviation.

2D Structure Characterization
-----------------------------
.. csv-table:: 2D Structure Features
   :file: data/structure2d_table.csv
   :header-rows: 1
   :align: center

.. [#] .→5, (→6, )→7
.. [#] for paired nucleotides

Structural Loops
^^^^^^^^^^^^^^^^

The RNA structure is mainly composed by loops, cataloged as follows:

.. list-table::
   :header-rows: 1

   * - Loop Kind
     - Description
   * - Hairpin loops
     - Single-stranded regions closed by base pairs; energy depends on loop
       length and closing nucleotides
   * - Interior loops
     - Unpaired regions flanked by two base pairs; includes asymmetry between
       both sides of the loop
   * - Multi-loops
     - Junctions between three or more helices; energy approximated with a
       linear model of length and branch count
   * - External loops
     - Unpaired regions on the exterior of the structure

For each loop, the following data is collected:

#. Loop kind
#. Start position
#. End position
#. Energy in hundredths of kcal/mol

This allow to calculate the global and partial MFE, and count for each loop
kind.

Recoded 2D Structure
^^^^^^^^^^^^^^^^^^^^

2D structures are recoded with the following table:

.. list-table::
   :header-rows: 1

   * - Symbol
     - Code
   * - .
     - 5
   * - (
     - 6
   * - )
     - 7
   * - Padding
     - 0

Adjacency Matrix
^^^^^^^^^^^^^^^^

Dot-bracket notation uniquely defines base pairs: each ( at position :math:`i`
is paired with exactly one ) at some position :math:`j>i`. From identified
pairs, construct adjacency matrix :math:`A\in\{0,1\}^{LxL}`:

.. math::
   A_{ij} = A_{ji} = 1 \iff \text{position } i \text{ is paired with position } j

.. rubric:: Normalization
   
To achieve fixed-size, length-independent representation:

#. Resize matrix :math:`A` to target dimensions :math:`N\times N` using bilinear
   interpolation.
#. Normalize values to range [0,1] via min-max scaling.

Interaction Characterization
----------------------------
.. csv-table:: Subsequence Features
   :file: data/subsequence_table.csv
   :header-rows: 1
   :align: center

.. csv-table:: 2D Substructure Features
   :file: data/substructure2d_table.csv
   :header-rows: 1
   :align: center

Subsequence Extraction
^^^^^^^^^^^^^^^^^^^^^^

All sequence descriptors defined for complete molecules are recalculated on
these subsequences:

.. list-table::
   :header-rows: 1

   * - Descriptor
     - Method
     - Notes
   * - k-mers
     - Frequency counts on subsequence only
     - Same dimensions (84 or 340) but represents only the interaction region
   * - CTD
     - Composition, Transition, Distribution calculated on subsequence
     - Captures physicochemical properties of contact region
   * - Doc2Vec
     - Use pre-trained models with model.infer_vector()
     - Projects trinucleotides to 256-dimensional space without retraining
   * - Role2Vec
     - Direct assignment from full molecule embeddings
     - Similarity graph is built at molecular level, not subsequence level
   * - Numeric encoding + resize
     - Same encoding scheme (A→1, C→2, G→3, U→4)
     - Resized to fixed length < complete molecule length

The binding site is represented by vectors of the same dimensions as complete
molecule descriptors, but calculated exclusively on the lncRNA–miRNA interaction
region. This allows consistent input dimensions for downstream modeling while
focusing on the biologically relevant interaction interface.

2D Substructure Extractor
^^^^^^^^^^^^^^^^^^^^^^^^^

The analysis is restricted to structural motifs whose position indices fall
within the [start, end] interval identified by miRanda. Only motifs overlapping
this window are included in computation:

.. rubric:: Binding Site MFE

The structural energy of the binding site (binding_site_energy) is defined as:

.. math::
   binding\_site\_energy = \frac{\sum_{i\in \text{binding region}}E_{loop,i}}{100}