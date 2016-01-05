treebest
========

TreeBeST: Tree Building guided by Species Tree (Ensembl Compara modifications)

This repository holds the necessary changes of [Heng Li's version](https://github.com/lh3/treebest) to run the latest Ensembl Compara pipeline.

You can find more documentation on SourceForge: (http://treesoft.sourceforge.net/treebest.shtml)

The main new features are:
* new `-s` option in `treebest sdi`, to allow a user-defined species tree. This change is from [Albert Vilella](https://sites.google.com/site/avilella/)
* new `T` node-tag in the NHX output: a bit-field listing the input trees that support the node. This is populated by the _mmerge_ algorithm
* new `-I` option in `treebest nj`, to carry on the `T` tags from the input tree
* new `-Z` option in `treebest best`, to redefine the PhyML variable `MIN_DIFF_LK`. It prevents PhyML from crashing during its computation
* new `-X` option in `treebest best`, to give a higher weight to the likelihood that comes from the reconciliation with the species tree (default 1)
* Species-intersection scores are now also reported as floating-point values under the `DCS` node-tag. The value is between 0 and 1, and displayed with 4 decimals.
* new `-I` option in `treebest best`, to start from the input tree instead of building one

Other changes include:
* bugfixes / tweaks when processing the filtered alignments (TreeBeST includes a Clustal-score-based MSA-filtering step)
* bugfixes / tweaks when merging the trees
* using `double` instead of `float` for floating-point values

## Branches and tags

There is a single branch (master) where all the development goes. The version number stated in the source code (1.9.2) is not maintained.

Tags are used to point at specific versions:
* lh3: when Heng Li handed-over the code to Ensembl Compara
* albert: changes made by Albert Vilella
* ensembl\_production\_XX: the version used for the production of Ensembl version XX
