# SynExtend 1.19.3
* `ExoLabel` will no longer crash when given a network lacking a trailing newline
* Various internal improvements and code refinements

# SynExtend 1.19.2
* Lots of bug fixes for `ExoLabel`
* `ExoLabel` now reports disk consumption during execution

# SynExtend 1.19.1
* `ExoLabel` is even faster due to in-place external sort for faster file I/O
* Other quality of life improvements to `ExoLabel`

# SynExtend 1.19.0
* First development version of Bioconductor 3.21

# SynExtend 1.18.0
* Official Bioconductor 3.20 release

# SynExtend 1.17.7
* `ExoLabel` is much much faster and does a better job cleaning up when aborted early
* `ExoLabel` now has fewer arguments
* Updates to man pages

# SynExtend 1.17.6
* Lots of internal improvements to `ExoLabel` to increase computational speed and decrease disk
usage.
* `ExoLabel` will no longer crash if given relative paths.
* Adds more internal error checking to prevent some rare bugs.
* Updates man pages to reflect new changes.
* Updates `EstimateExoLabel` to reflect new changes.

# SynExtend 1.17.5
* `ExoLabel` will no longer brick R during sorts on large files.
* `ExoLabel` reports more progress during some lengthy processing sections when `verbose=TRUE`
* Known issue: "Copying source file" step is still non-interruptable, will be fixed in a later update

# SynExtend 1.17.4
* `ExoLabel` now allows an `inflation` argument to control application of inflation

# SynExtend 1.17.3
* `predict.EvoWeaver` now supports returning p-values separately from raw score for some algorithms.
* OpenMP implementation for `EvoWeaver` algorithms that support it has been fixed

# SynExtend 1.17.2
* `RandForest` function added to train random forest models
* Associated man pages for `RandForest` and `DecisionTree` objects
* New methods for `DecisionTree` objects to plot and coerce to `dendrogram`
* Small bugfix to `subset.dendrogram`

# SynExtend 1.17.1
* Major updates to `EvoWeaver`:
  * `predict.EvoWeaver` now returns a `data.frame` by default
  * `Method` arguments are updated to match their names in the associated EvoWeaver manuscript
  * Above changes have propagated to documentation files
  * New Phylogenetic Profiling methods with improved accuracy
  * New meta-methods `PhylogeneticProfiling`, `PhylogeneticStructure`, `GeneOrganization`, `SequenceLevel` for `predict.EvoWeaver`
  * New pre-trained Ensemble models have been included
* Updates to `ExoLabel` for better status printing

# SynExtend 1.16.0
* Bioc 3.19 release

# SynExtend 1.15.3
* Addition of `FastLabelOOM` function to find communities in graphs/networks on disk space.

# SynExtend 1.15.2
* Addition of `PrepareSeqs` function, beginning the process of deprecating `PairSummaries` in favor of more cohesive and user friendly functions.

# SynExtend 1.15.1
* Fixes bug in JRF distance causing scores to be higher than expected.

# SynExtend 1.13.8
* Updates to all EvoWeaver documentation files
* Fixed small bug in `PhyloDistance` causing `Method='JRF'` to return similarity rather than the distance
* Fixed small bug in `TreeDistance.EvoWeaver` resulting in an inconsistent calculation of score when using `TreeMethods='JRF'`

# SynExtend 1.13.7
* Small fixes

# SynExtend 1.13.6
* `ProtWeaver` and `ProtWeb` have been renamed to `EvoWeaver` and `EvoWeb`, respectively
* New sequence level method for `EvoWeaver`
* Various small internal updates to `EvoWeaver`

# SynExtend 1.13.5
* Minor changes to `SelectByK` and vignette

# SynExtend 1.13.4
* New predictor `PAPV.ProtWeaver` to calculate p-values for presence/absence profiles.
* `ContextTree` now uses `MirrorTree` with species tree correction and p/a overlap correction
* Updates to documentation

# SynExtend 1.13.3
* `predict.ProtWeaver` now supports multiple algorithms at once (ex. `predict(ew, Method=c("Jaccard", "Hamming"))`)
* Documentation for `ProtWeaver` and associated methods has been updated to match recent updates.

# SynExtend 1.13.2

* `FastQFromSRR` function added as a convenience wrapper for the SRAtoolkit function `fastq-dump`.

# SynExtend 1.13.1

* `SuperTree` now works directly with `dist` objects, providing better performance and scaling
* Updates to `simMat` objects
  * No longer throw a warning when initialized in RStudio
  * Formatting is cleaner and supports larger object names
* Updates to `NVDC.ProtWeaver`
  * Now supports amino acid sequences using the `DNAseqs=FALSE` argument
  * Now calculates a p-value-weighted score
* Adds `MakeBlastDb` function to create a BLAST database from R, plus associated documentation updates
* Smaller fixes to some `ProtWeaver` methods
  * `predict.ProtWeaver` no longer returns using `invisible` (this was annoying and unneccessary)
  * APC correction for `MutualInformation.ProtWeaver` removed to allow for parallelization
  * `MirrorTree.ProtWeaver` now works correctly with `MTCorrection="speciestree"`
  * `CorrGL.ProtWeaver` now uses Fisher's Exact Test for p-values rather than the R value of spearman correlation
* Many internal performance improvements
  * `ProtWeaver` almost entirely uses `dist` objects rather than `matrix`, saving significantly on memory
  * faster `Cophenetic` function implemented internally
  * Copied internal `.Call('cophenetic')` from `DECIPHER` to `SynExtend` to avoid potential namespace issues
  * Small fixes to remove some notes from `BiocCheck::BiocCheck()`
* Variety of small updates to pass `BiocCheck`

# SynExtend 1.12.0

* Official Bioconductor 3.17 release (even with SynExtend 1.11.8)

# SynExtend 1.11.8

* Fixes various small bugs in `MoransI`
* Adds some multiprocessing support (more will be added in the future)
* Slight rework to species trees and their interaction with `ProtWeaver` objects
  - `ProtWeaver` has new attribute `speciesTree`, can be initialized with a `dendrogram object`
  - New method `SpeciesTree` to get species tree from a `ProtWeaver` object (or compute one, if it doesn't exist)
* Various internal improvements for Bioc style consistency
* Various documentation updates

# SynExtend 1.11.7

* Adds new optimized `dendrapply` implementation (overloads `stats::dendrapply`)
* Adds`HungarianAlgorithm` for optimal solving of the linear assignment problem (O(n^3) complexity)
* Adds new C code for fast computation of Pearson's R and p-value
* Adds new `Ancestral.ProtWeaver` algorithm for calculating coevolution from correlated residue changes
* Supporting code and documentation for `Ancestral.ProtWeaver`
* Other new internal methods
* Various updates and optimizations to internal methods and documentations
* Updates `GRF` method to be called `CI` (for **C**lustering **I**nformation Distance)
* `Method="CI"` in `PhyloDistance` now calculates an approximate p-value using simulated data from Smith (2020)

# SynExtend 1.11.6

* Adds new Residue method `NVDT` using gene sequence **N**atural **V**ector with **D**inucleotide and **T**rinucleotide frequency
* Adds some new C methods to speed up calculation of NVDT
* Fixes `.Call()` not using `PACKAGE="SynExtend"`
* Updates to documentation

# SynExtend 1.11.5

* Adds new colocalization algorithm `ColocMoran`, uses `Coloc` with `MoransI` to correct for phylogenetic signal
* Adds new colocalization algorithm `TranscripMI`, uses mutual information of transcriptional direction
* Adds new corrections/checks to allow for transcriptional direction to be in labels
* Various bugfixes to support new four number labelling scheme
* Various documentation updates
* Adds new function `MoransI` to calculate Moran's *I* for a set of spatially distributed signals

# SynExtend 1.11.4

* Internal code refactor
* `ShuffleC` now supports reproducibility using R's `set.seed`
* `ShuffleC` now support sampling with replacement, performance is around 2.25x faster than `sample`

# SynExtend 1.11.3

* Internal bugfixes for JRF Distance--previous commit was incorrectly calculating values
* Adds new `TreeDistance` predictor for `ProtWeaver`, incorporating all tree distance metrics; these metrics are bundled due to some backend optimizations that improve performance
* Bugfixes for `PhyloDistance`
* Adds Random Projection for `MirrorTree` predictor to solve memory problems and increase accuracy
* New internal random number generator using xorshift, significantly faster than `sample()`
* `HammingGL` changed to `CorrGL`, now uses Pearson's R weighted by p-value
* Refactors internal predictors to reduce size of codebase and remove redundancies
* Internal `ShuffleC` function to replicate `sample` functionality with 2-6x speedup
* Method `GainLoss` now uses bootstrapping to estimate a p-value
* Updates to documentation files

# SynExtend 1.11.2

* Adds KF Distance for trees
* Adds Jaccard Robinson Foulds Distance for trees
* Reworks tree distances into `PhyloDistance` function
* Numerous new documentation pages
* Updates internal functions to use `rapply` instead of `dendrapply` to avoid stack overflow issues due to R recursion

# SynExtend 1.11.1

* Minor bugfix to RF distance
* updates gitignore for workflows

# SynExtend 1.10.1

* Memory leak bugfix

# SynExtend 1.9.21

* Adds new `RFDist` function to calculate Robinson-Foulds Distance
* Adds normalization for `GeneralizedRF` to make the distance between 0 and 1
* Minor bugfixes
* Documentation for new functions

# SynExtend 1.9.20

* Adds new `GeneralizedRF` function to calculated information-theoretic Generalized Robinson-Foulds distance between two dendrograms.
* Documentation for new function
* New ProtWeaver predictor based off of `GeneralizedRF` metric
* New internal C source code for `GeneralizedRF`

# SynExtend 1.9.19

* Adds new `DPhyloStatistic` function to calculate the D-statistic for a binary state against a phylogeny following Fritz and Purvis (2009).
* Documentation for new function
* new internal C source code for `DPhyloStatistic`
* new internal C source code for random utility functions, currently only has functions to generate random numbers

# SynExtend 1.9.18

* Various internal improvements to presence/absence profile methods

# SynExtend 1.9.17

* Adds new prediction algorithm `GainLoss`
* Adds new internal C implementation of dendrograms, significantly faster than R dendrograms
* `ProtWeaver` methods `Behdenna` and `GainLoss` can now infer a species tree when possible
* Updates `Jaccard` and `Hamming` methods to use C implementations for distance calculation
* Adds `HammingGL` method to calculate Hamming distance of gain/loss events
* Minor bugfixes to `ProtWeaver` methods relating to subsetting
* Updates to various `man` pages

# SynExtend 1.9.16

* Removes `flatdendrapply`, function was already included in SynExtend
* minor bugfixes to `ProtWeaver`

# SynExtend 1.9.15

* Edits to `SelectByK`, function can work as intended, but is still too conservative at false positive removal.

# SynExtend 1.9.14

* Adds new function `flatdendrapply` for more options to apply functions to dendrograms. Function is used in `SuperTree`.
* Adds new function `SuperTree` to construct a species tree from a set of gene trees.
* Adds new dataset `SuperTreeEx` for `SuperTree` and `flatdendrapply` examples.

# SynExtend 1.9.13

* `SelectByK` function argument `ClusterSelect` switched to `ClusterScalar`. Cluster number selection now performed by fitting sum of total within cluster sum of squares to a right hyperbola and taking the ceiling of the half-max. Scalar allows a user to pick different tolerances to select more, or less clusters. Plotting behavior updated.

# SynExtend 1.9.12

* `simMat` class now supports empty indexing (`s[]`)
* `simMat` class now supports logical accession (`s[c(T,F,T),]`)

# SynExtend 1.9.11

* Added the function `SelectByK` that allows for quick removal of false positive predicted pairs based on a relatively simple k-means approach. Function is currently designed for use on the single genome-to-genome pairwise comparison, and not on an all-vs-all many genomes scale, though it may provide acceptable results on that scale.

# SynExtend 1.9.10

* New `simMat` class for `dist`-like similarity matrices that can be manipulated like base matrices
* Major update to `ProtWeaver` internals
* All internal calls use `simMat` objects whenever possible to decrease memory footprint
  * Note `ContextTree` and `ProfDCA` require matrices internally
* `ProtWeb` objects now inherit from `simMat`
  * `ProtWeb.show` and `ProtWeb.print` now display predictions in a more natural way
  * `GetProtWebData()` deprecated; `ProtWeb` now inherits `as.matrix.simMat` and `as.data.frame.simMat`
* New documentation pages for `simMat` class
* `GetProtWebData` documentation page reworked into `ProtWeb` documentation file.
* Fixes new bug in `Method='Hamming'` introduced in SynExtend 1.9.9

# SynExtend 1.9.9

* Fixes minor bug in `Method='Hamming'`
* Moves some code around

# SynExtend 1.9.8

* Major refactor to file structure of `ProtWeaver` to make individual files more manageable
* Adds new documentation files for individual prediction streams of `predict.ProtWeaver`

# SynExtend 1.9.7

* `BlockReconciliation` now returns a an object of class `PairSummaries`.

# SynExtend 1.9.6

* Fixes an error where warnings were mistakenly output to the user

# SynExtend 1.9.5

* Moves platform-specific files in `src/` (originally added by mistake)

# SynExtend 1.9.4

* Lots of bugfixes to `ResidueMI.ProtWeaver`
* `predict.ProtWeaver` now correctly labels rows/columns with gene names, not numbers
* `predict.ProtWeaver` now correctly handles `Subset` arguments
* `predict.ProtWeaver(..., Subset=3)` will correctly predict for all pairs involving gene `3` (or for any gene `x`, as long as `Subset` is a length 1 character or integer vector).

# SynExtend 1.9.2

* Adds residue MI method to `ProtWeaver`
* Various bugfixes for `ProtWeaver`

# SynExtend 1.7.14

* Various improvements for `GenRearrScen`, improves consistency and output formatting
* Major bugfix for `ProtWeaver` methods using dendrogram objects
* `ProtWeaver` now correctly guards against non-bifurcating dendrograms in methods that expect it

# SynExtend 1.7.13

* Introduces new `ProtWeaver` class to predict functional association of genes from COGs or gene trees. This implements many algorithms commonly used in the literature, such as MirrorTree and Inverse Potts Models.
* `predict(ProtWeaverObject)` returns a `ProtWeb` class with information on predicted associations.
* Adds `BlastSeqs` to run BLAST queries on sequences stored as an `XStringSet` or `FASTA` file.

# SynExtend 1.7.12

* Updates to `ExtractBy` function. Methods and inputs simplified and adjusted, and significant improvements to speed.

# SynExtend 1.7.11

* Updated `NucleotideOverlaps` to now correctly registers hits in genes with a large degree of overlap with the immediately preceding gene.
* Fixed aberrant behavior in `BlockExpansion` where contigs with zero features could cause an error in expansion attempts.

# SynExtend 1.7.10

* `BlockReconciliation` now allows for setting either block size or mean PID for reconciliation precedence.

# SynExtend 1.7.9

* Added retention thresholds to `BlockReconciliation`.

# SynExtend 1.7.8

* `BlockExpansion` cases corrected for zero added rows.

# SynExtend 1.7.7

* Improvements to `BlockExpansion` and `BlockReconciliation` functions.

# SynExtend 1.7.5

* Began integration of `DECIPHER`'s `ScoreAlignment` function.

# SynExtend 1.7.4

* Fixed a bug in `PairSummaries` function.

# SynExtend 1.7.3

* Added `BlockExpansion` function.

# SynExtend 1.7.2

* Adjustment in how `PairSummaries` handles default translation tables and GFF derived gene calls.

# SynExtend 1.7.1

* Large changes under the hood to `PairSummaries`.
* Failure to accurately assign neighbors in some cases should now be fixed.
* Extraction of genomic features is now faster.
* `OffSetsAllowed` argument now defaults to `FALSE`. This argument may be dropped in the future in favor of a more complex function post-summary.
* Small edits to `SequenceSimilarity`

# SynExtend 1.5.4

* Added the function `SubSetPairs` that allows for easy trimming of predicted pairs based on conflicting predictions and / or prediction statistics.
* Added the function `EstimageGenomeRearrangements` that generates rearrangement scenarios of large scale genomic events using the double cut and join model.

# SynExtend 1.5.3

* Added the function `SequenceSimilarity` and made improvements to runtime in `DisjointSet`.

# SynExtend 1.4.1

* Fixed a small bug in consensus scores in `PairSummaries` where features facing on different strands had their score computed incorrectly.

# SynExtend 1.3.15

* Changes to concensus score in `PairSummaries`.

# SynExtend 1.3.14

* Major changes to the `PairSummaries` function and minor changes to `NucleotideOverlaps`, `ExtractBy`, and `FindSets`. Adjustments to the model that `PairSummaries` calls on to predict PIDs.

# SynExtend 1.3.13

* `ExtractBy` function has been added. Allows extraction of feature sequences into `XStingSet`s organized by the a `PairSummaries` object or the single linkage clusters implied by pairings within the `PairSummaries` objects.
* `DisjointSet` function added to extract single linkage clusters from a `PairSummaries` object.

# SynExtend 1.3.12

* `PairSummaries` now computes 4-mer distance between predicted pairs.

# SynExtend 1.3.11

* `PairSummaries` now returns a column titled Adjacent that provides the number of directly adjacent neighbor pairs to a predicted pair. Gap filling code adjusted.
* The function `FindSets` has been added and performs single linkage clustering on a pairs list as represented by vectors of integers using the Union-Find algorithm. Long term this function will have a larger wrapper function for user ease of access but will remain exposed.

# SynExtend 1.3.10

* `NucleotideOverlap` now passes it's GeneCalls object forward, allowing `PairSummaries` to forego inclusion of that object as an argument.

# SynExtend 1.3.9

* Minor vignette and suggested package changes.

# SynExtend 1.3.8

* `PairSummaries` now allows users to fill in specific matching gaps in blocks of predicted pairs with the arguments `AllowGaps` and `OffSetsAllowed`.

# SynExtend 1.3.7

* Adjustments to progress bars in both `PairSummaries` and `NucleotideOverlap`.
* PID prediction models in `PairSummaries` adjusted.

# SynExtend 1.3.6

* Contig name matching has been implemented. Scheme expects users to follow NCBI contig naming and gff formats, accepting contig names from gffs directly, and removing the first whitespace and everything thereafter from FASTA headers. Contig name matching can be disabled if users wish, using the argument `AcceptContigNames`, but ensuring that the correct contigs in GeneCalls objects are matched to the appropriate contigs in Synteny objects are then the user's responsibility.

# SynExtend 1.3.5

* `PairSummaries` now translates sequences based on `transl_table` attributes provided by gene calls
* `PairSummaries` now uses a generic model for predicting PID
* `gffToDataFrame` now parses out the `transl_table` attribute

# SynExtend 1.3.2

* Minor changes to `NucleotideOverlap`
* Major changes to `PairSummaries` - can now take in objects of class `Genes` build by the DECIPHER function `FindGenes()`

# SynExtend 0.99.30

* Vignette and help files edited for clarity

# SynExtend 0.99.1

* `SynExtend` submitted to Bioconductor
* Added function `gffToDataFrame`
* Added function `NucleotideOverlap`
* Added function `PairSummaries`
