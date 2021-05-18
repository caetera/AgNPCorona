## File description

**coronaAnalysis.py** - Main script that performs quantification and generates figures and tables for the paper

**coronaBatcher.py** - Batch script that executes external tools (msconvert, MSGF+, Percolator) and prepares the input for KNIME workflow

**coronaHelpers.py** - Helper data handling functions used by other scripts

**coronaPlots.py** - Helper graphing functions used by main script

**msgf.mods** - Modification file used by MSGF+

**multiplier.csv** - Direction coefficients for alpha-, beta-, and turn propensities

**msgf-percolator-lfq.knwf** - KNIME workflow that executes OpenMS tools for label-free quantification

**parsePPD.py** - Script to convert Plasma Proteome Database to FASTA format and extract concentration data

**gravy.py** - functions to calculate GRAVY and amino acid content

**shrinkTest.py** - multithreaded permutation test for shrinkage measure

## Publication

V. Gorshkov, J. A. Bubis, E. M. Solovyeva, M. V. Gorshkov, and F. Kjeldsen *Protein corona formed on silver nanoparticles in blood plasma is highly selective and resistant to physicochemical changes of the solution* ***Environ. Sci.: Nano***, **2019**, 6, 1089-1098; DOI: [10.1039/C8EN01054D](https://pubs.rsc.org/fa/content/articlehtml/2019/en/c8en01054d)
