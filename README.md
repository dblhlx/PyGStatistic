
__Update__: Instead of the python module (scipy.stats.chi2_contingency), [the general formula](https://en.wikipedia.org/wiki/G-test) is used for G-statistic calculation in the new version of the script. The function based on this formula takes four one-dimensional numpy arrays or four Pandas series as input, which increases the speed of the script by a few hundred times.

# PyGStatistic
BSA-Seq data analysis - Python implementation of the G-statistic method

This is the Python implementation of [the G-statistic method](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255) for BSA-Seq data analysis. The script can be tested by issuing the command `python PyGStatistic.py` in a terminal. The [snp_final.tsv.bz2 file](https://github.com/dblhlx/PyBSASeq/blob/master/snp_final.tsv.bz2) needs to be uncompressed to the working directory before testing. Options as in [PyBSASeq](https://github.com/dblhlx/PyBSASeq) are needed to use the script for your own SNP datasets.
