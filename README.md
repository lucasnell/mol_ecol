`mol_ecol`
==========

Code for Molecular Ecology paper, Spring 2017
-----------


Author: Lucas Nell

## File and folder notes

__All files starting with `dv_`__ are files explaining how all steps were developed.
For a given step, the `.md` or `.pdf` file is the one you should examine to understand
the code's overall rationale, and the `.R` file includes all the raw R code, which is 
easier for testing.

__All files starting with `wr_`__ are "working R" files intended to be `source`d and 
the objects within used directly.

__The folder `bg_data`__ contains background data on fragment sizes in GBS and on 
next-gen sequencing simulators. References for these files are contained in the
files [`dv_size_filter.md`](./dv_size_filter.md) and 
[`dv_simulators.md`](./dv_simulators.md), respectively.

__The file `genome_data`__ (present in some scripts, but not in this repo) 
refers to a [symbolic link](https://kb.iu.edu/d/abbe) I made on my computer.
It links to the Box folder on my computer that contains all fasta files and other large,
genome-related data files.
I did it this way to reduce the size of this repo.

