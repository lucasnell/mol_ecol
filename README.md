`mol_ecol`
==========

Code for Molecular Ecology paper, Spring 2017
-----------


Author: Lucas Nell

## File and folder notes

__The folder `dv_files`__ contains files explaining how all steps were developed.
For a given step, the `.md` or `.pdf` file is the one you should examine to understand
the code's overall rationale, and the `.R` file includes all the raw R code, which is 
easier for testing.

__The folder `wr_files`__ contains "working R" files intended to be `source`d and 
the objects within used directly.
Because they are intended to be `source`d from the parent directory, you need to adjust
all file paths if these files are located in the working directory from which you'll
be `source`-ing them.

__The folder `bg_data`__ contains background data on fragment sizes in GBS and on 
next-gen sequencing simulators. References for these files are contained in the
files [`./dv_files/size_filter.md`](./dv_files/size_filter.md) and 
[`./dv_files/simulators.md`](./dv_files/simulators.md), respectively.

__The file `genome_data`__ (present in file paths in some scripts, but it's not in this 
repo) refers to a [symbolic link](https://kb.iu.edu/d/abbe) I made on my computer.
It links to the Box folder on my computer that contains all fasta files and other large,
genome-related data files.
I did it this way to reduce the size of this repo.

