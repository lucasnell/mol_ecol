
# This R file removes the many "N"s that are present on some reads, which would, under
# normal circumstances, be the second adapter sequence and would be trimmed anyway.

library(Rcpp)



Sys.setenv("PKG_LIBS" = "-lboost_iostreams")
sourceCpp('./wr_files/trim_fq.cpp')


# This took 25 minutes
fns <- paste0('/Volumes/750gb/fastq/', c('apeki_raw.fq', 'apeki.fq'))
cpp_trim_fq(fns[1], fns[2])

# If the input file was gzipped:
# cpp_trim_fqgz(fns[1], fns[2])

