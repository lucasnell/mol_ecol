
# I first need to source the "working R" files. They will load necessary packages and
# create functions to do the simulating.
source('./wr_files/preamble.R')
source('./wr_files/digest.R')
source('./wr_files/size_filter.R')
source('./wr_files/variants.R')
source('./wr_files/prep_seqs.R')


# Reference genome fasta file
ref_fa <- read_fasta('./genome_data/aphid_genome.fa.gz')

# ==============
# Digest
# ==============
# This is somewhat hidden (in the .digest_env environment), but it's the enzymes I 
# chose to simulate:
.digest_env$chosen_enz

# Now I'm digesting the reference genome by these restriction enzymes
digested_fa <- lapply(.digest_env$chosen_enz, function(enz) digest_genome(ref_fa, enz))

# To reduce memory consumption, I'm removing the ref_fa object and running garbage 
# collection.
# (I'll be doing this throughout for objects I no longer need.)
rm(ref_fa); invisible(gc())


# ==============
# Size filter
# ==============

set.seed(105)
filtered_fa <- lapply(digested_fa, function(.dna_ss) size_filter(.dna_ss))

rm(digested_fa); invisible(gc())


# ==============
# Add variants
# ==============

set.seed(195)
variant_fa <- lapply(filtered_fa, function(.dna_ss) make_variants(.dna_ss))

rm(filtered_fa); invisible(gc())



# ==============
# Prep sequences
# ==============

# Since make_variants returns a list of DNAStringSets for each input DNAStringSet 
# (1 for each variant), I now have to use nested lapply calls.
t0 <- Sys.time()
prepped_fa <- lapply(variant_fa, 
                     function(.dna_ss_list) {
                         lapply(.dna_ss_list, function(.dna_ss) prep_seqs(.dna_ss))
                     })
t1 <- Sys.time()
prepped_fa2 <- lapply(variant_fa, 
                      function(.dna_ss_list) {
                          lapply(.dna_ss_list, function(.dna_ss) prep_seqs2(.dna_ss))
                      })
t2 <- Sys.time()
t1 - t0; t2 - t1; rm(t0, t1, t2)

identical(prepped_fa, prepped_fa2)

# rm(variant_fa); invisible(gc())
# originally took 49 seconds





# ==============
# Write to new fasta files
# ==============

# The number of new files will be the number of samples times the number of restriction
# enzymes used.
# For my case, it was 30 files.
# I output these files to a flash drive because I wasn't sure how large they would be.


for (i in 1:3) {
    write_fastas(prepped_fa[[i]], 
                 file_names = paste0('/Volumes/64gb/', .digest_env$chosen_enz[i], '_s', 
                                     1:10, '.fa.gz'))
}; rm(i)


