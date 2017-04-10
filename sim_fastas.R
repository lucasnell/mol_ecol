
# Create fasta files of genomic variants, digested by restriction enzymes, filtered by 
# size, and prepared for sequencing simulation


# I first need to source the "working R" files. They will load necessary packages and
# create functions to do the simulating.
source('./wr_files/preamble.R')
wr_files <- list.files('./wr_files', '*.R', full.names = TRUE)
wr_files <- wr_files[!grepl('preamble.R', wr_files)]
for (f in wr_files) source(f)
rm(f, wr_files)


t0 <- Sys.time()
# Reference genome fasta file
ref_fa <- read_fasta('./genome_data/aphid_genome.fa.gz')


# ==============
# Digest
# ==============
# This is somewhat hidden (in the .wr_env environment), but it's the enzymes I 
# chose to simulate:
.wr_env$chosen_enz

# Now I'm digesting the reference genome by these restriction enzymes
digested_fa <- lapply(.wr_env$chosen_enz, function(enz) digest_genome(ref_fa, enz))
Sys.time() - t0

# To reduce memory consumption, I'm removing the ref_fa object and running garbage 
# collection.
# (I'll be doing this throughout for objects I no longer need.)
rm(ref_fa); invisible(gc())


# ==============
# Size filter
# ==============

set.seed(105)
filtered_fa <- lapply(digested_fa, function(dna_in) size_filter(dna_in))
Sys.time() - t0

rm(digested_fa); invisible(gc())


# ==============
# Add variants
# ==============

set.seed(195)
variant_fa <- lapply(filtered_fa, function(dna_in) make_variants(dna_in))
Sys.time() - t0

rm(filtered_fa); invisible(gc())



# ==============
# Prep sequences
# ==============

# (Since make_variants returns a list of dna_lists and prep_seqs uses dna_lists, I'm
# changing the name of the argument input to prep_seqs.)
prepped_fa <- lapply(variant_fa, function(dna_in_list) prep_seqs(dna_in_list))
Sys.time() - t0

rm(variant_fa); invisible(gc())


# ==============
# Write to new fasta files
# ==============

# The number of new files will be the number of restriction enzymes used.
# In my case, it was 3 files.
# I output these files to a flash drive because I wasn't sure how large they would be.
# All 3 (gzipped) files totalled 49.6 MB.

fasta_list <- lapply(1:length(prepped_fa), 
                     function(i) {
                         .dna(c(prepped_fa[[i]], recursive = TRUE))
                     })
fasta_names <- paste0('/Volumes/64gb/fasta/', .wr_env$chosen_enz, '.fa.gz')
write_fastas(fasta_list, file_names = fasta_names)



