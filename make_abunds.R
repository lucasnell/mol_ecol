# Make abundance files for different digestions

source('./wr_files/preamble.R')

enz <- c('ApeKI', 'BstBI', 'BspEI')
n_samps <- 10

fastas <- lapply(enz, function(e) read_fasta(paste0('/Volumes/64gb/fasta/', e, '.fa.gz')))

f <- fastas[[1]]

.names <- names(f)
length(f)

rep(1:3, each = 2)
