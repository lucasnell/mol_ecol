# This preamble loads packages and runs code necessary for >1 "working R" files

suppressPackageStartupMessages({
    library(ShortRead)
    library(magrittr)
    library(dplyr)
    library(purrr)
    library(ggplot2)
    library(Rcpp)
    library(gtools)
    library(parallel)
    library(RcppArmadillo)
})

# This is to hide certain functions and objects from view
.wr_env <- new.env()

# This sets the default ggplot theme
theme_set(theme_classic() %+replace% theme(strip.background = element_blank()))


# Constructing dna and dna_list objects
.dna <- function(x) {
    out_dna <- as.character(x)
    class(out_dna) <- 'dna'
    if (is.null(names(x))) {
        names(out_dna) <- paste0('seq_', 1:length(x))
    } else {
        names(out_dna) <- names(x)
    }
    return(out_dna)
}
.dna_list <- function(in_list) {
    out_list <- in_list
    class(out_list) <- 'dna_list'
    return(out_list)
}



# Printing only limited information from a dna type object
print.dna <- function(dna, all = FALSE) {
    wid <- options('width')$width
    len <- length(dna)
    max_print <- ifelse(all, len, 10)
    if (len <= max_print) {
        print_inds <- 1:len
    } else {
        print_inds <- c(1:5, NA, (len - 4):len)
    }
    print_mat <- sapply(print_inds, 
                        function(i) {
                            if (is.na(i)) {
                                return(rep('...', 3))
                            }
                            .seq <- dna[i]
                            .width <- nchar(.seq)
                            .name <- names(.seq)
                            return(paste(c(.width, .seq, .name)))
                        })
    tryCatch(wid_chars <- max(c(nchar(print_mat[1,]), 5)),
             error = function(e) stop(print_mat))
    tryCatch(name_chars <- max(c(nchar(print_mat[3,]), 5)),
             error = function(e) stop(print_mat))
    seq_chars <- min(c(wid - wid_chars - name_chars - 2, max(nchar(print_mat[2,]))))
    half_seq_len1 <- floor((seq_chars - 3) / 2)
    half_seq_len2 <- (seq_chars - 3) - half_seq_len1
    print_mat[2,] <- ifelse(nchar(print_mat[2,]) > seq_chars, 
                            paste0(substr(print_mat[2,], 1, half_seq_len1), '...',
                                   substr(print_mat[2,], 
                                          nchar(print_mat[2,]) - half_seq_len2 + 1, 
                                          nchar(print_mat[2,]))), 
                            print_mat[2,])
    cat('Number of sequences =', prettyNum(len, big.mark = ','), '\n\n')
    cat(paste0(c('width', rep(' ', wid_chars - 5)), collapse = ''), 
        paste0(c('seq', rep(' ', seq_chars - 3)), collapse = ''), 
        'names\n')
    for (i in 1:ncol(print_mat)) {
        cat(paste0(c(print_mat[1,i], rep(' ', wid_chars - nchar(print_mat[1,i]))), 
                   collapse = ''),
            paste0(c(print_mat[2,i], 
                     rep(' ', seq_chars - nchar(print_mat[2,i]))), 
                   collapse = ''),
            paste0(c(print_mat[3,i], '\n'), collapse = ''))
    }
}


# Printing only limited information from a dna_list type object
print.dna_list <- function(dna_list) {
    wid <- options('width')$width
    len <- length(dna_list)
    cat('dna_list with', len, 'samples\n\nSummary of sequence lengths:\n')
    widths <- lapply(1:len, 
                     function(i) {
                         x <- nchar(dna_list[[i]])
                         matrix(round(c(i, min(x), max(x), mean(x), sd(x), length(x)), 2),
                                nrow = 1)
                     })
    wid_df <- as.data.frame(do.call(rbind, widths))
    names(wid_df) <- c('sample', 'min', 'max', 'mean', 'sd', 'n')
    print(wid_df, row.names = FALSE, right = FALSE)
}



# Generic methods to extend the 'c' and '['
c.dna <- function(..., recursive = FALSE)  {
    dots <- list(...)
    classes <- rep("dna", length(dots))
    res <- structure(unlist(dots, recursive = FALSE), class = classes)
    return (res)
}


`[.dna` <- function (x, i)  {
    y <- unclass(x)[i]
    y <- .dna(y)
    return (y)
}

`[<-.dna` <- function(x, i, value) {
    y <- unclass(x)
    y[i] <- value
    y <- .dna(y)
    return(y)
}


`[.dna_list` <- function (x, i)  {
    y <- unclass(x)[i]
    y <- .dna_list(y)
    return (y)
}
`[[.dna_list` <- function (x, i)  {
    y <- unclass(x)[[i]]
    y <- .dna(y)
    return (y)
}

`[<-.dna_list` <- function(x, i, value) {
    y <- unclass(x)
    y[[i]] <- value
    y <- .dna_list(y)
    return(y)
}
`[[<-.dna_list` <- function(x, i, value) {
    y <- unclass(x)
    y[[i]] <- value
    y <- .dna_list(y)
    return(y)
}







# This reads a fasta file into a dna object
read_fasta <- function(file_name) {
    fasta <- readFasta(file_name)
    dna_ss <- sread(fasta)
    dna_in <- .dna(dna_ss)
    return(dna_in)
}




# This function can write a dna_list, a list of dna objects, or a single dna object
# to fasta file(s).
# It will write one file per dna object.
write_fastas <- function(dna, file_names = NULL) {
    if (!any(c('dna_list', 'dna', 'list') %in% class(dna))) {
        stop('dna argument must be a dna, dna_list, or list class')
    }
    dna_ss <- function(dna) {
        dna_ss_out <- DNAStringSet(dna)
        names(dna_ss_out) <- names(dna)
        return(dna_ss_out)
    }
    if ('dna' %in% class(dna)) {
        write_list <- list(dna_ss(dna))
    } else {
        if (all(sapply(dna, function(x) 'dna' %in% class(x)))) {
            write_list <- lapply(dna, dna_ss)
        } else {
            stop('if inputing a list, it must only contain dna objects')
        }
    }
    
    if (is.null(file_names)) {
        file_names <- paste0('frags_', 1:length(write_list), '.fa.gz')
    } else if (length(file_names) != length(write_list)) {
        stop('dna and file_names must be same length')
    }

    for (i in 1:length(write_list)) {
        writeFasta(write_list[[i]], file = file_names[i], mode = 'w', 
                   compress = grepl('.gz', file_names[i]))
        cat(file_names[i], 'written... \n')
    }
    return(invisible(NULL))
}

