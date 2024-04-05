#!/usr/bin/env Rscript
##-----------------------------------------------------------------------------##
## Script for extracting consensus coverage for multiple samples               ##
## Requires bedtools, for all, and samtools, for sorted bams                   ##
##-----------------------------------------------------------------------------##

#Functions
#Calculates coverage of consensus
consensusCov <- function(bed.path, consGr, max.val) {

    if (tolower(file_ext(bed.path)) != "bam") {
        tgt = read.delim(bed.path, header = F, stringsAsFactor = F)
        tgtGr = GRanges(data.frame(chr = tgt[,1], start = tgt[,2] + 1,
                                   end = tgt[,3], stringsAsFactors = F))
    } else {
        bamF = BamFile(bed.path)
        ga = readGAlignments(bamF)
        tgtGr = as(ga, "GRanges")            
    }
    cvg = countOverlaps(consGr, tgtGr)

    return (pmin(cvg, max.val))
}

#load libraries
suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(parallel))

"Coverage of consensus regions for multiple samples. 

Usage: 
    consensus_coverage.R [-h -t threads -m max_val -i ids_file -n filename_trim] <consensus_path> <bed_path>... 

Options:
    -h --help                       Show this help page.            
    -i --identifiers ids_file       File with sample identifiers in the same order than
                                    the bed paths passed in the command line. Either 
                                    use this option or the -n/fn_trim.
    -m --max max_val                Maximum value for coverage metric [default: 1e100]
    -n --fn_trim filename_trim      Number of chars used from filename for identifying 
                                    samples [default: 1000].  
    -t --threads threads            Threads to use in the computation [default: 1].
    <consensus_path>                Path to BED file containing consensus regions.
    <bed_path>                      Path(s) to BAM/BED file(s) containing mapped features.
" -> doc

opts <- docopt(doc)

max.val <- as.double(opts$max)
threads.n <- as.integer(opts$threads)
trim.n <- as.integer(opts$fn_trim)
cons.path <- opts$consensus_path
bed.paths <- unlist(opts$bed_path)

#Sort samples according to IDs
if (!is.null(opts$identifiers) > 0) {
  ids = read.delim(opts$identifiers, header = F, stringsAsFactors = F)[,1]
} else {
  ids <- substr(basename(bed.paths),1, trim.n)
}
ix <- order(ids)
ids <- ids[ix]

cons = read.delim(cons.path, header = F, stringsAsFactors = F)

consGr = GRanges(data.frame(chr = cons[,1], start = cons[,2] + 1,
                            end = cons[,3], stringsAsFactors = F))
cat(sprintf("chr\tstart\tend\t%s\n", paste(ids, collapse = "\t")))

covs <- mclapply(ix, function(i) consensusCov(bed.paths[i], consGr, max.val), 
                  mc.cores = threads.n)

#print results
if (length(covs) > 0) {
  
  ncols <- length(bed.paths) 
  nrows <- length(covs[[1]])
  res <- matrix(0, ncol = ncols, nrow = nrows)
  
  for (i in 1:ncols)
    res[,i] <- covs[[i]]
  
  df = cbind(cons[,1:3], res)

  write.table(df, file = stdout(), quote = F, sep = '\t', row.names = F,
              col.names = F)
}
