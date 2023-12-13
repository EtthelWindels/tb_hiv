## -----------------------------------------------------------------------------
## Replace names of sequences with complete names for 
## BEAST analyses (GNUMBER/date/type) with 'type' defined by HIV status
## 2023-03-24 Etthel Windels
## -----------------------------------------------------------------------------


## Load libraries

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(seqinr)))

## ---------------------------

## Parser

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata  file")
parser$add_argument("--alignment", type="character", help="Alignment file")
parser$add_argument("--output_alignment", type = "character",
                    help = "Output file for updated alignment")
args <- parser$parse_args()

## ---------------------------

## Read arguments

ALIGNMENT <- args$alignment
METADATA <- args$metadata
OUTPUT_ALIGNMENT <- args$output_alignment

print(paste("metadata: ", METADATA))
print(paste("alignment: ", ALIGNMENT))
print(paste("output alignment: ", OUTPUT_ALIGNMENT))

## ---------------------------

## Read files

alignment <- seqinr::read.fasta(file = ALIGNMENT, seqtype="DNA", forceDNAtolower = F)
metadata <- read.table(file = METADATA, header=T,sep= '\t')

## ---------------------------

## Filter metadata

metadata <- 
  metadata %>%
  filter(G_NUMBER %in% names(alignment))
 
## ---------------------------

## Adjust alignment

# 1. Filter alignment to remove outgroup G00157
alignment <- 
  alignment[names(alignment) %in% metadata$G_NUMBER]

# 2. Adjust alignment names
gnumber <- names(alignment)
full_names <- metadata$new_id[match(gnumber, metadata$G_NUMBER)]
print(full_names)
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names

## ---------------------------

## Save output alignment

seqinr::write.fasta(alignment, names = names(alignment), file.out = OUTPUT_ALIGNMENT)


