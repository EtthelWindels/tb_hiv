## -----------------------------------------------------------------------------
## Replace names of sequences with complete names for 
## BEAST analyses (GNUMBER/date/type) with 'type' defined by CD4 cell count
## 2023-03-24 Etthel Windels
## -----------------------------------------------------------------------------


## Load libraries

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(seqinr)))
suppressMessages(suppressWarnings(require(stringr)))

## ---------------------------

## Parser

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata file")
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

## Adjust metadata

metadata$type <-  ifelse(metadata$CD_counts<350, 1, 0) # 0 = high CD4 cell count (HIV neg); 1 = low CD4 cell count (HIV pos)
metadata <- metadata %>%
  mutate(new_id = paste(G_NUMBER, Isolation_date, type, sep = "/"))

## ---------------------------

## Adjust alignment

names(alignment) <- str_split_fixed(names(alignment),pattern="/",3)[,1]

alignment <- 
  # 1. Filter alignment to remove sequences without CD4 cell count
  alignment[names(alignment) %in% metadata$G_NUMBER[!is.na(metadata$type)]]

# 2. Adjust alignment names
full_names <- metadata$new_id[match(names(alignment), metadata$G_NUMBER)]
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names

## ---------------------------

## Save output alignment

seqinr::write.fasta(alignment, names = names(alignment), file.out = OUTPUT_ALIGNMENT)


