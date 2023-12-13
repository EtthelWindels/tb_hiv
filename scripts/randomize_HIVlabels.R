## -----------------------------------------------------------------------------
## Randomize HIV status labels of sequences
## 2023-05-11 Etthel Windels
## -----------------------------------------------------------------------------


## Load libraries

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(seqinr)))
suppressMessages(suppressWarnings(require(stringr)))

## ---------------------------

## Parser

parser <- argparse::ArgumentParser()
parser$add_argument("--alignment", type="character", help="Alignment file")
parser$add_argument("--metadata", type="character", help="Metadata file")
parser$add_argument("--country", type="character", help="Sampling location")
parser$add_argument("--randomization_type", type="character", help="permutation/reassignment")

args <- parser$parse_args()

## ---------------------------

## Read arguments

ALIGNMENT <- args$alignment
METADATA <- args$metadata
COUNTRY <- args$country
RANDOMIZATION_TYPE <- args$randomization_type

print(paste("alignment: ", ALIGNMENT))
print(paste("metadata: ", METADATA))
print(paste("country: ", COUNTRY))
print(paste("randomization_type: ", RANDOMIZATION_TYPE))

## ---------------------------

## Read files

alignment <- seqinr::read.fasta(file = ALIGNMENT, seqtype="DNA", forceDNAtolower = F)
metadata <- read.table(file = METADATA, header=T,sep= '\t')

## ---------------------------

## Filter metadata

names(alignment) <- str_split_fixed(names(alignment),pattern="/",3)[,1]
metadata <- metadata[metadata$G_NUMBER %in% names(alignment),]

## ---------------------------

## Get HIV frequencies in general population
# average freq over sampling period (based on data from The World Bank)

if (COUNTRY=="South Africa"){
  freq = 0.17
} else if (COUNTRY=="Tanzania"){
  freq = 0.05
} else if (COUNTRY=="Uganda"){
  freq = 0.07
} else if (COUNTRY=="Malawi"){
  freq = 0.13
}

## ---------------------------

## Label randomization

for (i in 1:10){
  
  if (RANDOMIZATION_TYPE=="permutation"){
    
    set.seed(i)
    nrHIV <- length(metadata$HIV[metadata$HIV==1])
    metadata$HIV_permuted <- rep(0, dim(metadata)[1])
    HIVsample <- sample(metadata$G_NUMBER, nrHIV, replace=F)
    metadata$HIV_permuted[metadata$G_NUMBER %in% HIVsample] <- 1
    metadata$HIV_permuted <- as.factor(metadata$HIV_permuted)
    output_name <- paste0(str_split_fixed(ALIGNMENT,'\\.',2)[1],i,"_permuted.fasta")
    
  } else if (RANDOMIZATION_TYPE=="reassignment"){
    
    set.seed(i)
    nrHIV <- round(length(metadata$HIV)*freq)
    metadata$HIV_permuted <- rep(0, dim(metadata)[1])
    HIVsample <- sample(metadata$G_NUMBER, nrHIV, replace=F)
    metadata$HIV_permuted[metadata$G_NUMBER %in% HIVsample] <- 1
    metadata$HIV_permuted <- as.factor(metadata$HIV_permuted)
    output_name <- paste0(str_split_fixed(ALIGNMENT,'\\.',2)[1],i,"_reassigned.fasta")
  }
  
  # Adjust new_id column
  metadata <- metadata %>%
    mutate(new_id = paste(G_NUMBER, Isolation_date, HIV_permuted, sep = "/"))
  
  # Adjust alignment
  gnumber <- str_split_fixed(names(alignment),'/',3)[,1]
  full_names <- metadata$new_id[match(gnumber, metadata$G_NUMBER)]
  if (any(is.na(full_names))) {
    stop("Not all strains have full names in metadata.")
  }
  names(alignment) <- full_names
  
  # Save output alignment
  seqinr::write.fasta(alignment, names = names(alignment), file.out = output_name)
  
}
