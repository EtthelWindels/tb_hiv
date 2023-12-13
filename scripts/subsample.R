## -----------------------------------------------------------------------------
## Subsample alignments with >400 sequences,
## for computational feasibility of BEAST analyses
## 2023-03-24 Etthel Windels
## -----------------------------------------------------------------------------


## Load libraries

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(seqinr)))
suppressMessages(suppressWarnings(require(ape)))

## ---------------------------

## Parser

parser <- argparse::ArgumentParser()
parser$add_argument("--path_to_alignments", type="character", 
                    help="Path to alignment files, folder should only contain fasta or nexus files")
parser$add_argument("--output_path", type="character", help="Path to output alignment files")

args <- parser$parse_args()

## ---------------------------

## Read arguments

PATH_ALIGNMENT <- args$path_to_alignments
OUTPUT_PATH <- args$output_path

print(paste("path to alignments: ", PATH_ALIGNMENT))
print(paste("path to output: ", OUTPUT_PATH))

## ---------------------------

## Subsample large alignments

alignment_list <- list.files(path=PATH_ALIGNMENT)

for (i in alignment_list){
  alignment <- tryCatch({alignment <- seqinr::read.fasta(file = paste(PATH_ALIGNMENT,i,sep="/"), seqtype="DNA", forceDNAtolower = F)},
           error = function(e){
             alignment <- ape::read.nexus.data(file = paste(PATH_ALIGNMENT,i,sep="/"))
           })
  if (length(names(alignment))>400){
    id <- sample(names(alignment), 400, replace=F)  
    alignment_new <- alignment[names(alignment) %in% id]
    write.fasta(alignment_new,names=names(alignment_new),file.out=paste0(OUTPUT_PATH,"/",str_split_fixed(i,'\\.',2)[1],"_subsampled.fasta"))
  }
}