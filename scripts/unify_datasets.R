## -----------------------------------------------------------------------------
## Unify variable names and formats for all metadata used in TB-HIV analyses
## 2023-03-02 Etthel Windels
## -----------------------------------------------------------------------------


## Load libraries

suppressMessages(suppressWarnings(require(argparse)))
suppressMessages(suppressWarnings(require(dplyr)))
suppressMessages(suppressWarnings(require(tibble)))
suppressMessages(suppressWarnings(require(stringr)))
suppressMessages(suppressWarnings(require(lubridate)))

## ---------------------------

## Parser

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata  file")
parser$add_argument("--country", type="character", help="Sampling location")
parser$add_argument("--output_metadata", type = "character",
                    help = "Output file for updated metadata")
args <- parser$parse_args()

## ---------------------------

## Read arguments

METADATA <- args$metadata
COUNTRY <- args$country
OUTPUT_METADATA <- args$output_metadata

print(paste("metadata: ", METADATA))
print(paste("country: ", COUNTRY))
print(paste("output metadata: ", OUTPUT_METADATA))

## ---------------------------

## Read files

metadata <- read.table(file = METADATA, header=T,sep= '\t')

## ---------------------------

## Adjust metadata

if (COUNTRY=="South Africa"){
  
  # 1. Remove samples with unknown HIV status, isolation date, or lineage
  metadata <- metadata[!is.na(metadata$HIV) & !is.na(metadata$Isolation_date) & !is.na(metadata$Lineage),]
  metadata <- metadata[metadata$HIV!="Unknown",]
  new_metadata <- 
    metadata %>%
    # 2. Adjust date format
    mutate(Isolation_date = str_split_fixed(date_decimal(2008.00+Isolation_date),pattern=" ",2)[,1]) %>% 
    mutate(Isolation_date=paste(str_split_fixed(Isolation_date,pattern="-",3)[,3],str_split_fixed(Isolation_date,pattern="-",3)[,2],str_split_fixed(Isolation_date,pattern="-",3)[,1], sep='-')) %>% 
    # 3. Adjust lineage format
    mutate(Lineage = paste0("L",substr(Lineage, nchar(Lineage), nchar(Lineage))))
  # 4. Adjust HIV status format
  new_metadata$HIV <- ifelse(as.factor(new_metadata$HIV)=="Pos", 1, 0)
  # 5. Create new_id column
  new_metadata <- new_metadata %>%
    mutate(new_id = paste(G_NUMBER, Isolation_date, HIV, sep = "/"))  
  
} else if (COUNTRY=="Tanzania"){
  
  # 1. Remove samples with unknown HIV status, isolation date, or lineage
  metadata <- metadata[!is.na(metadata$HIV_status) & !is.na(metadata$ids_date_interview) & !is.na(metadata$Lineage),]
  metadata <- metadata[metadata$HIV_status!="",]
  new_metadata <- 
    metadata %>%
    # 2. Adjust date format
    mutate(ids_date_interview = str_split_fixed(ids_date_interview,' ',2)[,1]) %>% 
    mutate(Isolation_date=paste(str_split_fixed(ids_date_interview,pattern="\\.",3)[,1],str_split_fixed(ids_date_interview,pattern="\\.",3)[,2],str_split_fixed(ids_date_interview,pattern="\\.",3)[,3], sep='-')) 
  # 3. Adjust HIV status format
  new_metadata$HIV <- ifelse(as.factor(new_metadata$HIV_status)=="infected", 1, 0)
  # 4. Create new_id column
  new_metadata <- new_metadata %>%
    mutate(new_id = paste(G_NUMBER, Isolation_date, HIV, sep = "/"))  
  
} else if (COUNTRY=="Malawi"){
  
  # 1. Remove samples with unknown HIV status, isolation date, or lineage
  metadata <- metadata[!is.na(metadata$HIV) & !is.na(metadata$date) & !is.na(metadata$LINEAGE),]
  new_metadata <- 
    metadata %>%
    # 2. Adjust date format
    mutate(Isolation_date = ifelse(str_split_fixed(metadata$date,pattern="\\.",3)[,3] > 90, 
                         paste(str_split_fixed(date,pattern="\\.",3)[,1],str_split_fixed(date,pattern="\\.",3)[,2],paste("19",str_split_fixed(date,pattern="\\.",3)[,3],sep="") , sep='-'),
                         paste(str_split_fixed(date,pattern="\\.",3)[,1],str_split_fixed(date,pattern="\\.",3)[,2],paste("20",str_split_fixed(date,pattern="\\.",3)[,3],sep="") , sep='-')))  %>% 
  # 3. Rename lineage column
  rename(Lineage=LINEAGE) %>% 
  # 4. Create new_id column
  new_metadata <- new_metadata %>%
    mutate(new_id = paste(G_NUMBER, Isolation_date, HIV, sep = "/"))  
  
} else if (COUNTRY=="Uganda"){
  
  # 1. Remove samples with unknown HIV status, isolation date, or lineage
  metadata <- metadata[!is.na(metadata$HIV.status) & !is.na(metadata$Year.isolation) & !is.na(metadata$lineage),]
  new_metadata <- 
    metadata %>%
    # 2. Rename G_NUMBER column
    rename(G_NUMBER=g_number) %>% 
    # 3. Adjust date format
    mutate(Isolation_date=paste("01-01",Year.isolation, sep='-'))  %>% 
    # 4. Rename HIV status column
    rename(HIV=HIV.status) %>% 
    # 5. Rename lineage column
    rename(Lineage=lineage) %>% 
    # 6. Rename CD4 count column
    rename(CD_counts=CD4.cells.ul) %>% 
  # 7. Create new_id column
  new_metadata <- new_metadata %>%
    mutate(new_id = paste(G_NUMBER, Isolation_date, HIV, sep = "/"))  
  
}

## ---------------------------

## Save output dataset

write.table(new_metadata, file = OUTPUT_METADATA, row.names=F, quote=F, sep= '\t')
