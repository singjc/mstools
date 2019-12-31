# Description
***
This package contains several functions that can be used for chromatogram visualization and .OSW and .PQP data extraction, and filtering.

# Installation
***
```
library(devtools)
install_github("https://github.com/Roestlab/mstools")
```
# SQL Database File Filtering

The following examples will show examples of filtering **OSW**, **PQP**, and **SQMASS** database files, for a specific list of peptides. This can be useful if you wish to share small example files for a specific few peptides. <br/>
**WARNING:** These filtering methods will overwrite the database! Create a copy of the original database to use as the filtered file.

##  Filer .OSW database File
```
## Specify list of unmodified peptide sequences you wish to filter for. Data will be retained for these peptide sequences in the database
filter_sequences <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
## Specify absolute osw file location
file <- list.files("/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/osw/", pattern = "*.osw", full.names = T )
## Call filterOSWdb
mstools::filterOSWdb(file, filter_sequences)
```

## Filter .PQP database file
```
## Specify list of unmodified peptide sequences you wish to filter for. Data will be retained for these peptide sequences in the database
filter_sequences <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
## Specify absolute pqp file location
file <- list.files("/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/pqp/", pattern = "*.pqp", full.names = T )
## Call filterPQPdb
mstools::filterPQPdb(file, filter_sequences)
```

## Filter .SQMASS database file
```
## Specify list of unmodified peptide sequences you wish to filter for. Data will be retained for these peptide sequences in the database
filter_sequences <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
## Specify list of sqMass chromatogram file(s)
files <- list.files("/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/sqmass/", pattern = "*.sqMass", full.names = T )
## Loop over each file to filter for the selected peptides
for ( file in files ){
  cat("Working on file:", file, "\n")
  ## Call filterSQMASSdb
  mstools::filterSQMASSdb(file, filter_sequences)
}
```

# Extracting data from SQL Database Files into R data.tables

The following examples will show examples of extracting data from **OSW** and **PQP** database files, into R data.tables that can be used for visualization or analysis.

## Extract Data from PQP file
```
## Specify aboslute path to library file
in_pqp <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/mstools/inst/extdata/Synthetic_Dilution_Phosphoproteomics/pqp/psgs_phospho_optimized_decoys.pqp"
## Call getPepLibData_ and store results into pqp_dt.
## Optional: peptide_id parameter can be used to filter for a specific sequence.
pqp_dt <- mstools::getPepLibData_( libfile = in_pqp, peptide_id = "" )
```

## Extract Data from OSW file
```
## Specify absolute path to osw file
in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/mstools/inst/extdata/Synthetic_Dilution_Phosphoproteomics/osw/merged.merged.osw"
## Call getOSWData_ and store results into osw_dt
osw_dt <- mstools::getOSWData_( oswfile = in_osw )
```
**NOTE:** getOSWData_ has other parameters to filter for a specififc run, precursor, peptide, peak group rank, mscore filter, ipf filter, use mscore, use ipf score and decoy filter
```
getOSWData_(
  oswfile,
  run_name = "",
  precursor_id = "",
  peptide_id = "",
  mod_peptide_id = "",
  mod_residue_position = "",
  peak_group_rank_filter = FALSE,
  pep_list = "",
  mscore_filter = "",
  ipf_filter = "",
  ipf_score = FALSE,
  ms2_score = TRUE,
  decoy_filter = TRUE
)
```
