## ---- cache= TRUE--------------------------------------------------------
## Load the saas library
library(saas)
## Load some convenience functions
library(dplyr, quietly = TRUE, warn.conflicts = FALSE) 

## Location of the zipped data files
zip_file_path = system.file("extdata", "extdata.zip", package = "saas")
## Unzip and get the (temporary) location of the mzid file with the MS-GF+ search results.
mzid_file_path = unzip(zip_file_path, 'pyrococcus.mzid',exdir = tempdir())

## Read and parse the mzid file
data = parse_msgf_mzid(mzid_file_path)

glimpse(data)

## ---- cache= TRUE--------------------------------------------------------
## Unzip and get the (temporary) location of the file with fasta headers containing 
## the protein_ids from the protein subset of interest.
fasta_file_path = unzip(zip_file_path, 'transferase_activity_[GO:0016740].fasta', exdir = tempdir())

## Preprocess the data.
data_prep = preprocess(data, is_subset = fasta_file_path)

glimpse(data_prep)

## ---- cache= TRUE,fig.width=10, fig.height= 6, out.width= 700------------
diagnostics = plot_diag(data_prep, score_higher = FALSE)
diagnostics$all

## ---- cache= TRUE--------------------------------------------------------
data_results = calculate_fdr(data_prep, score_higher = FALSE)
glimpse(data_results)

## ---- cache= TRUE--------------------------------------------------------
count(data_results, subset, decoy)

## ---- cache= TRUE--------------------------------------------------------
results_1_FDR = filter(data_results, subset, !decoy, FDR_stable >= .01) 
glimpse(results_1_FDR)

## ------------------------------------------------------------------------
sessionInfo()

