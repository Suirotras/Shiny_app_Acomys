### load libraries
library(dplyr)
library(readr)
library(purrr)
library(magrittr)

### set working directory
setwd(file.path("~", "Documents", "Github", "Shiny_app_Acomys"))

### load aBSREL data

load(file = file.path("data_files", "aBSREL_analysis_pval_count.RData"))

### Load MEME data

# MEME data with substitution columns
load(file = file.path("data_files", "MEME_substitutions_small.RData"))

# EBF data calculated from the MEME data
load(file = file.path("data_files", "MEME_EBF_values_simple.RData"))

# rename data
MEME_data <- MEME_simple_nsites_subs_small
rm(MEME_simple_nsites_subs_small)

aBSREL_data <- results_tibble_count_cor
rm(results_tibble_count_cor)

# Combine EBF_data with MEME_data 
MEME_data %<>%
  left_join(MEME_EBF_simple)

### Split the data in multiple files on a per row basis

# Function that takes the complete MEME or aBSREL data object, and creates a
# new RData object per row
#
# args:
#   RData_object: An RData object to split into multiple files using the rows
#   column_name: the column name of the column that will be used for the
#   splitting procedure
#   RData_basename: The basename that the separate files should have
#   subdir_name: The subdirectory that is created to store the generated RData
#   files
split_rows <- function(RData_object, column_name = "genename", 
                       RData_basename = "MEME_data_",
                       subdir_name = file.path("~", "Documents", "Github",
                                               "Shiny_app_Acomys",
                                               "Shiny_Acomys", "data",
                                               "MEME_data")) {
  column_vector <- RData_object[[column_name]]
  
  dir.create(subdir_name)
  
  # Loop through the rows and save each one as a separate RData file
  for (row_id in seq_along(RData_object[[column_name]])) {
    
    RData_row <- RData_object[row_id,]
    RData_genename <- RData_row[[column_name]]
    
    save(RData_row, file = file.path(subdir_name, paste0(RData_basename, RData_genename, ".RData")))
    cat(paste0("Saved ", RData_genename, " row as a separate RData object\n\n"))
  }
}

# get separate files for MEME data rows
split_rows(MEME_data)

# simplify the aBSREL RData
aBSREL_data <- aBSREL_data %>%
  select(-all_of(c("omega_rates", "tree", "TRUE_count", "FALSE_count")))

save(aBSREL_data, file = file.path("~", "Documents", "Github",
                                   "Shiny_app_Acomys",
                                   "Shiny_Acomys", "data",
                                   "aBSREL_data.RData"))
