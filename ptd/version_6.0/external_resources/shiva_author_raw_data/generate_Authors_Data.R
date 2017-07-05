
# ################################################################################
#  CREATE AN INPUT TABLE BY MERGING OLD AND NEW DATA FROM AUTHORS
##################################################################################
# Script that takes the input that was provided by the authors and merges it in one
# single input




# Data import
# --------------------------------------------------------------------------------
# Import new datasets of 496 ( with complete profile)
new_data_raw <- xlsx::read.xlsx("../Supplementary_file/SHIVA_496patients.xls", sheetIndex = 1)
cna_data <- xlsx::read.xlsx("../Supplementary_file/SHIVA_496patients.xls", sheetIndex = 2)
# Import the old dataset
# v2 is the one with the modification provided by the shiva authors. 
old_data_raw <- xlsx::read.xlsx("../Supplementary_file/SHIVA_741patients_old_v2.xls", sheetIndex = 1)


# Data reformatting
# --------------------------------------------------------------------------------
# Conver to tibbles
new_data_raw <- tibble::as_tibble(new_data_raw)
cna_data <- tibble::as_tibble(cna_data)
old_data_raw <- tibble::as_tibble(old_data_raw)

# Rename cols
new_data <- new_data_raw %>% 
  dplyr::rename(Patient.ID = ID.SHIVA
                , SHIVA.validated.molecular.alteration = SHIVA.validated.molecular.alteration
                , Drug = Drug 
                , Patient.eligible.for.randomization = Patient.eligible.for.randomization 
                , Signaling.pathway = Signaling.pathway
                , Tumor.cell.content = Tumor.cell.content....
                , DNA.extracted = DNA.extracted
                , Cytoscan.HD = Cytoscan.HD 
                , NGS = NGS
                , IHC.for.HR.expression = IHC.for.HR.expression
  ) %>%
  dplyr::mutate(Patient.ID = as.character(Patient.ID))

old_data <- old_data_raw %>% 
  dplyr::rename(ID = NA.
                , Patient.ID = Patient.ID
                , Cellularity.of.tumor.biopsy = Cellularity.of.tumor.biopsiy
                , DNA.extracted_old = DNA.extracted
                , Localization = Localization
                , Patient.eligible.for.randomization_old = Patient.eligible.for.randomization
                , Mutated_genes = Mutated.gene
                , cp_old = cp
                , pp_old = pp
                , NM_old = NM
                , IHC_ER = X..Expression.RO # Is this ER? 
                , IHC_PR = X..Expression.RP
                , IHC_AR = X..Expression.RA
                , complete.profile_old = complete.profile
                , Reason.for.no.randomization = Reason.for.no.randomization
  ) %>%
  dplyr::select(-ID, -Cellularity.of.tumor.biopsy) %>%
  dplyr::mutate(Patient.ID = as.character(Patient.ID))

# join new and old data by patient IDS
merged_df <- dplyr::left_join(new_data, old_data, by = "Patient.ID")


# Some QC on the input, check if we have perfect correspondence
# --------------------------------------------------------------------------------
#between patient IDs with complete profile
if (nrow(new_data) != sum((new_data$Patient.ID %in% old_data[which(old_data$complete.profile_old == "yes"),]$Patient.ID))) print("we have a problem")

# check for correspondance between rendomization
if (!any(merged_df$Patient.eligible.for.randomization == toupper(merged_df$Patient.eligible.for.randomization_old))) print("we have a problem")

# check for correspondance between DNA extracted
if (!any(merged_df$DNA.extracted == toupper(merged_df$DNA.extracted_old))) print("we have a problem")

# Ok, now that we have run the checks, remove useless fields
merged_df <- merged_df %>% dplyr::select(-DNA.extracted_old
                                         , -Patient.eligible.for.randomization_old
                                         , -complete.profile_old)


# Convert turmor localization to Tumor ID
# --------------------------------------------------------------------------------
# PROBLEM: The authors provided the localization but not the tumor type. We need to convert this into a tumor_id 
# compatible with the TCGA id in use. In order to address this problem we have created a conversion table.

# Import conversion table
conversiontable <- read.xlsx2("../../Supplementary_file/SHIVA_conversiontable.xlsx", sheetIndex = "Sheet1")
# Create a mapping array for the conversion values
mapping <- tolower(conversiontable$TCGA.1)
names(mapping) <- conversiontable$Localization
# Add TCGA turmor ID information
merged_df <- merged_df %>% dplyr::mutate(TCGA_tumor_ID = plyr::mapvalues(Localization
                                                                         , from = names(mapping)
                                                                         , to = as.character(mapping)
)
) 

# Convert SHIVA.validated.molecular.alteration into 3 columns
# --------------------------------------------------------------------------------
# col1 = NGS_report: short list of NGS selected genes.
# col2 = CNA_report: short list of CNA selected genes.
# col3 = OR_report  Short list of selected genes for hormon therapy.
NGS_report <- sapply(merged_df$SHIVA.validated.molecular.alteration
                     , function(x) return(grepl("mutation", as.character(x)))) 
OR_report <- sapply(merged_df$SHIVA.validated.molecular.alteration
                    , function(x) return(grepl("expression", as.character(x))))
CNA_report <- sapply(merged_df$SHIVA.validated.molecular.alteration
                     , function(x) return(grepl("gain|loss|deletion|amplification|activation", as.character(x)))) 

# Merge info back to the main dataframe
merged_df <- data.frame(merged_df, OR_report, CNA_report, NGS_report)


# Preview dataframe
# --------------------------------------------------------------------------------
datatable(merged_df # ADD BUTTONS TO THE TABLE
          , extensions = 'Buttons'
          , options = list(
            dom = 'lBfrtip'
            , buttons = c('copy', 'csv', 'excel')
          )
          , caption = ""
)

# export
write.xlsx(merged_df, file = "../../Supplementary_file/SHIVA_authorTable.xlsx")