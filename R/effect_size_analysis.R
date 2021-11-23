
# load packages ----- 


library(tidyverse)
library(cancereffectsizeR) # v2.3.3
library(ces.refset.hg19)


# Input data -----


# combined_maf.txt created in the input_data folder. 
# 
# YCCLSB_and_VAH_data <- readr::read_tsv("output_data/supp_tbl_sequence_data.tsv")
# NCI_and_YG_data <- readr::read_tsv("output_data/maf_NCI_YG_data.tsv")
# 
# LUSC_alldata <- rbind(YCCLSB_and_VAH_data,NCI_and_YG_data)

LUSC_alldata <- cancereffectsizeR::preload_maf(maf = "output_data/combined_maf.txt",refset = ces.refset.hg19,keep_extra_columns = "source")


# print(paste("There are ",length(unique(LUSC_alldata$Tumor_Sample_Barcode)),"LUSC tumors in the dataset."))



# Selection intensity run


# start analysis
analysis <- cancereffectsizeR::CESAnalysis(refset = "ces.refset.hg19")

# load in our sequencing data
analysis <- cancereffectsizeR::load_maf(cesa = analysis, maf = LUSC_alldata)



 

# calculate trinucleotide signature weights
analysis <- 
  cancereffectsizeR::trinuc_mutation_rates(cesa = analysis,
                                           signature_set = "COSMIC_v3.1",
                                           signature_exclusions = 
                                             c(cancereffectsizeR::suggest_cosmic_signature_exclusions(cancer_type = "LUSC",
                                                                                                      treatment_naive = T,quiet = T),"SBS89"),
                                           cores = 10,signature_extractor = "deconstructSigs")


# calculate the gene-level mutation rates
analysis <- cancereffectsizeR::gene_mutation_rates(cesa = analysis, 
                                                   covariates = "lung")

# calculate selection 
analysis <-  cancereffectsizeR::ces_variant(cesa = analysis, 
                                            cores=10) 



# load in file from manuscript 
source("R/population_scaled_effect_per_tumor.R")

nrsi_analysis <- population_scaled_effect_per_tumor(ces_output = analysis)

source("R/nrsi_postprocess.R")

nrsi_analysis <- nrsi_remove_double_annotation(analysis = analysis, nrsi_data = nrsi_analysis)




# nrsi_analysis_NFE2L2 <- nrsi_analysis[grep(x = nrsi_analysis$variant, pattern =  "NFE2L2"),]

# save.image("output_data/LUSC_APOBEC_NFE2L2_analysis_output_image.RData")

save_cesa(cesa = analysis,file = "output_data/LUSC_APOBEC_NFE2L2_cesa.rds")
saveRDS(object = nrsi_analysis, file = "output_data/LUSC_APOBEC_NFE2L2_nrsi_analysis.rds")

