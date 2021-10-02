

# loading and combining mutation data 

# remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@*release")
# remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release")

library(tidyverse)
library(cancereffectsizeR)
library(ces.refset.hg19)



# load in TCGA NCI data ----
# data from https://portal.gdc.cancer.gov/files/95258183-63ea-4c97-ae29-1bae9ed06334

NCI_maf <- cancereffectsizeR::preload_maf(maf = here::here("input_data/GDC/gdc_download_20210601_174949.427814/95258183-63ea-4c97-ae29-1bae9ed06334/TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf.gz"), 
                                          refset = ces.refset.hg19,
                                          chain_file = here::here("input_data/hg38ToHg19.over.chain"))

NCI_maf <- NCI_maf %>% 
  mutate(source = "TCGA")



# load in Yale-Gilead data ---- 
# from https://academic.oup.com/jnci/article/110/11/1171/5144449#supplementary-data 


YG_MAF <- readxl::read_excel(path = "input_data/YG/djy168_supp/18-0439R1 SuppTables 1and3_081318.xlsx",sheet = 1) %>%
  filter(tumor_type == "LUSC")

# need to add in chromosome data for cancereffectsizeR

# BiocManager::install("org.Hs.eg.db")
# from https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/org.Hs.eg.db/html/org.Hs.egCHRLOC.html
library(org.Hs.eg.db)
gene_data <- org.Hs.egCHRLOC
mapped_genes <- mappedkeys(gene_data)
mapped_genes_list <- as.list(gene_data[mapped_genes])

YG_MAF$Chromosome <- NA

not_found <- NULL

for(row_ind in 1:nrow(YG_MAF)){
  this_entry <- mapped_genes_list[as.character(YG_MAF$Entrez_Gene_Id[row_ind])]
  if(!is.null(this_entry[[1]])){
  this_chrom <- names(this_entry[[1]])
  YG_MAF$Chromosome[row_ind] <- this_chrom[1]
  }else{
    not_found <- c(not_found,row_ind)
  }
}

YG_MAF <- cancereffectsizeR::preload_maf(maf = YG_MAF,refset = ces.refset.hg19)

YG_MAF <- YG_MAF %>%
  filter(!is.na(Chromosome))

# YG_MAF <- cancereffectsizeR::preload_maf(maf = here::here("input_data/YG/mutationsTN_63_Lung_Squamous_Cell_Carcinoma_(Yale___MD_Anderson).maf"), 
#                                          refset=ces.refset.hg19,
#                                          chr_col = "Chrom")

YG_MAF <- YG_MAF %>% 
  mutate(source = "YG")

# supplemental sequencing data associated with this manuscript
YCCLSB_and_VAH_data <- readr::read_tsv("output_data/supp_tbl_sequence_data.tsv")

all_mafs <- list(NCI_maf = NCI_maf, 
                 YG_MAF = YG_MAF, 
                 YCCLSB_and_VAH_data = YCCLSB_and_VAH_data)


combined_maf <- data.table::rbindlist(all_mafs, idcol = "tumor_source", fill = T)

data.table::fwrite(combined_maf, here::here("output_data/combined_maf.txt"),sep="\t")


combined_maf %>% 
  group_by(tumor_source) %>% 
  count(Unique_Patient_Identifier) %>% 
  count(tumor_source)



