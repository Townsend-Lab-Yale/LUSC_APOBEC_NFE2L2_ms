Analyses associated with 

**Cannataro et al. APOBEC mutagenesis and selection for NFE2L2 contribute to the origin of lung squamous-cell carcinoma (submitted for review)**

To reproduce the calculations within the manuscript, you need to download the raw data from their free open sources. They should be stored within this directory as called within `input_data/LUSC_data_combine_ms.R`

1. Download `TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf.gz` from `https://portal.gdc.cancer.gov/files/95258183-63ea-4c97-ae29-1bae9ed06334`
2. Download `18-0439R1 SuppTables 1and3_081318.xlsx` from `https://academic.oup.com/jnci/article/110/11/1171/5144449#supplementary-data`
3. Download the supplemental sequencing data associated with this manuscript 

Then you may run the analysis by

1. Running `input_data/LUSC_data_combine_ms.R`
2. Running `LUSC_APOBEC_NFE2L2_analysis.Rmd` 



