# Load required libraries 
required_packages <- c('Matrix', 'dplyr', 'Seurat', 'patchwork', 'SeuratData', 'ggplot2', 'RColorBrewer', 'ggalluvial', 'ggrepel')

invisible(lapply(required_packages, library, character.only = T))

# obtain directories

# wget "https://zenodo.org/records/14199536/files/scAtlas.rds.gz"
# gunzip scAtlas.rds.gz

# https://aacrjournals.org/clincancerres/article/31/4/756/751743/Human-Pancreatic-Cancer-Single-Cell-Atlas-Reveals
# https://github.com/PDAC-MULTIOMIC/PDAC_Atlas

data_dir <- "path where the panc_ca_atlas and required files can be loaded"

####################################################################################################

# save pdac metadata separately
pdac_ref <- readRDS(paste0(data_dir, "scAtlas.rds"))
pdac_metadata <- pdac_ref@meta.data
write.csv(pdac_metadata, paste0(data_dir, "scAtlas_metadata.csv"))

##############################################################################################

# read metadata
pdac_metadata <- read.csv(paste0(data_dir, "scAtlas_metadata.csv"), row.names = 1)
colnames(pdac_metadata)
# [1] "nCount_RNA"              "nFeature_RNA"            "percent.mt"              "seurat_clusters"         "Count"                   "Study..Citation..PMID."  "GSE.SRA..Study."         "Name"                    "If.metastatic..location"
# [10] "Clusters"                "Treatment"               "DiseaseState"            "TreatmentType"   

table(pdac_metadata$GSE.SRA..Study.)
# EGAS00001002543       GSE154778       GSE155698       GSE156405       GSE158356       GSE194247       GSE202051       GSE205013       GSE211644       GSE229413 phs001840.v1.p1     PRJCA001063 
#           76094           20364           23991           17809            2746           29622          147643          167366           40971           33309           18285          147907 

table(pdac_metadata$If.metastatic..location)
# Liver         Lung      Omentum   Peritoneum Vaginal apex 
# 53145         1225          138          459         1475 

##############################################################################################
# add conflict info
pdac_metadata$sample_check <- NA
pdac_metadata$tissue_check <- NA
pdac_metadata$treatment_check <- NA
pdac_metadata$overall_check <- NA
##############################################################################################
# EGAS
study_id <- "EGAS00001002543"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# https://www.clinicaltrials.gov/study/NCT02750657
# inclusion criteria: no prior treatment
# primary or liver metastatic samples
# identification due to missing id links not possible
# sample
pdac_metadata$sample_check[study_index] <- "?"
# tissue
pdac_metadata$tissue_check[study_index] <- "?"
# treatment
pdac_metadata$treatment_check[study_index] <- "?"
##############################################################################################
# gse154778
study_id <- "GSE154778"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# sample: file name correlated information used && T1 != P01,K16733 or MET01, Y00008
# tissue: P02,Y00006, Primary != Liver Metastasis & MET03, Y00014, Omentum Metastasis != Liver Metastasis)
# treatment: no info available
# sample
pdac_metadata$sample_check[study_index] <- "conflict"
name_index <- which(pdac_metadata$Name == "T1")
pdac_metadata$sample_check[intersect(study_index, name_index)] <- "missing"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
name_index <- which(gsub("_.*", "", pdac_metadata$Name) %in% c("Y00006", "Y00014"))
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict" 
# treatment
pdac_metadata$treatment_check[study_index] <- "?"
##############################################################################################
# gse155698
study_id <- "GSE155698"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# multiple srr run ids match one sample: 12 samples, 20 srr ids, claimed to match 20 samples
# tissue match
# treatment naïve
# sample
pdac_metadata$sample_check[study_index] <- "conflict"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# gse156405
study_id <- "GSE156405"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# SRR IDs used & old srr run id: SRR12467643 -> new: SRR12467650 
# tissue: SRR24437450 Primary != Peritoneal Metastasis & old "SRR12467643" -> deprecated: replaced by SRR12467650, https://www.ncbi.nlm.nih.gov/sra/?term=SRR12467643 Peritoneal Metastasis != Primary
# treatment (4 out of 9 patients were pre-treated)
# sample
pdac_metadata$sample_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name == "SRR12467643")
pdac_metadata$sample_check[intersect(study_index, name_index)] <- "conflict"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name %in% c("SRR12467643", "SRR24437450"))
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name %in% c("SRR12467642", "SRR12467647", "SRR12467649", "SRR12467643"))
pdac_metadata$treatment_check[intersect(study_index, name_index)] <- "conflict"
##############################################################################################
# gse158356
study_id <- "GSE158356"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# sample names used
# tissue info: LiverMetD != Omentum Metastasis
# no treatment info available?
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name == "LiverMetD")
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "?"
##############################################################################################
# gse194247
study_id <- "GSE194247"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# amplification round info used
# tissue: all samples included cd45neg cells only
# treatment naïve
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# gse202051
study_id <- "GSE202051"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# srr ids
# tissue info: SRR19039355 & SRR19039354 organoid cultures != Primary)
# treatment correct
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index]  <- "match"
name_index <- which(pdac_metadata$Name %in% c("SRR19039355","SRR19039354"))
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# gse205013
study_id <- "GSE205013"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# sample (gsm ids) & tissue match
# treatment: GSM6204119 with prior treatment in supplementary
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name == "GSM6204119")
pdac_metadata$treatment_check[intersect(study_index, name_index)] <- "conflict"
##############################################################################################
# GSE211644
study_id <- "GSE211644"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# samples: srr ids
# tissue: all CD3+ enriched or Infusion Products from CD3+ enriched samples
# treatment: collected before treatment, correct
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "conflict"
# treatment info
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# gse229413
study_id <- "GSE229413"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# samples: gsm id derived names, no samn available
# tissue & treatment correct
# treatment
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# phs
study_id <- "phs001840.v1.p1"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# srr ids used
# one CAF enriched sample (SRR9274540)
# treatment naïve
# sample
pdac_metadata$sample_check[study_index] <- "match"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
name_index <- which(pdac_metadata$Name == "SRR9274540")
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################
# peng
study_id <- "PRJCA001063"
study_index <- which(pdac_metadata$GSE.SRA..Study.== study_id)
# sample 23 twice, 15 missing
# tissue: healthy, not adjacent normal
# treatment naïve
# sample
pdac_metadata$sample_check[study_index] <- "match"
name_index <-  which(pdac_metadata$Name == "new_CRR034518")
pdac_metadata$sample_check[intersect(study_index, name_index)] <- "conflict"
# tissue
pdac_metadata$tissue_check[study_index] <- "match"
healthy_id <- paste0("new_CRR0345", c(20:30))
name_index <- which(pdac_metadata$Name %in% healthy_id)
pdac_metadata$tissue_check[intersect(study_index, name_index)] <- "conflict"
# treatment
pdac_metadata$treatment_check[study_index] <- "match"
##############################################################################################

# Add overall match status
pdac_metadata <- pdac_metadata %>%
  mutate(overall_check = case_when(
    grepl("conflict", sample_check) ~ "conflict (ID)",
    grepl("conflict", treatment_check) ~ "conflict (treatment)",
    grepl("conflict", tissue_check) ~ "conflict (tissue)",
    is.na(sample_check) | is.na(tissue_check) | is.na(treatment_check) |
      grepl("?", sample_check, fixed = TRUE) |
      grepl("?", tissue_check, fixed = TRUE) |
      grepl("?", treatment_check, fixed = TRUE) ~ "not_checked",
    TRUE ~ "match"))

##############################################################################################

# check tables
table(pdac_metadata$sample_check, useNA="ifany")
#     ? conflict    match  missing 
# 76094    50727   598922      364

table(pdac_metadata$tissue_check, useNA="ifany")
# ? conflict    match 
# 76094   133203   516810 
       
table(pdac_metadata$treatment_check, useNA="ifany")
#     ? conflict    match 
# 99204    10489   616414 

table(pdac_metadata$overall_check, useNA="ifany")
# conflict (ID)    conflict (tissue) conflict (treatment)                match          not_checked 
#         50727               130856                 8282               457176                79066 

##############################################################################################

# create plots
# define colors
col_list <- c(
  "?"                = "#CACACA",
  "missing"          = "#A9A9A9",
  "not_checked"      = "#6F6A6A",    
  "match"            = "#2DC32D",
  "conflict (ID)"    = "#FFC7F3",
  "conflict (tissue)"= "#FF9BBF",
  "conflict (treatment)"= "#FC73FF",
  "conflict"         = "#C159D0")

# sankey
# overall check across studies & cell types
pal_blue <- colorRampPalette(brewer.pal(5,"Blues"))
pal_red <- colorRampPalette(brewer.pal(5,"YlOrRd"))
df_alluv <- pdac_metadata %>%
  dplyr::count(GSE.SRA..Study., overall_check, Clusters)

# df_alluv$overall_check <- gsub("conflict.*$", "conflict",df_alluv$overall_check)

p <- ggplot(df_alluv,
             aes(axis1 = GSE.SRA..Study., axis2 =overall_check, axis3=Clusters, y = n)) +
  geom_alluvium(aes(fill =overall_check)) +
  geom_stratum(width = 1/8, fill = c(rev(pal_blue(length(unique(df_alluv$GSE.SRA..Study.)))), rev(col_list[match(names(table(df_alluv$overall_check)),names(col_list))]), rev(pal_red(length(unique(df_alluv$Clusters))))), color = "black") +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), nudge_x = 0.3) +
  scale_fill_manual(values =col_list[match(names(table(df_alluv$overall_check)),names(col_list))], labels = c(paste0(names(table(pdac_metadata$overall_check))[match(names(table(df_alluv$overall_check)),names(table(pdac_metadata$overall_check)))], 
                                                                                                               ": n=",table(pdac_metadata$overall_check)[match(names(table(df_alluv$overall_check)),names(table(pdac_metadata$overall_check)))], 
                                                                                                               " (", (round(prop.table(table(pdac_metadata$overall_check)), 3)*100)[match(names(table(df_alluv$overall_check)),names(table(pdac_metadata$overall_check)))], "%)"))) +
  labs(title = "Sankey Plots showing the Relation between Overall Check, Studies and Cell Types",
       x = "", y = "Sample counts") + scale_x_continuous(labels = c("Study_ID", "overall_check", "Cell_Type"), breaks = c(1,2,3))

ggsave(paste0(data_dir,"sankey_study_cell_type_overall_check.png"), plot=p, width = 45, height = 25, units = "cm")

##############################################################################################

sessionInfo()
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8    LC_MONETARY=English_Germany.utf8
# [4] LC_NUMERIC=C                     LC_TIME=English_Germany.utf8    
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.5                RColorBrewer_1.1-3           ggalluvial_0.12.5           
# [4] tidyr_1.3.0                  stringr_1.5.1                rentrez_1.2.3               
# [7] GEOquery_2.70.0              Biobase_2.62.0               SRAdb_1.64.0                
# [10] RCurl_1.98-1.14              graph_1.80.0                 BiocGenerics_0.48.1         
# [13] RSQLite_2.3.4                ggplot2_3.5.1                stxBrain.SeuratData_0.1.2   
# [16] pbmcsca.SeuratData_3.0.0     pbmcref.SeuratData_1.0.0     pancreasref.SeuratData_1.0.0
# [19] panc8.SeuratData_3.0.2       adiposeref.SeuratData_1.0.0  SeuratData_0.2.2.9001       
# [22] patchwork_1.2.0              Seurat_5.1.0                 SeuratObject_5.0.2          
# [25] sp_2.1-2                     dplyr_1.1.4                  Matrix_1.6-4                
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                      matrixStats_1.2.0             spatstat.sparse_3.0-3        
# [4] bitops_1.0-7                  devtools_2.4.5                httr_1.4.7                   
# [7] doParallel_1.0.17             profvis_0.3.8                 tools_4.3.0                  
# [10] sctransform_0.4.1             utf8_1.2.4                    R6_2.5.1                     
# [13] lazyeval_0.2.2                uwot_0.2.2                    withr_3.0.0                  
# [16] urlchecker_1.0.1              prettyunits_1.2.0             gridExtra_2.3                
# [19] parallelDist_0.2.6            progressr_0.14.0              textshaping_0.3.7            
# [22] argparse_2.2.2                cli_3.6.2                     pacman_0.5.1                 
# [25] formatR_1.14                  spatstat.explore_3.2-5        fastDummies_1.7.3            
# [28] sandwich_3.1-0                labeling_0.4.3                mvtnorm_1.2-4                
# [31] spatstat.data_3.0-3           readr_2.1.5                   ggridges_0.5.5               
# [34] pbapply_1.7-2                 copykat_1.1.0                 systemfonts_1.3.1            
# [37] Rsamtools_2.18.0              R.utils_2.12.3                parallelly_1.36.0            
# [40] sessioninfo_1.2.2             limma_3.58.1                  readxl_1.4.3                 
# [43] rstudioapi_0.15.0             FNN_1.1.4                     BiocIO_1.12.0                
# [46] generics_0.1.3                gtools_3.9.5                  ica_1.0-3                    
# [49] spatstat.random_3.2-2         futile.logger_1.4.3           fansi_1.0.6                  
# [52] S4Vectors_0.40.2              abind_1.4-5                   R.methodsS3_1.8.2            
# [55] lifecycle_1.0.4               edgeR_4.0.6                   infercnv_1.18.1              
# [58] multcomp_1.4-25               yaml_2.3.8                    SummarizedExperiment_1.32.0  
# [61] gplots_3.1.3                  SparseArray_1.2.3             BiocFileCache_2.10.1         
# [64] Rtsne_0.17                    grid_4.3.0                    blob_1.2.4                   
# [67] promises_1.2.1                ExperimentHub_2.10.0          crayon_1.5.2                 
# [70] miniUI_0.1.1.1                lattice_0.22-5                cowplot_1.1.3                
# [73] GenomicFeatures_1.54.1        KEGGREST_1.42.0               pillar_1.9.0                 
# [76] knitr_1.46                    GenomicRanges_1.54.1          rjson_0.2.21                 
# [79] future.apply_1.11.1           codetools_0.2-19              leiden_0.4.3.1               
# [82] glue_1.6.2                    data.table_1.15.4             remotes_2.4.2.1              
# [85] vctrs_0.6.5                   png_0.1-8                     spam_2.10-0                  
# [88] cellranger_1.1.0              gtable_0.3.5                  cachem_1.0.8                 
# [91] xfun_0.44                     S4Arrays_1.2.0                mime_0.12                    
# [94] dlm_1.1-6                     libcoin_1.0-10                coda_0.19-4                  
# [97] survival_3.8-3                SingleCellExperiment_1.24.0   iterators_1.0.14             
# [100] DOTr_0.0.0.9000               statmod_1.5.0                 TH.data_1.1-2                
# [103] interactiveDisplayBase_1.40.0 ellipsis_0.3.2                fitdistrplus_1.1-11          
# [106] ROCR_1.0-11                   nlme_3.1-164                  usethis_2.2.2                
# [109] bit64_4.0.5                   progress_1.2.3                filelock_1.0.3               
# [112] RcppAnnoy_0.0.22              GenomeInfoDb_1.38.5           irlba_2.3.5.1                
# [115] KernSmooth_2.23-22            colorspace_2.1-0              DBI_1.2.3                    
# [118] celldex_1.12.0                tidyselect_1.2.1              bit_4.0.5                    
# [121] compiler_4.3.0                curl_5.2.0                    xml2_1.3.6                   
# [124] DelayedArray_0.28.0           plotly_4.10.4                 rtracklayer_1.62.0           
# [127] scales_1.3.0                  caTools_1.18.2                lmtest_0.9-40                
# [130] rappdirs_0.3.3                digest_0.6.33                 goftest_1.2-3                
# [133] spatstat.utils_3.1-0          rmarkdown_2.27                XVector_0.42.0               
# [136] htmltools_0.5.8.1             pkgconfig_2.0.3               sparseMatrixStats_1.14.0     
# [139] MatrixGenerics_1.14.0         dbplyr_2.4.0                  fastmap_1.1.1                
# [142] rlang_1.1.2                   htmlwidgets_1.6.4             shiny_1.8.0                  
# [145] DelayedMatrixStats_1.24.0     farver_2.1.2                  zoo_1.8-12                   
# [148] jsonlite_1.8.8                BiocParallel_1.36.0           R.oo_1.25.0                  
# [151] magrittr_2.0.3                modeltools_0.2-23             GenomeInfoDbData_1.2.11      
# [154] dotCall64_1.1-1               munsell_0.5.1                 Rcpp_1.0.11                  
# [157] ape_5.7-1                     reticulate_1.40.0             stringi_1.8.3                
# [160] zlibbioc_1.48.0               MASS_7.3-60                   AnnotationHub_3.10.1         
# [163] plyr_1.8.9                    pkgbuild_1.4.3                parallel_4.3.0               
# [166] listenv_0.9.0                 deldir_2.0-2                  Biostrings_2.70.1            
# [169] splines_4.3.0                 tensor_1.5                    hms_1.1.3                    
# [172] locfit_1.5-9.8                fastcluster_1.2.3             igraph_2.0.3                 
# [175] spatstat.geom_3.2-7           RcppHNSW_0.5.0                reshape2_1.4.4               
# [178] biomaRt_2.58.0                stats4_4.3.0                  pkgload_1.3.3                
# [181] futile.options_1.0.1          BiocVersion_3.18.1            XML_3.99-0.16                
# [184] evaluate_0.23                 RcppParallel_5.1.7            lambda.r_1.2.4               
# [187] BiocManager_1.30.22           tzdb_0.4.0                    phyclust_0.1-34              
# [190] foreach_1.5.2                 httpuv_1.6.13                 RANN_2.6.1                   
# [193] purrr_1.0.2                   polyclip_1.10-6               future_1.33.1                
# [196] scattermore_1.2               coin_1.4-3                    xtable_1.8-4                 
# [199] restfulr_0.0.15               RSpectra_0.16-1               later_1.3.2                  
# [202] ragg_1.2.7                    rjags_4-15                    viridisLite_0.4.2            
# [205] tibble_3.2.1                  GenomicAlignments_1.38.1      memoise_2.0.1                
# [208] AnnotationDbi_1.64.1          IRanges_2.36.0                cluster_2.1.6                
# [211] globals_0.16.2               
# > sessionInfo()
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_Germany.utf8  LC_CTYPE=English_Germany.utf8    LC_MONETARY=English_Germany.utf8 LC_NUMERIC=C                     LC_TIME=English_Germany.utf8    
# 
# time zone: Europe/Berlin
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggrepel_0.9.5                RColorBrewer_1.1-3           ggalluvial_0.12.5            tidyr_1.3.0                  stringr_1.5.1                rentrez_1.2.3                GEOquery_2.70.0             
# [8] Biobase_2.62.0               SRAdb_1.64.0                 RCurl_1.98-1.14              graph_1.80.0                 BiocGenerics_0.48.1          RSQLite_2.3.4                ggplot2_3.5.1               
# [15] stxBrain.SeuratData_0.1.2    pbmcsca.SeuratData_3.0.0     pbmcref.SeuratData_1.0.0     pancreasref.SeuratData_1.0.0 panc8.SeuratData_3.0.2       adiposeref.SeuratData_1.0.0  SeuratData_0.2.2.9001       
# [22] patchwork_1.2.0              Seurat_5.1.0                 SeuratObject_5.0.2           sp_2.1-2                     dplyr_1.1.4                  Matrix_1.6-4                
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                      matrixStats_1.2.0             spatstat.sparse_3.0-3         bitops_1.0-7                  devtools_2.4.5                httr_1.4.7                    doParallel_1.0.17            
# [8] profvis_0.3.8                 tools_4.3.0                   sctransform_0.4.1             utf8_1.2.4                    R6_2.5.1                      lazyeval_0.2.2                uwot_0.2.2                   
# [15] withr_3.0.0                   urlchecker_1.0.1              prettyunits_1.2.0             gridExtra_2.3                 parallelDist_0.2.6            progressr_0.14.0              textshaping_0.3.7            
# [22] argparse_2.2.2                cli_3.6.2                     pacman_0.5.1                  formatR_1.14                  spatstat.explore_3.2-5        fastDummies_1.7.3             sandwich_3.1-0               
# [29] labeling_0.4.3                mvtnorm_1.2-4                 spatstat.data_3.0-3           readr_2.1.5                   ggridges_0.5.5                pbapply_1.7-2                 copykat_1.1.0                
# [36] systemfonts_1.3.1             Rsamtools_2.18.0              R.utils_2.12.3                parallelly_1.36.0             sessioninfo_1.2.2             limma_3.58.1                  readxl_1.4.3                 
# [43] rstudioapi_0.15.0             FNN_1.1.4                     BiocIO_1.12.0                 generics_0.1.3                gtools_3.9.5                  ica_1.0-3                     spatstat.random_3.2-2        
# [50] futile.logger_1.4.3           fansi_1.0.6                   S4Vectors_0.40.2              abind_1.4-5                   R.methodsS3_1.8.2             lifecycle_1.0.4               edgeR_4.0.6                  
# [57] infercnv_1.18.1               multcomp_1.4-25               yaml_2.3.8                    SummarizedExperiment_1.32.0   gplots_3.1.3                  SparseArray_1.2.3             BiocFileCache_2.10.1         
# [64] Rtsne_0.17                    grid_4.3.0                    blob_1.2.4                    promises_1.2.1                ExperimentHub_2.10.0          crayon_1.5.2                  miniUI_0.1.1.1               
# [71] lattice_0.22-5                cowplot_1.1.3                 GenomicFeatures_1.54.1        KEGGREST_1.42.0               pillar_1.9.0                  knitr_1.46                    GenomicRanges_1.54.1         
# [78] rjson_0.2.21                  future.apply_1.11.1           codetools_0.2-19              leiden_0.4.3.1                glue_1.6.2                    data.table_1.15.4             remotes_2.4.2.1              
# [85] vctrs_0.6.5                   png_0.1-8                     spam_2.10-0                   cellranger_1.1.0              gtable_0.3.5                  cachem_1.0.8                  xfun_0.44                    
# [92] S4Arrays_1.2.0                mime_0.12                     dlm_1.1-6                     libcoin_1.0-10                coda_0.19-4                   survival_3.8-3                SingleCellExperiment_1.24.0  
# [99] iterators_1.0.14              DOTr_0.0.0.9000               statmod_1.5.0                 TH.data_1.1-2                 interactiveDisplayBase_1.40.0 ellipsis_0.3.2                fitdistrplus_1.1-11          
# [106] ROCR_1.0-11                   nlme_3.1-164                  usethis_2.2.2                 bit64_4.0.5                   progress_1.2.3                filelock_1.0.3                RcppAnnoy_0.0.22             
# [113] GenomeInfoDb_1.38.5           irlba_2.3.5.1                 KernSmooth_2.23-22            colorspace_2.1-0              DBI_1.2.3                     celldex_1.12.0                tidyselect_1.2.1             
# [120] bit_4.0.5                     compiler_4.3.0                curl_5.2.0                    xml2_1.3.6                    DelayedArray_0.28.0           plotly_4.10.4                 rtracklayer_1.62.0           
# [127] scales_1.3.0                  caTools_1.18.2                lmtest_0.9-40                 rappdirs_0.3.3                digest_0.6.33                 goftest_1.2-3                 spatstat.utils_3.1-0         
# [134] rmarkdown_2.27                XVector_0.42.0                htmltools_0.5.8.1             pkgconfig_2.0.3               sparseMatrixStats_1.14.0      MatrixGenerics_1.14.0         dbplyr_2.4.0                 
# [141] fastmap_1.1.1                 rlang_1.1.2                   htmlwidgets_1.6.4             shiny_1.8.0                   DelayedMatrixStats_1.24.0     farver_2.1.2                  zoo_1.8-12                   
# [148] jsonlite_1.8.8                BiocParallel_1.36.0           R.oo_1.25.0                   magrittr_2.0.3                modeltools_0.2-23             GenomeInfoDbData_1.2.11       dotCall64_1.1-1              
# [155] munsell_0.5.1                 Rcpp_1.0.11                   ape_5.7-1                     reticulate_1.40.0             stringi_1.8.3                 zlibbioc_1.48.0               MASS_7.3-60                  
# [162] AnnotationHub_3.10.1          plyr_1.8.9                    pkgbuild_1.4.3                parallel_4.3.0                listenv_0.9.0                 deldir_2.0-2                  Biostrings_2.70.1            
# [169] splines_4.3.0                 tensor_1.5                    hms_1.1.3                     locfit_1.5-9.8                fastcluster_1.2.3             igraph_2.0.3                  spatstat.geom_3.2-7          
# [176] RcppHNSW_0.5.0                reshape2_1.4.4                biomaRt_2.58.0                stats4_4.3.0                  pkgload_1.3.3                 futile.options_1.0.1          BiocVersion_3.18.1           
# [183] XML_3.99-0.16                 evaluate_0.23                 RcppParallel_5.1.7            lambda.r_1.2.4                BiocManager_1.30.22           tzdb_0.4.0                    phyclust_0.1-34              
# [190] foreach_1.5.2                 httpuv_1.6.13                 RANN_2.6.1                    purrr_1.0.2                   polyclip_1.10-6               future_1.33.1                 scattermore_1.2              
# [197] coin_1.4-3                    xtable_1.8-4                  restfulr_0.0.15               RSpectra_0.16-1               later_1.3.2                   ragg_1.2.7                    rjags_4-15                   
# [204] viridisLite_0.4.2             tibble_3.2.1                  GenomicAlignments_1.38.1      memoise_2.0.1                 AnnotationDbi_1.64.1          IRanges_2.36.0                cluster_2.1.6                
# [211] globals_0.16.2               

