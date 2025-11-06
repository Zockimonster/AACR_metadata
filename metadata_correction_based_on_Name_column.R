
# Code used to correct the metadata, adjusting the metadata information to match the information 
# available for the corresponding sample identifier used (Name)

# https://aacrjournals.org/clincancerres/article/31/4/756/751743/Human-Pancreatic-Cancer-Single-Cell-Atlas-Reveals

# Load required libraries 
required_packages <- c('Matrix', 'dplyr','rentrez', 'SRAdb', 'GEOquery', 'stringr', 'tidyr')

invisible(lapply(required_packages, library, character.only = T))

# obtain directories
data_dir <- "path, where the following files can be found:
  - scAtlas.rds 
    # wget 'https://zenodo.org/records/14199536/files/scAtlas.rds.gz'
    # gunzip scAtlas.rds.gz
  - metadata_correction_functions.R (from https://github.com/Zockimonster/AACR_metadata repository)
  - optionally soft files, downloaded using wget to avoid repetitive downloading using getGEO()
    # file_name <- paste0('https://ftp.ncbi.nlm.nih.gov/geo/series/', to_geo_nnn(geo_nr), '/', geo_nr, '/soft/', geo_nr, '_family.soft.gz')
    # wget file_name
  - metadata excel file from https://ngdc.cncb.ac.cn/gsa/browse/CRA001160, called peng_CRR.xlsx"
source(paste0(data_dir, "metadata_correction_functions.R"))

####################################################################################################

# save pdac metadata separately
pdac_ref <- readRDS(paste0(data_dir, "scAtlas.rds"))
pdac_metadata <- pdac_ref@meta.data
write.csv(pdac_metadata, paste0(data_dir, "scAtlas_metadata.csv"))

####################################################################################################

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

# simplify column names
pdac_metadata <-  pdac_metadata %>% rename(
  "Study_Citation_PMID"="Study..Citation..PMID.",
  "Study_ID"="GSE.SRA..Study.",
  "raw_name"="Name" ,
  "pre_treatment"="TreatmentType",
  "tissue_type"="DiseaseState",
  "tissue_location"="If.metastatic..location")
# add columns
pdac_metadata$sample_id <- pdac_metadata$raw_name
pdac_metadata$sample_name <- pdac_metadata$raw_name
pdac_metadata$sc_type <- "scRNA"

# adjust tissue location
pdac_metadata$tissue_location[is.na(pdac_metadata$tissue_location)] <- "Pancreas"
table(pdac_metadata$tissue_location, pdac_metadata$tissue_type)
#              Adjacent normal  Donor Metastatic lesion Primary tumor
# Liver                      0      0             53145             0
# Lung                       0      0              1225             0
# Omentum                    0      0               138             0
# Pancreas               74035  33309                 0        562321
# Peritoneum                 0      0               459             0
# Vaginal apex               0      0              1475             0

#############################################################################################
#############################################################################################

# EGAS00001002543
# EGAD00010001811
# https://www.clinicaltrials.gov/study/NCT02750657
# inclusion criteria: no prior treatment
# primary or liver metastatic samples
# identification due to missing id links not possible
# numbers  #to_adj#
# https://www.nature.com/articles/s41588-019-0566-9
study_nr <- "EGAS00001002543"
# only Primary tumor?
table(pdac_metadata$tissue_location[pdac_metadata$Study_ID == study_nr], pdac_metadata$tissue_type[pdac_metadata$Study_ID == study_nr])
#         n Primary tumor
# Pancreas         76094

cbind(table(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr], pdac_metadata$Treatment[pdac_metadata$Study_ID == study_nr]),
      table(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr], pdac_metadata$pre_treatment[pdac_metadata$Study_ID == study_nr]))
#        Treatment Treatment naïve RT & chemotherapy
# 100070         0            6637                 0
# 85948          0            9404                 0
# 87235          0            5673                 0
# 87784          0            4320                 0
# 91412          0            7845                 0
# 91610          0            1590                 0
# 91706          0            6498                 0
# 94930          0            2514                 0
# 95092          0            7881                 0
# 95373          0            4412                 0
# 96460          0            1767                 0
# 97727      11536               0             11536
# G9903          0            6017                 0

# df to adjust data
df <- data.frame(matrix(nrow=length(unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr])), ncol=1))
colnames(df) <- "sample_name"
df$sample_name <- unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr])
df$tissue_location <- "Pancreas"
df$tissue_type <- "Primary tumor"
df$Treatment <- NA
df$Treatment[df$sample_name == "97727"] <- "Treatment"
df$Treatment[df$sample_name != "97727"] <- "Treatment naïve"

df$pre_treatment <- NA
df$pre_treatment[df$sample_name == "97727"] <- "RCT"
df$pre_treatment[df$sample_name != "97727"] <- "none"
df
#    sample_name tissue_location   tissue_type       Treatment pre_treatment
# 1       100070        Pancreas Primary tumor Treatment naïve          none
# 2        85948        Pancreas Primary tumor Treatment naïve          none
# 3        87235        Pancreas Primary tumor Treatment naïve          none
# 4        87784        Pancreas Primary tumor Treatment naïve          none
# 5        91412        Pancreas Primary tumor Treatment naïve          none
# 6        91610        Pancreas Primary tumor Treatment naïve          none
# 7        91706        Pancreas Primary tumor Treatment naïve          none
# 8        94930        Pancreas Primary tumor Treatment naïve          none
# 9        95092        Pancreas Primary tumor Treatment naïve          none
# 10       95373        Pancreas Primary tumor Treatment naïve          none
# 11       96460        Pancreas Primary tumor Treatment naïve          none
# 12       97727        Pancreas Primary tumor       Treatment           RCT
# 13       G9903        Pancreas Primary tumor Treatment naïve          none

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = study_nr,
  df = df,
  match_meta_col = "sample_name",
  match_df_col = "sample_name",
  clean_pattern = "")

# Updating Study: EGAS00001002543 
# ✅ Matched 76094 of 76094 samples
# Detected overlapping columns: sample_name, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# gse154778
# tissue info (Primary != Metastasis & Omentum Metastasis != Liver Metastasis)
# no treatment information in geo nor in run table or supplementary file (?)
# sample names
geo_nr <- "GSE154778"
table(pdac_metadata$Treatment[pdac_metadata$Study_ID == geo_nr])
# Treatment naïve 
# 20364    
table(pdac_metadata$pre_treatment[pdac_metadata$Study_ID == geo_nr], useNA="ifany")
# <NA> 
#   20364
table(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr], pdac_metadata$tissue_location[pdac_metadata$Study_ID == geo_nr])
#                           Liver Pancreas
# T1                            0      364  # t1 missing id
# T10                           0     2850
# T2                            0      872
# T3                            0      654
# T4                            0      842
# T5                            0     1261
# T6                            0     1772
# T8                            0     1953
# T9                            0     1238
# Y00006_HT77VBCXY_AACCGTAA    94        0  # Y00006 primary
# Y00013_HW7W7AFXX_TTGGCATA   913        0
# Y00014_HW7W7AFXX_TGAAGTAC    46        0  # Y00014 omentum
# Y00016_HW7W7AFXX_GTTGCAGC   209        0
# Y00019_H5HVVBGX5_GATCTCAG  3171        0
# Y00027_H5HVVBGX5_ACGACATT  4125        0

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=5))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("sample_name","file_name","SAMN", "tissue_location", "tissue_type")

# metastatic location
# https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00776-9
# Supplementary File Table S1

for (x in 1:length(names(gse@gsms))){
  disease_index <- grep("disease state: ", gse@gsms[[x]]@header$characteristics_ch1)
  disease_state <- gsub("disease state: ", "", gse@gsms[[x]]@header$characteristics_ch1[disease_index])
  disease_state <- ifelse(disease_state == "Pancreatic ductal adenocarcinoma", "PDAC", "different")
  sample_name <- gsub("\\:.*", "", gse@gsms[[x]]@header$title) 
  file_name <-gsub("^.*suppl/GSM.*_", "", gsub("_barcodes.*$", "", gse@gsms[[x]]@header$supplementary_file_1))
  
  tissue_index <- grep("tissue type: ", gse@gsms[[x]]@header$characteristics_ch1)
  tissue_type <- gsub("tissue type: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  if(tissue_type == "Primary" & disease_state == "PDAC"){
    tissue_location <- "Pancreas"
    tissue_type <- "Primary tumor"
  }else if (tissue_type == "Metastasis"& sample_name != "MET03"){
    tissue_location <- "Liver"
    tissue_type <- "Metastatic lesion"
  }else if (tissue_type == "Metastasis"& sample_name == "MET03"){
    tissue_location <- "Omentum"
    tissue_type <- "Metastatic lesion"}
  
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  df[x,]<- c(sample_name, file_name, samn_id, tissue_location, tissue_type)}

# add treatment
df$Treatment <- "?"
df$pre_treatment <- "?"
df <- df %>%
  dplyr::rename(
    sample_id = SAMN)
df
#               sample_name file_name    sample_id tissue_location       tissue_type Treatment pre_treatment
#   GSM4679532         P01    K16733 SAMN15585742        Pancreas     Primary tumor         ?             ?
#   GSM4679533         P02    Y00006 SAMN15585741        Pancreas     Primary tumor         ?             ?
#   GSM4679534         P03        T2 SAMN15585740        Pancreas     Primary tumor         ?             ?
#   GSM4679535         P04        T3 SAMN15585739        Pancreas     Primary tumor         ?             ?
#   GSM4679536         P05        T4 SAMN15585738        Pancreas     Primary tumor         ?             ?
#   GSM4679537         P06        T5 SAMN15585737        Pancreas     Primary tumor         ?             ?
#   GSM4679538         P07        T6 SAMN15585736        Pancreas     Primary tumor         ?             ?
#   GSM4679539         P08        T8 SAMN15585735        Pancreas     Primary tumor         ?             ?
#   GSM4679540         P09        T9 SAMN15585734        Pancreas     Primary tumor         ?             ?
#   GSM4679541         P10       T10 SAMN15585733        Pancreas     Primary tumor         ?             ?
#   GSM4679542       MET01    Y00008 SAMN15585732           Liver Metastatic lesion         ?             ?
#   GSM4679543       MET02    Y00013 SAMN15585731           Liver Metastatic lesion         ?             ?
#   GSM4679544       MET03    Y00014 SAMN15585730         Omentum Metastatic lesion         ?             ?
#   GSM4679545       MET04    Y00016 SAMN15585729           Liver Metastatic lesion         ?             ?
#   GSM4679546       MET05    Y00019 SAMN15585728           Liver Metastatic lesion         ?             ?
#   GSM4679547       MET06    Y00027 SAMN15585727           Liver Metastatic lesion         ?             ?

# in original metadata
#  Y00006 primary
unique(pdac_metadata$tissue_location[pdac_metadata$sample_name =="Y00006_HT77VBCXY_AACCGTAA"])
# [1] "Liver"

# Y00014 omentum
unique(pdac_metadata$tissue_location[pdac_metadata$sample_name =="Y00014_HW7W7AFXX_TGAAGTAC"])
# [1] "Liver"

# all treatment naive?
table(pdac_metadata$Treatment[pdac_metadata$Study_ID == geo_nr])
# Treatment naïve 
# 20364 

# also no treatment info in run info (SRA Run Info)
run_info <- get_sra_run_info(geo_nr, data_dir = data_dir, use_soft_file = TRUE)

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_name",
  match_df_col = "file_name",
  clean_pattern = "_.*$")

# Updating Study: GSE154778 
# ✅ Matched 20000 of 20364 samples
# ⚠️ 364 samples had no match.
# Detected overlapping columns: sample_name, sample_id, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_name 
# Updated column: sample_id 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment

##############################################################################################
# gse155698
# srr don't match to geo data, info from rentrez -> one samn id & one sample_name matches multiple srr run ids
geo_nr <- "GSE155698"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=6))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("SAMN", "sample_name", "tissue_type", "tissue_location", "Treatment", "pre_treatment")

for (x in names(gse@gsms)){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  sample_name <- gse@gsms[[x]]@header$title
  
  if (grepl("PDAC_TISSUE", sample_name)){
    tissue_type <- "Primary tumor"
    tissue_location <- "Pancreas"
  }else if (grepl("AdjNorm_TISSUE", sample_name)){
    tissue_type <- "Adjacent normal"
    tissue_location <- "Pancreas"}
  
  if (grepl("PDAC_PBMC", sample_name)){
    tissue_type <- "PBMC from PDAC patients"
    tissue_location <- "PBMC"
  }else if (grepl("Healthy_PBMC", sample_name)){
    tissue_type <- "PBMC from healthy patients"
    tissue_location <- "PBMC"}
  treatment_index <- grepl("treatment: ", gse@gsms[[x]]@header$characteristics_ch1)
  treatment_info <-ifelse(grepl("Treatment Naïve", gse@gsms[[x]]@header$characteristics_ch1[treatment_index]), "Treatment naïve", "Healthy control")
  if (treatment_info %in% c("Treatment naïve", "Healthy control")){
    pre_treatment <- "none"}
  
  df[x,] <- c(samn_id, sample_name, tissue_type, tissue_location, treatment_info, pre_treatment)}

df
#                    SAMN      sample_name                tissue_type tissue_location       Treatment pre_treatment
# GSM4710689 SAMN15731050    PDAC_TISSUE_1              Primary tumor        Pancreas Treatment naïve          none
# GSM4710690 SAMN15731047    PDAC_TISSUE_2              Primary tumor        Pancreas Treatment naïve          none
# GSM4710691 SAMN15731045    PDAC_TISSUE_3              Primary tumor        Pancreas Treatment naïve          none
# GSM4710692 SAMN15731044    PDAC_TISSUE_4              Primary tumor        Pancreas Treatment naïve          none
# GSM4710693 SAMN15731043    PDAC_TISSUE_5              Primary tumor        Pancreas Treatment naïve          none
# GSM4710694 SAMN15731049    PDAC_TISSUE_6              Primary tumor        Pancreas Treatment naïve          none
# GSM4710695 SAMN15731041    PDAC_TISSUE_7              Primary tumor        Pancreas Treatment naïve          none
# GSM4710696 SAMN15731040    PDAC_TISSUE_8              Primary tumor        Pancreas Treatment naïve          none
# GSM4710697 SAMN15731038    PDAC_TISSUE_9              Primary tumor        Pancreas Treatment naïve          none
# GSM4710698 SAMN15731037   PDAC_TISSUE_10              Primary tumor        Pancreas Treatment naïve          none
# GSM4710699 SAMN15731036  PDAC_TISSUE_11A              Primary tumor        Pancreas Treatment naïve          none
# GSM4710700 SAMN15731035  PDAC_TISSUE_11B              Primary tumor        Pancreas Treatment naïve          none
# GSM4710701 SAMN15731031   PDAC_TISSUE_12              Primary tumor        Pancreas Treatment naïve          none
# GSM4710702 SAMN15731030   PDAC_TISSUE_13              Primary tumor        Pancreas Treatment naïve          none
# GSM4710703 SAMN15731029   PDAC_TISSUE_14              Primary tumor        Pancreas Treatment naïve          none
# GSM4710704 SAMN15731034   PDAC_TISSUE_15              Primary tumor        Pancreas Treatment naïve          none
# GSM4710705 SAMN15731033   PDAC_TISSUE_16              Primary tumor        Pancreas Treatment naïve          none
# GSM4710706 SAMN15731032 AdjNorm_TISSUE_1            Adjacent normal        Pancreas Treatment naïve          none
# GSM4710707 SAMN15731028 AdjNorm_TISSUE_2            Adjacent normal        Pancreas Treatment naïve          none
# GSM4710708 SAMN15731027 AdjNorm_TISSUE_3            Adjacent normal        Pancreas Treatment naïve          none
# GSM4710709 SAMN15731026      PDAC_PBMC_1    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710710 SAMN15731025      PDAC_PBMC_2    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710711 SAMN15731024      PDAC_PBMC_3    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710712 SAMN15731023      PDAC_PBMC_4    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710713 SAMN15731022      PDAC_PBMC_5    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710714 SAMN15731021      PDAC_PBMC_6    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710715 SAMN15731065      PDAC_PBMC_7    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710716 SAMN15731064      PDAC_PBMC_8    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710717 SAMN15731063      PDAC_PBMC_9    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710718 SAMN15731062    PDAC_PBMC_10A    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710719 SAMN15731061    PDAC_PBMC_10B    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710720 SAMN15731060     PDAC_PBMC_11    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710721 SAMN15731059     PDAC_PBMC_12    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710722 SAMN15731058     PDAC_PBMC_13    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710723 SAMN15731057     PDAC_PBMC_14    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710724 SAMN15731056     PDAC_PBMC_15    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710725 SAMN15731055     PDAC_PBMC_16    PBMC from PDAC patients            PBMC Treatment naïve          none
# GSM4710726 SAMN15731054   Healthy_PBMC_1 PBMC from healthy patients            PBMC Healthy control          none
# GSM4710727 SAMN15731053   Healthy_PBMC_2 PBMC from healthy patients            PBMC Healthy control          none
# GSM4710728 SAMN15731052   Healthy_PBMC_3 PBMC from healthy patients            PBMC Healthy control          none
# GSM4710729 SAMN15731051   Healthy_PBMC_4 PBMC from healthy patients            PBMC Healthy control          none

# use rentrez to get info for srr used in pdac_metadata
pdac_srr_id <- unique(pdac_metadata$sample_id[pdac_metadata$Study_ID == geo_nr])
info_155698 <- sapply(pdac_srr_id, get_samn_and_info_geo155698)
write.csv(info_155698, paste0(data_dir, "rentrez_info_GSE155698.csv"))
info_155698 <- read.csv(paste0(data_dir, "rentrez_info_GSE155698.csv"), row.names = 1)
info_155698 <- as.data.frame(t(info_155698))
colnames(info_155698) <- c("tissue_type", "SAMN", "sample_name")
# add tissue from reference
info_155698$in_ref <- pdac_metadata$tissue_type[match(rownames(info_155698), pdac_metadata$sample_name)]
info_155698
#                 tissue_type         SAMN        sample_name          in_ref
# SRR12461530   Primary tumor SAMN15750836   PDAC_1246_tissue   Primary tumor
# SRR12461533   Primary tumor SAMN15750849   PDAC_1229_tissue   Primary tumor
# SRR12461536   Primary tumor SAMN15750842   PDAC_1261_tissue   Primary tumor
# SRR12461538   Primary tumor SAMN15750838   PDAC_1262_tissue   Primary tumor
# SRR12461543 Adjacent normal SAMN15750850   Norm_1258_tissue Adjacent normal
# SRR12461544 Adjacent normal SAMN15750845 Norm_19_732_tissue Adjacent normal
# SRR12461548   Primary tumor SAMN15750826   PDAC_1253_tissue   Primary tumor
# SRR12507499 Adjacent normal SAMN15750829   Norm_1196_tissue Adjacent normal
# SRR12507500 Adjacent normal SAMN15750829   Norm_1196_tissue Adjacent normal
# SRR12507930   Primary tumor SAMN15750841   PDAC_1141_tissue   Primary tumor
# SRR12507931   Primary tumor SAMN15750841   PDAC_1141_tissue   Primary tumor
# SRR12507934   Primary tumor SAMN15750815   PDAC_3210_tissue   Primary tumor
# SRR12507935   Primary tumor SAMN15750815   PDAC_3210_tissue   Primary tumor
# SRR12507938   Primary tumor SAMN15750840   PDAC_1238_tissue   Primary tumor
# SRR12507939   Primary tumor SAMN15750840   PDAC_1238_tissue   Primary tumor
# SRR12508148   Primary tumor SAMN15750842   PDAC_1261_tissue   Primary tumor
# SRR12532883   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor
# SRR12532884   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor
# SRR12532885   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor
# SRR12532886   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor

# no SAMN id from geo in samn ids from rentrez?
any(df$SAMN %in% info_155698$SAMN)
# [1] FALSE

# only 12 samples with 20 run ids in metadata
length(unique(info_155698$SAMN))
# [1] 12
length(unique(rownames(info_155698)))
# [1] 20

# adjust info_155698
info_155698$SRR <- rownames(info_155698)
info_155698 <- info_155698 %>%
  mutate(tissue_location = case_when(
    grepl("_tissue$", sample_name) ~"Pancreas",
    TRUE ~ "?"))
info_155698$Treatment <- "Treatment naïve"
info_155698$pre_treatment <- "none"
info_155698 <- info_155698 %>%
  dplyr::rename(
    sample_id = SAMN)
info_155698
#                 tissue_type    sample_id        sample_name          in_ref         SRR tissue_location       Treatment pre_treatment
# SRR12461530   Primary tumor SAMN15750836   PDAC_1246_tissue   Primary tumor SRR12461530        Pancreas Treatment naïve          none
# SRR12461533   Primary tumor SAMN15750849   PDAC_1229_tissue   Primary tumor SRR12461533        Pancreas Treatment naïve          none
# SRR12461536   Primary tumor SAMN15750842   PDAC_1261_tissue   Primary tumor SRR12461536        Pancreas Treatment naïve          none
# SRR12461538   Primary tumor SAMN15750838   PDAC_1262_tissue   Primary tumor SRR12461538        Pancreas Treatment naïve          none
# SRR12461543 Adjacent normal SAMN15750850   Norm_1258_tissue Adjacent normal SRR12461543        Pancreas Treatment naïve          none
# SRR12461544 Adjacent normal SAMN15750845 Norm_19_732_tissue Adjacent normal SRR12461544        Pancreas Treatment naïve          none
# SRR12461548   Primary tumor SAMN15750826   PDAC_1253_tissue   Primary tumor SRR12461548        Pancreas Treatment naïve          none
# SRR12507499 Adjacent normal SAMN15750829   Norm_1196_tissue Adjacent normal SRR12507499        Pancreas Treatment naïve          none
# SRR12507500 Adjacent normal SAMN15750829   Norm_1196_tissue Adjacent normal SRR12507500        Pancreas Treatment naïve          none
# SRR12507930   Primary tumor SAMN15750841   PDAC_1141_tissue   Primary tumor SRR12507930        Pancreas Treatment naïve          none
# SRR12507931   Primary tumor SAMN15750841   PDAC_1141_tissue   Primary tumor SRR12507931        Pancreas Treatment naïve          none
# SRR12507934   Primary tumor SAMN15750815   PDAC_3210_tissue   Primary tumor SRR12507934        Pancreas Treatment naïve          none
# SRR12507935   Primary tumor SAMN15750815   PDAC_3210_tissue   Primary tumor SRR12507935        Pancreas Treatment naïve          none
# SRR12507938   Primary tumor SAMN15750840   PDAC_1238_tissue   Primary tumor SRR12507938        Pancreas Treatment naïve          none
# SRR12507939   Primary tumor SAMN15750840   PDAC_1238_tissue   Primary tumor SRR12507939        Pancreas Treatment naïve          none
# SRR12508148   Primary tumor SAMN15750842   PDAC_1261_tissue   Primary tumor SRR12508148        Pancreas Treatment naïve          none
# SRR12532883   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor SRR12532883        Pancreas Treatment naïve          none
# SRR12532884   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor SRR12532884        Pancreas Treatment naïve          none
# SRR12532885   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor SRR12532885        Pancreas Treatment naïve          none
# SRR12532886   Primary tumor SAMN15750832   PDAC_1314_tissue   Primary tumor SRR12532886        Pancreas Treatment naïve          none

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = info_155698,
  match_meta_col = "sample_id",
  match_df_col = "SRR",
  clean_pattern = "")

# Updating Study: GSE155698 
# ✅ Matched 23991 of 23991 samples
# Detected overlapping columns: tissue_type, sample_id, sample_name, tissue_location, Treatment, pre_treatment 
# Updated column: tissue_type 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: Treatment 
# Updated column: pre_treatment 

# 3 healthy, 9 pdac samples
table(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr])
# Norm_1196_tissue   Norm_1258_tissue Norm_19_732_tissue   PDAC_1141_tissue   PDAC_1229_tissue   PDAC_1238_tissue   PDAC_1246_tissue   PDAC_1253_tissue   PDAC_1261_tissue   PDAC_1262_tissue   PDAC_1314_tissue   PDAC_3210_tissue 
#             1035               1155               2312               1599               2070               4597                458                746               3420                850               4492               1257 

##############################################################################################
# gse156405
# tissue (Primary != Peritoneal Metastasis & Peritoneal Metastasis != Primary)
# treatment (a couple of patients were pre-treated) 
# one old srr id
geo_nr <- "GSE156405"
table(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr], pdac_metadata$tissue_location[pdac_metadata$Study_ID == geo_nr])
#             Liver Lung Pancreas Peritoneum Vaginal apex
# SRR12467642     0    0     1304          0            0
# SRR12467643     0    0     2207          0            0  # SRR12467643, old: "SRR12467643" -> deprecated: replaced by SRR12467650, https://www.ncbi.nlm.nih.gov/sra/?term=SRR12467643  & peritoneal metastasis
# SRR12467644     0    0     3975          0            0
# SRR12467645     0    0     1961          0            0
# SRR12467646     0    0      564          0            0
# SRR12467647     0    0        0          0         1475
# SRR12467648  4639    0        0          0            0
# SRR12467649     0 1225        0          0            0
# SRR24437450     0    0        0        459            0  # SRR24437450: primary

# replace SRR12467643 by current name
pdac_metadata$sample_id[pdac_metadata$sample_id == "SRR12467643"] <- "SRR12467650"
pdac_metadata$sample_name[pdac_metadata$sample_name == "SRR12467643"] <- "SRR12467650"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

# Sample Acquisition
# Five primary pancreatic cancers and four metastatic lesions in liver (LiM), lung (LuM), peritoneum (PM), and vaginal apex (VM) were used. 
# P2-P5 primary samples came from treatment-naïve patients, while P1 was a second FNA from a patient who had been treated with gemcitabine/paclitaxel. 
# For metastatic samples, the most recent therapies prior to sample acquisition were 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab and capecitabine for VM, LuM, and PM, respectively. 
# LiM was a treatment-naïve sample. Histologic confirmation of PDAC was performed by a pathologist.

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=8))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("SRR","SAMN","sample_name","sample_description", "tissue_location", "tissue_type", "Treatment", "pre_treatment")

# metastatic location
# https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00776-9
# Supplementary File Table S1

for (x in 1:length(names(gse@gsms))){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  srr_index <- grep("SRA:", gse@gsms[[x]]@header$relation)
  srr_id <-  gsub("^.*term\\=", "", gse@gsms[[x]]@header$relation[srr_index])
  srr_id <- get_srr_id(srr_id)
  
  sample_description <-  gse@gsms[[x]]@header$title
  sample_name <-  gse@gsms[[x]]@header$description
  
  tissue_index <- grep("tissue: ", gse@gsms[[x]]@header$characteristics_ch1)
  tissue_location <- gsub("tissue: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  
  if (tissue_location == "Pancreas"){
    tissue_type <- "Primary tumor"
  } else{
    tissue_type <- "Metastatic lesion"}
  
  if (sample_name == "P1"){
    pre_treatment <- "gemcitabine/paclitaxel"
  }else if (sample_name %in% c("VM", "LuM", "PM")){
    pre_treatment <- "5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab and capecitabine"}
  if (sample_name %in% c("P1", "VM", "LuM", "PM")){
    treatment_info <- "Treatment"
  } else{
    treatment_info <- "Treatment naïve"
    pre_treatment <- "none"}
  
  df[x,]<- c(srr_id, samn_id, sample_name, sample_description, tissue_location, tissue_type, treatment_info, pre_treatment)}

# adjust df
df$pre_treatment <- gsub("5-fluorouracil/liposmal irinotecan", "FOLFIRI", df$pre_treatment)
df <- df %>%
  dplyr::rename(
    sample_id = SAMN)

df[df$SRR %in% unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr]),]
#                    SRR    sample_id sample_name          sample_description tissue_location       tissue_type       Treatment                                     pre_treatment
# GSM4730260 SRR12467642 SAMN15848242          P1 Fine needle aspiration - P1        Pancreas     Primary tumor       Treatment                            gemcitabine/paclitaxel
# GSM4730261 SRR24437450 SAMN15848241          P2 Fine needle aspiration - P2        Pancreas     Primary tumor Treatment naïve                                              none
# GSM4730262 SRR12467644 SAMN15848240          P3 Fine needle aspiration - P3        Pancreas     Primary tumor Treatment naïve                                              none
# GSM4730263 SRR12467645 SAMN15848239          P4 Fine needle aspiration - P4        Pancreas     Primary tumor Treatment naïve                                              none
# GSM4730264 SRR12467646 SAMN15848238          P5 Fine needle aspiration - P5        Pancreas     Primary tumor Treatment naïve                                              none
# GSM4730265 SRR12467647 SAMN15848237          VM          Vaginal metastasis    Vaginal apex Metastatic lesion       Treatment FOLFIRI, evofosfamide/ipilimumab and capecitabine
# GSM4730266 SRR12467648 SAMN15848236         LiM            Liver metastasis           Liver Metastatic lesion Treatment naïve                                              none
# GSM4730267 SRR12467649 SAMN15848235         LuM             Lung metastasis            Lung Metastatic lesion       Treatment FOLFIRI, evofosfamide/ipilimumab and capecitabine
# GSM4730268 SRR12467650 SAMN15848234          PM       Peritoneal metastasis      Peritoneum Metastatic lesion       Treatment FOLFIRI, evofosfamide/ipilimumab and capecitabine

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "SRR",
  clean_pattern = "_.*$")

# Updating Study: GSE156405 
# ✅ Matched 17809 of 17809 samples
# Detected overlapping columns: sample_id, sample_name, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# gse158356
# tissue info (LivMetD != Omentum Metastasis)
# no treatment info available?
# sample names  #to_adj#
geo_nr <- "GSE158356"
table(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr], pdac_metadata$tissue_location[pdac_metadata$Study_ID == geo_nr])
#           Liver Omentum
# LiverMetA  1252       0
# LiverMetB   455       0
# LiverMetC   184       0
# LiverMetD     0     138
# LiverMetE   717       0

table(pdac_metadata$Treatment[pdac_metadata$Study_ID == geo_nr])
# Treatment Treatment naïve 
# 1969             777  
table(pdac_metadata$pre_treatment[pdac_metadata$Study_ID == geo_nr], useNA="ifany")
# <NA> 
#   2746 

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

# Exclude mouse samples
for (x in names(gse@gsms)){
  if (gse@gsms[[x]]@header$organism_ch1 == "Homo sapiens"){
    gse@gsms[[x]] <- gse@gsms[[x]]
  }else{
    gse@gsms[[x]] <- NULL}}

length(gse@gsms)
# [1] 5

names(gse@gsms)
# [1] "GSM4798244" "GSM4798245" "GSM4798246" "GSM4798247" "GSM4798248"


# obtain geo sample characteristics
df <- data.frame(matrix(nrow = length(gse@gsms), ncol=4))
rownames(df) <- names(gse@gsms)
colnames(df) <- c("SAMN", "sample_name", "tissue_location", "tissue_type")
for (x in 1:length(names(gse@gsms))){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  tissue_index <- grep("tissue: ", gse@gsms[[x]]@header$characteristics_ch1)
  if (gse@gsms[[x]]@header$characteristics_ch1[tissue_index] == "tissue: human liver metastasis sample"){
    tissue_type <- "Metastatic lesion"
    tissue_location <- "Liver"}
  
  sample_name <- gse@gsms[[x]]@header$title
  
  df[x,]<- c(samn_id, sample_name, tissue_location, tissue_type)}

df$Treatment <- "?"
df$pre_treatment <- "?"
df
#                      SAMN sample_name tissue_location       tissue_type Treatment pre_treatment
#   GSM4798244 SAMN16237818 Liver_Met_A           Liver Metastatic lesion         ?             ?
#   GSM4798245 SAMN16237817 Liver_Met_B           Liver Metastatic lesion         ?             ?
#   GSM4798246 SAMN16237816 Liver_Met_C           Liver Metastatic lesion         ?             ?
#   GSM4798247 SAMN16237815 Liver_Met_D           Liver Metastatic lesion         ?             ?
#   GSM4798248 SAMN16237814 Liver_Met_E           Liver Metastatic lesion         ?             ?

df <- df %>%
  dplyr::rename(
    sample_id = SAMN)

# also no treatment info in run info (SRA Run Info)
run_info <- get_sra_run_info(geo_nr, data_dir=data_dir, use_soft_file=TRUE)
# select human samples

run_info[sapply("Homo sapiens",grepl,  run_info)]
# $GSM4798244
# [1] "  <Summary><Title>GSM4798244: Liver_Met_A; Homo sapiens; RNA-Seq</Title><Platform instrument_model=\"Illumina NovaSeq 6000\">ILLUMINA</Platform><Statistics total_runs=\"8\" total_spots=\"839088205\" total_bases=\"119989613315\" total_size=\"36107829829\" load_done=\"true\" cluster_name=\"public\"/></Summary><Submitter acc=\"SRA1129781\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/><Experiment acc=\"SRX9172133\" ver=\"1\" status=\"public\" name=\"GSM4798244: Liver_Met_A; Homo sapiens; RNA-Seq\"/><Study acc=\"SRP284880\" name=\"Pancreatic cancer is marked by complement-high blood monocytes and tumor-associated macrophages\"/><Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/><Sample acc=\"SRS7409916\" name=\"\"/><Instrument ILLUMINA=\"Illumina NovaSeq 6000\"/><Library_descriptor><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT><LIBRARY_CONSTRUCTION_PROTOCOL>For both human and mouse, samples were mechanically and enzymatically (Collagenase P; 1mg/mL in DMEM) digested under constant agitation. Single cell suspension was obtained by filtering through a 40 µm mesh. Dead cells were removed using MACS®Dead Cell Removal Kit (Miltenyi Biotec Inc.). Single-cell cDNA libraries were  prepared  and sequenced using  the  10x Genomics  Platform at the University of Michigan Advanced Genomics Core</LIBRARY_CONSTRUCTION_PROTOCOL></Library_descriptor><Bioproject>PRJNA664962</Bioproject><Biosample>SAMN16237818</Biosample>  "
# 
# $GSM4798245
# [1] "  <Summary><Title>GSM4798245: Liver_Met_B; Homo sapiens; RNA-Seq</Title><Platform instrument_model=\"Illumina NovaSeq 6000\">ILLUMINA</Platform><Statistics total_runs=\"4\" total_spots=\"1214818392\" total_bases=\"366875154384\" total_size=\"132777398027\" load_done=\"true\" cluster_name=\"public\"/></Summary><Submitter acc=\"SRA1129781\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/><Experiment acc=\"SRX9172134\" ver=\"1\" status=\"public\" name=\"GSM4798245: Liver_Met_B; Homo sapiens; RNA-Seq\"/><Study acc=\"SRP284880\" name=\"Pancreatic cancer is marked by complement-high blood monocytes and tumor-associated macrophages\"/><Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/><Sample acc=\"SRS7409917\" name=\"\"/><Instrument ILLUMINA=\"Illumina NovaSeq 6000\"/><Library_descriptor><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT><LIBRARY_CONSTRUCTION_PROTOCOL>For both human and mouse, samples were mechanically and enzymatically (Collagenase P; 1mg/mL in DMEM) digested under constant agitation. Single cell suspension was obtained by filtering through a 40 µm mesh. Dead cells were removed using MACS®Dead Cell Removal Kit (Miltenyi Biotec Inc.). Single-cell cDNA libraries were  prepared  and sequenced using  the  10x Genomics  Platform at the University of Michigan Advanced Genomics Core</LIBRARY_CONSTRUCTION_PROTOCOL></Library_descriptor><Bioproject>PRJNA664962</Bioproject><Biosample>SAMN16237817</Biosample>  "
# 
# $GSM4798246
# [1] "  <Summary><Title>GSM4798246: Liver_Met_C; Homo sapiens; RNA-Seq</Title><Platform instrument_model=\"Illumina NovaSeq 6000\">ILLUMINA</Platform><Statistics total_runs=\"4\" total_spots=\"931450989\" total_bases=\"281298198678\" total_size=\"100968942588\" load_done=\"true\" cluster_name=\"public\"/></Summary><Submitter acc=\"SRA1129781\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/><Experiment acc=\"SRX9172135\" ver=\"1\" status=\"public\" name=\"GSM4798246: Liver_Met_C; Homo sapiens; RNA-Seq\"/><Study acc=\"SRP284880\" name=\"Pancreatic cancer is marked by complement-high blood monocytes and tumor-associated macrophages\"/><Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/><Sample acc=\"SRS7409918\" name=\"\"/><Instrument ILLUMINA=\"Illumina NovaSeq 6000\"/><Library_descriptor><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT><LIBRARY_CONSTRUCTION_PROTOCOL>For both human and mouse, samples were mechanically and enzymatically (Collagenase P; 1mg/mL in DMEM) digested under constant agitation. Single cell suspension was obtained by filtering through a 40 µm mesh. Dead cells were removed using MACS®Dead Cell Removal Kit (Miltenyi Biotec Inc.). Single-cell cDNA libraries were  prepared  and sequenced using  the  10x Genomics  Platform at the University of Michigan Advanced Genomics Core</LIBRARY_CONSTRUCTION_PROTOCOL></Library_descriptor><Bioproject>PRJNA664962</Bioproject><Biosample>SAMN16237816</Biosample>  "
# 
# $GSM4798247
# [1] "  <Summary><Title>GSM4798247: Liver_Met_D; Homo sapiens; RNA-Seq</Title><Platform instrument_model=\"Illumina NovaSeq 6000\">ILLUMINA</Platform><Statistics total_runs=\"4\" total_spots=\"1106827162\" total_bases=\"334261802924\" total_size=\"121239256251\" load_done=\"true\" cluster_name=\"public\"/></Summary><Submitter acc=\"SRA1129781\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/><Experiment acc=\"SRX9172136\" ver=\"1\" status=\"public\" name=\"GSM4798247: Liver_Met_D; Homo sapiens; RNA-Seq\"/><Study acc=\"SRP284880\" name=\"Pancreatic cancer is marked by complement-high blood monocytes and tumor-associated macrophages\"/><Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/><Sample acc=\"SRS7409919\" name=\"\"/><Instrument ILLUMINA=\"Illumina NovaSeq 6000\"/><Library_descriptor><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT><LIBRARY_CONSTRUCTION_PROTOCOL>For both human and mouse, samples were mechanically and enzymatically (Collagenase P; 1mg/mL in DMEM) digested under constant agitation. Single cell suspension was obtained by filtering through a 40 µm mesh. Dead cells were removed using MACS®Dead Cell Removal Kit (Miltenyi Biotec Inc.). Single-cell cDNA libraries were  prepared  and sequenced using  the  10x Genomics  Platform at the University of Michigan Advanced Genomics Core</LIBRARY_CONSTRUCTION_PROTOCOL></Library_descriptor><Bioproject>PRJNA664962</Bioproject><Biosample>SAMN16237815</Biosample>  "
# 
# $GSM4798248
# [1] "  <Summary><Title>GSM4798248: Liver_Met_E; Homo sapiens; RNA-Seq</Title><Platform instrument_model=\"Illumina NovaSeq 6000\">ILLUMINA</Platform><Statistics total_runs=\"4\" total_spots=\"1039051907\" total_bases=\"143389163166\" total_size=\"46051995476\" load_done=\"true\" cluster_name=\"public\"/></Summary><Submitter acc=\"SRA1129781\" center_name=\"GEO\" contact_name=\"Gene Expression Omnibus (GEO), NCBI, NLM, NIH, htt\" lab_name=\"\"/><Experiment acc=\"SRX9172137\" ver=\"1\" status=\"public\" name=\"GSM4798248: Liver_Met_E; Homo sapiens; RNA-Seq\"/><Study acc=\"SRP284880\" name=\"Pancreatic cancer is marked by complement-high blood monocytes and tumor-associated macrophages\"/><Organism taxid=\"9606\" ScientificName=\"Homo sapiens\"/><Sample acc=\"SRS7409920\" name=\"\"/><Instrument ILLUMINA=\"Illumina NovaSeq 6000\"/><Library_descriptor><LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY><LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE><LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION><LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT><LIBRARY_CONSTRUCTION_PROTOCOL>For both human and mouse, samples were mechanically and enzymatically (Collagenase P; 1mg/mL in DMEM) digested under constant agitation. Single cell suspension was obtained by filtering through a 40 µm mesh. Dead cells were removed using MACS®Dead Cell Removal Kit (Miltenyi Biotec Inc.). Single-cell cDNA libraries were  prepared  and sequenced using  the  10x Genomics  Platform at the University of Michigan Advanced Genomics Core</LIBRARY_CONSTRUCTION_PROTOCOL></Library_descriptor><Bioproject>PRJNA664962</Bioproject><Biosample>SAMN16237814</Biosample>  "

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "sample_name",
  clean_pattern = "_")

# Updating Study: GSE158356 
# ✅ Matched 2746 of 2746 samples
# Detected overlapping columns: sample_id, sample_name, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# gse194247
# tissue info (all samples included cd45neg cells only)
# treatment correct
# amplification round
geo_nr <- "GSE194247"
table(pdac_metadata$Treatment[pdac_metadata$Study_ID == geo_nr])
# Treatment naïve 
# 29622   
table(pdac_metadata$pre_treatment[pdac_metadata$Study_ID == geo_nr], useNA="ifany")
# <NA> 
#   29622 

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

# Obtain GEO-Sample Characteristics
df <- data.frame(matrix(nrow = length(gse@gsms), ncol=6))
rownames(df) <- names(gse@gsms)
colnames(df) <- c("SAMN", "amplification_round", "tissue_location", "tissue_type", "Treatment", "pre_treatment")
for (x in 1:length(names(gse@gsms))){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  tissue_index <- grep("tissue: ", gse@gsms[[x]]@header$characteristics_ch1)
  tissue_location <- to_capital(gsub("tissue: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index]))
  
  cell_index <- grep("cell type: ", gse@gsms[[x]]@header$characteristics_ch1)
  cell_type <- gsub("cell type: ","", gse@gsms[[x]]@header$characteristics_ch1[cell_index])
  type_index <- grep("tissue type: ", gse@gsms[[x]]@header$characteristics_ch1)
  if (gse@gsms[[x]]@header$characteristics_ch1[type_index]== "tissue type: pancreatic cancer tumor tissue" & 
      gse@gsms[[x]]@header$characteristics_ch1[cell_index] == "cell type: CD45(-) cells"){
    tissue_type <- "Primary tumor, CD45neg cells"}
  
  amplification_round <- gsub("\\]$","",gsub("^5_","", gsub("^.*\\[", "", gse@gsms[[x]]@header$title)))
  
  treatment_info <- gse@gsms[[x]]@header$treatment_protocol_ch1
  
  if (treatment_info == "Untreated"){
    treatment_info <- "Treatment naïve"
    pre_treatment <- "none"}
  
  df[x,]<- c(samn_id, amplification_round, tissue_location, tissue_type, treatment_info, pre_treatment)}

df
#                    SAMN amplification_round tissue_location                  tissue_type       Treatment pre_treatment
# GSM5831620 SAMN25221948               GEX_4        Pancreas Primary tumor, CD45neg cells Treatment naïve          none
# GSM5831621 SAMN25221947               GEX_5        Pancreas Primary tumor, CD45neg cells Treatment naïve          none
# GSM5831622 SAMN25221946               GEX_6        Pancreas Primary tumor, CD45neg cells Treatment naïve          none
# GSM5831623 SAMN25221945               GEX_9        Pancreas Primary tumor, CD45neg cells Treatment naïve          none
# GSM5831624 SAMN25221944           GEX_45_MM        Pancreas Primary tumor, CD45neg cells Treatment naïve          none

# adjust df
df <- df %>%
  dplyr::rename(
    sample_id = SAMN,
    sample_name = amplification_round)

unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr])
# [1] "GEX_4"               "GEX_45_MM--ATTTGCTA" "GEX_5"               "GEX_6"               "GEX_9"   

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "sample_name",
  clean_pattern = "\\-.*$")

# Updating Study: GSE194247 
# ✅ Matched 29622 of 29622 samples
# Detected overlapping columns: sample_id, sample_name, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# gse202051
# tissue info (organoid cultures != Primary)
# treatment correct
# srr ids
geo_nr <- "GSE202051"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=8))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("SRR","SAMN","sample_name", "tissue_location", "tissue_type", "Treatment", "pre_treatment", "sc_type")

for (x in 1:length(names(gse@gsms))){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  srr_index <- grep("SRA:", gse@gsms[[x]]@header$relation)
  srr_id <-  gsub("^.*term\\=", "", gse@gsms[[x]]@header$relation[srr_index])
  srr_id <- get_srr_id(srr_id)
  
  tissue_index <- grep("tissue:", gse@gsms[[x]]@header$characteristics_ch1)  
  tissue_location <-  gsub("tissue: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  
  tissue_index <- grep("cell type:", gse@gsms[[x]]@header$characteristics_ch1)  
  tissue_type <-  gsub("cell type: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  
  treatment_index <- grep("treatment:", gse@gsms[[x]]@header$characteristics_ch1)  
  pre_treatment <-  gsub("treatment: ", "", gse@gsms[[x]]@header$characteristics_ch1[treatment_index])
  
  if (pre_treatment == "Untreated"){
    pre_treatment <- "none"
    treatment_info <- "Treatment naïve"
  }else{
    treatment_info <- "Treatment"}
  
  sample_description <-  gse@gsms[[x]]@header$title
  
  # adjust sample name
  sample_name <- gsub(",","", gsub(", .*$ | replicate.*$", "", sample_description))
  
  for (i in 1:3){
    look_for <- paste0("replicate ", i)
    if (grepl(look_for, sample_description)){
      sample_name <- paste0(sample_name,"_R", i)}}
  
  if (grepl("snRNA", sample_description)){
    sc_type <- "snRNA"}
  
  df[x,]<- c(srr_id, samn_id, sample_name, tissue_location, tissue_type, treatment_info, pre_treatment, sc_type)}

# adjust df
df <- df %>%
  mutate(tissue_type = case_when(
    grepl("PDAC", sample_name) ~ "Primary tumor",
    grepl("Normal", sample_name) ~ "Adjacent normal",
    grepl("^010", sample_name) ~ "Organoid cultures from Primary tumor samples"))

df <- df %>%
  dplyr::rename(
    sample_id = SAMN)

df[df$SRR %in% unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr]),]
#                    SRR    sample_id      sample_name tissue_location                                  tissue_type       Treatment            pre_treatment sc_type
# GSM6090425 SRR19039353 SAMN28041187      PDAC_U_1_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090426 SRR19039306 SAMN28041186      PDAC_U_2_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090427 SRR19039376 SAMN28041185      PDAC_U_3_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090428 SRR19039375 SAMN28041184      PDAC_U_4_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090429 SRR19039377 SAMN28041183      PDAC_U_5_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090430 SRR19039374 SAMN28041182      PDAC_U_6_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090431 SRR19039373 SAMN28041181      PDAC_U_7_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090432 SRR19039372 SAMN28041180      PDAC_U_8_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090433 SRR19039371 SAMN28041179      PDAC_U_8_R2        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090434 SRR19039370 SAMN28041178      PDAC_U_9_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090435 SRR19039369 SAMN28041177     PDAC_U_10_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090436 SRR19039367 SAMN28041176     PDAC_U_10_R2        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090437 SRR19039366 SAMN28041175     PDAC_U_11_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090438 SRR19039368 SAMN28041174     PDAC_U_12_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090439 SRR19039365 SAMN28041173     PDAC_U_13_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090440 SRR19039352 SAMN28041172     PDAC_U_14_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090441 SRR19039350 SAMN28041171     PDAC_U_15_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090442 SRR19039349 SAMN28041170     PDAC_U_16_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090443 SRR19039348 SAMN28041169     PDAC_U_17_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090444 SRR19039347 SAMN28041168     PDAC_U_17_R2        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090445 SRR19039346 SAMN28041167     PDAC_U_18_R1        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090446 SRR19039345 SAMN28041166     PDAC_U_18_R2        Pancreas                                Primary tumor Treatment naïve                     none   snRNA
# GSM6090447 SRR19039344 SAMN28041165      PDAC_T_1_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090448 SRR19039343 SAMN28041164      PDAC_T_1_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090449 SRR19039342 SAMN28041163      PDAC_T_2_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090452 SRR19039339 SAMN28041160      PDAC_T_4_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090453 SRR19039338 SAMN28041159      PDAC_T_5_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090454 SRR19039337 SAMN28041158      PDAC_T_6_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090455 SRR19039336 SAMN28041157      PDAC_T_6_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090456 SRR19039335 SAMN28041156      PDAC_T_7_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090457 SRR19039334 SAMN28041155      PDAC_T_8_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090458 SRR19039333 SAMN28041154      PDAC_T_8_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090459 SRR19039332 SAMN28041153      PDAC_T_8_R3        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090460 SRR19039331 SAMN28041152      PDAC_T_9_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090461 SRR19039330 SAMN28041151      PDAC_T_9_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090462 SRR19039329 SAMN28041150     PDAC_T_10_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090463 SRR19039328 SAMN28041149     PDAC_T_10_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090464 SRR19039327 SAMN28041148     PDAC_T_11_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090465 SRR19039326 SAMN28041147     PDAC_T_11_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090466 SRR19039325 SAMN28041146     PDAC_T_12_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090467 SRR19039324 SAMN28041145     PDAC_T_12_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090468 SRR19039323 SAMN28041144     PDAC_T_13_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090469 SRR19039322 SAMN28041143     PDAC_T_13_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090470 SRR19039321 SAMN28041142     PDAC_T_13_R3        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090471 SRR19039320 SAMN28041141     PDAC_T_14_R1        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090472 SRR19039319 SAMN28041140     PDAC_T_14_R2        Pancreas                                Primary tumor       Treatment                      CRT   snRNA
# GSM6090473 SRR19039318 SAMN28041139     PDAC_T_15_R1        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090474 SRR19039317 SAMN28041138     PDAC_T_15_R2        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090475 SRR19039316 SAMN28041137     PDAC_T_16_R1        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090476 SRR19039315 SAMN28041136     PDAC_T_17_R1        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090477 SRR19039351 SAMN28041135     PDAC_T_18_R1        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090478 SRR19039314 SAMN28041134     PDAC_T_19_R1        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090479 SRR19039313 SAMN28041133     PDAC_T_19_R2        Pancreas                                Primary tumor       Treatment                     CRTL   snRNA
# GSM6090480 SRR19039312 SAMN28041132     PDAC_T_20_R1        Pancreas                                Primary tumor       Treatment                     CRTN   snRNA
# GSM6090481 SRR19039311 SAMN28041131     PDAC_T_21_R1        Pancreas                                Primary tumor       Treatment                    CRTLN   snRNA
# GSM6090482 SRR19039310 SAMN28041130     PDAC_T_21_R2        Pancreas                                Primary tumor       Treatment                    CRTLN   snRNA
# GSM6090483 SRR19039309 SAMN28041129     PDAC_T_22_R1        Pancreas                                Primary tumor       Treatment FOLFIRINOX, Gy/cisplatin   snRNA
# GSM6090484 SRR19039308 SAMN28041128     PDAC_T_23_R1        Pancreas                                Primary tumor       Treatment FOLFIRINOX, Gy/cisplatin   snRNA
# GSM6090486 SRR19039305 SAMN28041126     PDAC_T_24_R2        Pancreas                                Primary tumor       Treatment    Gem/abraxane, Gy/cape   snRNA
# GSM6090487 SRR19039304 SAMN28041125     PDAC_T_25_R1        Pancreas                                Primary tumor       Treatment                  Gy/cape   snRNA
# GSM6090488 SRR19039364 SAMN28041124 Normal_MGHR13_R1        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090489 SRR19039363 SAMN28041123 Normal_MGHR13_R2        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090490 SRR19039362 SAMN28041122  Normal_2675N_R1        Pancreas                                       Normal       Treatment                      CRT   snRNA
# GSM6090491 SRR19039361 SAMN28041121   Normal_2548_R1        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090492 SRR19039360 SAMN28041120  Normal_2540N_R1        Pancreas                                       Normal       Treatment                      CRT   snRNA
# GSM6090493 SRR19039359 SAMN28041119   Normal_2517_R1        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090494 SRR19039358 SAMN28041118   Normal_2517_R2        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090495 SRR19039357 SAMN28041117   Normal_2517_R3        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090496 SRR19039356 SAMN28041116   Normal_010N_R1        Pancreas                                       Normal Treatment naïve                     none   snRNA
# GSM6090497 SRR19039355 SAMN28041115     010orgCRT_R1        Pancreas Organoid cultures from Primary tumor samples       Treatment                      CRT   snRNA
# GSM6090498 SRR19039354 SAMN28041114        010nuc_R1        Pancreas Organoid cultures from Primary tumor samples Treatment naïve                     none   snRNA

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "SRR")

# Updating Study: GSE202051 
# ✅ Matched 147643 of 147643 samples
# Detected overlapping columns: sample_id, sample_name, tissue_location, tissue_type, Treatment, pre_treatment, sc_type 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment 
# Updated column: sc_type 

##############################################################################################
# gse205013 conflict for P11
# gsm ids
geo_nr <- "GSE205013"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=6))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("SAMN","sample_name", "tissue_location", "tissue_type", "Treatment", "pre_treatment")

for (x in 1:length(names(gse@gsms))){
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  tissue_index <- grep("tissue:", gse@gsms[[x]]@header$characteristics_ch1)  
  tissue_location <-  gsub("tissue: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  
  if (tissue_location == "PDAC liver met"){
    tissue_location <- "Liver"
    tissue_type <- "Metastatic lesion"
  } else if (tissue_location == "Primary PDAC"){
    tissue_location <- "Pancreas"
    tissue_type <- "Primary tumor"}
  
  treatment_index <- grep("treatment:", gse@gsms[[x]]@header$characteristics_ch1)  
  pre_treatment <-  gsub("treatment: ", "", gse@gsms[[x]]@header$characteristics_ch1[treatment_index])
  
  if (pre_treatment == "Untreated"){
    pre_treatment <- "none"
    treatment_info <- "Treatment naïve"
  }else if (pre_treatment == "Treated"){
    treatment_info <- "Treatment"
    pre_treatment <- NA}
  
  sample_name <-  gse@gsms[[x]]@header$title
  
  
  df[x,]<- c(samn_id, sample_name, tissue_location, tissue_type, treatment_info, pre_treatment)}

df
#                    SAMN sample_name tissue_location       tissue_type       Treatment pre_treatment
# GSM6204109 SAMN28701976         P01           Liver Metastatic lesion Treatment naïve          none
# GSM6204110 SAMN28701975         P02           Liver Metastatic lesion Treatment naïve          none
# GSM6204111 SAMN28701974         P03        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204112 SAMN28701973         P04        Pancreas     Primary tumor Treatment naïve          none
# GSM6204113 SAMN28701972         P05        Pancreas     Primary tumor Treatment naïve          none
# GSM6204114 SAMN28701971         P06        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204115 SAMN28701970         P07        Pancreas     Primary tumor Treatment naïve          none
# GSM6204116 SAMN28701969         P08        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204117 SAMN28701968         P09        Pancreas     Primary tumor Treatment naïve          none
# GSM6204118 SAMN28701967         P10        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204119 SAMN28701966         P11           Liver Metastatic lesion Treatment naïve          none
# GSM6204120 SAMN28701965         P12        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204121 SAMN28701964         P13        Pancreas     Primary tumor Treatment naïve          none
# GSM6204122 SAMN28701963         P14        Pancreas     Primary tumor       Treatment          <NA>
# GSM6204123 SAMN28701962         P15        Pancreas     Primary tumor Treatment naïve          none
# GSM6204124 SAMN28701961         P16           Liver Metastatic lesion Treatment naïve          none
# GSM6204125 SAMN28701960         P17           Liver Metastatic lesion       Treatment          <NA>
# GSM6204126 SAMN28701959         P18           Liver Metastatic lesion Treatment naïve          none
# GSM6204127 SAMN28701958         P19        Pancreas     Primary tumor Treatment naïve          none
# GSM6204128 SAMN28701957         P20        Pancreas     Primary tumor Treatment naïve          none
# GSM6204129 SAMN28701956         P21           Liver Metastatic lesion Treatment naïve          none
# GSM6204130 SAMN28701955         P22        Pancreas     Primary tumor Treatment naïve          none
# GSM6204131 SAMN28701954         P23        Pancreas     Primary tumor Treatment naïve          none
# GSM6204132 SAMN28701953         P24           Liver Metastatic lesion Treatment naïve          none
# GSM6204133 SAMN28701952         P25           Liver Metastatic lesion Treatment naïve          none
# GSM6204134 SAMN28701951         P26        Pancreas     Primary tumor Treatment naïve          none
# GSM6204135 SAMN28701950         P27           Liver Metastatic lesion Treatment naïve          none

# adjust df
df <- df %>%
  dplyr::rename(
    sample_id = SAMN)
df$GSM <- rownames(df)

# add pre_treatment info
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-36296-4/MediaObjects/41467_2023_36296_MOESM1_ESM.pdf
# Supplementary Table 2. Clinical outcomes for treated patients. 
treatment_info <- "Patient 
Treatment 
Treatment Duration 
RECIST Criteria Response 

P03 
FFX 
7M 
PR

P06 
FFX, G/A 
17M, 3M 
PD 

P08 
G/A 
4M 
SD 

P10 
FFX/SBRT 
5M 
SD 

P11 
FFX 
6M 
PD 

P12 
FFX, G/A 
20M, 2M 
PD

P14 
FFX 
4M 
PR 

P17 
FFX 
3M 
Non-evaluable"

# convert into table
treatment_info <- convert_text_with_empty_lines_into_table(text = treatment_info, n_cols=4)
treatment_info
#   Patient Treatment Treatment Duration RECIST Criteria Response
# 1     P03       FFX                 7M                       PR
# 2     P06  FFX, G/A            17M, 3M                       PD
# 3     P08       G/A                 4M                       SD
# 4     P10  FFX/SBRT                 5M                       SD
# 5     P11       FFX                 6M                       PD
# 6     P12  FFX, G/A            20M, 2M                       PD
# 7     P14       FFX                 4M                       PR
# 8     P17       FFX                 3M            Non-evaluable

df$pre_treatment <- df$sample_name
df$pre_treatment <- ifelse(
  df$pre_treatment %in% treatment_info$Patient,
  treatment_info$Treatment[match(df$pre_treatment, treatment_info$Patient)],
  "none")

# check treatment
df <- df %>%
  mutate(treatment_check= case_when(
    Treatment == "Treatment"& pre_treatment !="none"  ~ "match",
    Treatment == "Treatment naïve" & pre_treatment == "none"  ~ "match",
    TRUE ~ "conflict"))
df
#               sample_id sample_name tissue_location       tissue_type       Treatment pre_treatment        GSM treatment_check
# GSM6204109 SAMN28701976         P01           Liver Metastatic lesion Treatment naïve          none GSM6204109           match
# GSM6204110 SAMN28701975         P02           Liver Metastatic lesion Treatment naïve          none GSM6204110           match
# GSM6204111 SAMN28701974         P03        Pancreas     Primary tumor       Treatment           FFX GSM6204111           match
# GSM6204112 SAMN28701973         P04        Pancreas     Primary tumor Treatment naïve          none GSM6204112           match
# GSM6204113 SAMN28701972         P05        Pancreas     Primary tumor Treatment naïve          none GSM6204113           match
# GSM6204114 SAMN28701971         P06        Pancreas     Primary tumor       Treatment      FFX, G/A GSM6204114           match
# GSM6204115 SAMN28701970         P07        Pancreas     Primary tumor Treatment naïve          none GSM6204115           match
# GSM6204116 SAMN28701969         P08        Pancreas     Primary tumor       Treatment           G/A GSM6204116           match
# GSM6204117 SAMN28701968         P09        Pancreas     Primary tumor Treatment naïve          none GSM6204117           match
# GSM6204118 SAMN28701967         P10        Pancreas     Primary tumor       Treatment      FFX/SBRT GSM6204118           match
# GSM6204119 SAMN28701966         P11           Liver Metastatic lesion Treatment naïve           FFX GSM6204119        conflict
# GSM6204120 SAMN28701965         P12        Pancreas     Primary tumor       Treatment      FFX, G/A GSM6204120           match
# GSM6204121 SAMN28701964         P13        Pancreas     Primary tumor Treatment naïve          none GSM6204121           match
# GSM6204122 SAMN28701963         P14        Pancreas     Primary tumor       Treatment           FFX GSM6204122           match
# GSM6204123 SAMN28701962         P15        Pancreas     Primary tumor Treatment naïve          none GSM6204123           match
# GSM6204124 SAMN28701961         P16           Liver Metastatic lesion Treatment naïve          none GSM6204124           match
# GSM6204125 SAMN28701960         P17           Liver Metastatic lesion       Treatment           FFX GSM6204125           match
# GSM6204126 SAMN28701959         P18           Liver Metastatic lesion Treatment naïve          none GSM6204126           match
# GSM6204127 SAMN28701958         P19        Pancreas     Primary tumor Treatment naïve          none GSM6204127           match
# GSM6204128 SAMN28701957         P20        Pancreas     Primary tumor Treatment naïve          none GSM6204128           match
# GSM6204129 SAMN28701956         P21           Liver Metastatic lesion Treatment naïve          none GSM6204129           match
# GSM6204130 SAMN28701955         P22        Pancreas     Primary tumor Treatment naïve          none GSM6204130           match
# GSM6204131 SAMN28701954         P23        Pancreas     Primary tumor Treatment naïve          none GSM6204131           match
# GSM6204132 SAMN28701953         P24           Liver Metastatic lesion Treatment naïve          none GSM6204132           match
# GSM6204133 SAMN28701952         P25           Liver Metastatic lesion Treatment naïve          none GSM6204133           match
# GSM6204134 SAMN28701951         P26        Pancreas     Primary tumor Treatment naïve          none GSM6204134           match
# GSM6204135 SAMN28701950         P27           Liver Metastatic lesion Treatment naïve          none GSM6204135           match

# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA843078&search=GSM6204119&o=acc_s%3Aa
# as well as soft file for GSE205013 -> P11 as untreated sample?

table(pdac_metadata$Treatment[pdac_metadata$sample_id == df$GSM[df$treatment_check == "conflict"]])
# Treatment naïve 
# 4278 
table(pdac_metadata$pre_treatment[pdac_metadata$sample_id == df$GSM[df$treatment_check == "conflict"]], useNA="ifany")
# <NA> 
#   4278 

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "GSM")

# Updating Study: GSE205013 
# ✅ Matched 167366 of 167366 samples
# Detected overlapping columns: sample_id, sample_name, tissue_location, tissue_type, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: tissue_type 
# Updated column: Treatment 
# Updated column: pre_treatment

##############################################################################################
# GSE211644
# mda1 only, mda2 == peng data
# tissue (CD3+ enriched != Primary, Infusion Product from CD3+ enriched != Primary, U -> Adjacent normal) 
# treatment likely after sample acquisition (geo -> samples are collected before treatment)
# https://pmc.ncbi.nlm.nih.gov/articles/instance/9547957/bin/NIHMS1826134-supplement-2.pdf
# srr ids
geo_nr <- "GSE211644"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=6))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("sample_name", "tissue_type", "SAMN", "SRR", "Treatment", "pre_treatment")

for (x in names(gse@gsms)){
  sample_name <- gse@gsms[[x]]@header$title
  if (grepl("_IP$", sample_name)){
    tissue_type <- "Infusion Product, Tumor-derived"
  }else if (grepl("_T", sample_name)){
    tissue_type <- "Primary tumor, CD3pos cells"
  }else if (grepl("_U", sample_name)){
    tissue_type <- "Adjacent normal, CD3pos cells"}
  
  samn_index <- grep("BioSample:", gse@gsms[[x]]@header$relation)
  samn_id <-  gsub("^.*biosample/", "", gse@gsms[[x]]@header$relation[samn_index])
  
  srr_index <- grep("SRA:", gse@gsms[[x]]@header$relation)
  srr_id <-  gsub("^.*term\\=", "", gse@gsms[[x]]@header$relation[srr_index])
  srr_id <- get_srr_id(srr_id)
  
  treatment_info <- gse@gsms[[x]]@header$treatment_protocol_ch1
  treatment_info <- ifelse(treatment_info == "samples are collected before treatment", "Treatment naïve", "?")
  pre_treatment <- ifelse(treatment_info ==  "Treatment naïve", "none", NA)
  df[x,] <- c(sample_name, tissue_type, samn_id,srr_id,treatment_info,pre_treatment)}

# tissue
df$tissue_location <- "Pancreas"

df[df$SRR %in% unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr]),]
#            sample_name                     tissue_type         SAMN         SRR       Treatment pre_treatment tissue_location
# GSM6481941    MDA1_U05         Adjacent normal, CD3pos cells SAMN25871451 SRR18045264 Treatment naïve          none        Pancreas
# GSM6481942    MDA1_T05           Primary tumor, CD3pos cells SAMN25871450 SRR18045265 Treatment naïve          none        Pancreas
# GSM6481943    MDA1_T02           Primary tumor, CD3pos cells SAMN25871447 SRR18045268 Treatment naïve          none        Pancreas
# GSM6481944    MDA1_T03           Primary tumor, CD3pos cells SAMN25871448 SRR18045267 Treatment naïve          none        Pancreas
# GSM6481945    MDA1_T04           Primary tumor, CD3pos cells SAMN25871449 SRR18045266 Treatment naïve          none        Pancreas
# GSM6481946    MDA1_T06           Primary tumor, CD3pos cells SAMN25871452 SRR18045263 Treatment naïve          none        Pancreas
# GSM6481947    MDA1_T07           Primary tumor, CD3pos cells SAMN25871453 SRR18045262 Treatment naïve          none        Pancreas
# GSM6481948    MDA1_T01           Primary tumor, CD3pos cells SAMN25948361 SRR18025515 Treatment naïve          none        Pancreas
# GSM6541649 MDA1_T01_IP       Infusion Product, Tumor-derived SAMN30431574 SRR21152734 Treatment naïve          none        Pancreas
# GSM6541650 MDA1_T02_IP       Infusion Product, Tumor-derived SAMN30431575 SRR21152733 Treatment naïve          none        Pancreas
# GSM6541651 MDA1_T03_IP       Infusion Product, Tumor-derived SAMN30431576 SRR21152732 Treatment naïve          none        Pancreas
# GSM6541652 MDA1_T04_IP       Infusion Product, Tumor-derived SAMN30431577 SRR21152731 Treatment naïve          none        Pancreas
# GSM6541653 MDA1_T05_IP       Infusion Product, Tumor-derived SAMN30431578 SRR21152730 Treatment naïve          none        Pancreas
# GSM6541654 MDA1_T07_IP       Infusion Product, Tumor-derived SAMN30431579 SRR21152729 Treatment naïve          none        Pancreas

# adjust df
df <- df %>%
  dplyr::rename(
    sample_id = SAMN)

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "SRR")

# Updating Study: GSE211644 
# ✅ Matched 40971 of 40971 samples
# Detected overlapping columns: sample_name, tissue_type, sample_id, Treatment, pre_treatment, tissue_location 
# Updated column: sample_name 
# Updated column: tissue_type 
# Updated column: sample_id 
# Updated column: Treatment 
# Updated column: pre_treatment 
# Updated column: tissue_location 

##############################################################################################
# gse229413 correct
# brain dead donors with healthy pancreas tissue
# but: Treatment due to distinct health issues would be possible?
# gsm ids -> no samn ids available?
geo_nr <- "GSE229413"

# get geo info
gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
# wget file_name
gse <- getGEO(filename = paste0(data_dir, geo_nr,"_family.soft.gz"))

unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr])
# [1] "GSM7162998_3829-EC"  "GSM7162999_3861-EC"  "GSM7163000_3862-EC"  "GSM7163001_4160-EC"  "GSM7163002_4161-EC"  "GSM7163003_4346-EC"  "GSM7163004_4347-EC"  "GSM7163007_4637-EC1" "GSM7163008_4637-EC2" "GSM7163009_4741-EC1" "GSM7163010_4741-EC2"
samples <- names(gse@gsms)[names(gse@gsms) %in% gsub("_.*$","", unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == geo_nr]))]
samples
# [1] "GSM7162998" "GSM7162999" "GSM7163000" "GSM7163001" "GSM7163002" "GSM7163003" "GSM7163004" "GSM7163007" "GSM7163008" "GSM7163009" "GSM7163010"
df <- data.frame(matrix(nrow = length(samples), ncol=4))
rownames(df) <- samples
colnames(df) <- c("tissue_type", "sample_name", "patient_id", "tissue_location")

for (x in samples){
  tissue_type <- ifelse(grepl("Donor", gse@gsms[[x]]@header$title), "Healthy Donor", "?")
  sample_id <- gse@gsms[[x]]@header$title
  pat_id <- sub("^(Donor_[0-9]+).*", "\\1", sample_id)
  tissue_location <-  ifelse(grepl("Pancreas", gse@gsms[[x]]@header$source_name_ch1), "Pancreas", "?")
  df[x,] <- c( tissue_type, sample_id, pat_id, tissue_location)}
df
#              tissue_type           sample_name patient_id tissue_location
# GSM7162998 Healthy Donor Donor_1_head_and_tail    Donor_1        Pancreas
# GSM7162999 Healthy Donor          Donor_2_head    Donor_2        Pancreas
# GSM7163000 Healthy Donor          Donor_2_tail    Donor_2        Pancreas
# GSM7163001 Healthy Donor          Donor_4_head    Donor_4        Pancreas
# GSM7163002 Healthy Donor          Donor_4_tail    Donor_4        Pancreas
# GSM7163003 Healthy Donor          Donor_5_head    Donor_5        Pancreas
# GSM7163004 Healthy Donor          Donor_5_tail    Donor_5        Pancreas
# GSM7163007 Healthy Donor          Donor_8_head    Donor_8        Pancreas
# GSM7163008 Healthy Donor          Donor_8_tail    Donor_8        Pancreas
# GSM7163009 Healthy Donor          Donor_9_head    Donor_9        Pancreas
# GSM7163010 Healthy Donor          Donor_9_tail    Donor_9        Pancreas

# adjust df
df$GSM <- rownames(df)

# check treatment
unique(pdac_metadata$Treatment[pdac_metadata$Study_ID == geo_nr])
# [1] "N_A"
unique(pdac_metadata$pre_treatment[pdac_metadata$Study_ID == geo_nr])
# [1] NA

#################################
# in case of no prior treatment
df$Treatment <- "Treatment naïve"
df$pre_treatment <- "none"

df
#              tissue_type           sample_name patient_id tissue_location        GSM       Treatment pre_treatment
# GSM7162998 Healthy Donor Donor_1_head_and_tail    Donor_1        Pancreas GSM7162998 Treatment naïve          none
# GSM7162999 Healthy Donor          Donor_2_head    Donor_2        Pancreas GSM7162999 Treatment naïve          none
# GSM7163000 Healthy Donor          Donor_2_tail    Donor_2        Pancreas GSM7163000 Treatment naïve          none
# GSM7163001 Healthy Donor          Donor_4_head    Donor_4        Pancreas GSM7163001 Treatment naïve          none
# GSM7163002 Healthy Donor          Donor_4_tail    Donor_4        Pancreas GSM7163002 Treatment naïve          none
# GSM7163003 Healthy Donor          Donor_5_head    Donor_5        Pancreas GSM7163003 Treatment naïve          none
# GSM7163004 Healthy Donor          Donor_5_tail    Donor_5        Pancreas GSM7163004 Treatment naïve          none
# GSM7163007 Healthy Donor          Donor_8_head    Donor_8        Pancreas GSM7163007 Treatment naïve          none
# GSM7163008 Healthy Donor          Donor_8_tail    Donor_8        Pancreas GSM7163008 Treatment naïve          none
# GSM7163009 Healthy Donor          Donor_9_head    Donor_9        Pancreas GSM7163009 Treatment naïve          none
# GSM7163010 Healthy Donor          Donor_9_tail    Donor_9        Pancreas GSM7163010 Treatment naïve          none
#################################

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = geo_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "GSM",
  clean_pattern = "_.*$")

# Updating Study: GSE229413 
# ✅ Matched 33309 of 33309 samples
# Detected overlapping columns: tissue_type, sample_name, tissue_location, Treatment, pre_treatment 
# Updated column: tissue_type 
# Updated column: sample_name 
# Updated column: tissue_location 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# phs
# srr ids
study_nr <- "phs001840.v1.p1"

# For one tumor sample, viable, CD45-negative, CD31-negative, and EpCAM-negative cells were also isolated to enrich for CAFs.
# https://aacrjournals.org/cancerdiscovery/article/9/8/1102/42174/Cross-Species-Single-Cell-Analysis-of-Pancreatic
# 6 untreated pdac patients
# Using this approach, we were able to isolate a fibroblast-enriched fraction from sample hT137 (hT137 FE).
# hT97 hT99 hT103 hT137 hN149 hT149 hN150 hT150

# mapping approach generated with help of chat gpt
phs001840_map <- data.frame(
  raw_name = c("DT17007","DT17012","DT17024","DT17039", "DT17040",
               "DT17042","DT17043","DT17044","DT17045"),
  study_sample_id = c("hT97","hT99","hT103","hT137", "hT137 FE",
                      "hN149","hT149","hN150","hT150"),
  enrichment = c(rep("Standard",4),"CAF-enriched",
                 rep("Standard",4)))

phs001840_map
# raw_name study_sample_id   enrichment
# 1  DT17007            hT97     Standard
# 2  DT17012            hT99     Standard
# 3  DT17024           hT103     Standard
# 4  DT17039           hT137     Standard
# 5  DT17040        hT137 FE CAF-enriched
# 6  DT17042           hN149     Standard
# 7  DT17043           hT149     Standard
# 8  DT17044           hN150     Standard
# 9  DT17045           hT150     Standard

table(pdac_metadata$raw_name[pdac_metadata$Study_ID == study_nr])
# SRR9274536 SRR9274537 SRR9274538 SRR9274539 SRR9274540 SRR9274541 SRR9274542 SRR9274543 SRR9274544 
# 2229       5298       1150       1192        865        643       1677       3538       1693

srr_id <- unique(pdac_metadata$raw_name[pdac_metadata$Study_ID == study_nr])

#  obtain corresponding samn  id & tissue information
phs_info <- sapply(srr_id, get_samn_and_info_phs) 
phs_info
#       SRR9274536       SRR9274537       SRR9274538       SRR9274539       SRR9274540       SRR9274541        SRR9274542       SRR9274543        SRR9274544      
# [1,] "Primary tumor"  "Primary tumor"  "Primary tumor"  "Primary tumor"  "Primary tumor"  "Adjacent normal" "Primary tumor"  "Adjacent normal" "Primary tumor" 
# [2,] "SAMN11970422"   "SAMN11970417"   "SAMN11970415"   "SAMN11970423"   "SAMN11970416"   "SAMN11970418"    "SAMN11970419"   "SAMN11970421"    "SAMN11970420"  
# [3,] "Sample DT17007" "Sample DT17012" "Sample DT17024" "Sample DT17039" "Sample DT17040" "Sample DT17042"  "Sample DT17043" "Sample DT17044"  "Sample DT17045"

write.csv(phs_info, paste0(data_dir, "samn_and_info_",study_nr, ".csv"))

df <- read.csv(paste0(data_dir, "samn_and_info_",study_nr, ".csv"), row.names = 1)
df <- as.data.frame(t(df))
df$SRR <- rownames(df)

df <- df %>%
  dplyr::rename(
    tissue_type = 1,
    sample_id = 2,
    sample_name = 3)
df$sample_name <- gsub("Sample ", "", df$sample_name)
df$Treatment <- "Treatment naïve"
df$pre_treatment <- "none"
df$sample_name <- paste0(df$sample_name, ", ", phs001840_map$study_sample_id[match(df$sample_name, phs001840_map$raw_name)])
df$tissue_type[grep("FE", df$sample_name)] <- paste0(df$tissue_type[grep("FE", df$sample_name)], ", CAF-enriched")
df
#                            tissue_type    sample_id       sample_name        SRR       Treatment pre_treatment
# SRR9274536               Primary tumor SAMN11970422     DT17007, hT97 SRR9274536 Treatment naïve          none
# SRR9274537               Primary tumor SAMN11970417     DT17012, hT99 SRR9274537 Treatment naïve          none
# SRR9274538               Primary tumor SAMN11970415    DT17024, hT103 SRR9274538 Treatment naïve          none
# SRR9274539               Primary tumor SAMN11970423    DT17039, hT137 SRR9274539 Treatment naïve          none
# SRR9274540 Primary tumor, CAF-enriched SAMN11970416 DT17040, hT137 FE SRR9274540 Treatment naïve          none
# SRR9274541             Adjacent normal SAMN11970418    DT17042, hN149 SRR9274541 Treatment naïve          none
# SRR9274542               Primary tumor SAMN11970419    DT17043, hT149 SRR9274542 Treatment naïve          none
# SRR9274543             Adjacent normal SAMN11970421    DT17044, hN150 SRR9274543 Treatment naïve          none
# SRR9274544               Primary tumor SAMN11970420    DT17045, hT150 SRR9274544 Treatment naïve          none

# update metadata based on df
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = study_nr,
  df = df,
  match_meta_col = "sample_id",
  match_df_col = "SRR")

# Updating Study: phs001840.v1.p1 
# ✅ Matched 18285 of 18285 samples
# Detected overlapping columns: tissue_type, sample_id, sample_name, Treatment, pre_treatment 
# Updated column: tissue_type 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: Treatment 
# Updated column: pre_treatment 

##############################################################################################
# peng
# https://pubmed.ncbi.nlm.nih.gov/31273297/
# Adjacent normal -> normal (control pancreas samples, healthy donors)
# crr ids
study_nr <- "PRJCA001063"
table(pdac_metadata$raw_name[pdac_metadata$Study_ID == study_nr])

# new_CRR034499 new_CRR034500 new_CRR034501 new_CRR034503 new_CRR034504 new_CRR034505 new_CRR034506 new_CRR034507 new_CRR034509 new_CRR034511 new_CRR034512 new_CRR034513 new_CRR034516 new_CRR034517 new_CRR034518 new_CRR034519 new_CRR034520 
# 2215          2615          3650          2466          5840          1278          5097          3334          3271          3651          5998          5051          2085          4148          4529          6370          7361 
# new_CRR034521 new_CRR034522 new_CRR034523 new_CRR034524 new_CRR034525 new_CRR034526 new_CRR034527 new_CRR034528 new_CRR034529 new_CRR034530 new_CRR241798 new_CRR241799 new_CRR241800 new_CRR241801 new_CRR241802 new_CRR241803 new_CRR241804 
# 5269          2343          5235          3062          5181          4464          4502          7648          3474          6153          4045          2730          2335          5423          8333          4577          2034 
# new_CRR241805 
# 2140
# check that only samples in peng start with new, and if so, remove new in sample_name column
if (all(table(pdac_metadata$raw_name[pdac_metadata$Study_ID == study_nr]) == table(pdac_metadata$raw_name[grepl("^new_", pdac_metadata$sample_name)]))){
  pdac_metadata$sample_name <- gsub("^new_", "", pdac_metadata$sample_name)}

table(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr])
# CRR034499 CRR034500 CRR034501 CRR034503 CRR034504 CRR034505 CRR034506 CRR034507 CRR034509 CRR034511 CRR034512 CRR034513 CRR034516 CRR034517 CRR034518 CRR034519 CRR034520 CRR034521 CRR034522 CRR034523 CRR034524 CRR034525 CRR034526 CRR034527 CRR034528 
# 2215      2615      3650      2466      5840      1278      5097      3334      3271      3651      5998      5051      2085      4148      4529      6370      7361      5269      2343      5235      3062      5181      4464      4502      7648 
# CRR034529 CRR034530 CRR241798 CRR241799 CRR241800 CRR241801 CRR241802 CRR241803 CRR241804 CRR241805 
# 3474      6153      4045      2730      2335      5423      8333      4577      2034      2140 


# metadata file including crr codes
# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160
df <- as.data.frame(readxl::read_xlsx(paste0(data_dir, "peng_CRR.xlsx")))
head(df)
#   ID       CRR Run title BioProject accession Experiment accession
# 1  1 CRR241805        T1          PRJCA001063            CRX030762
# 2  2 CRR241798        T2          PRJCA001063            CRX030763
# 3  3 CRR241799        T3          PRJCA001063            CRX030764
# 4  4 CRR034499        T4          PRJCA001063            CRX030765
# 5  5 CRR034500        T5          PRJCA001063            CRX030766
# 6  6 CRR034501        T6          PRJCA001063            CRX030767

df <- df %>%
  mutate(tissue_type = case_when(
    grepl("^T", `Run title`) ~ "Primary tumor",
    grepl("^N", `Run title`) ~ "Healthy Donor"))
df$tissue_location <- "Pancreas"
df$Treatment <- "Treatment naïve"
df$pre_treatment <- "none"
df <- df %>%
  dplyr::rename(
    sample_name = `Run title`,
    sample_id = CRR)

# check mismatch
peng_aacr <-  unique(pdac_metadata$sample_name[pdac_metadata$Study_ID == study_nr])
peng_df <- unique(df$sample_id)
setdiff(peng_aacr, peng_df)
# [1] "CRR034518" # old CRR for T23
setdiff(peng_df, peng_aacr)
# [1] "CRR034510" # T15 in peng_df

# check replacement ids
# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160
# Note: 8 runs in this dataset were replaced on 2021-01-14:
replacements <- "Experiment accession	Run alias	Old run accession	New run accession
CRX030762	T1	CRR034496	CRR241805
CRX030763	T2	CRR034497	CRR241798
CRX030764	T3	CRR034498	CRR241799
CRX030768	T7	CRR034502	CRR241800
CRX030774	T13	CRR034508	CRR241801
CRX030780	T19	CRR034514	CRR241802
CRX030781	T20	CRR034515	CRR241804
CRX030784	T23	CRR034518	CRR241803"

replacements <- read.delim(text = replacements, header=T, sep = "\t")
colnames(replacements) <-c("CRX", "sample_id", "old_CRR", "new_CRR")
replacements
#         CRX sample_id   old_CRR   new_CRR
# 1 CRX030762        T1 CRR034496 CRR241805
# 2 CRX030763        T2 CRR034497 CRR241798
# 3 CRX030764        T3 CRR034498 CRR241799
# 4 CRX030768        T7 CRR034502 CRR241800
# 5 CRX030774       T13 CRR034508 CRR241801
# 6 CRX030780       T19 CRR034514 CRR241802
# 7 CRX030781       T20 CRR034515 CRR241804
# 8 CRX030784       T23 CRR034518 CRR241803

table(pdac_metadata$sample_name[pdac_metadata$sample_name %in% c(setdiff(peng_df, peng_aacr), setdiff(peng_aacr, peng_df), replacements$new_CRR[replacements$old_CRR == setdiff(peng_aacr, peng_df)])])
# CRR034518 CRR241803 
# 4529      4577
# For T23: old and new CRR in sample_names
# No sample_name corresponding to T15

# add info about old CRR IDs
df$old_CRR <- df$sample_id
df$old_CRR <- ifelse(
  df$old_CRR %in% replacements$new_CRR,
  replacements$old_CRR[match(df$old_CRR, replacements$new_CRR)],
  df$old_CRR)
df
#    ID sample_id sample_name BioProject accession Experiment accession   tissue_type tissue_location       Treatment pre_treatment   old_CRR
# 1   1 CRR241805          T1          PRJCA001063            CRX030762 Primary tumor        Pancreas Treatment naïve          none CRR034496
# 2   2 CRR241798          T2          PRJCA001063            CRX030763 Primary tumor        Pancreas Treatment naïve          none CRR034497
# 3   3 CRR241799          T3          PRJCA001063            CRX030764 Primary tumor        Pancreas Treatment naïve          none CRR034498
# 4   4 CRR034499          T4          PRJCA001063            CRX030765 Primary tumor        Pancreas Treatment naïve          none CRR034499
# 5   5 CRR034500          T5          PRJCA001063            CRX030766 Primary tumor        Pancreas Treatment naïve          none CRR034500
# 6   6 CRR034501          T6          PRJCA001063            CRX030767 Primary tumor        Pancreas Treatment naïve          none CRR034501
# 7   7 CRR241800          T7          PRJCA001063            CRX030768 Primary tumor        Pancreas Treatment naïve          none CRR034502
# 8   8 CRR034503          T8          PRJCA001063            CRX030769 Primary tumor        Pancreas Treatment naïve          none CRR034503
# 9   9 CRR034504          T9          PRJCA001063            CRX030770 Primary tumor        Pancreas Treatment naïve          none CRR034504
# 10 10 CRR034505         T10          PRJCA001063            CRX030771 Primary tumor        Pancreas Treatment naïve          none CRR034505
# 11 11 CRR034506         T11          PRJCA001063            CRX030772 Primary tumor        Pancreas Treatment naïve          none CRR034506
# 12 12 CRR034507         T12          PRJCA001063            CRX030773 Primary tumor        Pancreas Treatment naïve          none CRR034507
# 13 13 CRR241801         T13          PRJCA001063            CRX030774 Primary tumor        Pancreas Treatment naïve          none CRR034508
# 14 14 CRR034509         T14          PRJCA001063            CRX030775 Primary tumor        Pancreas Treatment naïve          none CRR034509
# 15 15 CRR034510         T15          PRJCA001063            CRX030776 Primary tumor        Pancreas Treatment naïve          none CRR034510
# 16 16 CRR034511         T16          PRJCA001063            CRX030777 Primary tumor        Pancreas Treatment naïve          none CRR034511
# 17 17 CRR034512         T17          PRJCA001063            CRX030778 Primary tumor        Pancreas Treatment naïve          none CRR034512
# 18 18 CRR034513         T18          PRJCA001063            CRX030779 Primary tumor        Pancreas Treatment naïve          none CRR034513
# 19 19 CRR241802         T19          PRJCA001063            CRX030780 Primary tumor        Pancreas Treatment naïve          none CRR034514
# 20 20 CRR241804         T20          PRJCA001063            CRX030781 Primary tumor        Pancreas Treatment naïve          none CRR034515
# 21 21 CRR034516         T21          PRJCA001063            CRX030782 Primary tumor        Pancreas Treatment naïve          none CRR034516
# 22 22 CRR034517         T22          PRJCA001063            CRX030783 Primary tumor        Pancreas Treatment naïve          none CRR034517
# 23 23 CRR241803         T23          PRJCA001063            CRX030784 Primary tumor        Pancreas Treatment naïve          none CRR034518
# 24 24 CRR034519         T24          PRJCA001063            CRX030785 Primary tumor        Pancreas Treatment naïve          none CRR034519
# 25 25 CRR034520          N1          PRJCA001063            CRX030786 Healthy Donor        Pancreas Treatment naïve          none CRR034520
# 26 26 CRR034521          N2          PRJCA001063            CRX030787 Healthy Donor        Pancreas Treatment naïve          none CRR034521
# 27 27 CRR034522          N3          PRJCA001063            CRX030788 Healthy Donor        Pancreas Treatment naïve          none CRR034522
# 28 28 CRR034523          N4          PRJCA001063            CRX030789 Healthy Donor        Pancreas Treatment naïve          none CRR034523
# 29 29 CRR034524          N5          PRJCA001063            CRX030790 Healthy Donor        Pancreas Treatment naïve          none CRR034524
# 30 30 CRR034525          N6          PRJCA001063            CRX030791 Healthy Donor        Pancreas Treatment naïve          none CRR034525
# 31 31 CRR034526          N7          PRJCA001063            CRX030792 Healthy Donor        Pancreas Treatment naïve          none CRR034526
# 32 32 CRR034527          N8          PRJCA001063            CRX030793 Healthy Donor        Pancreas Treatment naïve          none CRR034527
# 33 33 CRR034528          N9          PRJCA001063            CRX030794 Healthy Donor        Pancreas Treatment naïve          none CRR034528
# 34 34 CRR034529         N10          PRJCA001063            CRX030795 Healthy Donor        Pancreas Treatment naïve          none CRR034529
# 35 35 CRR034530         N11          PRJCA001063            CRX030796 Healthy Donor        Pancreas Treatment naïve          none CRR034530

# update metadata based on df
# using the new CRR IDs
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = study_nr,
  df = df,
  match_meta_col = "sample_name",
  match_df_col = "sample_id")
# Updating Study: PRJCA001063 
# ✅ Matched 143378 of 147907 samples
# ⚠️ 4529 samples had no match.
# Detected overlapping columns: sample_id, sample_name, tissue_type, tissue_location, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_type 
# Updated column: tissue_location 
# Updated column: Treatment 
# Updated column: pre_treatment

# using the old CRR IDs
pdac_metadata <- update_metadata_from_df(
  meta = pdac_metadata,
  study_id = study_nr,
  df = df,
  match_meta_col = "sample_name",
  match_df_col = "old_CRR")
# Updating Study: PRJCA001063 
# ✅ Matched 4529 of 147907 samples
# ⚠️ 143378 samples had no match.
# Detected overlapping columns: sample_id, sample_name, tissue_type, tissue_location, Treatment, pre_treatment 
# Updated column: sample_id 
# Updated column: sample_name 
# Updated column: tissue_type 
# Updated column: tissue_location 
# Updated column: Treatment 
# Updated column: pre_treatment 


##############################################################################################

# Check seq_type
table(pdac_metadata$sc_type, pdac_metadata$Study_ID, useNA="ifany")
#       EGAS00001002543 GSE154778 GSE155698 GSE156405 GSE158356 GSE194247 GSE202051 GSE205013 GSE211644 GSE229413 phs001840.v1.p1 PRJCA001063
# scRNA           76094     20364     23991     17809      2746     29622         0    167366     40971     33309           18285      147907
# snRNA               0         0         0         0         0         0    147643         0         0         0               0           0

# Check tissue type
table(pdac_metadata$tissue_type, pdac_metadata$Study_ID, useNA="ifany")
#                                              EGAS00001002543 GSE154778 GSE155698 GSE156405 GSE158356 GSE194247 GSE202051 GSE205013 GSE211644 GSE229413 phs001840.v1.p1 PRJCA001063
# Adjacent normal                                            0         0      4502         0         0         0     10660         0         0         0            4181           0
# Adjacent normal, CD3pos cells                              0         0         0         0         0         0         0         0       473         0               0           0
# Healthy Donor                                              0         0         0         0         0         0         0         0         0     33309               0       54692
# Infusion Product, Tumor-derived                            0         0         0         0         0         0         0         0     37360         0               0           0
# Metastatic lesion                                          0      8464         0      9546      2746         0         0     37340         0         0               0           0
# Organoid cultures from Primary tumor samples               0         0         0         0         0         0      4109         0         0         0               0           0
# Primary tumor                                          76094     11900     19489      8263         0         0    132874    130026         0         0           13239       93215
# Primary tumor, CAF-enriched                                0         0         0         0         0         0         0         0         0         0             865           0
# Primary tumor, CD3pos cells                                0         0         0         0         0         0         0         0      3138         0               0           0
# Primary tumor, CD45neg cells                               0         0         0         0         0     29622         0         0         0         0               0           0

# Check tissue location
table(pdac_metadata$tissue_location, pdac_metadata$Study_ID, useNA="ifany")
#              EGAS00001002543 GSE154778 GSE155698 GSE156405 GSE158356 GSE194247 GSE202051 GSE205013 GSE211644 GSE229413 phs001840.v1.p1 PRJCA001063
# Liver                      0      8418         0      4639      2746         0         0     37340         0         0               0           0
# Lung                       0         0         0      1225         0         0         0         0         0         0               0           0
# Omentum                    0        46         0         0         0         0         0         0         0         0               0           0
# Pancreas               76094     11900     23991      8263         0     29622    147643    130026     40971     33309           18285      147907
# Peritoneum                 0         0         0      2207         0         0         0         0         0         0               0           0
# Vaginal apex               0         0         0      1475         0         0         0         0         0         0               0           0


# Check treatment
table(pdac_metadata$Treatment, pdac_metadata$Study_ID, useNA="ifany")
#                 EGAS00001002543 GSE154778 GSE155698 GSE156405 GSE158356 GSE194247 GSE202051 GSE205013 GSE211644 GSE229413 phs001840.v1.p1 PRJCA001063
# ?                             0     20000         0         0      2746         0         0         0         0         0               0           0
# Treatment                 11536         0         0      6211         0         0     80811     56781         0         0               0           0
# Treatment naïve           64558       364     23991     11598         0     29622     66832    110585     40971     33309           18285      147907

# Check type of prior treatment
table(pdac_metadata$pre_treatment, pdac_metadata$Study_ID, useNA="ifany")
#                                                   EGAS00001002543 GSE154778 GSE155698 GSE156405 GSE158356 GSE194247 GSE202051 GSE205013 GSE211644 GSE229413 phs001840.v1.p1 PRJCA001063
# ?                                                               0     20000         0         0      2746         0         0         0         0         0               0           0
# CRT                                                             0         0         0         0         0         0     48550         0         0         0               0           0
# CRTL                                                            0         0         0         0         0         0     13215         0         0         0               0           0
# CRTLN                                                           0         0         0         0         0         0      7215         0         0         0               0           0
# CRTN                                                            0         0         0         0         0         0      3139         0         0         0               0           0
# FFX                                                             0         0         0         0         0         0         0     45317         0         0               0           0
# FFX, G/A                                                        0         0         0         0         0         0         0      5226         0         0               0           0
# FFX/SBRT                                                        0         0         0         0         0         0         0      4793         0         0               0           0
# FOLFIRI, evofosfamide/ipilimumab and capecitabine               0         0         0      4907         0         0         0         0         0         0               0           0
# FOLFIRINOX, Gy/cisplatin                                        0         0         0         0         0         0      6412         0         0         0               0           0
# G/A                                                             0         0         0         0         0         0         0      5723         0         0               0           0
# Gem/abraxane, Gy/cape                                           0         0         0         0         0         0      1140         0         0         0               0           0
# gemcitabine/paclitaxel                                          0         0         0      1304         0         0         0         0         0         0               0           0
# Gy/cape                                                         0         0         0         0         0         0      1140         0         0         0               0           0
# none                                                        64558         0     23991     11598         0     29622     66832    106307     40971     33309           18285      147907
# RCT                                                         11536         0         0         0         0         0         0         0         0         0               0           0
# <NA>                                                            0       364         0         0         0         0         0         0         0         0               0           0

# save corrected metadata
write.csv(pdac_metadata, paste0(data_dir, "updated_metadata_samn.csv"), row.names = TRUE)
