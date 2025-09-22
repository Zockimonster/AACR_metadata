
# wget "https://zenodo.org/records/14199536/files/scAtlas.rds.gz"
# gunzip scAtlas.rds.gz

# https://aacrjournals.org/clincancerres/article/31/4/756/751743/Human-Pancreatic-Cancer-Single-Cell-Atlas-Reveals
# https://github.com/PDAC-MULTIOMIC/PDAC_Atlas

# Load required libraries 
required_packages <- c('Matrix', 'dplyr', 'Seurat', 'patchwork', 'SeuratData', 'ggplot2', 'SRAdb', 'GEOquery')

invisible(lapply(required_packages, library, character.only = T))

# load wp
# obtain directories
# ref_dir with a folder panc_ca_atlas, where the scAtlas file can be loaded

# obtain SRAdb sqlite file
# wget https://gbnci.cancer.gov/backup/SRAmetadb.sqlite.gz
# gunzip SRAmetadb.sqlite.gz

# load SRA Metadata info (for detailed code below visualizations)
#sra_dbname <- paste0(ref_dir, "panc_ca_atlas/SRAmetadb.sqlite")	
#sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)

# pdac data
pdac_ref <- readRDS(paste0(ref_dir, "panc_ca_atlas/scAtlas.rds"))
# to enable correction
pdac_ref@meta.data$DiseaseState <- as.character(pdac_ref@meta.data$DiseaseState)

dim(pdac_ref)
# [1]  36601 726107

##############################################################################################

# check treatment info
table(pdac_ref@meta.data$GSE.SRA..Study., pdac_ref@meta.data$Treatment, useNA = "ifany")
#                     N_A Treatment Treatment naïve Unknown
# EGAS00001002543      0     11536           64558       0 (?)
# GSE154778            0         0           20364       0 (?)
# GSE155698            0         0           23991       0 - tn
# GSE156405            0         0           17809       0 - X
# -> TO BE CORRECTED: 
# p1 had been treated with gemcitabine/abraxane, 
# all metastatic samples apart from liver metastasis had been treated with folfiri, evofosfamide/ipilimumab, and capecitabine
# GSE158356            0      1969             777       0 (?, no treatment info?)
# GSE194247            0         0           29622       0 - tn
# GSE202051            0     80811           66832       0 - ok
# GSE205013            0     56781          110585       0 - ok
# GSE211644            0         0            3611   37360 - X
# -> TO BE ADJUSTED: cd3+ positive cells & IP
# Treatment naive were mainly samples from neoadjuvant treated patients (m2 - m7, m1 unknown)
# Unknown, were infusion product enriched cd3+ positive cells
# GSE229413        33309         0               0       0 -healthy
# phs001840.v1.p1      0         0           18285       0 -tn
# PRJCA001063          0         0          147907       0 -tn

##############################################################################################

# check tissue type
table(pdac_ref@meta.data$GSE.SRA..Study., pdac_ref@meta.data$DiseaseState, useNA = "ifany")
#                  Donor Adjacent normal Primary tumor Metastatic lesion
# EGAS00001002543      0               0         76094                 0 -?
# GSE154778            0               0         11806              8558 - X & sample ID issues
# -> ambiguous sample ids hindering validation & correction
# GSE155698            0            4502         19489                 0 - srr id issue
# GSE156405            0               0         10011              7798 - primary <-> peritoneal metastasis
# GSE158356            0               0             0              2746 - LivMetD: omentum -> liver
# GSE194247            0               0         29622                 0 -> CD45 neg
# GSE202051            0           10660        136983                 0 - X
# snRNA, 4109 primary tumor cells were patient-derived organoid cultures
# GSE205013            0               0        130026             37340 - ok
# GSE211644            0               0         40971                 0 - X
# cd3 enriched cells & ex vivo cultured cd3+ enriched cells in MDA1 samples 
# (also treatment info needs to be adjusted: m1 unknown, m2-m7 treated)
# GSE229413        33309               0             0                 0 - ok
# phs001840.v1.p1      0            4181         14104                 0 - ok
# PRJCA001063          0           54692         93215                 0 - ok


##############################################################################################
# add conflict info
pdac_ref@meta.data$treatment_check <- NA
pdac_ref@meta.data$tissue_check <- NA
pdac_ref@meta.data$overall_check <- NA
##############################################################################################
# EGAS
# https://www.clinicaltrials.gov/study/NCT02750657
# inclusion criteria: no prior treatment
# primary or liver metastatic samples
# identification due to missing id links not possible
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "EGAS00001002543"] <- "?"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "EGAS00001002543"] <- "?"
##############################################################################################
# gse154778
# tissue info (Primary != Metastasis & Omentum Metastasis != Liver Metastasis)
df <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE154778_df.csv"))
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778"] <- df$tissue_match[match(gsub("_.*", "", pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778"]), df$sample_ID)] 
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778" & is.na(pdac_ref@meta.data$tissue_check)] <- "missing"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778"] <- "?"
##############################################################################################
# gse155698
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE155698"] <- "match"
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE155698"] <- "match"
##############################################################################################
# gse156405
# treatment (a couple of patients were pre-treated) & tissue (Primary != Peritoneal Metastasis & Peritoneal Metastasis != Primary)
# replace SRR12467643 by current name
incl_df <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE156405_incl_df.csv"))
pdac_ref@meta.data$Name[pdac_ref@meta.data$Name == "SRR12467643"] <- "SRR12467650"
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"] <- incl_df$tissue_match[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"], incl_df$SRR)] 
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"] <- incl_df$treatment_match[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"], incl_df$SRR)] 
##############################################################################################
# gse158356
# tissue info (LivMetD != Omentum Metastasis)
# no treatment info available?
df <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE158356_df.csv"))
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"] <- df$tissue_match[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"], gsub("_", "", df$sample_ID))] 
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"] <- "?"
##############################################################################################
# gse194247
# tissue info (all samples included cd45neg cells only)
df <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE194247_df.csv"))
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"] <- "conflict"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"] <- "match"
##############################################################################################
# gse202051
# tissue info (oragnoid cultures != Primary)
in_ref <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE202051_df.csv"))
organoid_cultures <- in_ref[grepl("^010nuc", in_ref$sample_id) | grepl("^010orgCRT", in_ref$sample_id),c("experiment", "run", "geo_id", "sample_id", "treatment")]
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE202051"]  <- "match"
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$Name %in% organoid_cultures$run] <- "conflict"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE202051"] <- "match"
##############################################################################################
# gse205013
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE205013"] <- "match"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE205013"] <- "match"
##############################################################################################
# GSE211644
# treatment (mainly pre-treated) & tissue (CD3+ enriched != Primary, Infusion Product from CD3+ enriched != Primary) 
incl_data <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE211644_incl_df.csv"))
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"] <- incl_data$tissue_match[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"], incl_data$SRR)]
# treatment info
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"] <- incl_data$treatment_match[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"], incl_data$SRR)]
##############################################################################################
# gse229413
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"] <- "match"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"] <- "match"
##############################################################################################
# phs
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "phs001840.v1.p1"] <- "match"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "phs001840.v1.p1"] <- "match"
##############################################################################################
# peng
# tissue
pdac_ref@meta.data$tissue_check[pdac_ref@meta.data$GSE.SRA..Study.== "PRJCA001063"] <- "match"
# treatment
pdac_ref@meta.data$treatment_check[pdac_ref@meta.data$GSE.SRA..Study.== "PRJCA001063"] <- "match"
##############################################################################################

# Add match status
pdac_ref@meta.data <- pdac_ref@meta.data %>%
  mutate(overall_check = case_when(
    grepl("GSE155698", GSE.SRA..Study.) ~ "conflict (ID)",  # SRR
    grepl("GSE154778", GSE.SRA..Study.) & grepl("missing", tissue_check)  ~ "conflict (ID)",  # sample_ID
    grepl("EGAS00001002543", GSE.SRA..Study.) ~ "?",
    grepl("conflict", treatment_check) | grepl("?", treatment_check, fixed = TRUE) 
    ~ "conflict (treatment)",
    grepl("conflict", tissue_check) 
    ~ "conflict (tissue)",
    is.na(tissue_check) | is.na(treatment_check)
    ~ "not_checked",
    TRUE ~ "match"))

write.csv(pdac_ref@meta.data, paste0(ref_dir,"panc_ca_atlas/pdac_ref_metadata.csv"))

# Load metadata
meta <- read.csv(
  paste0(ref_dir,"panc_ca_atlas/pdac_ref_metadata.csv"),
  row.names = 1)
# Create a minimal dummy counts matrix
dummy_counts <- matrix(0, nrow = 2, ncol = nrow(meta))
rownames(dummy_counts) <- c("Gene A", "Gene B")
colnames(dummy_counts) <- rownames(meta)
# Initialize a minimal Seurat object, to assign meta data to the related slot
pdac_ref <- CreateSeuratObject(dummy_counts, meta.data = meta)

# tissue
table(pdac_ref@meta.data$tissue_check, useNA="ifany")
# ?       conflict    match     missing
# 76094   77646       572003    364

prop.table(table(pdac_ref@meta.data$tissue_check, useNA="ifany"))
# ?             conflict        match         missing
# 0.1047972269  0.1069346529    0.7877668167  0.0005013035

# treatment
table(pdac_ref@meta.data$treatment_check, useNA="ifany")
# ?       conflict    match
# 99204   42892       584011

prop.table(table(pdac_ref@meta.data$treatment_check, useNA="ifany"))
# ?           conflict      match
# 0.13662449  0.05907118    0.80430432

# overall
table(pdac_ref@meta.data$overall_check, useNA="ifany")
# ?           conflict (ID)    conflict (tissue)  conflict (treatment)      match
# 76094       24355            38480              65638                     521540

prop.table(table(pdac_ref@meta.data$overall_check, useNA="ifany"))
# ?             conflict (ID)    conflict (tissue)    conflict (treatment)      match
# 0.10479723    0.03354189       0.05299494           0.09039715                0.71826880


# cross-table
table(pdac_ref@meta.data$tissue_check,pdac_ref@meta.data$treatment_check, useNA="ifany")
#               ? conflict  match
# ?         76094        0      0
# conflict    278    38888  38480
# match     22468     4004 545531
# missing     364        0      0

prop.table(table(pdac_ref@meta.data$tissue_check,pdac_ref@meta.data$treatment_check, useNA="ifany"))
#                     ?     conflict        match
# ?        0.1047972269 0.0000000000 0.0000000000
# conflict 0.0003828637 0.0535568449 0.0529949443
# match    0.0309430979 0.0055143388 0.7513093800
# missing  0.0005013035 0.0000000000 0.0000000000


##############################################################################################
##############################################################################################

# create plots
# define colors
col_list <- c(
  "?"                = "#A9A9A9",
  "missing"          = "#6F6A6A",    
  "match"            = "#2DC32D",
  "conflict (ID)"    = "#FFC7F3",
  "conflict (tissue)"= "#FF9BBF",
  "conflict (treatment)"= "#FC73FF",
  "conflict"         = "#C159D0")


# umaps
library(patchwork)
p1 <- DimPlot(pdac_ref, group.by = "Clusters", label = TRUE, label.size=3, repel = T, pt.size = 0.5) + ggtitle("Cell Type")
p2 <-  DimPlot(pdac_ref, group.by = "GSE.SRA..Study.",label = TRUE, label.size=3, repel = T, pt.size = 0.5) + ggtitle("Study ID")
p3 <- DimPlot(pdac_ref, group.by = "DiseaseState", label = TRUE, label.size=3, repel = T, pt.size = 0.5) + ggtitle("Tissue Type")
p4 <- DimPlot(pdac_ref, group.by = "Treatment", label = TRUE, label.size=3, repel = T, pt.size = 0.5)

p5 <- DimPlot(pdac_ref, group.by = "seurat_clusters", label = TRUE, label.size=3, repel = T, pt.size = 0.5) + ggtitle("Clusters")
p6 <- DimPlot(pdac_ref, group.by = "overall_check", label = TRUE, label.size=3, repel = T, pt.size = 0.5, cols=col_list) + ggtitle("Overall Check (Tissue, Treatment, ID)")
p7 <- DimPlot(pdac_ref, group.by = "tissue_check", label = TRUE, label.size=3, repel = T, pt.size = 0.5, cols=col_list) + ggtitle("Tissue Check (Tissue Type, Metastasis Location)")
p8 <-  DimPlot(pdac_ref, group.by = "treatment_check",label = TRUE, label.size=3, repel = T, pt.size = 0.5, cols=col_list) + ggtitle("Treatment Check (treated vs. untrated)")
plot_patched <- wrap_plots(p1,p2,p3,p4,p5,p6,p7,p8, ncol=4, byrow=T) + plot_annotation("UMAP for Pancreatic Cancer Atlas") + theme(plot.title = element_text(size = 18))
ggsave(paste0(ref_dir, "panc_ca_atlas/umap_check.png"), plot = plot_patched, width = 70, height =25, units = "cm")


# barplots
library(RColorBrewer)
pal_set <- colorRampPalette(brewer.pal(10,"Paired"))
png(paste0(ref_dir, "panc_ca_atlas/prop_celltypes_for_overall_check.png"), 
    width = 35, height = 20, units = "cm", res = 600)
par(mar=c(4, 10, 4, 2))
barplot(prop.table(table(pdac_ref@meta.data$Clusters,pdac_ref@meta.data$overall_check)), col=pal_set(14), legend.text=T, horiz=T, args.legend=list(x="bottomright", bty = "n"), las=1, main="Cell Type Proportions grouped by Overall Check")
dev.off()

prop.table(table(pdac_ref@meta.data$Clusters,gsub("conflict.*$", "conflict",pdac_ref@meta.data$overall_check)))
barplot(prop.table(table(gsub("conflict.*$", "conflict",pdac_ref@meta.data$overall_check), pdac_ref@meta.data$Clusters)), col=c("?"="#A9A9A9", "conflict"="#C159D0", "match"="#2DC32D"), legend.text=T, horiz=T, args.legend=list(x="bottomright", bty = "n"), las=1, main="Cell Type Proportions grouped by Overall Check")

barplot(prop.table(table(gsub("conflict.*$", "conflict",pdac_ref@meta.data$overall_check), pdac_ref@meta.data$GSE.SRA..Study.)), col=c("?"="#A9A9A9", "conflict"="#C159D0", "match"="#2DC32D"), legend.text=T, horiz=T, args.legend=list(x="bottomright", bty = "n"), las=1, main="Overall Check across Studies")
barplot(prop.table(table(pdac_ref@meta.data$GSE.SRA..Study.,gsub("conflict.*$", "conflict",pdac_ref@meta.data$overall_check))), col=pal_set(14), legend.text=T, horiz=T, args.legend=list(x="bottomright", bty = "n"), las=1, main="Study ID Proportions grouped by Overall Check")


# heatmap
library(pheatmap)
png(paste0(ref_dir, "panc_ca_atlas/prop_celltypes_per_study.png"), 
    width = 35, height = 20, units = "cm", res = 600)

tab <- prop.table(table(pdac_ref@meta.data$GSE.SRA..Study.,
                        pdac_ref@meta.data$overall_check), 1)
p1 <- pheatmap(tab, cluster_rows=TRUE, cluster_cols=FALSE,
         color=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         main="Overall Check per Study", legend_breaks = 0:1, legend = FALSE)[[4]]

tab <- prop.table(table(pdac_ref@meta.data$Clusters,
                        pdac_ref@meta.data$overall_check), 1)
p2 <- pheatmap(tab, cluster_rows=TRUE, cluster_cols=FALSE,
         color=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         main="Overall Check per Cell Type", legend_breaks = 0:1, legend=FALSE)[[4]]

tab <- prop.table(table(pdac_ref@meta.data$seurat_clusters,
                        pdac_ref@meta.data$overall_check), 1)
p3 <- pheatmap(tab, cluster_rows=TRUE, cluster_cols=FALSE,
         color=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
         main="Overall Check per integrated Leiden Clusters")[[4]]

tab <- prop.table(table(pdac_ref@meta.data$GSE.SRA..Study.,
                        pdac_ref@meta.data$Clusters), 1)
p4 <- pheatmap(tab, cluster_rows=TRUE, cluster_cols=FALSE,
               color=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
               main="Cell Type per Study", legend_breaks = 0:1, legend = FALSE)[[4]]

tab <- prop.table(table(pdac_ref@meta.data$GSE.SRA..Study.,
                        pdac_ref@meta.data$seurat_clusters), 1)
p5 <- pheatmap(tab, cluster_rows=TRUE, cluster_cols=FALSE,
               color=colorRampPalette(brewer.pal(9, "YlOrRd"))(100),
               main="Integrated Leiden Clusters per Study", legend_breaks = 0:1, legend = FALSE)[[4]]
layout <- "AABBCC
DDDEEE"
wrap_plots(p1,p2,p3, p4,p5, guides="collect", design = layout)  + plot_annotation("Heatmaps showing the proportions of:")


# sankey
library(ggalluvial)
# study <-> overall_check
df_alluv <- pdac_ref@meta.data %>%
  dplyr::count(GSE.SRA..Study., overall_check)

p1 <- ggplot(df_alluv,
       aes(axis1 = GSE.SRA..Study., axis2 =overall_check, y = n)) +
  geom_alluvium(aes(fill = GSE.SRA..Study.)) +
  geom_stratum(width = 1/8, fill = c(rev(pal_set(length(unique(df_alluv$GSE.SRA..Study.)))), rev(col_list[match(unique(df_alluv$overall_check),names(col_list))])), color = "black") +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), nudge_x = 0.2) +
  scale_fill_manual(values =pal_set(length(unique(df_alluv$GSE.SRA..Study.)))) +
  theme_minimal(base_size = 14) +
  labs(title = "Study ID and Overall Check",
       x = "", y = "Sample count")

# cell type <-> overall_check
df_alluv <- pdac_ref@meta.data %>%
  dplyr::count(Clusters, overall_check)

p2 <- ggplot(df_alluv,
       aes(axis1 = Clusters, axis2 =overall_check, y = n)) +
  geom_alluvium(aes(fill =Clusters)) +
  geom_stratum(width = 1/8, fill = c(rev(pal_set(length(unique(df_alluv$Clusters)))), rev(col_list[match(unique(df_alluv$overall_check),names(col_list))])), color = "black") +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), nudge_x = 0.2) +
  scale_fill_manual(values =pal_set(length(unique(df_alluv$Clusters)))) +
  theme_minimal(base_size = 14) +
  labs(title = "Cell Type and Overall Check",
       x = "", y = "Sample count")

p_patched <- wrap_plots(p1,p2,ncol=2)  + plot_annotation("Sankey Plots showing the Relation between")

# overall check across studies & cell types
pal_blue <- colorRampPalette(brewer.pal(5,"Blues"))
pal_red <- colorRampPalette(brewer.pal(5,"YlOrRd"))
df_alluv <- pdac_ref@meta.data %>%
  dplyr::count(GSE.SRA..Study., overall_check, Clusters)

# df_alluv$overall_check <- gsub("conflict.*$", "conflict",df_alluv$overall_check)

p <- ggplot(df_alluv,
             aes(axis1 = GSE.SRA..Study., axis2 =overall_check, axis3=Clusters, y = n)) +
  geom_alluvium(aes(fill =overall_check)) +
  geom_stratum(width = 1/8, fill = c(rev(pal_blue(length(unique(df_alluv$GSE.SRA..Study.)))), rev(col_list[match(unique(df_alluv$overall_check),names(col_list))]), rev(pal_red(length(unique(df_alluv$Clusters))))), color = "black") +
  geom_text_repel(stat = "stratum", aes(label = after_stat(stratum)), nudge_x = 0.2) +
  scale_fill_manual(values =col_list[match(unique(df_alluv$overall_check),names(col_list))], labels = c(paste0(names(table(pdac_ref@meta.data$overall_check))[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], ": n=",table(pdac_ref@meta.data$overall_check)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], " (", (round(prop.table(table(pdac_ref@meta.data$overall_check)), 3)*100)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], "%)"))) +
  labs(title = "Sankey Plots showing the Relation between Overall Check, Studies and Cell Types",
       x = "", y = "Sample count") + scale_x_continuous(labels = c("Study_ID", "overall_check", "Cell_Type"), breaks = c(1,2,3))

ggsave(paste0(ref_dir,"panc_ca_atlas/sankey_study_cell_type_overall_check.png"), plot=p, width = 40, height = 20, units = "cm")

# code to identify correct counts and percentages in legend
names(table(pdac_ref@meta.data$overall_check))[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))]
table(pdac_ref@meta.data$overall_check)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))]
(round(prop.table(table(pdac_ref@meta.data$overall_check)), 3)*100)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))]
paste0(names(table(pdac_ref@meta.data$overall_check))[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], ": n=",table(pdac_ref@meta.data$overall_check)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], " (", (round(prop.table(table(pdac_ref@meta.data$overall_check)), 3)*100)[match(unique(df_alluv$overall_check),names(table(pdac_ref@meta.data$overall_check)))], "%)")

##############################################################################################
##############################################################################################
##############################################################################################

# EGAS
table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "EGAS00001002543"])
# 100070  85948  87235  87784  91412  91610  91706  94930  95092  95373  96460
# 6637   9404   5673   4320   7845   1590   6498   2514   7881   4412   1767
# 97727  G9903
# 11536   6017

table(pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study. == "EGAS00001002543"])


table(pdac_ref@meta.data$Study..Citation..PMID.[pdac_ref@meta.data$GSE.SRA..Study. == "EGAS00001002543"])
# Chan-Seng-Yue, M., Kim, J.C., Wilson, G.W. et al. Transcription phenotypes of pancreatic cancer are driven by genomic events during tumor evolution. Nat Genet 52, 231–240 (2020). https://doi.org/10.1038/s41588-019-0566-9 PMID: 31932696
# 76094

##############################################################################################

# GSE154778
# tissue adjustments
# get sample names in pdac_ref
pid <- unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778"])
pid
# [1] "T1"                        "T10"
# [3] "T2"                        "T3"
# [5] "T4"                        "T5"
# [7] "T6"                        "T8"
# [9] "T9"                        "Y00006_HT77VBCXY_AACCGTAA"
# [11] "Y00013_HW7W7AFXX_TTGGCATA" "Y00014_HW7W7AFXX_TGAAGTAC"
# [13] "Y00016_HW7W7AFXX_GTTGCAGC" "Y00019_H5HVVBGX5_GATCTCAG"
# [15] "Y00027_H5HVVBGX5_ACGACATT"

tissue_pid <- pdac_ref@meta.data[pdac_ref@meta.data$Name %in% pid & pdac_ref@meta.data$GSE.SRA..Study.== "GSE154778", c("DiseaseState", "If.metastatic..location", "Name")]
str(tissue_pid)
# 'data.frame':   20364 obs. of  3 variables:
#   $ DiseaseState           : Factor w/ 4 levels "Donor","Adjacent normal",..: 3 3 3 3 3 3 3 3 3 3 ...
# $ If.metastatic..location: chr  "NA" "NA" "NA" "NA" ...
# $ Name                   : chr  "T1" "T1" "T1" "T1" ...

# simplify
pid <- gsub("_.*", "",pid)
pid
# [1] "T1"     "T10"    "T2"     "T3"     "T4"     "T5"     "T6"     "T8"
# [9] "T9"     "Y00006" "Y00013" "Y00014" "Y00016" "Y00019" "Y00027"

tissue_pid$simplified <- gsub("_.*", "",tissue_pid$Name)
tissue_pid$DiseaseState <- as.character(tissue_pid$DiseaseState)
colnames(tissue_pid)[colnames(tissue_pid) == "If.metastatic..location"] <- "metastatic_location"
str(tissue_pid)
# 'data.frame':   20364 obs. of  4 variables:
#   $ DiseaseState       : chr  "Primary tumor" "Primary tumor" "Primary tumor" "Primary tumor" ...
# $ metastatic_location: chr  "NA" "NA" "NA" "NA" ...
# $ Name               : chr  "T1" "T1" "T1" "T1" ...
# $ simplified         : chr  "T1" "T1" "T1" "T1" ...

# load gse info
geo_nr <- "GSE154778"
# gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154778/soft/GSE154778_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE154778_family.soft.gz"))

# Obtain GEO-Sample Characteristics
df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=6))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("tissue_type", "patient_ID", "sample_ID", "in_ref", "ref_tissue", "ref_metastatic_location")

sid <- c()
for (x in 1:length(names(gse@gsms))){
  disease_index <- grep("disease state: ", gse@gsms[[x]]@header$characteristics_ch1)
  disease_state <- gsub("disease state: ", "", gse@gsms[[x]]@header$characteristics_ch1[disease_index])
  disease_state <- ifelse(disease_state == "Pancreatic ductal adenocarcinoma", "PDAC", "different")
  
  tissue_index <- grep("tissue type: ", gse@gsms[[x]]@header$characteristics_ch1)
  tissue_type <- gsub("tissue type: ", "", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  if(tissue_type == "Primary"){
    tissue_type <- "Pancreas"
  }else if (tissue_type == "Metastasis"){
    tissue_type <- "Metastasis"
  }else{
    tissue_type <- "different"}
  
  tissue <- paste0(disease_state, "_", tissue_type)
  patient_id <- gsub("\\:.*", "", gse@gsms[[x]]@header$title) 
  
  sample_id <-gsub("^.*suppl/GSM.*_", "", gsub("_barcodes.*$", "", gse@gsms[[x]]@header$supplementary_file_1))
  sid <- c(sid, sample_id)
  tissue_in_pid <- unique(tissue_pid$DiseaseState[tissue_pid$simplified %in% sample_id])
  if(length(tissue_in_pid)==0){
    tissue_in_pid <- NA}
  met_loc_in_pid <- unique(tissue_pid$metastatic_location[tissue_pid$simplified %in% sample_id])
  if(length(met_loc_in_pid)==0){
    met_loc_in_pid <- NA}  
  df[x,]<- c(tissue, patient_id, sample_id, sample_id %in% pid, tissue_in_pid, met_loc_in_pid)}

# add tissue location
# https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00776-9
# Supplementary File Table S1
df <- df %>%
  mutate(tissue_location = case_when(
    patient_ID == "MET03" ~ "Omentum",
    grepl("PDAC_Pancreas", tissue_type)  ~ "Pancreas",
    grepl("PDAC_Metastasis", tissue_type) & grepl("MET03", tissue_type)==FALSE  ~ "Liver"))

# Add match status
df <- df %>%
  mutate(tissue_match = case_when(
    in_ref == "FALSE" ~ "missing",
    grepl("PDAC_Pancreas", tissue_type) & grepl("Primary tumor", ref_tissue) ~ "match",
    grepl("PDAC_Metastasis", tissue_type) & grepl("Metastatic lesion", ref_tissue)  & 
      tissue_location == ref_metastatic_location   ~ "match",
    TRUE ~ "conflict"))


df
#                tissue_type patient_ID sample_ID in_ref        ref_tissue ref_metastatic_location tissue_location tissue_match
# GSM4679532   PDAC_Pancreas        P01    K16733  FALSE              <NA>                    <NA>        Pancreas      missing
# GSM4679533   PDAC_Pancreas        P02    Y00006   TRUE Metastatic lesion                   Liver        Pancreas     conflict
# GSM4679534   PDAC_Pancreas        P03        T2   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679535   PDAC_Pancreas        P04        T3   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679536   PDAC_Pancreas        P05        T4   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679537   PDAC_Pancreas        P06        T5   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679538   PDAC_Pancreas        P07        T6   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679539   PDAC_Pancreas        P08        T8   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679540   PDAC_Pancreas        P09        T9   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679541   PDAC_Pancreas        P10       T10   TRUE     Primary tumor                      NA        Pancreas        match
# GSM4679542 PDAC_Metastasis      MET01    Y00008  FALSE              <NA>                    <NA>           Liver      missing
# GSM4679543 PDAC_Metastasis      MET02    Y00013   TRUE Metastatic lesion                   Liver           Liver        match
# GSM4679544 PDAC_Metastasis      MET03    Y00014   TRUE Metastatic lesion                   Liver         Omentum     conflict
# GSM4679545 PDAC_Metastasis      MET04    Y00016   TRUE Metastatic lesion                   Liver           Liver        match
# GSM4679546 PDAC_Metastasis      MET05    Y00019   TRUE Metastatic lesion                   Liver           Liver        match
# GSM4679547 PDAC_Metastasis      MET06    Y00027   TRUE Metastatic lesion                   Liver           Liver        match


sid
# [1] "K16733" "Y00006" "T2"     "T3"     "T4"     "T5"     "T6"     "T8"
# [9] "T9"     "T10"    "Y00008" "Y00013" "Y00014" "Y00016" "Y00019" "Y00027"

pid
# [1] "T1"     "T10"    "T2"     "T3"     "T4"     "T5"     "T6"     "T8"
# [9] "T9"     "Y00006" "Y00013" "Y00014" "Y00016" "Y00019" "Y00027"

setdiff(pid, sid)
# [1] "T1"
setdiff(sid, pid)
# [1] "K16733" "Y00008"

any(c(unlist(lapply("treat", grepl, tolower(gse@gsms))),  unlist(lapply("treat", grepl, tolower(gse@header))), unlist(lapply("treat", grepl, tolower(gse@gpls)))))
# [1] FALSE

write.csv(df, paste0(ref_dir, "panc_ca_atlas/GSE154778_df.csv"))

# add conflict info
df[df$tissue_match == "conflict",]
#                tissue_type patient_ID sample_ID in_ref        ref_tissue   ref_metastatic_location tissue_location tissue_match
# GSM4679533   PDAC_Pancreas        P02    Y00006   TRUE Metastatic lesion                     Liver        Pancreas     conflict
# GSM4679544 PDAC_Metastasis      MET03    Y00014   TRUE Metastatic lesion                     Liver         Omentum     conflict

pdac_ref@meta.data$DiseaseState[gsub("_.*", "", pdac_ref@meta.data$Name)  == "Y00006"] <- "conflict"
pdac_ref@meta.data$"If.metastatic..location"[gsub("_.*", "", pdac_ref@meta.data$Name)  == "Y00006"] <- "conflict"
pdac_ref@meta.data$"If.metastatic..location"[gsub("_.*", "", pdac_ref@meta.data$Name)  == "Y00014"] <- "conflict"


##############################################################################################

# GSE155698
# load gse info
geo_nr <- "GSE155698"
gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155698/soft/GSE155698_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE155698_family.soft.gz"))

length(gse@gsms)
# [1] 41

# Obtain GEO-Sample Characteristics
df <- data.frame(matrix(nrow = length(names(gse@gsms)), ncol=6))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("disease_state", "tissue_type", "patient_ID", "sample_ID", "treatment_info", "SAMN")
for (x in names(gse@gsms)){
  disease_state <- ifelse(grepl("PDAC", gse@gsms[[x]]@header$title), "PDAC", "Healthy")
  tissue_type <- ifelse(grepl("TISSUE", gse@gsms[[x]]@header$title), "Pancreas", "PBMC")
  patient_id <- paste0(substr(disease_state,1,1), gsub(".*_", "", gse@gsms[[x]]@header$title))
  sample_id <- gse@gsms[[x]]@header$title
  treatment_info <-ifelse(grepl("Treatment Naïve", gse@gsms[[x]]@header$characteristics_ch1), "Treatment Naïve", "Healthy Control")
  samn <- gsub("^.*SAMN", "SAMN", gse@gsms[[x]]@header$relation)
  df[x,] <- c(disease_state,tissue_type, patient_id, sample_id, treatment_info, samn)}

write.csv(df, paste0(ref_dir, "panc_ca_atlas/GSE155698_df.csv"))

# use srr ids to obtain samn ids
unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study == "GSE155698"])
# [1] "SRR12461530" "SRR12461533" "SRR12461536" "SRR12461538" "SRR12461543"
# [6] "SRR12461544" "SRR12461548" "SRR12507499" "SRR12507500" "SRR12507930"
# [11] "SRR12507931" "SRR12507934" "SRR12507935" "SRR12507938" "SRR12507939"
# [16] "SRR12508148" "SRR12532883" "SRR12532884" "SRR12532885" "SRR12532886"


pdac_srr_id <- unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study == "GSE155698"])

#  obtain corresponding samn  id & tissue information
get_samn_id <- function(srr) {
  exp_search <- entrez_search(db = "sra", term = srr)
  if (length(exp_search$ids) > 0) {
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession <- exp_summary$expxml
    if (grepl("<Summary><Title>PDAC_Tissue:", run_accession)){
      tissue_type <- "Primary Tumor"
    } else if (grepl("<Summary><Title>Normal_Tissue:", run_accession)){
      tissue_type <- "Adjacent Normal"}
    samn_id <- gsub("^.*SAMN","SAMN", run_accession)
    samn_id <- gsub("<.*$","", samn_id)
    return(list(tissue_type, samn_id))
  } else {
    return(NA)}}
srr_samn <- sapply(pdac_srr_id, get_samn_id)
write.csv(srr_samn, paste0(ref_dir, "panc_ca_atlas/srr_samn_GSE155698.csv"))

srr_samn <- read.csv(paste0(ref_dir, "panc_ca_atlas/srr_samn_GSE155698.csv"), row.names = 1)
srr_samn <- t(srr_samn)
colnames(srr_samn) <- c("tissue_type", "SAMN")

any(srr_samn$SAMN %in% df$SAMN)
# [1] FALSE

freq <- table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study == "GSE155698"], pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study == "GSE155698"])
freq
#             Adjacent normal Primary tumor
# SRR12461530               0           458
# SRR12461533               0          2070
# SRR12461536               0          1383
# SRR12461538               0           850
# SRR12461543            1155             0
# SRR12461544            2312             0
# SRR12461548               0           746
# SRR12507499             514             0
# SRR12507500             521             0
# SRR12507930               0           800
# SRR12507931               0           799
# SRR12507934               0           632
# SRR12507935               0           625
# SRR12507938               0          2280
# SRR12507939               0          2317
# SRR12508148               0          2037
# SRR12532883               0          1131
# SRR12532884               0          1126
# SRR12532885               0          1115
# SRR12532886               0          1120

# Add tissue type
df <- data.frame(matrix(nrow=nrow(freq), ncol=1))
rownames(df) <- rownames(freq)
colnames(df) <- "tissue_type"
for (row in 1:nrow(freq)){
  if (freq[row,"Adjacent normal"] == 0){
    df$tissue_type[row] <- "Primary Tumor"
  }else{
    df$tissue_type[row] <- "Adjacent Normal"}}

merge(srr_samn, df, by=0, suffixes=c("rentrez", "_in_ref"))
#      Row.names tissue_typerentrez         SAMN tissue_type_in_ref
# 1  SRR12461530      Primary Tumor SAMN15750836      Primary Tumor
# 2  SRR12461533      Primary Tumor SAMN15750849      Primary Tumor
# 3  SRR12461536      Primary Tumor SAMN15750842      Primary Tumor
# 4  SRR12461538      Primary Tumor SAMN15750838      Primary Tumor
# 5  SRR12461543    Adjacent Normal SAMN15750850    Adjacent Normal
# 6  SRR12461544    Adjacent Normal SAMN15750845    Adjacent Normal
# 7  SRR12461548      Primary Tumor SAMN15750826      Primary Tumor
# 8  SRR12507499    Adjacent Normal SAMN15750829    Adjacent Normal
# 9  SRR12507500    Adjacent Normal SAMN15750829    Adjacent Normal
# 10 SRR12507930      Primary Tumor SAMN15750841      Primary Tumor
# 11 SRR12507931      Primary Tumor SAMN15750841      Primary Tumor
# 12 SRR12507934      Primary Tumor SAMN15750815      Primary Tumor
# 13 SRR12507935      Primary Tumor SAMN15750815      Primary Tumor
# 14 SRR12507938      Primary Tumor SAMN15750840      Primary Tumor
# 15 SRR12507939      Primary Tumor SAMN15750840      Primary Tumor
# 16 SRR12508148      Primary Tumor SAMN15750842      Primary Tumor
# 17 SRR12532883      Primary Tumor SAMN15750832      Primary Tumor
# 18 SRR12532884      Primary Tumor SAMN15750832      Primary Tumor
# 19 SRR12532885      Primary Tumor SAMN15750832      Primary Tumor
# 20 SRR12532886      Primary Tumor SAMN15750832      Primary Tumor



##############################################################################################

# GSE156405
# -> TO BE CORRECTED: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8563410/
# Materials and Methods
# Sample Acquisition
# Five primary pancreatic cancers and four metastatic lesions in liver (LiM), lung (LuM), peritoneum (PM), and vaginal apex (VM) were used. 
# P2-P5 primary samples came from treatment-naïve patients, while P1 was a second FNA from a patient who had been treated with gemcitabine/paclitaxel. 
# For metastatic samples, the most recent therapies prior to sample acquisition were 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab and capecitabine for VM, LuM, and PM, respectively. 
# LiM was a treatment-naïve sample. Histologic confirmation of PDAC was performed by a pathologist.

# load gse info
geo_nr <- "GSE156405"
# gse <- getGEO(geo_nr, GSEMatrix = F)
# alternative on the cluster, if ssl refuses to connect
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156405/soft/GSE156405_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE156405_family.soft.gz"))

# P2–P5 primary samples came from treatment-naïve patients, while P1 was a second FNA from a patient who had been treated with gemcitabine/paclitaxel. 
# For metastatic samples, the most recent therapies prior to sample acquisition were 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab, and capecitabine for VM, LuM, and PM, respectively. 
# LiM was a treatment-naïve sample. Histologic confirmation of PDAC was performed by a pathologist.

table(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"])
# Treatment naïve
# 17809

table(pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"], pdac_ref@meta.data$If.metastatic..location[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"])
#                   Liver  Lung    NA Peritoneum Vaginal apex
# Metastatic lesion  4639  1225     0        459         1475
# Primary tumor         0     0 10011          0            0

table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"], pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"])
#             Metastatic lesion Primary tumor
# SRR12467642                 0          1304
# SRR12467643                 0          2207
# SRR12467644                 0          3975
# SRR12467645                 0          1961
# SRR12467646                 0           564
# SRR12467647              1475             0
# SRR12467648              4639             0
# SRR12467649              1225             0
# SRR24437450               459             0

table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"], pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"])
#             Treatment naïve
# SRR12467642            1304
# SRR12467643            2207 -> deprecated: replaced by SRR12467650, https://www.ncbi.nlm.nih.gov/sra/?term=SRR12467643
# SRR12467644            3975
# SRR12467645            1961
# SRR12467646             564
# SRR12467647            1475
# SRR12467648            4639
# SRR12467649            1225
# SRR24437450             459

# replace SRR12467643 by current name
pdac_ref@meta.data$Name[pdac_ref@meta.data$Name == "SRR12467643"] <- "SRR12467650"

#  obtain corresponding sra id
sra_id <- grep("SRA:", gse@header$relation)
sra_id <-  gsub("^.*=", "", gse@header$relation[sra_id])
sra_id
# [1] "SRP277925"

# Obtain Sample Info & SRX IDs
sample_id <- list()
sample_description <- list()
srx_id <- list()
for (x in 1:length(names(gse@gsms))){
  sample_id[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$title
  sample_description[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$description
  srx_index <- grep("term=SRX", gse@gsms[[x]]@header$relation)
  srx_id[[names(gse@gsms)[x]]] <-  gsub("^.*=", "", gse@gsms[[x]]@header$relation[srx_index])}

# obtain sra & srx ids for SRP ID
conversion <- sraConvert("SRP277925", sra_con = sra_con)

# obtain corresponding srr_id
srr_id <- srx_id
for (x in 1:length(names(gse@gsms))){
  index <- which(conversion$experiment == srr_id[[names(gse@gsms)[x]]])
  srr_id[[names(gse@gsms)[x]]] <-  conversion$run[index]}


# data in reference
incl_data <- srr_id[match(names(table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE156405"])), srr_id)]

incl_df <- data.frame(matrix(nrow=length(incl_data), ncol=4))
colnames(incl_df) <-c("SRR", "GEO", "sample_id", "sample_description")
for (x in 1:length(incl_data)){
  incl_df[x,] <- c(unname(incl_data[x]), names(incl_data)[x], unname(sample_id[names(incl_data)[x]]), unname(sample_description[names(incl_data)[x]]))}

# P2–P5 primary samples came from treatment-naïve patients, while P1 was a second FNA from a patient who had been treated with gemcitabine/paclitaxel. 
# For metastatic samples, the most recent therapies prior to sample acquisition were 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab, and capecitabine for VM, LuM, and PM, respectively. 
# LiM was a treatment-naïve sample. Histologic confirmation of PDAC was performed by a pathologist.

incl_df
#           SRR        GEO                   sample_id sample_description
# 1 SRR12467642 GSM4730260 Fine needle aspiration - P1                 P1 -> gemcitabine/paclitaxel
# 2 SRR12467644 GSM4730262 Fine needle aspiration - P3                 P3 
# 3 SRR12467645 GSM4730263 Fine needle aspiration - P4                 P4
# 4 SRR12467646 GSM4730264 Fine needle aspiration - P5                 P5
# 5 SRR12467647 GSM4730265          Vaginal metastasis                 VM -> 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab, and capecitabine
# 6 SRR12467648 GSM4730266            Liver metastasis                LiM 
# 7 SRR12467649 GSM4730267             Lung metastasis                LuM -> 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab, and capecitabine
# 8 SRR12467650 GSM4730268       Peritoneal metastasis                 PM -> 5-fluorouracil/liposmal irinotecan, evofosfamide/ipilimumab, and capecitabine
# 9 SRR24437450 GSM4730261 Fine needle aspiration - P2                 P2

# treatment
incl_df$pre_treatment <- NA
incl_df[1,"pre_treatment"] <- "Gemcitabine, Abraxane"
incl_df[c(2:4,6,9),"pre_treatment"] <- "Treatment naïve"
incl_df[c(5,7,8),"pre_treatment"] <- "FOLFIRI, Evofosfamide/Ipilimumab, Capecitabine"

incl_df$treatment_ref <- "Treatment naïve"

# Add match status
incl_df <- incl_df %>%
  mutate(treatment_match = case_when(
    grepl("Treatment naïve",pre_treatment ) & grepl("Treatment naïve", treatment_ref) ~ "match",
    TRUE ~ "conflict"))

# tissue
incl_df$tissue_ref <- pdac_ref@meta.data$DiseaseState[match(incl_df$SRR, pdac_ref@meta.data$Name)]
incl_df$met_loc <-  pdac_ref@meta.data$"If.metastatic..location"[match(incl_df$SRR, pdac_ref@meta.data$Name)]

# Add match status
incl_df <- incl_df %>%
  mutate(tissue_match = case_when(
    grepl("Fine needle aspiration", sample_id) & grepl("Primary tumor", tissue_ref) ~ "match",
    grepl("metastasis", sample_id) & grepl("Metastatic lesion", tissue_ref) & 
    grepl("Vaginal", sample_id) & grepl("Vaginal", met_loc) |
    grepl("Liver", sample_id) & grepl("Liver", met_loc) |
    grepl("Lung", sample_id) & grepl("Lung", met_loc) |
    grepl("Peritone", sample_id) & grepl("Peritone", met_loc) 
    ~ "match",
    TRUE ~ "conflict"))

write.csv(incl_df, paste0(ref_dir, "panc_ca_atlas/GSE156405_incl_df.csv"), row.names = F)

# adjust prior treatment info for GSE156405
table(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$Name %in% incl_df$SRR])
# Treatment naïve
# 17809

pdac_ref@meta.data$Treatment[pdac_ref@meta.data$Name %in% incl_df$SRR[incl_df$pre_treatment != "Treatment naïve"]] <- "Treatment"
table(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study. == "GSE156405"])
# Treatment Treatment naïve
# 6211           11598

# add pre-treatment info
table(pdac_ref@meta.data$TreatmentType[pdac_ref@meta.data$Name %in% incl_df$SRR], useNA="ifany")
#  <NA>
# 17809

incl_df <- incl_df[incl_df$pre_treatment != "Treatment naïve",]

all(pdac_ref@meta.data$Name[pdac_ref@meta.data$Name %in% incl_df$SRR] == incl_df$SRR[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$Name %in% incl_df$SRR], incl_df$SRR)])
# [1] TRUE

pdac_ref@meta.data$TreatmentType[pdac_ref@meta.data$Name %in% incl_df$SRR] <- incl_df$pre_treatment[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$Name %in% incl_df$SRR], incl_df$SRR)]
table(pdac_ref@meta.data$TreatmentType[pdac_ref@meta.data$Name %in% incl_df$SRR], useNA="ifany")
# FOLFIRI, Evofosfamide/Ipilimumab, Capecitabine
# 4907
# Gemcitabine, Abraxane
# 1304

# adjust tissue type
pdac_ref@meta.data[pdac_ref@meta.data$Name %in% incl_df$SRR[incl_df$tissue_ref == "Primary tumor" & incl_df$tissue_match == "conflict"], c("DiseaseState", "If.metastatic..location")] <- "conflict"
pdac_ref@meta.data[pdac_ref@meta.data$Name %in% incl_df$SRR[incl_df$tissue_ref == "Metastatic lesion" & incl_df$tissue_match == "conflict"], c("DiseaseState", "If.metastatic..location")] <- "conflict"


##############################################################################################

# GSE211644
# adjust treatment & tissue info
table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"], pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"])
#             Treatment naïve Unknown
# SRR18025515             392       0
# SRR18045262             701       0
# SRR18045263              82       0
# SRR18045264             473       0
# SRR18045265            1065       0
# SRR18045266             742       0
# SRR18045267              97       0
# SRR18045268              59       0
# SRR21152729               0    6226
# SRR21152730               0    6474
# SRR21152731               0    6842
# SRR21152732               0    6935
# SRR21152733               0    6593
# SRR21152734               0    4290

# load gse info
geo_nr <- "GSE211644"  # mda1_t2-t7 were treated before sample acquisition, no info available for mda_t1, ip = infusion product, ex vivo grown samples (with 37360 ex vivo cultured cells)
gse <- getGEO(geo_nr, GSEMatrix = F)
# alternative on the cluster, if ssl refuses to connect
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211644/soft/GSE211644_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE211644_family.soft.gz"))

gse@header$overall_design
# TCR enrichment for cohort MDA1
# [1] "we profiled 80,000 T cells from 57 PDAC, 22 uninvolved/normal samples, and cultured TIL using single-cell transcriptomic and T-cell receptor analysis."
# To further increase the data sample size, CD3+ TIL scRNA-seq data from another cohort at the institution (MDA2 data set; 26 primary PDAC and 10 uninvolved) 
# and from a previous study by Peng and colleagues (PUMCH data set; 24 PDAC and 11 normal pancreas) were combined (39)

# Obtain Sample Info & SRX IDs
sample_id <- list()
sample_description <- list()
srx_id <- list()
for (x in 1:length(names(gse@gsms))){
  sample_id[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$title
  sample_description[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$description
  srx_index <- grep("term=SRX", gse@gsms[[x]]@header$relation)
  srx_id[[names(gse@gsms)[x]]] <-  gsub("^.*=", "", gse@gsms[[x]]@header$relation[srx_index])}

# obtain corresponding sra number
library(rentrez)
get_srr_id <- function(srx) {
  exp_search <- entrez_search(db = "sra", term = srx, config = httr::config(connecttimeout = 100000))
  if (length(exp_search$ids) > 0) {
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession <- exp_summary$runs
    run_accession <- strsplit(run_accession, ",")[[1]]
    run_id <- gsub("^.*SRR","SRR", run_accession)
    run_id <- gsub("total_spots.*$","", run_id)
    run_id <- gsub("\" ", "", run_id)
    return(run_id)
  } else {
    return(NA)}}

srr_geo <- sapply(srx_id, get_srr_id)
write.csv(srr_geo, paste0(ref_dir, "panc_ca_atlas/GSE211644_geo_srr.csv"))
srr_geo <- read.csv( paste0(ref_dir, "panc_ca_atlas/GSE211644_geo_srr.csv"), col.names = c("GEO", "SRR"))

# data in reference
incl_data <- srr_geo[match(unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE211644"]), srr_geo$SRR),]
incl_data
# GEO         SRR
# 44 GSM6481948 SRR18025515
# 43 GSM6481947 SRR18045262
# 42 GSM6481946 SRR18045263
# 37 GSM6481941 SRR18045264
# 38 GSM6481942 SRR18045265
# 41 GSM6481945 SRR18045266
# 40 GSM6481944 SRR18045267
# 39 GSM6481943 SRR18045268
# 50 GSM6541654 SRR21152729
# 49 GSM6541653 SRR21152730
# 48 GSM6541652 SRR21152731
# 47 GSM6541651 SRR21152732
# 46 GSM6541650 SRR21152733
# 45 GSM6541649 SRR21152734


incl_data$sample_id <- unlist(unname(sample_id))[match(incl_data$GEO,names(sample_id))]
incl_data
#            SRR        GEO   sample_id
# 1  SRR18025515 GSM6481948    MDA1_T01 -> NA
# 2  SRR18045262 GSM6481947    MDA1_T07 -> T | FOLFIRINOX
# 3  SRR18045263 GSM6481946    MDA1_T06 -> T | Gemcitabine, Abraxane, Capecitabine
# 4  SRR18045264 GSM6481941    MDA1_U05 -> T | FOLFIRINOX
# 5  SRR18045265 GSM6481942    MDA1_T05 -> T | FOLFIRINOX
# 6  SRR18045266 GSM6481945    MDA1_T04 -> T | Gemcitabine
# 7  SRR18045267 GSM6481944    MDA1_T03 -> T | Gemcitabine, Abraxane
# 8  SRR18045268 GSM6481943    MDA1_T02 -> T | Pembrolizumab, Gemcitabine, Abraxane
# 9  SRR21152729 GSM6541654 MDA1_T07_IP -> T | FOLFIRINOX
# 10 SRR21152730 GSM6541653 MDA1_T05_IP -> T | FOLFIRINOX
# 11 SRR21152731 GSM6541652 MDA1_T04_IP -> T | Gemcitabine
# 12 SRR21152732 GSM6541651 MDA1_T03_IP -> T | Gemcitabine, Abraxane
# 13 SRR21152733 GSM6541650 MDA1_T02_IP -> T | Pembrolizumab, Gemcitabine, Abraxane
# 14 SRR21152734 GSM6541649 MDA1_T01_IP -> NA

# add & check treatment
incl_data$pre_treatment <- NA
incl_data[c(2,4,5,9,10),"pre_treatment"] <- "FOLFIRINOX"
incl_data[3,"pre_treatment"] <- "Gemcitabine, Abraxane, Capecitabine"
incl_data[c(6,11),"pre_treatment"] <- "Gemcitabine"
incl_data[c(7,12),"pre_treatment"] <- "Gemcitabine, Abraxane"
incl_data[c(8,13),"pre_treatment"] <- "Pembrolizumab, Gemcitabine, Abraxane"
incl_data[c(1,14),"pre_treatment"] <- "NA"

incl_data$treatment_ref <- pdac_ref@meta.data$Treatment[match(incl_data$SRR, pdac_ref@meta.data$Name)]

# Add match status
incl_data <- incl_data %>%
  mutate(treatment_match = case_when(
    grepl("NA",pre_treatment ) & grepl("Unknown", treatment_ref) ~ "match",
    TRUE ~ "conflict"))

# tissue
incl_data <- incl_data %>%
  mutate(tissue = case_when(
    grepl("_IP$", sample_id) ~ "Infusion Product, Tumor-derived",
    TRUE ~ "Primary tumor, CD3pos"))

incl_data$tissue_ref <- pdac_ref@meta.data$DiseaseState[match(incl_data$SRR, pdac_ref@meta.data$Name)]

# Add match status
incl_data <- incl_data %>%
  mutate(tissue_match = case_when(
    grepl("_IP$", sample_id) & grepl("Primary tumor", tissue_ref) ~ "conflict",
    grepl("CD3pos", tissue) & grepl("Primary tumor", tissue_ref) ~ "conflict",
    TRUE ~ "match"))
incl_data
#           GEO         SRR   sample_id                        pre_treatment   treatment_ref treatment_match                          tissue    tissue_ref tissue_match
# 44 GSM6481948 SRR18025515    MDA1_T01                                   NA Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 43 GSM6481947 SRR18045262    MDA1_T07                           FOLFIRINOX Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 42 GSM6481946 SRR18045263    MDA1_T06  Gemcitabine, Abraxane, Capecitabine Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 37 GSM6481941 SRR18045264    MDA1_U05                           FOLFIRINOX Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 38 GSM6481942 SRR18045265    MDA1_T05                           FOLFIRINOX Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 41 GSM6481945 SRR18045266    MDA1_T04                          Gemcitabine Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 40 GSM6481944 SRR18045267    MDA1_T03                Gemcitabine, Abraxane Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 39 GSM6481943 SRR18045268    MDA1_T02 Pembrolizumab, Gemcitabine, Abraxane Treatment naïve        conflict           Primary tumor, CD3pos Primary tumor     conflict
# 50 GSM6541654 SRR21152729 MDA1_T07_IP                           FOLFIRINOX         Unknown        conflict Infusion Product, Tumor-derived Primary tumor     conflict
# 49 GSM6541653 SRR21152730 MDA1_T05_IP                           FOLFIRINOX         Unknown        conflict Infusion Product, Tumor-derived Primary tumor     conflict
# 48 GSM6541652 SRR21152731 MDA1_T04_IP                          Gemcitabine         Unknown        conflict Infusion Product, Tumor-derived Primary tumor     conflict
# 47 GSM6541651 SRR21152732 MDA1_T03_IP                Gemcitabine, Abraxane         Unknown        conflict Infusion Product, Tumor-derived Primary tumor     conflict
# 46 GSM6541650 SRR21152733 MDA1_T02_IP Pembrolizumab, Gemcitabine, Abraxane         Unknown        conflict Infusion Product, Tumor-derived Primary tumor     conflict
# 45 GSM6541649 SRR21152734 MDA1_T01_IP                                   NA         Unknown           match Infusion Product, Tumor-derived Primary tumor     conflict

write.csv(incl_data, paste0(ref_dir, "panc_ca_atlas/GSE211644_incl_df.csv"), row.names = F)


# correct data
# tissue
pdac_ref@meta.data$DiseaseState[match(incl_data$SRR, pdac_ref@meta.data$Name)] <- incl_data$tissue
# treatment info
pdac_ref@meta.data$Treatment[match(incl_data$SRR, pdac_ref@meta.data$Name)] <-  ifelse(is.na(incl_data$pre_treatment), "Unknown", "Treatment")
pdac_ref@meta.data$TreatmentType[match(incl_data$SRR, pdac_ref@meta.data$Name)] <- incl_data$pre_treatment

##############################################################################################
##############################################################################################

# tissue adjustments for gse202051
# load gse info
geo_nr <- "GSE202051"  # change srr in geo_id, 4109 primary tumor samples were patient-derived organoid cultures
gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202051/soft/GSE202051_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE202051_family.soft.gz"))

length(gse@gsms)
# [1] 74

# Here we construct a high-resolution molecular landscape of the multicellular subtypes and spatial communities 
# that compose PAC using single-nucleus RNA-seq and whole-transcriptome digital spatial profiling (SP) of 
# 43 primary PDAC tumor specimens that either received neoadjuvant therapy or were treatment-naïve. 
# 18 specimens received no treatment prior to resection; 
# 25 specimens received neoadjuvant chemoradiation therapy (e.g. CRT = FOLFIRINOX), losartan (L), and/or nivolumab (N) prior to resection; 
# 2 specimens were non-malignant

# Obtain Sample ID & treatment info & corresponding srx number
sample_id <- list()
treatment_info <- list()
srx_id <- list()
for (x in 1:length(names(gse@gsms))){
  sample_id[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$title
  index <- grep("treatment:", gse@gsms[[x]]@header$characteristics_ch1)  
  treatment_info[[names(gse@gsms)[x]]] <-  gsub("treatment: ", "", gse@gsms[[x]]@header$characteristics_ch1[index])
  srx_index <- grep("term=SRX", gse@gsms[[x]]@header$relation)
  srx_id[[names(gse@gsms)[x]]] <-  gsub("^.*=", "", gse@gsms[[x]]@header$relation[srx_index])}

length(unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE202051"]))
# [1] 71

in_ref <- sraConvert(unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE202051"]), sra_con = sra_con)
str(in_ref)
# 'data.frame':   71 obs. of  5 variables:
#   $ run       : chr  "SRR19039304" "SRR19039305" "SRR19039306" "SRR19039308" ...
# $ submission: chr  "SRA1413638" "SRA1413638" "SRA1413638" "SRA1413638" ...
# $ study     : chr  "SRP373187" "SRP373187" "SRP373187" "SRP373187" ...
# $ sample    : chr  "SRS12855924" "SRS12855922" "SRS12855877" "SRS12855920" ...
# $ experiment: chr  "SRX15110829" "SRX15110828" "SRX15110782" "SRX15110826" ...

# add geo id, sample info & treatment info
in_ref$geo_id <- names(srx_id)[match(in_ref$experiment, unlist(unname(srx_id)))]
in_ref$sample_id <- unlist(unname(sample_id))[match(in_ref$geo_id, unlist(names(gse@gsms)))]
in_ref$treatment <- unlist(unname(treatment_info))[match(in_ref$geo_id, unlist(names(gse@gsms)))]
str(in_ref[,c("experiment", "run","submission", "geo_id", "sample_id", "treatment")])
# 'data.frame':   71 obs. of  6 variables:
#   $ experiment: chr  "SRX15110829" "SRX15110828" "SRX15110782" "SRX15110826" ...
# $ run       : chr  "SRR19039304" "SRR19039305" "SRR19039306" "SRR19039308" ...
# $ submission: chr  "SRA1413638" "SRA1413638" "SRA1413638" "SRA1413638" ...
# $ geo_id    : chr  "GSM6090487" "GSM6090486" "GSM6090426" "GSM6090484" ...
# $ sample_id : chr  "PDAC_T_25, replicate 1, snRNAseq" "PDAC_T_24, replicate 2, snRNAseq" "PDAC_U_2, replicate 1, snRNAseq" "PDAC_T_23, replicate 1, snRNAseq" ...
# $ treatment : chr  "Gy/cape" "Gem/abraxane, Gy/cape" "Untreated" "FOLFIRINOX, Gy/cisplatin" ...

# mark organoid cultures
organoid_cultures <- in_ref[grepl("^010nuc", in_ref$sample_id) | grepl("^010orgCRT", in_ref$sample_id),c("experiment", "run", "geo_id", "sample_id", "treatment")]
organoid_cultures
#     experiment         run     geo_id                        sample_id  treatment
# 48 SRX15110780 SRR19039354 GSM6090498    010nuc, replicate 1, snRNAseq  Untreated
# 49 SRX15110779 SRR19039355 GSM6090497 010orgCRT, replicate 1, snRNAseq        CRT 

pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$Name %in% organoid_cultures$run] <- "Primary tumor, Organoid Cultures"

# add & check treatment info from pdac_ref
table(pdac_ref@meta.data$Treatment[match(in_ref$run, pdac_ref@meta.data$Name)])
# Treatment Treatment naïve
# 41              30

in_ref$ref_treatment <- pdac_ref@meta.data$Treatment[match(in_ref$run, pdac_ref@meta.data$Name)]

table(in_ref$treatment, in_ref$ref_treatment)
#                          Treatment Treatment naïve
# CRT                             27               0
# CRTL                             7               0
# CRTLN                            2               0
# CRTN                             1               0
# FOLFIRINOX, Gy/cisplatin         2               0
# Gem/abraxane, Gy/cape            1               0
# Gy/cape                          1               0
# Untreated                        0              30

write.csv(in_ref, paste0(ref_dir, "panc_ca_atlas/GSE202051_df.csv"), row.names = FALSE)

##############################################################################################

# tissue adjustments for gse194247
# load gse info
geo_nr <- "GSE194247"
gse <- getGEO(geo_nr, GSEMatrix = F)
names(gse@gsms)
# [1] "GSM5831620" "GSM5831621" "GSM5831622" "GSM5831623" "GSM5831624"
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194247/soft/GSE194247_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE194247_family.soft.gz"))

length(gse@gsms)
# [1] 5

# Obtain GEO-Sample Characteristics
df <- data.frame(matrix(nrow = length(gse@gsms), ncol=5))
rownames(df) <- names(gse@gsms)
colnames(df) <- c("tissue", "cell type", "tissue type", "amplification round", "Treatment Protocol")
for (x in 1:length(names(gse@gsms))){
  sample_characteristics <- c(gse@gsms[[x]]@header$characteristics_ch1,gsub("\\]$","",gsub("^5_","", gsub("^.*\\[", "", gse@gsms[[x]]@header$title))), gse@gsms[[x]]@header$treatment_protocol_ch1)
  for (col in 1:ncol(df)){
    df[x,col] <- gsub(paste0(colnames(df)[col], ": "), "", sample_characteristics[col])}}  

df
#              tissue     cell type                    tissue type  amplification round  Treatment Protocol
# GSM5831620 pancreas CD45(-) cells pancreatic cancer tumor tissue                GEX_4           Untreated
# GSM5831621 pancreas CD45(-) cells pancreatic cancer tumor tissue                GEX_5           Untreated
# GSM5831622 pancreas CD45(-) cells pancreatic cancer tumor tissue                GEX_6           Untreated
# GSM5831623 pancreas CD45(-) cells pancreatic cancer tumor tissue                GEX_9           Untreated
# GSM5831624 pancreas CD45(-) cells pancreatic cancer tumor tissue            GEX_45_MM           Untreated


write.csv(df, paste0(ref_dir, "panc_ca_atlas/GSE194247_df.csv"))

cbind(table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"], pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"]), table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"], pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"]))
#                      Treatment naïve Primary tumor
# GEX_4                          8688          8688
# GEX_45_MM--ATTTGCTA            2127          2127
# GEX_5                         10034         10034
# GEX_6                          3710          3710
# GEX_9                          5063          5063

# adjust tissue type
pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE194247"] <- "Primary tumor, CD45neg"

#############################################################################################################

# gse158356
# correct LiverMetD: Omentum -> Liver
# Treatment Info???

# load gse info
geo_nr <- "GSE158356"  # change livermet names for geo_ids, no treatment info available
gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158356/soft/GSE158356_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE158356_family.soft.gz"))

length(gse@gsms)
# [1] 15

# This 10X Genomics single cell RNA sequencing repository includes raw and processed data files from biomaterial scaffolds from control (n=2) and tumor-bearingmice using an orthotopic model of pancreatic cancer (PDA) (n=1) and the iKras* p53* (n=1) genetically engineered model of PDA, mouse orthotopic pancreatic tumors (n=2), mouse normal pancreas samples (n=2), mouse iKras* p53* tumors (n=2), and human liver metastasis samples (n=5).

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
df <- data.frame(matrix(nrow = length(gse@gsms), ncol=8))
rownames(df) <- names(gse@gsms)
colnames(df) <- c("disease_state", "tissue", "sample_ID", "treatment", "tissue_ref", "met_loc_ref", "treatment_ref", "treatment_type_ref")
for (x in 1:length(names(gse@gsms))){
  disease_index <- grep("condition: ", gse@gsms[[x]]@header$characteristics_ch1)
  disease_state <- gsub("condition\\: pancreatic cancer \\(PDA\\)", "PDAC", gse@gsms[[x]]@header$characteristics_ch1[disease_index])
  
  tissue_index <- grep("tissue: ", gse@gsms[[x]]@header$characteristics_ch1)
  tissue_type <- gsub("tissue: human liver metastasis sample", "Liver Metastasis", gse@gsms[[x]]@header$characteristics_ch1[tissue_index])
  
  sample_id <- gse@gsms[[x]]@header$title
  
  gsub("_", "", sample_id)  
  tissue_in_pdac <- unique(pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$Name == gsub("_", "", sample_id)])
  if(length(tissue_in_pdac)==0){
    tissue_in_pdac <- NA}
  met_loc_in_pdac <- unique(pdac_ref@meta.data$"If.metastatic..location"[pdac_ref@meta.data$Name == gsub("_", "", sample_id)])
  if(length(met_loc_in_pdac)==0){
    met_loc_in_pdac <- NA}  
  treatment_in_pdac <- unique(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$Name == gsub("_", "", sample_id)])
  if(length(treatment_in_pdac)==0){
    treatment_in_pdac <- NA}  
  treatment_info_in_pdac <- unique(pdac_ref@meta.data$TreatmentType[pdac_ref@meta.data$Name == gsub("_", "", sample_id)])
  if(length(treatment_info_in_pdac)==0){
    treatment_info_in_pdac <- NA}  
  
  if (gse@gsms[[x]]@header$treatment_protocol_ch1 == "iKras* and iKras* p53*mice were administered doxycycline chow (BioServ, F3949) to induce expression of KrasG12D for 72 hours, followed by two days of 8 intraperitoneal injections of caerulein (Sigma, 75 ?g/kg) to induce pancreatitis. Doxycycline chow was continuously administered until experiment endpoint."){
    treatment_info <- "?"}
  
  df[x,]<- c(disease_state, tissue_type, sample_id, treatment_info, tissue_in_pdac, met_loc_in_pdac, treatment_in_pdac, treatment_info_in_pdac)}


cbind(table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"], pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"]),
      table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"], pdac_ref@meta.data$"If.metastatic..location"[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"]))
#           Metastatic lesion Liver Omentum
# LiverMetA              1252  1252       0
# LiverMetB               455   455       0
# LiverMetC               184   184       0
# LiverMetD               138     0     138
# LiverMetE               717   717       0


# Add match status
df <- df %>%
  mutate(tissue_match = case_when(
    grepl("Liver Metastasis", tissue) & grepl("Metastatic lesion", tissue_ref)  & 
      grepl("Liver", met_loc_ref) ~ "match",
    TRUE ~ "conflict"))

df
#              disease_state           tissue   sample_ID treatment          tissue_ref met_loc_ref   treatment_ref treatment_type_ref    tissue_match
#   GSM4798244          PDAC Liver Metastasis Liver_Met_A         ?   Metastatic lesion       Liver       Treatment               <NA>           match
#   GSM4798245          PDAC Liver Metastasis Liver_Met_B         ?   Metastatic lesion       Liver Treatment naïve               <NA>           match
#   GSM4798246          PDAC Liver Metastasis Liver_Met_C         ?   Metastatic lesion       Liver Treatment naïve               <NA>           match
#   GSM4798247          PDAC Liver Metastasis Liver_Met_D         ?   Metastatic lesion     Omentum Treatment naïve               <NA>        conflict
#   GSM4798248          PDAC Liver Metastasis Liver_Met_E         ?   Metastatic lesion       Liver       Treatment               <NA>           match

write.csv(df, paste0(ref_dir, "panc_ca_atlas/GSE158356_df.csv"))

pdac_ref@meta.data$"If.metastatic..location"[pdac_ref@meta.data$Name  == "LiverMetD"] <- "conflict"


# ???????????????????????????????????????????????????????????????????????????????????????
pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE158356"] <- "?"
##############################################################################################
##############################################################################################

# gse205013 correct
# load gse info
geo_nr <- "GSE205013"  # no changes required
gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205013/soft/GSE205013_family.soft.gz 
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE205013_family.soft.gz"))

length(gse@gsms)
# [1] 27

# Obtain GEO-Sample Characteristics
df <- data.frame(matrix(nrow = length(gse@gsms), ncol=2))
rownames(df) <- unlist(names(gse@gsms))
colnames(df) <- c("tissue", "treatment")
for (x in 1:length(gse@gsms)){
  sample_info <- gse@gsms[[x]]@header$characteristics_ch1
  for (col in 1:ncol(df)){
    df[x,col] <- gsub(paste0(colnames(df)[col], ": "), "", sample_info[col])}}

write.csv(df, paste0(ref_dir, "panc_ca_atlas/GSE205013_df.csv"))
df <- read.csv(paste0(ref_dir, "panc_ca_atlas/GSE205013_df.csv"), row.names = 1)

unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"])
# [1] "GSM6204109" "GSM6204110" "GSM6204111" "GSM6204112" "GSM6204113"
# [6] "GSM6204114" "GSM6204115" "GSM6204116" "GSM6204117" "GSM6204118"
# [11] "GSM6204119" "GSM6204120" "GSM6204121" "GSM6204122" "GSM6204123"
# [16] "GSM6204124" "GSM6204125" "GSM6204126" "GSM6204127" "GSM6204128"
# [21] "GSM6204129" "GSM6204130" "GSM6204131" "GSM6204132" "GSM6204133"
# [26] "GSM6204134" "GSM6204135"

# all samples from study in pdac_ref
length(unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"])) == nrow(df) &
  length(intersect(unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"]), rownames(df))) == nrow(df)
# [1] TRUE

# check treatment info
df$treatment[df$treatment == "Untreated"] <- "Treatment naïve"
df$treatment[df$treatment == "Treated"] <- "Treatment"

all(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"]== df$treatment[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"], rownames(df))])
# [1] TRUE

# check tissue info
df$tissue_type <- NA
df$tissue_type[df$tissue == "PDAC liver met"] <- "Metastatic lesion"
df$tissue_type[df$tissue == "Primary PDAC"] <- "Primary tumor"

all(pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"]== df$tissue_type[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"], rownames(df))])
# [1] TRUE

# check metastatic info
df$met_loc <- NA
df$met_loc[df$tissue == "PDAC liver met"] <- "Liver"
df$met_loc[df$tissue == "Primary PDAC"] <- "NA"

all(pdac_ref@meta.data$If.metastatic..location[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"]== df$met_loc[match(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE205013"], rownames(df))])
# [1] TRUE

#############################################################################################################
#############################################################################################################

# load gse info
geo_nr <- "GSE229413"
# gse <- getGEO(geo_nr, GSEMatrix = F)
# gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/soft/GSE229413_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE229413_family.soft.gz"))

length(gse@gsms)
# [1] 32

gse@header$summary
#  We compared this dataset of normal pancreata to single cell sequencing of tumor samples 
#  previously published in Steele, et al, Nature Cancer 2020, (raw fastq in dbGap phs002071) realigned to GRCh38 reference genome (CellRanger 6.0).
#  This GEO series contains the raw and filtered feature matrices for the donor pancreata as well as the realigned tumor samples.

# https://www.nature.com/articles/s43018-020-00121-4
# GSE155698
# already in PDAC Atlas reference -> including realigned samples from GSE229413 would result in duplicated samples
# therefore exclude realigned samples from GSE155698

unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"])
# [1] "GSM7162998_3829-EC"  "GSM7162999_3861-EC"  "GSM7163000_3862-EC"
# [4] "GSM7163001_4160-EC"  "GSM7163002_4161-EC"  "GSM7163003_4346-EC"
# [7] "GSM7163004_4347-EC"  "GSM7163007_4637-EC1" "GSM7163008_4637-EC2"
# [10] "GSM7163009_4741-EC1" "GSM7163010_4741-EC2"

# Obtain GEO-Sample Characteristics
sample_type <- list()
for (x in 1:length(names(gse@gsms))){
  sample_type[[names(gse@gsms)[x]]]<-  gse@gsms[[x]]@header$title}

# only healthy samples were included
sample_type[gsub("_.*$","", unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"]))]
# $GSM7162998
# [1] "Donor_1_head_and_tail"
# 
# $GSM7162999
# [1] "Donor_2_head"
# 
# $GSM7163000
# [1] "Donor_2_tail"
# 
# $GSM7163001
# [1] "Donor_4_head"
# 
# $GSM7163002
# [1] "Donor_4_tail"
# 
# $GSM7163003
# [1] "Donor_5_head"
# 
# $GSM7163004
# [1] "Donor_5_tail"
# 
# $GSM7163005
# [1] "Donor_6_head"
# 
# $GSM7163006
# [1] "Donor_6_tail"
# 
# $GSM7163007
# [1] "Donor_8_head"
# 
# $GSM7163008
# [1] "Donor_8_tail"
# 
# $GSM7163009
# [1] "Donor_9_head"
# 
# $GSM7163010
# [1] "Donor_9_tail"

unique(pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"])
# [1] "Donor"

unique(pdac_ref@meta.data$Treatment[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"])
# [1] "N_A"


#############################################################################################################

geo_nr <- "GSE229413"
# gse <- getGEO(geo_nr, GSEMatrix = F)
# gse <- getGEO(geo_nr, GSEMatrix = F)
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229413/soft/GSE229413_family.soft.gz
gse <- getGEO(filename = paste0(ref_dir, "/panc_ca_atlas/GSE229413_family.soft.gz"))

length(gse@gsms)
# [1] 32


unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study.== "GSE229413"])
# [1] "GSM7162998_3829-EC"  "GSM7162999_3861-EC"  "GSM7163000_3862-EC"
# [4] "GSM7163001_4160-EC"  "GSM7163002_4161-EC"  "GSM7163003_4346-EC"
# [7] "GSM7163004_4347-EC"  "GSM7163007_4637-EC1" "GSM7163008_4637-EC2"
# [10] "GSM7163009_4741-EC1" "GSM7163010_4741-EC2"


unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "GSE155698"])
# [1] "SRR12461530" "SRR12461533" "SRR12461536" "SRR12461538" "SRR12461543"
# [6] "SRR12461544" "SRR12461548" "SRR12507499" "SRR12507500" "SRR12507930"
# [11] "SRR12507931" "SRR12507934" "SRR12507935" "SRR12507938" "SRR12507939"
# [16] "SRR12508148" "SRR12532883" "SRR12532884" "SRR12532885" "SRR12532886"

# no related SAMN or SRR ID available

#############################################################################################################

unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "phs001840.v1.p1"])
# [1] "SRR9274536" "SRR9274537" "SRR9274538" "SRR9274539" "SRR9274540"
# [6] "SRR9274541" "SRR9274542" "SRR9274543" "SRR9274544"

srr_id <- unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "phs001840.v1.p1"])

#  obtain corresponding samn  id & tissue information
get_samn_id <- function(srr) {
  exp_search <- entrez_search(db = "sra", term = srr)
  if (length(exp_search$ids) > 0) {
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession <- exp_summary$expxml
    if (grepl("tissue adjacent to pancreatic ductal adenocarcinoma", run_accession)){
      tissue_type <- "Adjacent Normal"
    } else if (grepl("single-cell RNA-seq of human tissue: pancreatic ductal adenocarcinoma", run_accession)){
      tissue_type <- "Primary tumor"}
    samn_id <- gsub("^.*SAMN","SAMN", run_accession)
    samn_id <- gsub("</Biosample>  ","", samn_id)
    return(list(tissue_type, samn_id))
  } else {
    return(NA)}}

srr_id <- c("SRR9274536", "SRR9274537", "SRR9274538", "SRR9274539", "SRR9274540", "SRR9274541", "SRR9274542", "SRR9274543", "SRR9274544")
srr_samn <- sapply(srr_id, get_samn_id)
srr_samn
#          SRR9274536      SRR9274537      SRR9274538      SRR9274539      SRR9274540      SRR9274541        SRR9274542      SRR9274543        SRR9274544     
# [1,] "Primary tumor" "Primary tumor" "Primary tumor" "Primary tumor" "Primary tumor" "Adjacent Normal" "Primary tumor" "Adjacent Normal" "Primary tumor"
# [2,] "SAMN11970422"  "SAMN11970417"  "SAMN11970415"  "SAMN11970423"  "SAMN11970416"  "SAMN11970418"    "SAMN11970419"  "SAMN11970421"    "SAMN11970420" 

write.csv(srr_samn, paste0(ref_dir, "panc_ca_atlas/srr_samn_phs001840v1p1.csv"))

srr_samn <- read.csv(paste0(ref_dir, "panc_ca_atlas/srr_samn_phs001840v1p1.csv"), row.names = 1)
srr_samn <- t(srr_samn)
colnames(srr_samn) <- c("tissue_type", "SAMN")

freq <- table(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study == "phs001840.v1.p1"], pdac_ref@meta.data$DiseaseState[pdac_ref@meta.data$GSE.SRA..Study == "phs001840.v1.p1"])
freq
#             Adjacent normal Primary tumor
# SRR9274536               0          2229
# SRR9274537               0          5298
# SRR9274538               0          1150
# SRR9274539               0          1192
# SRR9274540               0           865
# SRR9274541             643             0
# SRR9274542               0          1677
# SRR9274543            3538             0
# SRR9274544               0          1693


# Add tissue type
df <- data.frame(matrix(nrow=nrow(freq), ncol=1))
rownames(df) <- rownames(freq)
colnames(df) <- "tissue_type"
for (row in 1:nrow(freq)){
  if (freq[row,"Adjacent normal"] == 0){
    df$tissue_type[row] <- "Primary Tumor"
  }else{
    df$tissue_type[row] <- "Adjacent Normal"}}

merge(srr_samn, df, by=0, suffixes=c("rentrez", "_in_ref"))
#    Row.names tissue_typerentrez         SAMN tissue_type_in_ref
# 1 SRR9274536      Primary tumor SAMN11970422      Primary Tumor
# 2 SRR9274537      Primary tumor SAMN11970417      Primary Tumor
# 3 SRR9274538      Primary tumor SAMN11970415      Primary Tumor
# 4 SRR9274539      Primary tumor SAMN11970423      Primary Tumor
# 5 SRR9274540      Primary tumor SAMN11970416      Primary Tumor
# 6 SRR9274541    Adjacent Normal SAMN11970418    Adjacent Normal
# 7 SRR9274542      Primary tumor SAMN11970419      Primary Tumor
# 8 SRR9274543    Adjacent Normal SAMN11970421    Adjacent Normal
# 9 SRR9274544      Primary tumor SAMN11970420      Primary Tumor


#############################################################################################################
#############################################################################################################

# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160
# metadata file including crr codes
peng_crr <- readxl::read_xlsx(paste0(ref_dir, "panc_ca_atlas/peng_CRR.xlsx"))
head(peng_crr)
# # A tibble: 6 × 5
# ID CRR       `Run title` `BioProject accession` `Experiment accession`
# <dbl> <chr>     <chr>       <chr>                  <chr>
# 1     1 CRR241805 T1          PRJCA001063            CRX030762
# 2     2 CRR241798 T2          PRJCA001063            CRX030763
# 3     3 CRR241799 T3          PRJCA001063            CRX030764
# 4     4 CRR034499 T4          PRJCA001063            CRX030765
# 5     5 CRR034500 T5          PRJCA001063            CRX030766
# 6     6 CRR034501 T6          PRJCA001063            CRX030767

unique(pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "PRJCA001063"])
# [1] "new_CRR034499" "new_CRR034500" "new_CRR034501" "new_CRR034503"
# [5] "new_CRR034504" "new_CRR034505" "new_CRR034506" "new_CRR034507"
# [9] "new_CRR034509" "new_CRR034511" "new_CRR034512" "new_CRR034513"
# [13] "new_CRR034516" "new_CRR034517" "new_CRR034518" "new_CRR034519"
# [17] "new_CRR034520" "new_CRR034521" "new_CRR034522" "new_CRR034523"
# [21] "new_CRR034524" "new_CRR034525" "new_CRR034526" "new_CRR034527"
# [25] "new_CRR034528" "new_CRR034529" "new_CRR034530" "new_CRR241798"
# [29] "new_CRR241799" "new_CRR241800" "new_CRR241801" "new_CRR241802"
# [33] "new_CRR241803" "new_CRR241804" "new_CRR241805"

table(peng_crr$"Run title"[match(gsub("new_", "", pdac_ref@meta.data$Name[pdac_ref@meta.data$GSE.SRA..Study. == "PRJCA001063"]), peng_crr$CRR)])
# N1  N10  N11   N2   N3   N4   N5   N6   N7   N8   N9   T1  T10  T11  T12  T13
# 7361 3474 6153 5269 2343 5235 3062 5181 4464 4502 7648 2140 1278 5097 3334 5423
# T14  T16  T17  T18  T19   T2  T20  T21  T22  T23  T24   T3   T4   T5   T6   T7
# 3271 3651 5998 5051 8333 4045 2034 2085 4148 4577 6370 2730 2215 2615 3650 2335
# T8   T9
# 2466 5840

peng_study <- pdac_ref@meta.data[pdac_ref@meta.data$GSE.SRA..Study. == "PRJCA001063", c("Name", "DiseaseState", "Treatment", "TreatmentType", "If.metastatic..location")]
str(peng_study)
# 'data.frame':   147907 obs. of  5 variables:
# $ Name                   : chr  "new_CRR034499" "new_CRR034499" "new_CRR034499" "new_CRR034499" ...
# $ DiseaseState           : chr  "Primary tumor" "Primary tumor" "Primary tumor" "Primary tumor" ...
# $ Treatment              : chr  "Treatment naïve" "Treatment naïve" "Treatment naïve" "Treatment naïve" ...
# $ TreatmentType          : chr  "NA" "NA" "NA" "NA" ...
# $ If.metastatic..location: chr  "NA" "NA" "NA" "NA" ...

# all treatment naïve
table(peng_study$Treatment, useNA = "ifany")
# Treatment naïve
# 147907

# all treatmenttype & metastatic locations == "NA"
all(peng_study$TreatmentType == "NA") && all(peng_study$"If.metastatic..location"=="NA")
# [1] TRUE

# remove "new_" in Names
peng_study$Name <- gsub("new_", "", peng_study$Name)

peng_study$sample_ID <- peng_crr$"Run title"[match(peng_study$Name, peng_crr$CRR)]

table(peng_study$DiseaseState, peng_study$sample_ID)
#                   N1  N10  N11   N2   N3   N4   N5   N6   N7   N8   N9   T1
# Adjacent normal 7361 3474 6153 5269 2343 5235 3062 5181 4464 4502 7648    0
# Primary tumor      0    0    0    0    0    0    0    0    0    0    0 2140
# 
#                  T10  T11  T12  T13  T14  T16  T17  T18  T19   T2  T20  T21
# Adjacent normal    0    0    0    0    0    0    0    0    0    0    0    0
# Primary tumor   1278 5097 3334 5423 3271 3651 5998 5051 8333 4045 2034 2085
# 
#                  T22  T23  T24   T3   T4   T5   T6   T7   T8   T9
# Adjacent normal    0    0    0    0    0    0    0    0    0    0
# Primary tumor   4148 4577 6370 2730 2215 2615 3650 2335 2466 5840

#############################################################################################################