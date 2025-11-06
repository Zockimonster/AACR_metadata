# AACR_metadata

Check of scAtlas metadata (https://doi.org/10.1158/1078-0432.CCR-24-2183, https://aacrjournals.org/clincancerres/article/31/4/756/751743/Human-Pancreatic-Cancer-Single-Cell-Atlas-Reveals)

aacr_scRNA_pancreas_rev.R: Coarse summary of conflicts found in the scAtlas metadata.

<img width="5314" height="2952" alt="sankey_study_cell_type_overall_check" src="https://github.com/user-attachments/assets/8993218e-0e66-48c9-96fc-9d0d22137587" />

metadata_correction_based_on_Name_column.R: Code used to check and correct the metadata based on the Name column. To run the code, the metadata_correction_functions.R script, the metadata file from https://ngdc.cncb.ac.cn/gsa/browse/CRA001160 (peng_CRR.xlsx) as well as sc_Atlas.rds are required. The metadata will be extracted in the code, and updated based on information from SRA Run Selector (rentrez, SRAdb), NCBI GEO soft files (GEOquery), https://ngdc.cncb.ac.cn/gsa/browse/CRA001160 (peng_CRR.xlsx), as well as information extracted from the related publications for studies used. Further information can be found as comment in the metadata_correction_based_on_Name_column.R script.

At the end of the code, some tables, summarizing the adjusted information across all studies, based on the corrected metadata are returned, and the corrected metadata can be saved separately in the data_dir as updated_metadata_samn.csv file, where it can be loaded using 

adjusted_metadata <- read.csv(paste0(data_dir, "updated_metadata_samn.csv"), row.names = 1)

