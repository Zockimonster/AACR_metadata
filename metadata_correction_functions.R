
# load functions, used to run code in metadata_correction_based_on_Name

#############################################################################################

#  obtain corresponding samn id & tissue information
# in gse155698
get_samn_and_info_geo155698 <- function(srr) {
  exp_search <- entrez_search(db = "sra", term = srr)
  if (length(exp_search$ids) > 0) {
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession <- exp_summary$expxml
    if (grepl("<Summary><Title>PDAC_Tissue:", run_accession)){
      sample_nr <- gsub("^  <Summary><Title>PDAC_Tissue: Sample ", "", run_accession)
      sample_nr <- gsub("</Title.*$", "", sample_nr)      
      tissue_type <- "Primary tumor"
    } else if (grepl("<Summary><Title>Normal_Tissue:", run_accession)){
      sample_nr <- gsub("^  <Summary><Title>Normal_Tissue: Sample ", "", run_accession)
      sample_nr <- gsub("</Title.*$", "", sample_nr)  
      tissue_type <- "Adjacent normal"}
    samn_id <- gsub("^.*SAMN","SAMN", run_accession)
    samn_id <- gsub("<.*$","", samn_id)
    return(list(tissue_type, samn_id, sample_nr))
  } else {
    return(NA)}}

# in phs
get_samn_and_info_phs <- function(srr) {
  exp_search <- entrez_search(db = "sra", term = srr)
  if (length(exp_search$ids) > 0) {
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession <- exp_summary$expxml
    if (grepl("tissue adjacent to pancreatic ductal adenocarcinoma", run_accession)){
      sample_name <- gsub("  <Summary><Title>single-cell RNA-seq of human tissue: tissue adjacent to pancreatic ductal adenocarcinoma: ","", run_accession)
      sample_name <- gsub("<.*","", sample_name)
      tissue_type <- "Adjacent normal"
    } else if (grepl("single-cell RNA-seq of human tissue: pancreatic ductal adenocarcinoma", run_accession)){
      sample_name <- gsub("  <Summary><Title>single-cell RNA-seq of human tissue: pancreatic ductal adenocarcinoma: ","", run_accession)
      sample_name <- gsub("<.*","", sample_name)
      tissue_type <- "Primary tumor"}
    samn_id <- gsub("^.*SAMN","SAMN", run_accession)
    samn_id <- gsub("</Biosample>  ","", samn_id)
    return(list(tissue_type, samn_id, sample_name))
  } else {
    return(NA)}}

# obtain corresponding srr number (run id) based on srx id
get_srr_id <- function(srx) {
  exp_search <- entrez_search(db = "sra", term = srx)
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

# get sra run info in general based on geo soft file
get_sra_run_info <- function(geo_nr, data_dir=NULL, use_soft_file=FALSE){
  if (use_soft_file == TRUE && is.null(data_dir)){
    stop("Please download the soft.gz file to the data_dir first.")}
  if (use_soft_file == TRUE && !is.null(data_dir)){
    # alternative on the cluster, if ssl refuses to connect
    # file_name <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", to_geo_nnn(geo_nr), "/", geo_nr, "/soft/", geo_nr, "_family.soft.gz")
    # wget file_name
    gse <- getGEO(filename = paste0(data_dir,geo_nr,"_family.soft.gz"))
  } else{
    gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)}
  #  obtain corresponding sra run info
  run_accession <- list()
  for (x in 1:length(names(gse@gsms))){
    srr_index <- grep("SRA:", gse@gsms[[x]]@header$relation)
    srr<-  gsub("^.*term\\=", "", gse@gsms[[x]]@header$relation[srr_index])
    exp_search <- entrez_search(db = "sra", term = srr)
    exp_summary <- entrez_summary(db = "sra", id = exp_search$ids)
    run_accession[[names(gse@gsms)[x]]] <- exp_summary$expxml}
  
  return(run_accession)}

#############################################################################################

to_geo_nnn <- function(x){
    x <- paste0(substr(x, 1, nchar(x)-3), "nnn")
   return(x)}

to_capital <- function(x){
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  return(x)}

convert_text_with_empty_lines_into_table <-  function(text, n_cols){
  # obtain individual entries, previously separated by an empty line
  lines <- str_split(text, "\\n+")[[1]] |> str_trim()
  # extract header
  header <- lines[1:n_cols]
  # get data lines
  data_lines <- lines[-(1:n_cols)]
  n_rows <- length(data_lines) / n_cols
  # reshape into data frame
  text_df <- matrix(data_lines, ncol = n_cols, byrow = TRUE) |>
    as.data.frame(stringsAsFactors = FALSE)
  colnames(text_df) <- header
  # clean up
  text_df <- text_df |> mutate(across(everything(), str_squish))
  return(text_df)}

#############################################################################################

update_metadata_from_df <- function(
    meta,               # pdac_metadata
    study_id,           # study_id used in pdac_metadata
    df,                 # df to use for adjustments
    match_meta_col,     # column name in pdac_meta used for matching
    match_df_col,       # column in df used for matching
    update_cols = NULL, # if NULL, update all overlapping columns
    clean_pattern = ""){
  
  # Clean matching keys
  meta_clean <- gsub(clean_pattern, "", meta[[match_meta_col]])
  df_clean <- gsub(clean_pattern, "", df[[match_df_col]])
  
  # Restrict to relevant study
  idx <- which(meta$Study_ID == study_id)
  match_idx <- match(meta_clean[idx], df_clean)
  matched <- !is.na(match_idx)
  replace_rows <- idx[matched]
  
  # Reporting
  cat("\nUpdating Study:", study_id, "\n")
  cat("✅ Matched", sum(matched), "of", length(idx), "samples\n")
  if (sum(!matched) > 0)
    cat("⚠️", sum(!matched), "samples had no match.\n")
  
  # Detect overlapping columns
  if (is.null(update_cols)) {
    update_cols <- intersect(names(df), names(meta))
    cat("Detected overlapping columns:", paste(update_cols, collapse=", "), "\n")}
  
  # Update specified columns
  for (col in update_cols) {
    if (col %in% names(df)) {
      new_values <- df[[col]][match_idx[matched]]
      meta[[col]][replace_rows] <- new_values
      cat("Updated column:", col, "\n")
    }else{
      cat("Column", col, "not found.\n")}}
  
  return(meta)}

