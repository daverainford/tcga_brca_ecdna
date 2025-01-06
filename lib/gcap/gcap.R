suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(gcap))
suppressPackageStartupMessages(library(argparse))


# ================================
# Function Definitions
# ================================

# Function to parse ASCAT3 segment data from files
parse_ascat = function(ascat_dir) {
  log_info("Parsing ASCAT3 segment data.")
  files = list.files(path = ascat_dir, full.names = TRUE)
  
  # Read each file and append its data to a list
  ascat_df_list = lapply(files, function(file) {
    filename = basename(file)
    sample_name = strsplit(filename, "_")[[1]][1]
    df = read.table(file, header = TRUE)
    df$sample_name = sample_name
    return(df)
  })
  
  # Combine all dataframes into a single dataframe
  ascat_df = bind_rows(ascat_df_list)
  
  # Clean and select specific columns
  ascat_df_clean = data.frame(
    cbind(
      ascat_df$Chromosome, 
      ascat_df$Start, 
      ascat_df$End, 
      ascat_df$Copy_Number, 
      ascat_df$Minor_Copy_Number, 
      ascat_df$sample_name
    )
  )
  
  # Rename columns for consistency
  colnames(ascat_df_clean) = c("chromosome", "start", "end", "total_cn", "minor_cn", "sample")
  
  # Ensure columns are numeric
  ascat_df_clean$start <- as.numeric(ascat_df_clean$start)
  ascat_df_clean$end <- as.numeric(ascat_df_clean$end)
  ascat_df_clean$total_cn <- as.numeric(ascat_df_clean$total_cn)
  ascat_df_clean$minor_cn <- as.numeric(ascat_df_clean$minor_cn)
  
  return(ascat_df_clean)
}

# Function to run GCAP on the parsed ASCAT data
run_gcap = function(ascat_df) {
  log_info("Running GCAP workflow.")
  
  # Run the GCAP workflow
  gcap_result = gcap.ASCNworkflow(
    ascat_df,
    outdir = tempdir(), 
    tightness = 2L,
    model = "XGB11", 
    genome_build = "hg38"
  )
  
  # Extract sample-level and gene-level results
  gcap_sample_summary = gcap_result$sample_summary   
  gcap_gene_data = gcap_result$data                
  
  # Combine results into a list
  gcap_res_clean = list(
    gcap_res = gcap_sample_summary,
    gcap_genes = gcap_gene_data
  )
  
  return(gcap_res_clean)
}

# Function to clean GCAP results and write outputs
clean_gcap_res = function(gcap_res, output_dir) {
  log_info("Cleaning GCAP results and writing outputs.")
  
  # Extract GCAP results
  gcap_summary = gcap_res$gcap_res
  gcap_genes = gcap_res$gcap_genes
  
  # Map Ensembl gene IDs to gene symbols
  gcap_genes$gene_symbol <- mapIds(
    org.Hs.eg.db,
    keys = gcap_genes$gene_id,          
    column = "SYMBOL",          
    keytype = "ENSEMBL",       
    multiVals = "first"         
  )
  
  # Write GCAP results to output files
  write.csv(gcap_summary, file.path(output_dir, "gcap_result_summary.csv"), row.names = FALSE)
  write.csv(gcap_genes, file.path(output_dir, "gcap_gene-level_summary.csv"), row.names = FALSE)
  
  # Create and update model matrix
  model_matrix = data.frame(sample = gcap_summary$sample, ecdna_status = gcap_summary$class)
  model_matrix <- model_matrix %>%
    mutate(
      ecdna_status = case_when(
        ecdna_status == "circular" ~ "positive",     
        ecdna_status %in% c("noncircular", "nofocal") ~ "negative",  
        TRUE ~ ecdna_status                         
      )
    )
  
  # Return model_matrix if function called in another module
  return(model_matrix)
}

# ================================
# Main Execution Block
# ================================

main = function(ascat_dir, output_dir) {
  # Parse ASCAT3 segment data
  ascat_df = parse_ascat(ascat_dir)
  
  # Run GCAP workflow
  gcap_res = run_gcap(ascat_df)
  
  # Clean and write GCAP results
  model_matrix = clean_gcap_res(gcap_res, output_dir)
}

# Ensure the script is being executed as the main script
if (sys.nframe() == 0) {
  # Parse command-line arguments
  parser <- ArgumentParser()
  parser$add_argument("--ascat_dir", required=TRUE, help="Path to the directory containing ASCAT3 data.")
  parser$add_argument("--output_dir", required=TRUE, help="Path to the GCAP output directory.")
  args <- parser$parse_args()
  
  # Call the main function with command-line arguments
  main(args$ascat_dir, args$output_dir)
}
