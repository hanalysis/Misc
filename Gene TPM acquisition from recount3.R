library(recount3)
library(tidyverse)
library(RCurl)
library(readr)
library(R.utils)

# Read in project info

project_info <- read_csv("project_info.csv")


# Read in gene info 

gene_info <- read_csv("recount3_genes.csv")

# Loop over project info to get all the samples

  # Get number of projects

  project_info_SRA <- project_info %>%
    filter(file_source == "sra")
  
  x = length(project_info_SRA$project)

  # Gene IDs we want 
  gene_IDs <- gene_info$Ensemble_ID
  
    # KB for genes
    
    kb <- gene_info$KB
    
    error_list <- c()
    
  
  # Empty matrix for answers 
  gene_matrix <- matrix(,nrow = 8, ncol = 1)
  

  # Specify file path for back up CSV
  
  # Skipped proj number record
  
  skipped <- c()
  
  # Loop over each 
    
    # Going to 614 because this was missed on most recent re-run
    for (h in 1:x){
    
    # Get specific project information
    specific_proj <- project_info_SRA[h,]
    
      # Proj suffix
    
      suff <- str_sub(specific_proj[1],-2,-1)
    
    # Check if project URL exists
    URL <- paste("https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/sra/metadata/", suff, "/", specific_proj[1], "/sra.recount_pred.", specific_proj[1], ".MD.gz", sep = "")
    
    if(url.exists(URL) == FALSE) {
      
      skipped <- c(skipped, project_info[1])
      next()
      
    }
    
    # Make RSE for project
    tryCatch({
      withTimeout({
        rse <- create_rse(
        project_info = specific_proj, 
        type = "gene",
        annotation = "gencode_v29",
        recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release"
        )}, timeout = 360)
    }, TimeoutException = function(ex) {
      message(paste("Timeout, skipping": specific_proj[1]))
      next()
    })
    
    # Get all raw counts
    raw_counts <- as.data.frame((assays(rse)$raw_counts))
    
    # Get raw counts
    raw_count_OI <- raw_counts%>% 
      filter(row.names(raw_counts) %in% gene_IDs)
    
    gene_matrix <- cbind(gene_matrix, raw_count_OI)
    
    print(paste("Project", specific_proj[1], "completed"))
    
  
    }
  
  gene_matrix <- gene_matrix[, colnames(gene_matrix) != "gene_matrix"]
  
  write_csv(gene_matrix, "Gene_matrix_to_615.csv")
    
  # TPM
  
  total_counts <- rowSums(gene_matrix)
  
  pm_scaling_factor <- total_counts/ 1000000
  
  scaled_count_OI <- gene_matrix / pm_scaling_factor
  
  TPM <- scaled_count_OI / kb
  
  write.csv(TPM, "TPM_up_to_615.csv")
  
  write.table(skipped, "skipped_2.txt")
  