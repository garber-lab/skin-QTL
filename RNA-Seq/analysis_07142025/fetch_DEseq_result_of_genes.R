fetch_DEseq_result_of_genes <- function(ct, this_condition) {
  
  library(data.table)
  library(magrittr)
  
  metagene_dict <- data.table::fread("/pi/manuel.garber-umw/human/skin/eQTLs/literature/metaidname.txt") %>% 
    dplyr::select(-id) %>% distinct()
  metagene_dict_collapsed <- metagene_dict %>% group_by(meta) %>% 
    summarize(source=paste(name,collapse=','))
  
  convert_meta_genes <- function(gene_vector) {
    # Identify meta genes in the input vector
    meta_genes <- gene_vector[grepl("meta", gene_vector)]
    non_meta_genes <- gene_vector[!grepl("meta", gene_vector)]
    
    # Convert meta genes to actual gene names using metagene_dict
    meta_to_gene <- metagene_dict$name[metagene_dict$meta %in% meta_genes]
    
    # Combine non-meta genes and converted meta genes, ensuring uniqueness
    converted_genes <- unique(c(non_meta_genes, meta_to_gene))
    
    return(converted_genes)
  }
  
  load("/pi/manuel.garber-umw/human/skin/eQTLs/RNA-Seq/analysis_07142025/all_DESeq2_results_post_outlier_exclusion.RData")
  vars <- ls(pattern = "^res")
  
  target_var <- vars[grepl(tolower(ct), vars, ignore.case = TRUE) & 
                       grepl(tolower(this_condition), vars, ignore.case = TRUE)]
  
  deg_obj <- get(target_var) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    dplyr::mutate(celltype = ct,
                  condition = this_condition) %>%
    dplyr::select(c(gene, celltype, condition, log2FoldChange, padj))
  
  
  deg_dt <- as.data.table(deg_obj)
  meta_dt <- as.data.table(metagene_dict)
  
  deg_obj_metagene_converted <- merge(
    deg_dt,
    meta_dt,
    by.x = "gene",
    by.y = "meta",
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  deg_obj_metagene_converted[
    !is.na(name),
    gene := name
  ]
  
  deg_obj_metagene_converted[, name := NULL]
  
  return(deg_obj_metagene_converted)
  
}