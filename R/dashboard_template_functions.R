# ---------------------------------------------------------------------------
# Title: dashboard_template_functions.R
# Description: This script contains helper functions used in 
# templates/dashboardTemplate.Rmd
# ---------------------------------------------------------------------------

#' This function gets the database to synapse id mapping table,
#' maps the provided database_name to its synapse id and returns it
#'
#' @param database_name (str) database name in database 
#' to synapse id mapping table
#' @param database_synid_mappingid (str) synapse id of the database 
#' to synapse id mapping table
#'
#' @return (str) synapse id of the mapped database name
get_syn_id_from_mapped_database <- function(database_name, database_synid_mappingid){
    database_synid_mapping = synTableQuery(sprintf('select * from %s',
                                               database_synid_mappingid))
    database_synid_mappingdf = as.data.frame(database_synid_mapping)
    table_synid = database_synid_mappingdf$Id[database_synid_mappingdf$Database == database_name]
    return(table_synid)
}

#' This function creates a table of failed annotation counts by grouped columns
#' @param maf_data (data.frame) input maf data frame
#' @param group_by_cols (str vector) list of columns to create counts by
#' @param counts_col_name (str) name to give to the counts column
#'
#' @return (data.frame) counts table
get_failed_annotation_table_counts <- function(maf_data, group_by_cols, counts_col_name){
    table_counts <- table(maf_data[maf_data$Annotation_Status == "FAILED", group_by_cols])

    if (nrow(table_counts) == 0){
        counts_table <- data.frame(matrix(ncol = length(group_by_cols) + 1, nrow = 0))
    } else{
        counts_table <- as.data.frame(table_counts)
    }
    colnames(counts_table) <- c(group_by_cols, counts_col_name)
    counts_table <- counts_table[do.call(order, counts_table[group_by_cols]), ]
    return(counts_table)
}