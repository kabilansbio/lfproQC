#' Best combination of normalization and imputation method
#'
#' @description This function will provide the best combinations of normalization
#' and imputation methods for the user given dataset based on the intragroup
#' variation evaluation parameters called PCV, PEV, PMAD and NRMSE.
#'
#' @param data_input Label-free proteomics expression data as a dataframe
#' @param groups Group information about the input data
#' @param data_type A character string specifying the type of data being used. Use "Peptide" if your dataset contains peptide information, where the first column represents peptide IDs and the second column represents protein IDs. Use "Protein" if your dataset consists of protein data, where the first column represents protein IDs. The function will handle the data accordingly based on this parameter.
#' @param aggr_method A character string specifying the method for aggregating peptide data to corresponding protein values. Use "sum" to aggregate by the sum of peptide intensities, "mean" to aggregate by the mean of peptide intensities, or "median" to aggregate by the median of peptide intensities. This parameter is only applicable when "data_type" is set to "Peptide".
#'
#' @details Label-free LC-MS proteomics expression data is often affected by heterogeneity and missing values. 
#' Normalization and missing value imputation are the commonly used techniques to solve these issues and make the dataset suitable for further downstream analysis. 
#' This function provides the best combination of normalization and imputation methods for the dataset, choosing from the three normalization methods (vsn, loess, and rlr) and three imputation methods (knn, lls, svd). 
#' The intragroup variation evaluation measures named pooled co-efficient of variance (PCV), pooled estimate of variance (PEV) and pooled median absolute deviation (PMAD) are used for selecting the best combination of normalization and imputation method for the given dataset.
#' It will return the best combinations based on each evaluation parameters of PCV, PEV, and PMAD. 
#'
#' Along with this, the user can get all three normalized datasets, nine combinations of normalized and missing values imputed datasets, and the PCV, PEV, and PMAD result values. The user can also obtain the Normalized Root Mean Square Error (NRMSE) values, 
#' calculated by comparing the normalized and imputed dataset with the original dataset across all nine combinations.
#' These NRMSE values provide insight into the accuracy of the imputation and normalization processes.
#'
#' @returns
#'  This function gives the list  which consist of following results.
#'
#' `Best Combinations`  The best combinations based on each PCV, PEV and PMAD
#' for the given dataset.
#'
#' `PCV Result` Values of groupwise PCV, overall PCV, PCV mean, PCV median and
#'              PCV standard deviation for all combinations.
#'
#' `PEV Result` Values of groupwise PEV, overall PEV, PEV mean, PEV median and
#'              PEV standard deviation for all combinations.
#'
#' `PMAD Result` Values of groupwise PMAD, overall PMAD, PMAD mean, PMAD median
#'               and PMAD standard deviation for all combinations.
#'
#' `NRMSE Result` NRMSE values calculated for the normalized and imputed dataset 
#'                to the original dataset. 
#'                
#' `vsn_data` The `vsn` normalized dataset
#'
#' `loess_data` The `loess` normalized dataset
#'
#' `rlr_data` The `rlr` normalized dataset
#'
#' `vsn_knn_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `knn` method.
#'              
#' `vsn_lls_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `lls` method.
#'              
#' `vsn_svd_data` The dataset normalized by `vsn` method and missing values imputed
#'              by `svd` method.
#'
#' `loess_knn_data` The dataset normalized by `loess` method and missing values imputed
#'              by `knn` method.
#'
#' `loess_lls_data` The dataset normalized by `loess` method and missing values imputed
#'              by `lls` method.
#'              
#' `loess_svd_data` The dataset normalized by `loess` method and missing values imputed
#'              by `svd` method.
#'              
#' `rlr_knn_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `knn` method.
#'
#' `rlr_lls_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `lls` method.
#'
#' `rlr_svd_data` The dataset normalized by `rlr` method and missing values imputed
#'              by `svd` method.
#'
#' @author
#'  Dr Sudhir Srivastava ("Sudhir.Srivastava@icar.gov.in")
#'
#'  Kabilan S ("kabilan151414@gmail.com")
#'
#' @export
#'
#' @examples
#' \donttest{
#' result <- best_combination(yeast_data, yeast_groups, data_type = "Protein")
#' result$`Best combinations`
#' result$`PCV Result`
#' result$`PMAD Result`
#' result$`rlr_knn_data`
#' }

#Main function for finding out the top three combinations
best_combination <- function (data_input, groups, data_type, aggr_method){
  
  sink("output.txt")
  
  Type <- Group <- name <- value <- . <- value_Mean <- median <-
    value_Median <- sd <- value_SD <- PCV_mean <- PCV_median <- PCV_sd <-
    PEV_mean <- PEV_median <- PEV_sd <- PMAD_mean <- PMAD_median <- PMAD_sd <-
    rowid <- original_order <- group <- mass <- all_na <- where <- NULL
  
  # Adding peptide to protein aggregation functionality
  aggregate_peptide_to_protein <- function(data, method) {
    # Ensure the data has at least three columns: Peptide ID, Protein ID, and expression values
    if (ncol(data) < 3) {
      stop("Peptide data must have at least three columns (Peptide ID, Protein ID, and expression values).")
    }
    
    peptide_column <- names(data)[1]  # First column is assumed to be Peptide ID
    protein_column <- names(data)[2]  # Second column is assumed to be Protein ID
    
    # Check if the aggregation method is valid
    if (!method %in% c("sum", "mean", "median")) {
      stop("Invalid aggregation method. Choose either 'sum', 'mean', or 'median'.")
    }
    
    # Convert column names to symbols for dplyr functions
    protein_col_sym <- dplyr::sym(protein_column)
    
    # Perform aggregation based on the Protein ID column
    aggregated_data <- data %>%
      dplyr::group_by(!!protein_col_sym) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), ~ match.fun(method)(., na.rm = TRUE)))  # Aggregate numeric columns
    
    # Rename the Protein ID column back to original after aggregation
    colnames(aggregated_data)[1] <- protein_column
    
    return(aggregated_data)
  }
  
  # If peptide data is chosen, perform aggregation, otherwise proceed with the protein data
  if (data_type == "Peptide" && !is.null(aggr_method)) {
    data_input <- aggregate_peptide_to_protein(data_input, aggr_method)
  }
  
  #Converting all zeros to NAs
  data_input[data_input == 0] <- NA
  
  complete_data_fn <- function(data, groups) {
    # Rename the columns in the groups data
    colnames(groups) <- c("name", "group")
    
    # Add a new column to preserve the original order of rows
    data <- data %>% dplyr::mutate(original_order = dplyr::row_number())
    
    # Rename the first column of data
    colnames(data)[1] <- "rowid"
    
    # Reshape the data into long format
    long_data <- data %>%
      tidyr::pivot_longer(-c(rowid, original_order), names_to = "name", values_to = "mass")
    
    # Merge with groups to add group information
    long_data <- long_data %>%
      dplyr::left_join(groups, by = "name")
    
    # Identify rows where any group has all missing values
    group_summary <- long_data %>%
      dplyr::group_by(rowid, group) %>%
      dplyr::summarise(all_na = all(is.na(mass)), .groups = 'drop') %>%
      dplyr::ungroup()
    
    # Identify rows to remove (where any group has all missing values)
    rows_to_remove <- group_summary %>%
      dplyr::group_by(rowid) %>%
      dplyr::summarise(remove = any(all_na)) %>%
      dplyr::filter(remove) %>%
      dplyr::pull(rowid) %>%
      unique()
    
    # Filter out rows with any completely missing group
    filtered_data <- long_data %>%
      dplyr::filter(!(rowid %in% rows_to_remove)) %>%
      dplyr::select(-group)
    
    # Reshape back to wide format and reorder based on the original order
    com_data <- filtered_data %>%
      tidyr::pivot_wider(names_from = name, values_from = mass) %>%
      dplyr::arrange(original_order) %>%
      dplyr::select(-original_order)
    
    return(com_data)
  } 
  
  com_data <- complete_data_fn (data_input, groups)
  
  #Giving original name to the first column
  com_data1 <- com_data
  colnames(com_data1)[1] <- colnames(data_input)[1]
  
  #Extracting the first column contains ID information
  com_data_ID <- com_data1[,1]
  
  #Removing the ID column and selecting remaining data
  com_data2 <- com_data1[,-1]
  
  #Grouping of dataframe according to sample groups
  grouping_data <- function(df, groups) {  # df = dataframe, groups = groups dataframe
    # Rename columns in groups dataframe for consistency
    colnames(groups) <- c("name", "group")
    
    # Create a list to store grouped columns
    grouped_data <- list()
    
    # Loop through each unique group
    for (group_name in unique(groups$group)) {
      # Get the names of the columns that belong to the current group
      group_columns <- groups$name[groups$group == group_name]
      
      # Extract the columns from df that match the current group columns
      grouped_data[[group_name]] <- df[, group_columns, drop = FALSE]
    }
    
    return(grouped_data)
  }
  
  # Function to replace NaN with NA in list elements
  replace_NaN_with_NA <- function(data_list) {
    lapply(data_list, function(df) {
      df[] <- lapply(df, function(col) {
        if (is.numeric(col)) {
          col[is.nan(col)] <- NA
        }
        return(col)
      })
      return(df)
    })
  }
  
  #Grouping of dataframe as a triplicate groups
  group_data <- grouping_data(com_data2, groups) 
  group_data <- replace_NaN_with_NA(group_data)
  
  #VSN Normalization function
  VSN_Norm <- function(dat) {
    dat<-as.data.frame(dat)
    vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
    colnames(vsnNormed) <- colnames(dat)
    row.names(vsnNormed) <- rownames(dat)
    return(as.matrix(vsnNormed))
  }
  
  #VSN normalized data
  vsn.dat <- do.call("cbind", lapply(group_data, VSN_Norm))
  
  #Loess normalization
  LOESS_Norm <- function(dat) {
    newdata<- as.data.frame(log2(dat))
    cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(newdata), method="fast")
    colnames(cycLoessNormed) <- colnames(newdata)
    row.names(cycLoessNormed) <- rownames(newdata)
    return(as.matrix(cycLoessNormed))
  }
  
  
  #Loess normalized data
  loess.dat <- do.call("cbind", lapply(group_data, LOESS_Norm))
  
  #RLR normalization
  RLR_Norm <- function (dat) {
    log2Matrix <- log2(dat)
    log2Matrix <- as.matrix(log2Matrix)
    
    # Extract the numeric component from the list
    log2Matrix <- unlist(log2Matrix)
    log2Matrix[is.infinite(log2Matrix)] <- 0
    
    sampleLog2Median <- matrixStats::rowMedians(log2Matrix, na.rm = TRUE)
    calculateRLMForCol <- function(colIndex, sampleLog2Median,
                                   log2Matrix) {
      lrFit <- MASS::rlm(as.matrix(log2Matrix[, colIndex]) ~
                           sampleLog2Median, na.action = stats::na.exclude, maxit = 20)
      coeffs <- lrFit$coefficients
      coefIntercept <- coeffs[1]
      coefSlope <- coeffs[2]
      globalFittedRLRCol <- (log2Matrix[, colIndex] - coefIntercept)/coefSlope
      globalFittedRLRCol
    }
    
    data <- log2Matrix
    globalFittedRLR <- vapply(seq_len(ncol(data)), calculateRLMForCol,
                              rep(0, nrow(data)), sampleLog2Median = sampleLog2Median,
                              log2Matrix = data)
    colnames(globalFittedRLR) <- colnames(dat)
    globalFittedRLR[globalFittedRLR < 0] <- 0
    globalFittedRLR
  }
  
  #RLR normalized data
  rlr.dat <- do.call("cbind", lapply(group_data, RLR_Norm))
  
  #Grouping of normalized datasets
  vsn_group_data <- grouping_data(vsn.dat, groups)
  
  loess_group_data <- grouping_data(loess.dat, groups)
  
  rlr_group_data <- grouping_data(rlr.dat, groups)
  
  #Imputation of normalized datasets
  #KNN imputation
  KNN_Imputation <- function (dat)
  {
    resultkNN <- VIM::kNN(dat, numFun = laeken::weightedMean, weightDist = TRUE,
                          imp_var = FALSE, k= 10)
    return(resultkNN)
  }
  
  #LLS imputation
  LLS_Imputation <- function (dat)
  {
    resultLLS <- pcaMethods::llsImpute(dat, k=2, correlation = "pearson", allVariables = TRUE)
    dataSet.imputed <- resultLLS@completeObs
    return(dataSet.imputed)
  }
  
  #SVD imputation
  SVD_Imputation <- function (dat)
  {
    resultSVD <- pcaMethods::pca(dat, method = "svdImpute", nPcs = 2)
    dataSet.imputed <- resultSVD@completeObs
    return(dataSet.imputed)
  }
  
  #VSN and KNN
  vsn.knn.dat <- do.call("cbind", lapply(vsn_group_data, KNN_Imputation))
  
  #VSN and LLS
  vsn.lls.dat <- do.call("cbind", lapply(vsn_group_data, LLS_Imputation))
  
  #VSN and SVD
  vsn.svd.dat <- do.call("cbind", lapply(vsn_group_data, SVD_Imputation))
  
  #LOESS and KNN
  loess.knn.dat <- do.call("cbind", lapply(loess_group_data, KNN_Imputation))
  
  #LOESS and LLS
  loess.lls.dat <- do.call("cbind", lapply(loess_group_data, LLS_Imputation))
  
  #LOESS and SVD
  loess.svd.dat <- do.call("cbind", lapply(loess_group_data, SVD_Imputation))
  
  #RLR and KNN
  rlr.knn.dat <- do.call("cbind", lapply(rlr_group_data, KNN_Imputation))
  
  #RLR and LLS
  rlr.lls.dat <- do.call("cbind", lapply(rlr_group_data, LLS_Imputation))
  
  #RLR and SVD
  rlr.svd.dat <- do.call("cbind", lapply(rlr_group_data, SVD_Imputation))
  
  #Transposing the sample to wide format
  new_sample <- as.data.frame(t(groups))
  names(new_sample) <- new_sample[1,]
  sample <- new_sample[-1,]
  
  #Changing all data files column names
  new_colnames <- colnames(sample)
  colnames(vsn.knn.dat) <- new_colnames
  colnames(vsn.lls.dat) <- new_colnames
  colnames(vsn.svd.dat) <- new_colnames
  colnames(loess.knn.dat) <- new_colnames
  colnames(loess.lls.dat) <- new_colnames
  colnames(loess.svd.dat) <- new_colnames
  colnames(rlr.knn.dat) <- new_colnames
  colnames(rlr.lls.dat) <- new_colnames
  colnames(rlr.svd.dat) <- new_colnames
  
  #Group wise PCV
  Group_data_PCV = function(data, groups){
    PCV = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      CVs = apply(data, 1, sd, na.rm = FALSE)/
        rowMeans(tempData, na.rm = FALSE)
      PCV_mean[group] = mean(CVs, na.rm = T)
      PCV_median[group] = median(CVs, na.rm = T)
      PCV_sd[group] = sd(CVs, na.rm = T)
    }
    Group_data_PCV_list = list("Group_data_PCV_mean" = as.data.frame(PCV_mean),
                               "Group_data_PCV_median" = as.data.frame(PCV_median),
                               "Group_data_PCV_sd" = as.data.frame(PCV_sd))
    return(Group_data_PCV_list)
  }
  
  #PCV calculation_groupwise
  test1 <- Group_data_PCV(vsn.knn.dat, sample)
  vsn_knn_PCV_mean <- cbind(Type ="vsn_knn", test1$Group_data_PCV_mean)
  vsn_knn_PCV_median <- cbind(Type ="vsn_knn", test1$Group_data_PCV_median)
  vsn_knn_PCV_sd <- cbind(Type ="vsn_knn", test1$Group_data_PCV_sd)
  
  test2 <- Group_data_PCV(vsn.lls.dat, sample)
  vsn_lls_PCV_mean <- cbind(Type ="vsn_lls", test2$Group_data_PCV_mean)
  vsn_lls_PCV_median <- cbind(Type ="vsn_lls", test2$Group_data_PCV_median)
  vsn_lls_PCV_sd <- cbind(Type ="vsn_lls", test2$Group_data_PCV_sd)
  
  test3 <- Group_data_PCV(vsn.svd.dat, sample)
  vsn_svd_PCV_mean <- cbind(Type ="vsn_svd", test3$Group_data_PCV_mean)
  vsn_svd_PCV_median <- cbind(Type ="vsn_svd", test3$Group_data_PCV_median)
  vsn_svd_PCV_sd <- cbind(Type ="vsn_svd", test3$Group_data_PCV_sd)
  
  test4 <- Group_data_PCV(loess.knn.dat, sample)
  loess_knn_PCV_mean <- cbind(Type ="loess_knn", test4$Group_data_PCV_mean)
  loess_knn_PCV_median <- cbind(Type ="loess_knn", test4$Group_data_PCV_median)
  loess_knn_PCV_sd <- cbind(Type ="loess_knn", test4$Group_data_PCV_sd)
  
  test5 <- Group_data_PCV(loess.lls.dat, sample)
  loess_lls_PCV_mean <- cbind(Type ="loess_lls", test5$Group_data_PCV_mean)
  loess_lls_PCV_median <- cbind(Type ="loess_lls", test5$Group_data_PCV_median)
  loess_lls_PCV_sd <- cbind(Type ="loess_lls", test5$Group_data_PCV_sd)
  
  test6 <- Group_data_PCV(loess.svd.dat, sample)
  loess_svd_PCV_mean <- cbind(Type ="loess_svd", test6$Group_data_PCV_mean)
  loess_svd_PCV_median <- cbind(Type ="loess_svd", test6$Group_data_PCV_median)
  loess_svd_PCV_sd <- cbind(Type ="loess_svd", test6$Group_data_PCV_sd)
  
  test7 <- Group_data_PCV(rlr.knn.dat, sample)
  rlr_knn_PCV_mean <- cbind(Type ="rlr_knn", test7$Group_data_PCV_mean)
  rlr_knn_PCV_median <- cbind(Type ="rlr_knn", test7$Group_data_PCV_median)
  rlr_knn_PCV_sd <- cbind(Type ="rlr_knn", test7$Group_data_PCV_sd)
  
  test8 <- Group_data_PCV(rlr.lls.dat, sample)
  rlr_lls_PCV_mean <- cbind(Type ="rlr_lls", test8$Group_data_PCV_mean)
  rlr_lls_PCV_median <- cbind(Type ="rlr_lls", test8$Group_data_PCV_median)
  rlr_lls_PCV_sd <- cbind(Type ="rlr_lls", test8$Group_data_PCV_sd)
  
  test9 <- Group_data_PCV(rlr.svd.dat, sample)
  rlr_svd_PCV_mean <- cbind(Type ="rlr_svd", test9$Group_data_PCV_mean)
  rlr_svd_PCV_median <- cbind(Type ="rlr_svd", test9$Group_data_PCV_median)
  rlr_svd_PCV_sd <- cbind(Type ="rlr_svd", test9$Group_data_PCV_sd)
  
  ###PCV_mean
  #Combining all the above results
  total_pcv_mean <- plyr::rbind.fill(vsn_knn_PCV_mean, vsn_lls_PCV_mean, vsn_svd_PCV_mean, 
                                     loess_knn_PCV_mean, loess_lls_PCV_mean, loess_svd_PCV_mean, 
                                     rlr_knn_PCV_mean, rlr_lls_PCV_mean, rlr_svd_PCV_mean)
  
  #Separating the results groupwise
  total_Group_data_PCV_mean <- total_pcv_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_mean)
  
  #Extract the top combination in each group
  total_Group_data_PCV_mean2<-
    total_Group_data_PCV_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  grouping_result <- function (data){
    result2 <- as.data.frame(data |>
                               dplyr::mutate(row = row_number()) |>
                               tidyr::pivot_longer(-row, values_transform = as.character) |>
                               dplyr::mutate(pair_num = (row_number() + 1) %/% 2, 
                                      type = if_else(row_number() %% 2 == 1, "val", "grp"), .by = row) |>
                               dplyr::select(-name) |>
                               tidyr::pivot_wider(names_from = type, values_from = value) |>
                               dplyr::summarize(vals = paste0(val, collapse = ", "),
                                         .by = c(pair_num, grp)) |>
                               dplyr::mutate(row = row_number(), .by = pair_num) |>
                               tidyr::pivot_wider(names_from = pair_num, values_from = c(vals, grp), names_vary = "slowest") |>
                               dplyr::select(-row) |>
                               `colnames<-`(colnames(data)))
    
    result <- result2[1,]
    return(result)
  } 
  
  final_Group_data_PCV_mean1 <- subset(total_Group_data_PCV_mean2, select = -row)
  
  final_Group_data_PCV_mean <- grouping_result(final_Group_data_PCV_mean1)
  
  ###PCV_median
  #Combining all the above results
  total_pcv_median <- plyr::rbind.fill(vsn_knn_PCV_median, vsn_lls_PCV_median, vsn_svd_PCV_median, 
                                       loess_knn_PCV_median, loess_lls_PCV_median, loess_svd_PCV_median, 
                                       rlr_knn_PCV_median, rlr_lls_PCV_median, rlr_svd_PCV_median)
  
  #Separating the results groupwise
  total_Group_data_PCV_median <- total_pcv_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_median)
  
  #Extract the top combination in each group
  total_Group_data_PCV_median2<-
    total_Group_data_PCV_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result 
  final_Group_data_PCV_median1 <- subset(total_Group_data_PCV_median2, select = -row)
  
  final_Group_data_PCV_median <- grouping_result (final_Group_data_PCV_median1)
  
  ###PCV_sd
  #Combining all the above results
  total_pcv_sd <- plyr::rbind.fill(vsn_knn_PCV_sd, vsn_lls_PCV_sd, vsn_svd_PCV_sd, 
                                   loess_knn_PCV_sd, loess_lls_PCV_sd, loess_svd_PCV_sd, 
                                   rlr_knn_PCV_sd, rlr_lls_PCV_sd, rlr_svd_PCV_sd)
  
  #Separating the results groupwise
  total_Group_data_PCV_sd <- total_pcv_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PCV_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PCV_sd)
  
  #Extract the top combination in each group
  total_Group_data_PCV_sd2<-
    total_Group_data_PCV_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PCV_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PCV_sd1 <- subset(total_Group_data_PCV_sd2, select = -row)
  
  final_Group_data_PCV_sd <- grouping_result (final_Group_data_PCV_sd1)
  
  #Overall data PCV
  Total_data_PCV = function(data){
    CVs = apply(data, 1, sd, na.rm = FALSE)/
      rowMeans(data, na.rm = FALSE)
    PCV_mean = mean(CVs, na.rm = T)
    PCV_median =  median(CVs, na.rm = T)
    PCV_sd = sd(CVs, na.rm = T)
    PCV_list = list("PCV_mean" = as.data.frame(PCV_mean),
                    "PCV_median" = as.data.frame(PCV_median),
                    "PCV_sd" = as.data.frame(PCV_sd))
    return(PCV_list)
  }
  
  
  #Overall PCV_mean estimation
  test1 <- Total_data_PCV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PCV_mean)
  
  test2 <- Total_data_PCV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PCV_mean)
  
  test3 <- Total_data_PCV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PCV_mean)
  
  test4 <- Total_data_PCV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PCV_mean)
  
  test5 <- Total_data_PCV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PCV_mean)
  
  test6 <- Total_data_PCV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PCV_mean)
  
  test7 <- Total_data_PCV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PCV_mean)
  
  test8 <- Total_data_PCV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PCV_mean)
  
  test9 <- Total_data_PCV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PCV_mean)
  
  #Combining all the above results
  total_pcv_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_mean1 <-total_pcv_overall_mean2%>%dplyr::slice_min(PCV_mean, n=1, with_ties = TRUE)
  total_pcv_overall_mean <- grouping_result(total_pcv_overall_mean1)
  original_cols <- c("Overall_Type.PCV_mean", "Overall_value.PCV_mean")
  colnames(total_pcv_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PCV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PCV_median)
  
  test2 <- Total_data_PCV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PCV_median)
  
  test3 <- Total_data_PCV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PCV_median)
  
  test4 <- Total_data_PCV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PCV_median)
  
  test5 <- Total_data_PCV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PCV_median)
  
  test6 <- Total_data_PCV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PCV_median)
  
  test7 <- Total_data_PCV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PCV_median)
  
  test8 <- Total_data_PCV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PCV_median)
  
  test9 <- Total_data_PCV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PCV_median)
  
  #Combining all the above results
  total_pcv_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_median1 <-total_pcv_overall_median2%>%dplyr::slice_min(PCV_median, n=1, with_ties = TRUE)
  total_pcv_overall_median <- grouping_result(total_pcv_overall_median1)
  original_cols <- c("Overall_Type.PCV_median", "Overall_value.PCV_median")
  colnames(total_pcv_overall_median) <- original_cols
  
  #Overall PCV_sd estimation
  test1 <- Total_data_PCV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PCV_sd)
  
  test2 <- Total_data_PCV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PCV_sd)
  
  test3 <- Total_data_PCV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PCV_sd)
  
  test4 <- Total_data_PCV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PCV_sd)
  
  test5 <- Total_data_PCV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PCV_sd)
  
  test6 <- Total_data_PCV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PCV_sd)
  
  test7 <- Total_data_PCV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PCV_sd)
  
  test8 <- Total_data_PCV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PCV_sd)
  
  test9 <- Total_data_PCV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PCV_sd)
  
  #Combining all the above results
  total_pcv_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pcv_overall_sd1 <-total_pcv_overall_sd2%>%dplyr::slice_min(PCV_sd, n=1, with_ties = TRUE)
  total_pcv_overall_sd <- grouping_result(total_pcv_overall_sd1)
  original_cols <- c("Overall_Type.PCV_sd", "Overall_value.PCV_sd")
  colnames(total_pcv_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PCV<- cbind(final_Group_data_PCV_mean, final_Group_data_PCV_median,
                     final_Group_data_PCV_sd, total_pcv_overall_mean,
                     total_pcv_overall_median, total_pcv_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PCV), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PCV_names1 <- result_PCV[, -col_indices]
  n <- matrix(t(result_PCV_names1), ncol=1)
  # Split the values by comma and remove extra spaces
  split_data <- strsplit(n, ",\\s*")
  
  # Flatten the list into a single vector
  separated_data <- unlist(split_data)
  
  # Get the frequency of each element in the dataframe
  freq <- table(separated_data)
  
  # Find the most occurring element
  PCV_best_combination <- names(freq)[which.max(freq)]
  
  #Groupwise PEV estimation
  Group_data_PEV = function(data1, groups){
    data <- as.matrix(data1)
    PEV = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      rowNonNACnt = rowSums(!is.na(tempData)) - 1
      EV = rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
      PEV_mean[group] = mean(EV, na.rm = T)
      PEV_median[group] = median(EV, na.rm = T)
      PEV_sd[group] = sd(EV, na.rm = T)
    }
    Group_data_PEV_list = list("Group_data_PEV_mean" = as.data.frame(PEV_mean),
                               "Group_data_PEV_median" = as.data.frame(PEV_median),
                               "Group_data_PEV_sd" = as.data.frame(PEV_sd))
    return(Group_data_PEV_list)
  }
  
  
  #PEV calculation_groupwise
  test1 <- Group_data_PEV(vsn.knn.dat, sample)
  vsn_knn_PEV_mean <- cbind(Type ="vsn_knn", test1$Group_data_PEV_mean)
  vsn_knn_PEV_median <- cbind(Type ="vsn_knn", test1$Group_data_PEV_median)
  vsn_knn_PEV_sd <- cbind(Type ="vsn_knn", test1$Group_data_PEV_sd)
  
  test2 <- Group_data_PEV(vsn.lls.dat, sample)
  vsn_lls_PEV_mean <- cbind(Type ="vsn_lls", test2$Group_data_PEV_mean)
  vsn_lls_PEV_median <- cbind(Type ="vsn_lls", test2$Group_data_PEV_median)
  vsn_lls_PEV_sd <- cbind(Type ="vsn_lls", test2$Group_data_PEV_sd)
  
  test3 <- Group_data_PEV(vsn.svd.dat, sample)
  vsn_svd_PEV_mean <- cbind(Type ="vsn_svd", test3$Group_data_PEV_mean)
  vsn_svd_PEV_median <- cbind(Type ="vsn_svd", test3$Group_data_PEV_median)
  vsn_svd_PEV_sd <- cbind(Type ="vsn_svd", test3$Group_data_PEV_sd)
  
  test4 <- Group_data_PEV(loess.knn.dat, sample)
  loess_knn_PEV_mean <- cbind(Type ="loess_knn", test4$Group_data_PEV_mean)
  loess_knn_PEV_median <- cbind(Type ="loess_knn", test4$Group_data_PEV_median)
  loess_knn_PEV_sd <- cbind(Type ="loess_knn", test4$Group_data_PEV_sd)
  
  test5 <- Group_data_PEV(loess.lls.dat, sample)
  loess_lls_PEV_mean <- cbind(Type ="loess_lls", test5$Group_data_PEV_mean)
  loess_lls_PEV_median <- cbind(Type ="loess_lls", test5$Group_data_PEV_median)
  loess_lls_PEV_sd <- cbind(Type ="loess_lls", test5$Group_data_PEV_sd)
  
  test6 <- Group_data_PEV(loess.svd.dat, sample)
  loess_svd_PEV_mean <- cbind(Type ="loess_svd", test6$Group_data_PEV_mean)
  loess_svd_PEV_median <- cbind(Type ="loess_svd", test6$Group_data_PEV_median)
  loess_svd_PEV_sd <- cbind(Type ="loess_svd", test6$Group_data_PEV_sd)
  
  test7 <- Group_data_PEV(rlr.knn.dat, sample)
  rlr_knn_PEV_mean <- cbind(Type ="rlr_knn", test7$Group_data_PEV_mean)
  rlr_knn_PEV_median <- cbind(Type ="rlr_knn", test7$Group_data_PEV_median)
  rlr_knn_PEV_sd <- cbind(Type ="rlr_knn", test7$Group_data_PEV_sd)
  
  test8 <- Group_data_PEV(rlr.lls.dat, sample)
  rlr_lls_PEV_mean <- cbind(Type ="rlr_lls", test8$Group_data_PEV_mean)
  rlr_lls_PEV_median <- cbind(Type ="rlr_lls", test8$Group_data_PEV_median)
  rlr_lls_PEV_sd <- cbind(Type ="rlr_lls", test8$Group_data_PEV_sd)
  
  test9 <- Group_data_PEV(rlr.svd.dat, sample)
  rlr_svd_PEV_mean <- cbind(Type ="rlr_svd", test9$Group_data_PEV_mean)
  rlr_svd_PEV_median <- cbind(Type ="rlr_svd", test9$Group_data_PEV_median)
  rlr_svd_PEV_sd <- cbind(Type ="rlr_svd", test9$Group_data_PEV_sd)
  
  ###PEV_mean
  #Combining all the above results
  total_pev_mean <- plyr::rbind.fill(vsn_knn_PEV_mean, vsn_lls_PEV_mean, vsn_svd_PEV_mean, 
                                     loess_knn_PEV_mean, loess_lls_PEV_mean, loess_svd_PEV_mean, 
                                     rlr_knn_PEV_mean, rlr_lls_PEV_mean, rlr_svd_PEV_mean)
  
  #Separating the results groupwise
  total_Group_data_PEV_mean <- total_pev_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_mean)
  
  #Extract the top combination in each group
  total_Group_data_PEV_mean2<-
    total_Group_data_PEV_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PEV_mean1 <- subset(total_Group_data_PEV_mean2, select = -row)
  final_Group_data_PEV_mean <- grouping_result(final_Group_data_PEV_mean1)
  
  ###PEV_median
  #Combining all the above results
  total_pev_median <- plyr::rbind.fill(vsn_knn_PEV_median, vsn_lls_PEV_median, vsn_svd_PEV_median, 
                                       loess_knn_PEV_median, loess_lls_PEV_median, loess_svd_PEV_median, 
                                       rlr_knn_PEV_median, rlr_lls_PEV_median, rlr_svd_PEV_median)
  
  #Separating the results groupwise
  total_Group_data_PEV_median <- total_pev_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_median)
  
  #Extract the top combination in each group
  total_Group_data_PEV_median2<-
    total_Group_data_PEV_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PEV_median1 <- subset(total_Group_data_PEV_median2, select = -row)
  final_Group_data_PEV_median <- grouping_result(final_Group_data_PEV_median1)
  
  ###PEV_sd
  #Combining all the above results
  total_pev_sd <- plyr::rbind.fill(vsn_knn_PEV_sd, vsn_lls_PEV_sd, vsn_svd_PEV_sd, 
                                   loess_knn_PEV_sd, loess_lls_PEV_sd, loess_svd_PEV_sd, 
                                   rlr_knn_PEV_sd, rlr_lls_PEV_sd, rlr_svd_PEV_sd)
  
  #Separating the results groupwise
  total_Group_data_PEV_sd <- total_pev_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PEV_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PEV_sd)
  
  #Extract the top combination in each group
  total_Group_data_PEV_sd2<-
    total_Group_data_PEV_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PEV_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PEV_sd1 <- subset(total_Group_data_PEV_sd2, select = -row)
  final_Group_data_PEV_sd <- grouping_result(final_Group_data_PEV_sd1)
  
  #Overall PEV function
  Total_data_PEV = function(data1){
    data <- as.matrix(data1)
    tempData = data[]
    rowNonNACnt = rowSums(!is.na(tempData)) - 1
    EVs = rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
    PEV_mean = mean(EVs, na.rm = FALSE)
    PEV_median = median(EVs, na.rm = FALSE)
    PEV_sd = sd(EVs, na.rm = FALSE)
    PEV_list = list("PEV_mean" = as.data.frame(PEV_mean),
                    "PEV_median" = as.data.frame(PEV_median),
                    "PEV_sd" = as.data.frame(PEV_sd))
    return(PEV_list)
  }
  
  #Overall PEV_mean estimation
  test1 <- Total_data_PEV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PEV_mean)
  
  test2 <- Total_data_PEV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PEV_mean)
  
  test3 <- Total_data_PEV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PEV_mean)
  
  test4 <- Total_data_PEV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PEV_mean)
  
  test5 <- Total_data_PEV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PEV_mean)
  
  test6 <- Total_data_PEV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PEV_mean)
  
  test7 <- Total_data_PEV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PEV_mean)
  
  test8 <- Total_data_PEV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PEV_mean)
  
  test9 <- Total_data_PEV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PEV_mean)
  
  #Combining all the above results
  total_pev_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_mean1 <-total_pev_overall_mean2%>%dplyr::slice_min(PEV_mean, n=1, with_ties = TRUE)
  total_pev_overall_mean <- grouping_result(total_pev_overall_mean1)
  original_cols <- c("Overall_Type.PEV_mean", "Overall_value.PEV_mean")
  colnames(total_pev_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PEV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PEV_median)
  
  test2 <- Total_data_PEV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PEV_median)
  
  test3 <- Total_data_PEV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PEV_median)
  
  test4 <- Total_data_PEV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PEV_median)
  
  test5 <- Total_data_PEV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PEV_median)
  
  test6 <- Total_data_PEV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PEV_median)
  
  test7 <- Total_data_PEV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PEV_median)
  
  test8 <- Total_data_PEV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PEV_median)
  
  test9 <- Total_data_PEV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PEV_median)
  
  #Combining all the above results
  total_pev_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_median1 <-total_pev_overall_median2%>%dplyr::slice_min(PEV_median, n=1, with_ties = TRUE)
  total_pev_overall_median <- grouping_result(total_pev_overall_median1)
  original_cols <- c("Overall_Type.PEV_median", "Overall_value.PEV_median")
  colnames(total_pev_overall_median) <- original_cols
  
  #Overall PEV_sd estimation
  test1 <- Total_data_PEV(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PEV_sd)
  
  test2 <- Total_data_PEV(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PEV_sd)
  
  test3 <- Total_data_PEV(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PEV_sd)
  
  test4 <- Total_data_PEV(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PEV_sd)
  
  test5 <- Total_data_PEV(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PEV_sd)
  
  test6 <- Total_data_PEV(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PEV_sd)
  
  test7 <- Total_data_PEV(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PEV_sd)
  
  test8 <- Total_data_PEV(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PEV_sd)
  
  test9 <- Total_data_PEV(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PEV_sd)
  
  #Combining all the above results
  total_pev_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pev_overall_sd1 <-total_pev_overall_sd2%>%dplyr::slice_min(PEV_sd, n=1, with_ties = TRUE)
  total_pev_overall_sd <- grouping_result(total_pev_overall_sd1)
  original_cols <- c("Overall_Type.PEV_sd", "Overall_value.PEV_sd")
  colnames(total_pev_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PEV<- cbind(final_Group_data_PEV_mean, final_Group_data_PEV_median,
                     final_Group_data_PEV_sd, total_pev_overall_mean,
                     total_pev_overall_median, total_pev_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PEV), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PEV_names1 <- result_PEV[, -col_indices]
  n <- matrix(t(result_PEV_names1), ncol=1)
  
  # Split the values by comma and remove extra spaces
  split_data <- strsplit(n, ",\\s*")
  
  # Flatten the list into a single vector
  separated_data <- unlist(split_data)
  
  # Get the frequency of each element in the dataframe
  freq <- table(separated_data)
  
  # Find the most occurring element
  PEV_best_combination <- names(freq)[which.max(freq)]
  
  #Groupwise PMAD function
  Group_data_PMAD = function(data1, groups){
    data <- as.matrix(data1)
    PMAD = NULL
    for(group in unique(groups)){
      tempData = data[,groups %in% group]
      MAD = matrixStats::rowMads(tempData, na.rm = FALSE)
      PMAD_mean[group] = mean(MAD, na.rm = T)
      PMAD_median[group] = median(MAD, na.rm = T)
      PMAD_sd[group] = sd(MAD, na.rm = T)
    }
    Group_data_PMAD_list = list("Group_data_PMAD_mean" = as.data.frame(PMAD_mean),
                                "Group_data_PMAD_median" = as.data.frame(PMAD_median),
                                "Group_data_PMAD_sd" = as.data.frame(PMAD_sd))
    return(Group_data_PMAD_list)
  }
  
  #PMAD calculation_groupwise
  test1 <- Group_data_PMAD(vsn.knn.dat, sample)
  vsn_knn_PMAD_mean <- cbind(Type ="vsn_knn", test1$Group_data_PMAD_mean)
  vsn_knn_PMAD_median <- cbind(Type ="vsn_knn", test1$Group_data_PMAD_median)
  vsn_knn_PMAD_sd <- cbind(Type ="vsn_knn", test1$Group_data_PMAD_sd)
  
  test2 <- Group_data_PMAD(vsn.lls.dat, sample)
  vsn_lls_PMAD_mean <- cbind(Type ="vsn_lls", test2$Group_data_PMAD_mean)
  vsn_lls_PMAD_median <- cbind(Type ="vsn_lls", test2$Group_data_PMAD_median)
  vsn_lls_PMAD_sd <- cbind(Type ="vsn_lls", test2$Group_data_PMAD_sd)
  
  test3 <- Group_data_PMAD(vsn.svd.dat, sample)
  vsn_svd_PMAD_mean <- cbind(Type ="vsn_svd", test3$Group_data_PMAD_mean)
  vsn_svd_PMAD_median <- cbind(Type ="vsn_svd", test3$Group_data_PMAD_median)
  vsn_svd_PMAD_sd <- cbind(Type ="vsn_svd", test3$Group_data_PMAD_sd)
  
  test4 <- Group_data_PMAD(loess.knn.dat, sample)
  loess_knn_PMAD_mean <- cbind(Type ="loess_knn", test4$Group_data_PMAD_mean)
  loess_knn_PMAD_median <- cbind(Type ="loess_knn", test4$Group_data_PMAD_median)
  loess_knn_PMAD_sd <- cbind(Type ="loess_knn", test4$Group_data_PMAD_sd)
  
  test5 <- Group_data_PMAD(loess.lls.dat, sample)
  loess_lls_PMAD_mean <- cbind(Type ="loess_lls", test5$Group_data_PMAD_mean)
  loess_lls_PMAD_median <- cbind(Type ="loess_lls", test5$Group_data_PMAD_median)
  loess_lls_PMAD_sd <- cbind(Type ="loess_lls", test5$Group_data_PMAD_sd)
  
  test6 <- Group_data_PMAD(loess.svd.dat, sample)
  loess_svd_PMAD_mean <- cbind(Type ="loess_svd", test6$Group_data_PMAD_mean)
  loess_svd_PMAD_median <- cbind(Type ="loess_svd", test6$Group_data_PMAD_median)
  loess_svd_PMAD_sd <- cbind(Type ="loess_svd", test6$Group_data_PMAD_sd)
  
  test7 <- Group_data_PMAD(rlr.knn.dat, sample)
  rlr_knn_PMAD_mean <- cbind(Type ="rlr_knn", test7$Group_data_PMAD_mean)
  rlr_knn_PMAD_median <- cbind(Type ="rlr_knn", test7$Group_data_PMAD_median)
  rlr_knn_PMAD_sd <- cbind(Type ="rlr_knn", test7$Group_data_PMAD_sd)
  
  test8 <- Group_data_PMAD(rlr.lls.dat, sample)
  rlr_lls_PMAD_mean <- cbind(Type ="rlr_lls", test8$Group_data_PMAD_mean)
  rlr_lls_PMAD_median <- cbind(Type ="rlr_lls", test8$Group_data_PMAD_median)
  rlr_lls_PMAD_sd <- cbind(Type ="rlr_lls", test8$Group_data_PMAD_sd)
  
  test9 <- Group_data_PMAD(rlr.svd.dat, sample)
  rlr_svd_PMAD_mean <- cbind(Type ="rlr_svd", test9$Group_data_PMAD_mean)
  rlr_svd_PMAD_median <- cbind(Type ="rlr_svd", test9$Group_data_PMAD_median)
  rlr_svd_PMAD_sd <- cbind(Type ="rlr_svd", test9$Group_data_PMAD_sd)
  
  ###PMAD_mean
  #Combining all the above results
  total_pmad_mean <- plyr::rbind.fill(vsn_knn_PMAD_mean, vsn_lls_PMAD_mean, vsn_svd_PMAD_mean, 
                                      loess_knn_PMAD_mean, loess_lls_PMAD_mean, loess_svd_PMAD_mean, 
                                      rlr_knn_PMAD_mean, rlr_lls_PMAD_mean, rlr_svd_PMAD_mean)
  
  #Separating the results groupwise
  total_Group_data_PMAD_mean <- total_pmad_mean %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_mean_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_mean)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_mean2<-
    total_Group_data_PMAD_mean%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_mean_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PMAD_mean1 <- subset(total_Group_data_PMAD_mean2, select = -row)
  final_Group_data_PMAD_mean <- grouping_result(final_Group_data_PMAD_mean1)
  
  ###PMAD_median
  #Combining all the above results
  total_pmad_median <- plyr::rbind.fill(vsn_knn_PMAD_median, vsn_lls_PMAD_median, vsn_svd_PMAD_median, 
                                        loess_knn_PMAD_median, loess_lls_PMAD_median, loess_svd_PMAD_median, 
                                        rlr_knn_PMAD_median, rlr_lls_PMAD_median, rlr_svd_PMAD_median)
  
  #Separating the results groupwise
  total_Group_data_PMAD_median <- total_pmad_median %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_median_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_median)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_median2<-
    total_Group_data_PMAD_median%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_median_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result 
  final_Group_data_PMAD_median1 <- subset(total_Group_data_PMAD_median2, select = -row)
  final_Group_data_PMAD_median <- grouping_result(final_Group_data_PMAD_median1)
  
  ###PMAD_sd
  #Combining all the above results
  total_pmad_sd <- plyr::rbind.fill(vsn_knn_PMAD_sd, vsn_lls_PMAD_sd, vsn_svd_PMAD_sd, 
                                    loess_knn_PMAD_sd, loess_lls_PMAD_sd, loess_svd_PMAD_sd, 
                                    rlr_knn_PMAD_sd, rlr_lls_PMAD_sd, rlr_svd_PMAD_sd)
  
  #Separating the results groupwise
  total_Group_data_PMAD_sd <- total_pmad_sd %>%
    dplyr::group_by(Type) %>%
    dplyr::mutate(Group = paste0("PMAD_sd_Group", 1:dplyr::n())) %>%
    tidyr::pivot_wider(names_from = Group, values_from = PMAD_sd)
  
  #Extract the top combination in each group
  total_Group_data_PMAD_sd2<-
    total_Group_data_PMAD_sd%>%
    tidyr::pivot_longer(-Type)%>%
    dplyr::group_by(name)%>%
    dplyr::slice_min(value, n=1)%>%
    dplyr::mutate(row = dplyr::row_number())%>%
    dplyr::ungroup()%>%
    tidyr::pivot_wider(names_from = name, names_prefix = "PMAD_sd_", names_sep = ".",values_from = c(Type, value),
                       names_vary = "slowest")%>%
    stats::setNames(nm = sub("(.*)_(.*)", "\\2_\\1", names(.)))
  
  #Final result
  final_Group_data_PMAD_sd1 <- subset(total_Group_data_PMAD_sd2, select = -row)
  final_Group_data_PMAD_sd <- grouping_result(final_Group_data_PMAD_sd1)
  
  #Overall PMAD function
  Total_data_PMAD = function(data1){
    data <- as.matrix(data1)
    tempData = data[]
    MADs = matrixStats::rowMads(tempData, na.rm = FALSE)
    PMAD_mean = mean(MADs, na.rm = T)
    PMAD_median = median(MADs, na.rm = T)
    PMAD_sd = sd(MADs, na.rm = T)
    PMAD_list = list("PMAD_mean" = as.data.frame(PMAD_mean),
                     "PMAD_median" = as.data.frame(PMAD_median),
                     "PMAD_sd" = as.data.frame(PMAD_sd))
    return(PMAD_list)
  }
  
  #Overall PMAD_mean estimation
  test1 <- Total_data_PMAD(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PMAD_mean)
  
  test2 <- Total_data_PMAD(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PMAD_mean)
  
  test3 <- Total_data_PMAD(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PMAD_mean)
  
  test4 <- Total_data_PMAD(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PMAD_mean)
  
  test5 <- Total_data_PMAD(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PMAD_mean)
  
  test6 <- Total_data_PMAD(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PMAD_mean)
  
  test7 <- Total_data_PMAD(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PMAD_mean)
  
  test8 <- Total_data_PMAD(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PMAD_mean)
  
  test9 <- Total_data_PMAD(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PMAD_mean)
  
  
  #Combining all the above results
  total_pmad_overall_mean2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_mean1 <-total_pmad_overall_mean2%>%dplyr::slice_min(PMAD_mean, n=1, with_ties = TRUE)
  total_pmad_overall_mean <- grouping_result(total_pmad_overall_mean1)
  original_cols <- c("Overall_Type.PMAD_mean", "Overall_value.PMAD_mean")
  colnames(total_pmad_overall_mean) <- original_cols
  
  #Overall PCV_median estimation
  test1 <- Total_data_PMAD(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PMAD_median)
  
  test2 <- Total_data_PMAD(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PMAD_median)
  
  test3 <- Total_data_PMAD(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PMAD_median)
  
  test4 <- Total_data_PMAD(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PMAD_median)
  
  test5 <- Total_data_PMAD(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PMAD_median)
  
  test6 <- Total_data_PMAD(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PMAD_median)
  
  test7 <- Total_data_PMAD(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PMAD_median)
  
  test8 <- Total_data_PMAD(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PMAD_median)
  
  test9 <- Total_data_PMAD(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PMAD_median)
  
  #Combining all the above results
  total_pmad_overall_median2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_median1 <-total_pmad_overall_median2%>%dplyr::slice_min(PMAD_median, n=1, with_ties = TRUE)
  total_pmad_overall_median <- grouping_result(total_pmad_overall_median1)
  original_cols <- c("Overall_Type.PMAD_median", "Overall_value.PMAD_median")
  colnames(total_pmad_overall_median) <- original_cols
  
  #Overall PMAD_sd estimation
  test1 <- Total_data_PMAD(vsn.knn.dat)
  data1 <- cbind(Type ="vsn_knn", test1$PMAD_sd)
  
  test2 <- Total_data_PMAD(vsn.lls.dat)
  data2 <- cbind(Type ="vsn_lls", test2$PMAD_sd)
  
  test3 <- Total_data_PMAD(vsn.svd.dat)
  data3 <- cbind(Type ="vsn_svd", test3$PMAD_sd)
  
  test4 <- Total_data_PMAD(loess.knn.dat)
  data4 <- cbind(Type ="loess_knn", test4$PMAD_sd)
  
  test5 <- Total_data_PMAD(loess.lls.dat)
  data5 <- cbind(Type ="loess_lls", test5$PMAD_sd)
  
  test6 <- Total_data_PMAD(loess.svd.dat)
  data6 <- cbind(Type ="loess_svd", test6$PMAD_sd)
  
  test7 <- Total_data_PMAD(rlr.knn.dat)
  data7 <- cbind(Type ="rlr_knn", test7$PMAD_sd)
  
  test8 <- Total_data_PMAD(rlr.lls.dat)
  data8 <- cbind(Type ="rlr_lls", test8$PMAD_sd)
  
  test9 <- Total_data_PMAD(rlr.svd.dat)
  data9 <- cbind(Type ="rlr_svd", test9$PMAD_sd)
  
  #Combining all the above results
  total_pmad_overall_sd2 <- as.data.frame(plyr::rbind.fill(data1, data2, data3, data4, data5, data6, data7, data8, data9))
  
  #Extract the top combination in overall
  total_pmad_overall_sd1 <-total_pmad_overall_sd2%>%dplyr::slice_min(PMAD_sd, n=1, with_ties = TRUE)
  total_pmad_overall_sd <- grouping_result(total_pmad_overall_sd1)
  original_cols <- c("Overall_Type.PMAD_sd", "Overall_value.PMAD_sd")
  colnames(total_pmad_overall_sd) <- original_cols
  
  #Combining the groupwise, overall_ mean, overall_median and overall_SD results
  result_PMAD<- cbind(final_Group_data_PMAD_mean, final_Group_data_PMAD_median,
                      final_Group_data_PMAD_sd, total_pmad_overall_mean,
                      total_pmad_overall_median, total_pmad_overall_sd)
  
  #Processing
  #Removing the even number columns to extract the combination names
  # Get the column indices
  col_indices <- seq(2, ncol(result_PMAD), by = 2)  # Get even-numbered column indices
  
  # Remove even-numbered columns
  result_PMAD_names1 <- result_PMAD[, -col_indices]
  n <- matrix(t(result_PMAD_names1), ncol=1)
  # Split the values by comma and remove extra spaces
  split_data <- strsplit(n, ",\\s*")
  
  # Flatten the list into a single vector
  separated_data <- unlist(split_data)
  
  # Get the frequency of each element in the dataframe
  freq <- table(separated_data)
  
  # Find the most occurring element
  PMAD_best_combination <- names(freq)[which.max(freq)]
  
  #Finding the best combination
  Best_combinations <- cbind(PCV_best_combination, PEV_best_combination, PMAD_best_combination)
  
  #Adding names to table
  Combinations <- c("vsn_knn", "vsn_lls", "vsn_svd",
                    "loess_knn", "loess_lls", "loess_svd",
                    "rlr_knn", "rlr_lls", "rlr_svd")
  #Extracting PCV values for all combinations
  pcv_group_mean <- total_Group_data_PCV_mean
  pcv_group_median <- total_Group_data_PCV_median
  pcv_group_sd <- total_Group_data_PCV_sd
  
  pcv_overall_mean <- as.data.frame(total_pcv_overall_mean2[,-1])
  colnames(pcv_overall_mean) <- "Overall_PCV_mean"
  pcv_overall_median <- as.data.frame(total_pcv_overall_median2[,-1])
  colnames(pcv_overall_median) <- "Overall_PCV_median"
  pcv_overall_sd <- as.data.frame(total_pcv_overall_sd2[,-1])
  colnames(pcv_overall_sd) <- "Overall_PCV_sd"
  PCV_table2 <- suppressMessages(cbind(pcv_group_mean, pcv_group_median, pcv_group_sd,
                                       pcv_overall_mean, pcv_overall_median, pcv_overall_sd))
  PCV_table1 <- PCV_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PCV_table <- cbind(Combinations, PCV_table1)
  
  #Extracting PEV values for all combinations
  pev_group_mean <- total_Group_data_PEV_mean
  pev_group_median <- total_Group_data_PEV_median
  pev_group_sd <- total_Group_data_PEV_sd
  
  pev_overall_mean <- as.data.frame(total_pev_overall_mean2[,-1])
  colnames(pev_overall_mean) <- "Overall_PEV_mean"
  pev_overall_median <- as.data.frame(total_pev_overall_median2[,-1])
  colnames(pev_overall_median) <- "Overall_PEV_median"
  pev_overall_sd <- as.data.frame(total_pev_overall_sd2[,-1])
  colnames(pev_overall_sd) <- "Overall_PEV_sd"
  PEV_table2 <- suppressMessages(cbind(pev_group_mean, pev_group_median, pev_group_sd,
                                       pev_overall_mean, pev_overall_median, pev_overall_sd))
  PEV_table1 <- PEV_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PEV_table <- cbind(Combinations, PEV_table1)
  
  #Extracting PMAD values for all combinations
  pmad_group_mean <- total_Group_data_PMAD_mean
  pmad_group_median <- total_Group_data_PMAD_median
  pmad_group_sd <- total_Group_data_PMAD_sd
  
  pmad_overall_mean <- as.data.frame(total_pmad_overall_mean2[,-1])
  colnames(pmad_overall_mean) <- "Overall_PMAD_mean"
  pmad_overall_median <- as.data.frame(total_pmad_overall_median2[,-1])
  colnames(pmad_overall_median) <- "Overall_PMAD_median"
  pmad_overall_sd <- as.data.frame(total_pmad_overall_sd2[,-1])
  colnames(pmad_overall_sd) <- "Overall_PMAD_sd"
  PMAD_table2 <- suppressMessages(cbind(pmad_group_mean, pmad_group_median, pmad_group_sd,
                                        pmad_overall_mean, pmad_overall_median, pmad_overall_sd))
  PMAD_table1 <- PMAD_table2 %>% dplyr::select(-tidyselect::contains("Type"))
  PMAD_table <- cbind(Combinations, PMAD_table1)
  
  # Apply imputation methods and handle missing values
  vsn.knn.dat[is.na(vsn.knn.dat)] <- loess.knn.dat[is.na(rlr.knn.dat)] <- rlr.knn.dat[is.na(rlr.knn.dat)] <- 0
  vsn.lls.dat[is.na(vsn.lls.dat)] <- loess.lls.dat[is.na(rlr.lls.dat)] <- rlr.lls.dat[is.na(rlr.lls.dat)] <- 0
  vsn.svd.dat[is.na(vsn.svd.dat)] <- loess.svd.dat[is.na(rlr.svd.dat)] <- rlr.svd.dat[is.na(rlr.svd.dat)] <- 0
  com_data2[is.na(com_data2)] <- 0
  
  #Calculating NRMSE
  nrmse <- function(ximp, xtrue) {
    # Convert both inputs to numeric vectors
    ximp_vector <- as.numeric(as.matrix(ximp))
    xtrue_vector <- as.numeric(as.matrix(xtrue))
    
    # Calculate NRMSE
    sqrt(mean((ximp_vector - xtrue_vector)^2, na.rm = TRUE) / stats::var(xtrue_vector, na.rm = TRUE))
  }  
  nrmse_vsn_knn <- nrmse(2^vsn.knn.dat, com_data2)
  nrmse_vsn_lls <- nrmse(2^vsn.lls.dat, com_data2)
  nrmse_vsn_svd <- nrmse(2^vsn.svd.dat, com_data2)
  nrmse_loess_knn <- nrmse(2^loess.knn.dat, com_data2)
  nrmse_loess_lls <- nrmse(2^loess.lls.dat, com_data2)
  nrmse_loess_svd <- nrmse(2^loess.svd.dat, com_data2)
  nrmse_rlr_knn <- nrmse(2^rlr.knn.dat, com_data2)
  nrmse_rlr_lls <- nrmse(2^rlr.lls.dat, com_data2)
  nrmse_rlr_svd <- nrmse(2^rlr.svd.dat, com_data2)
  
  nrmse <- c(nrmse_vsn_knn, nrmse_vsn_lls, nrmse_vsn_svd,
             nrmse_loess_knn, nrmse_loess_lls, nrmse_loess_svd,
             nrmse_rlr_knn, nrmse_rlr_lls, nrmse_rlr_svd)
  
  nrmse_1 <- as.data.frame(round(nrmse, digits = 5))
  
  nrmse_result <- cbind(Combinations, nrmse_1)
  colnames(nrmse_result) <- c("Combinations", "NRMSE")  
  
  result_list <- list("Best combinations" = as.data.frame(Best_combinations), "PCV Result" = PCV_table, "PEV Result" = PEV_table, "PMAD Result" =PMAD_table, "NRMSE Result" = nrmse_result,
                      "vsn_data" = cbind(com_data_ID, vsn.dat), "loess_data" = cbind(com_data_ID, loess.dat), "rlr_data" = cbind(com_data_ID, rlr.dat),
                      "vsn_knn_data" = cbind(com_data_ID,vsn.knn.dat),  "vsn_lls_data" = cbind(com_data_ID,vsn.lls.dat),"vsn_svd_data" = cbind(com_data_ID,vsn.svd.dat),
                      "loess_knn_data" = cbind(com_data_ID, as.data.frame(loess.knn.dat)), "loess_lls_data" =  cbind(com_data_ID, as.data.frame(loess.lls.dat)), "loess_svd_data" = cbind(com_data_ID, as.data.frame(loess.svd.dat)),
                      "rlr_knn_data" =  cbind(com_data_ID, as.data.frame(rlr.knn.dat)), "rlr_lls_data" =  cbind(com_data_ID, as.data.frame(rlr.lls.dat)),"rlr_svd_data" =  cbind(com_data_ID, as.data.frame(rlr.svd.dat)))
  sink()
  return (result_list)
}
