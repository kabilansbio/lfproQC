% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/best_combination.R
\name{best_combination}
\alias{best_combination}
\title{Best combination of normalization and imputation method}
\usage{
best_combination(data_input, groups, data_type, aggr_method)
}
\arguments{
\item{data_input}{Label-free proteomics expression data as a dataframe}

\item{groups}{Group information about the input data}

\item{data_type}{A character string specifying the type of data being used. Use "Peptide" if your dataset contains peptide information, where the first column represents peptide IDs and the second column represents protein IDs. Use "Protein" if your dataset consists of protein data, where the first column represents protein IDs. The function will handle the data accordingly based on this parameter.}

\item{aggr_method}{A character string specifying the method for aggregating peptide data to corresponding protein values. Use "sum" to aggregate by the sum of peptide intensities, "mean" to aggregate by the mean of peptide intensities, or "median" to aggregate by the median of peptide intensities. This parameter is only applicable when "data_type" is set to "Peptide".}
}
\value{
This function gives the list  which consist of following results.

`Best Combinations`  The best combinations based on each PCV, PEV and PMAD
for the given dataset.

`PCV Result` Values of groupwise PCV, overall PCV, PCV mean, PCV median and
             PCV standard deviation for all combinations.

`PEV Result` Values of groupwise PEV, overall PEV, PEV mean, PEV median and
             PEV standard deviation for all combinations.

`PMAD Result` Values of groupwise PMAD, overall PMAD, PMAD mean, PMAD median
              and PMAD standard deviation for all combinations.

`NRMSE Result` NRMSE values calculated for the normalized and imputed dataset 
               to the original dataset. 

`rollup_protein` The aggregated protein values for the peptide dataset are 
                 based on either the sum, mean, or median.
               
`vsn_data` The `vsn` normalized dataset

`loess_data` The `loess` normalized dataset

`rlr_data` The `rlr` normalized dataset

`vsn_knn_data` The dataset normalized by `vsn` method and missing values imputed
             by `knn` method.
             
`vsn_lls_data` The dataset normalized by `vsn` method and missing values imputed
             by `lls` method.
             
`vsn_svd_data` The dataset normalized by `vsn` method and missing values imputed
             by `svd` method.

`loess_knn_data` The dataset normalized by `loess` method and missing values imputed
             by `knn` method.

`loess_lls_data` The dataset normalized by `loess` method and missing values imputed
             by `lls` method.
             
`loess_svd_data` The dataset normalized by `loess` method and missing values imputed
             by `svd` method.
             
`rlr_knn_data` The dataset normalized by `rlr` method and missing values imputed
             by `knn` method.

`rlr_lls_data` The dataset normalized by `rlr` method and missing values imputed
             by `lls` method.

`rlr_svd_data` The dataset normalized by `rlr` method and missing values imputed
             by `svd` method.
}
\description{
This function will provide the best combinations of normalization
and imputation methods for the user given dataset based on the intragroup
variation evaluation parameters called PCV, PEV, PMAD and NRMSE.
}
\details{
Label-free LC-MS proteomics expression data is often affected by heterogeneity and missing values. 
Normalization and missing value imputation are the commonly used techniques to solve these issues and make the dataset suitable for further downstream analysis. 
This function provides the best combination of normalization and imputation methods for the dataset, choosing from the three normalization methods (vsn, loess, and rlr) and three imputation methods (knn, lls, svd). 
The intragroup variation evaluation measures named pooled co-efficient of variance (PCV), pooled estimate of variance (PEV) and pooled median absolute deviation (PMAD) are used for selecting the best combination of normalization and imputation method for the given dataset.
It will return the best combinations based on each evaluation parameters of PCV, PEV, and PMAD. 

Along with this, the user can get all three normalized datasets, nine combinations of normalized and missing values imputed datasets, and the PCV, PEV, and PMAD result values. The user can also obtain the Normalized Root Mean Square Error (NRMSE) values, 
calculated by comparing the normalized and imputed dataset with the original dataset across all nine combinations.
These NRMSE values provide insight into the accuracy of the imputation and normalization processes.
}
\examples{
\donttest{
result <- best_combination(yeast_data, yeast_groups, data_type = "Protein")
result$`Best combinations`
result$`PCV Result`
result$`PMAD Result`
result$`rlr_knn_data`
}
}
\author{
Dr Sudhir Srivastava ("Sudhir.Srivastava@icar.gov.in")

 Kabilan S ("kabilan151414@gmail.com")
}
