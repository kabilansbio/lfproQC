## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(lfproQC)

## ----results='markup'---------------------------------------------------------
yeast <- best_combination(yeast_data, yeast_groups)

## -----------------------------------------------------------------------------
yeast$`PCV Result`

## -----------------------------------------------------------------------------
yeast$`PEV Result`

## -----------------------------------------------------------------------------
yeast$`PMAD Result`

## -----------------------------------------------------------------------------
yeast$`Best combinations`

## ----out.width = "400px"------------------------------------------------------
Boxplot_data(yeast$knn_rlr_data) 

## ----out.width = "400px"------------------------------------------------------
Densityplot_data(yeast$knn_rlr_data)

## ----out.width = "400px"------------------------------------------------------
Corrplot_data(yeast$knn_rlr_data)

## ----out.width = "400px"------------------------------------------------------
MDSplot_data(yeast$knn_rlr_data)

## ----out.width = "400px"------------------------------------------------------
QQplot_data(yeast$knn_rlr_data)

## ----results='hide'-----------------------------------------------------------
top_table_yeast <- top_table_fn(yeast$knn_rlr_data, yeast_groups, 2, 1)

## ----out.width = "400px"------------------------------------------------------
de_yeast_MA <- MAplot_DE_fn(top_table_yeast,-1,1,0.05)
de_yeast_MA$`MA Plot`

## ----out.width = "400px"------------------------------------------------------
de_yeast_volcano <- volcanoplot_DE_fn (top_table_yeast,-1,1,0.05)
de_yeast_volcano$`Volcano Plot`

## ----results= 'hide'----------------------------------------------------------
de_yeast_MA$`Result `

## ----results='hide'-----------------------------------------------------------
de_yeast_MA$`Up-regulated`

## ----results='hide'-----------------------------------------------------------
de_yeast_MA$`Down-regulated`

## ----results='hide'-----------------------------------------------------------
de_yeast_MA$`Significant`

## ----results='hide'-----------------------------------------------------------
de_yeast_MA$`Non-significant`

## ----setup, include=FALSE-----------------------------------------------------
library(knitr)

## ----echo=FALSE, out.width = "800px"------------------------------------------
knitr::include_graphics("images/flow1.jpg")

## -----------------------------------------------------------------------------
sessionInfo()

