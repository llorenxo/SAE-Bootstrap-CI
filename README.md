# Bootstrap Confidence Intervals for Unit-Level Small Area Estimation

This repository contains the source code of the Shiny application accompanying
the paper:

**Confidence Intervals for Unit-Level Small Area Models: A Bootstrap Comparison**

## Online application

The application is available at:

https://sae-boot-conf-int.shinyapps.io/shinyapp/

## Description

The application implements confidence interval procedures for unit-level
empirical best prediction in Small Area Estimation, including:

- Normal
- Percentile
- Reverse Percentile
- BCa (area-level jackknife)

The application also provides bootstrap diagnostics and supports alternative
response transformations.


```r
renv::restore()
shiny::runApp()
