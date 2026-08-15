# Bootstrap Confidence Intervals for Unit-Level Small Area Estimation

This repository contains the source code of the Shiny application accompanying the paper:

**From Variability to Confidence: Exploring Bootstrap Confidence Interval Methods for Unit-Level Small Area Models**

by Eleonora Tizianel, Lorenzo Mori, and Maria Rosaria Ferrante  
Department of Statistical Sciences, University of Bologna.

## Online application

The application is publicly available at:

https://sae-boot-conf-int.shinyapps.io/shinyapp/

## Overview

The application provides an interactive implementation of confidence interval procedures for unit-level Empirical Best Prediction (EBP) in Small Area Estimation.

It builds on the `emdi` framework and extends its bootstrap routines to retain the full bootstrap distribution of the target indicators, allowing alternative confidence interval constructions and bootstrap diagnostics to be computed.

## Confidence interval methods

The application currently implements:

- Normal
- Percentile
- Reverse Percentile
- BCa with area-level jackknife

The area-level BCa implementation is used because it provides a computationally more convenient alternative to the unit-level jackknife while retaining similar performance in the simulation study.

## Target indicators

The application supports estimation of:

- Mean
- Headcount Ratio (HCR)
- Gini coefficient

## Response transformations

The following transformations are available:

- No transformation
- Logarithmic
- Box-Cox
- Dual Power
- Log-shift
- Logit
- Probit
- Arc-Sine

The Logit, Probit, and Arc-Sine transformations extend the standard workflow to responses defined on the unit interval.

## Main features

The application allows users to:

- upload population and sample datasets;
- specify a unit-level EBP model;
- choose the response transformation and target indicator;
- compute alternative confidence intervals;
- compare confidence interval methods across domains;
- inspect domain-specific bootstrap distributions and diagnostics;
- visualize confidence intervals;
- download numerical results and reproducible reports.

## Repository contents

- `app.R` — Shiny user interface, server functions, confidence interval computation, and diagnostic tools.
- `Bootstrap_emdi_new_function.R` — modified `emdi` routines, additional transformations, and bootstrap functions.
- `pop.rds` — example population dataset.
- `smp.rds` — example sample dataset.
- `renv.lock` — R package versions required to reproduce the software environment.
- `CITATION.cff` — citation information for the software.

## Running the application locally

Clone or download this repository and open the repository directory in R.

Install `renv`, if necessary:

```r
install.packages("renv")
```

Restore the package environment:

```r
renv::restore()
```

Then launch the application:

```r
shiny::runApp()
```

Alternatively, open `app.R` in RStudio and click **Run App**.

## Example data

The repository includes example population (`pop.rds`) and sample (`smp.rds`) datasets.

They can be loaded directly from the application by selecting:

**Use example data included in the app**

This allows the application to be tested without uploading external data.

## Citation

If you use this software in your research, please cite the repository using the citation information provided in `CITATION.cff`.

Once the accompanying paper is published, the full bibliographic reference will also be added here.

## Authors

- **Eleonora Tizianel**  
  Department of Statistical Sciences, University of Bologna

- **Lorenzo Mori**  
  Department of Statistical Sciences, University of Bologna

- **Maria Rosaria Ferrante**  
  Department of Statistical Sciences, University of Bologna

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.
