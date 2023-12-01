## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}
source('functions.R')
library(tidyverse)
proj_dir <- getwd()
data_dir <- file.path(proj_dir, 'data')
results_dir <- file.path(proj_dir, 'output')
fig_dir <- file.path(proj_dir, 'analysis')
