if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tidyverse")

#clean workspace
rm(list=ls())

setwd("C:/Users/friederike.ehrhart/Documents/EJP-RD/Workshop BYOOD 17.11.22 Nijmegen")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

#import data
Metabolomics = read.table(file = 'Metabolomics.txt', sep = '\t', header = TRUE)
Proteomics = read.table(file = 'Proteomics.txt', sep = '\t', header = TRUE)

#shove data in appropriate format
# create metadata file

#example data
data <- make_example_data(
  n_views = 2, 
  n_samples = 200, 
  n_features = 1000,
  n_factors = 10
)[[1]]

lapply(data,dim)

#Create the MOFA object:
MOFAobject <- create_mofa(data)

plot_data_overview(MOFAobject)

#Define data groups
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

#Define model options
model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

#Define train options
train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

#Prepare the MOFA object
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#Train the MOFA model
outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

slotNames(MOFAobject.trained)

