# Information ----
# Analyzing test data from breathalyzer
# Author: Diego Rodríguez Esperante
# Date of creation: 12/11/2024
# Last edited: 12/11/2024

# Loading ----
## Loading packages ----
require(dplyr)
require(ggfortify)
require(reshape2)
require(ggplot2)
require(rstudioapi)

# Open files ----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_raw = read.csv("dev4891_imps_0_sample_00000.csv", sep = ";")

# Experimental constants ----
chunk_duration = 20 # Each chunk represents a 5-second frequency sweep
conc_duration = 300 # Every 5 minutes concentration increases
conc_increase = c(0, 0.1, 0.5, 1, 10, 50, 100) # Concentrations increase to these values in nM

# Cleanup ----
data = data_raw

# Cleanup rows inbetween chunks
data$timestamp = as.numeric(data$timestamp)
data = data[!is.na(data$timestamp),]
data$chunk = as.numeric(data$chunk)
data$size = as.numeric(data$size)

# List of dataframes, one for each chunk
datalist = split(data, data$chunk)
names(datalist) = paste0(rep("chunk_", length(datalist)), names(datalist))
for (i in 1:length(datalist)){
  # From the source data
  nms = datalist[[i]]$fieldname
  chnk = datalist[[i]]$chunk[1]
  tmstmp = datalist[[i]]$timestamp[1]
  sz = datalist[[i]]$size[1]
  
  # From experimental constants
  time = chnk*chunk_duration
  conc = conc_increase[floor(time/conc_duration)+1]
  
  data_t = as.data.frame(t(datalist[[i]][,c(5:ncol(data))]))
  data_t = cbind(rep(chnk, nrow(data_t)), rep(tmstmp, nrow(data_t)), rep(sz, nrow(data_t)), data_t)
  colnames(data_t) = c("chunk", "timestamp", "size", nms)
  
  datalist[[i]] = data_t
}

data_join = datalist[[1]]
for (i in 2:length(datalist)){
  data_join = rbind(data_join, datalist[[2]])
}
rownames(data_join) = NULL

rm(data_t, i, nms, sz, tmstmp, chnk)

# Replicate Joe's results ----
