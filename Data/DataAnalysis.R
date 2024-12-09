# Information ----
# Analyzing test data from breathalyzer
# Author: Diego Rodríguez Esperante
# Date of creation: 12/11/2024
# Last edited: 29/11/2024

# Loading ----
## Loading packages ----
require(dplyr)
require(ggfortify)
require(reshape2)
require(ggplot2)
require(stringr)
require(gridExtra)
require(cowplot)
require(rstudioapi)
require(NMF)

# Open files ----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

analyte = "Acetaldehyde"
file_analyze = choose.files(paste0(c(getwd(), "Data/TMPS", analyte, "Choose"), collapse = "/"))

data_raw = read.csv(file_analyze, sep = ";")

# str_match : sample_\d{5}
# Options ----
opts = list(deep_clean = TRUE,
            joe_plots = FALSE,
            make_plots = FALSE,
            save_plots = FALSE)

# Experimental constants ----
chunk_duration = 20 # Each chunk represents a 5-second frequency sweep
conc_duration = 300 # Every 5 minutes concentration increases
conc_increase = c(0, 0.1, 0.5, 1, 10, 50, 100) # Concentrations increase to these values in nM

# Support functions ----
extract_freq = function(dataframe, target_phase = -45){
  phase_deg = dataframe$phasez[dataframe$chunk == 0] * 180/pi
  phase_close = order(abs(phase_deg-target_phase))
  
  freqs = dataframe$frequency[dataframe$chunk == 0]
  target_freq = freqs[phase_close[1]]
  
  data_target = dataframe[dataframe$frequency == target_freq,]
  
  return(data_target)
}

cutie_layer = function(){
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = 'bottom',
        legend.key.width = unit(3, 'lines'))
}

customggsave = function(plot, upscale=1.5, save_path = '', name = NULL){
  save_path = paste0('./Plots', save_path)
  if (is.null(name)){
      name = deparse(substitute(plot))
  }
  ggsave(paste0(name,".png"),
         plot = plot,
         device = 'png',
         width = round(1920*upscale),
         height = round(1080*upscale),
         units = 'px',
         path = save_path)
}

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
  data_t = cbind(rep(chnk, nrow(data_t)), rep(tmstmp, nrow(data_t)), rep(sz, nrow(data_t)), 
                 rep(time, nrow(data_t)), rep(conc, nrow(data_t)), 
                 data_t)
  colnames(data_t) = c("chunk", "timestamp", "size", "time", "concentration", nms)
  
  datalist[[i]] = data_t
}

data_join = datalist[[1]]
for (i in 2:length(datalist)){
  data_join = rbind(data_join, datalist[[i]])
}
rownames(data_join) = NULL

if(opts$deep_clean){
  col_nms = colnames(data_join)
  # Whitelist specific columns
  white_cols = c("frequency")
  
  # Delete specific columns
  manu_cols = c("timestamp", "tc", "tcmeas", "count", "flags", "nexttimestamp", "settimestamp",
               "param0", "param0pwr", "param0stddev", "param1", "param1pwr", "param1stddev")
  
  # Delete power columns for redundancy
  pwr_cols = col_nms[grepl('pwr', col_nms)]
  
  # Delete columns whose standard deviation is 0
  std_cols = col_nms[grepl('stddev', col_nms)]
  std_del = c()
  for (i in 1:length(std_cols)){
    if ((sum(data_join[,std_cols[i]]==0)/(nrow(data_join)))>=.9){
      std_del = c(std_del,std_cols[i])
    }
  }
  
  std_cols = col_nms[grepl(paste(gsub("stddev", "", std_del), collapse = "|"), col_nms)]
  
  # Join together
  del_cols = c(manu_cols, pwr_cols, std_cols)
  del_cols = del_cols[!duplicated(del_cols)]
  
  # Remove whitelisted columns
  for (i in 1:length(white_cols)){
    del_cols = del_cols[!(del_cols == white_cols[i])]
  }
  
  data_join = data_join %>% select(-one_of(del_cols))
  
  rm(col_nms, white_cols, manu_cols, pwr_cols, std_cols, std_del, del_cols)
}

# Reorder columns
col_order = c("chunk", "size", "time", "concentration", "frequency", "settling", "bandwidth", "grid")
col_order = c(col_order, colnames(data_join)[!(colnames(data_join) %in% col_order)])

data_join = data_join[,col_order]

rm(data_t, i, nms, sz, tmstmp, chnk, time, conc, col_order)

# Replicate Joe's results ----
if(opts$joe_plots){
  # Plotting time against absz, only the frequency associated with -45 degree phase on the first chunk
  
  # Extract phases of first chunk
  data_target = extract_freq(data_join, -45)
  
  # Find first local maxima
  data_diff = diff(data_target$absz, differences = 2) > 0
  i = 1
  while (data_diff[i] == data_diff[1]) {
    i = i+1
  }
  data_effsz = cbind(Time = data_target$time, EffectSize = data_target$absz/(data_target$absz[i-1]))
  
  # Plotting
  ggplot(data_target, aes(x = time, y = absz)) +
    geom_line() + 
    geom_vline(xintercept = conc_duration*(1:floor(max(data_target$time)/conc_duration)))
  
  ggplot(data_effsz, aes(x = Time, y = EffectSize)) +
    geom_line()
}
# Other analyses ----
## PCA ----
pca_opts = list(targeting = "untargeted",
                fields = "reduced",
                scaling = "unscaled",
                color = "absz")

# Frequency targeting
if(pca_opts$targeting == "targeted"){
  data_pca = extract_freq(data_join, -45)
}else{
  data_pca = data_join  
}

# Z-score normalization just in case
if(pca_opts$scaling == "scaled"){
  for (i in 7:ncol(data_pca)){
    data_pca[,i] = (data_pca[,i] - mean(data_pca[,i])) / sd(data_pca[,i])
  }
}

# Take in extraneous fields left behind
# ignored_fields = c("chunk", "size", "time", "concentration", "frequency")
if(pca_opts$fields == "reduced"){
  pca_fields = c("absz", "abszstddev", "imagz", "imagzstddev", "phasez", "phasezstddev", "realz", "realzstddev")
  pca = prcomp(data_pca[,pca_fields])
  
  rm(pca_fields)
}else{
  pca = prcomp(data_pca[,-c(1:6)])
}

if(opts$make_plots){
  pca_loadings = signif(-1*pca$rotation, digits = 3)
  pca_obj = pca$x
  plottype = c("PCA", paste(c("Target:", "Fields:", "Scaling:", "Color:"), pca_opts))
  
  pca_plot = ggplot(pca_obj, aes(x = PC1, y= PC2, color = data_pca[, pca_opts$color])) +
    geom_point(size = 2) +
    labs(x = paste0('PC1', ' (', (summary(pca)$importance[2,1]), ')'),
         y = paste0('PC2', ' (', (summary(pca)$importance[2,2]), ')'),
         color = pca_opts$color) +
    scale_color_continuous(low = "magenta", high = "cyan") +
    ggtitle(paste(plottype, collapse = ", ")) +
    theme_bw() + cutie_layer()
  
  pca_plot_tab = plot_grid(pca_plot, tableGrob(pca_loadings[,c("PC1","PC2")]), rel_widths = c(3,1)) +
    theme_bw()
  if(opts$save_plots){
    customggsave(pca_plot_tab, upscale = 2, name = paste(c("PCA", pca_opts), collapse = "_"))
  }
}

rm(plottype)

## Correlation plot ----


## Non-negative matrix factorization ----


## Other features ----


# Save data ----
