# Information ----
# Analyzing test data from breathalyzer
# Author: Diego Rodríguez Esperante
# Date of creation: 12/11/2024
# Last edited: 29/11/2024
rm(list=ls())
# Loading ----
## Loading packages ----
require(dplyr)
require(utils)
require(ggfortify)
require(reshape2)
require(ggplot2)
require(stringr)
require(gridExtra)
require(cowplot)
require(rstudioapi)
require(NMF)

# str_match : sample_\d{5}
# Options ----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

pca_opts = list(
  targeting = -45,
  fields = "reduced",
  scaling = "scaled",
  color = "absz",
  pca_fields = c("absz", "imagz", "phasez", "realz")
)

opts_multi = list(all_files = TRUE,
                  field = "absz",
                  scaling = "unscaled")

# Experimental constants ----
chunk_duration = 20 # Each chunk represents a 20-second frequency sweep
conc_duration = 300 # Every 5 minutes concentration increases
conc_increase = c(0, 0.1, 0.5, 1, 10, 50, 100) # Concentrations increase to these values in nM

# Support functions ----
# Opening data
clean_data = function(data, deep_clean = T) {
  data$timestamp = as.numeric(data$timestamp)
  data = data[!is.na(data$timestamp), ]
  data$chunk = as.numeric(data$chunk)
  data$size = as.numeric(data$size)
  
  # List of dataframes, one for each chunk
  datalist = split(data, data$chunk)
  names(datalist) = paste0(rep("chunk_", length(datalist)), names(datalist))
  for (i in 1:length(datalist)) {
    # From the source data
    nms = datalist[[i]]$fieldname
    chnk = datalist[[i]]$chunk[1]
    tmstmp = datalist[[i]]$timestamp[1]
    sz = datalist[[i]]$size[1]
    
    # From experimental constants
    time = chnk * chunk_duration
    conc = conc_increase[floor(time / conc_duration) + 1]
    
    data_t = as.data.frame(t(datalist[[i]][, c(5:ncol(data))]))
    data_t = cbind(
      rep(chnk, nrow(data_t)),
      rep(tmstmp, nrow(data_t)),
      rep(sz, nrow(data_t)),
      rep(time, nrow(data_t)),
      rep(conc, nrow(data_t)),
      data_t
    )
    colnames(data_t) = c("chunk", "timestamp", "size", "time", "concentration", nms)
    
    datalist[[i]] = data_t
  }
  
  data_join = datalist[[1]]
  if (length(datalist)>1){
    for (i in 2:length(datalist)) {
      data_join = rbind(data_join, datalist[[i]])
    }
  }
  rownames(data_join) = NULL
  
  if (deep_clean) {
    col_nms = colnames(data_join)
    # Whitelist specific columns
    white_cols = c("frequency")
    
    # Delete specific columns
    manu_cols = c(
      "settling",
      "bandwidth",
      "grid",
      "timestamp",
      "tc",
      "tcmeas",
      "count",
      "flags",
      "nexttimestamp",
      "settimestamp",
      "param0",
      "param0pwr",
      "param0stddev",
      "param1",
      "param1pwr",
      "param1stddev"
    )
    
    # Delete power columns for redundancy
    pwr_cols = col_nms[grepl('pwr', col_nms)]
    
    # Delete columns whose standard deviation is 0
    std_cols = col_nms[grepl('stddev', col_nms)]
    std_del = c()
    for (k in 1:length(std_cols)) {
      if (((sum(data_join[, std_cols[k]] == 0, na.rm = T) / (nrow(data_join))) >= .9)){
        std_del = c(std_del, std_cols[k])
      }
    }
    
    std_cols = col_nms[grepl(paste(gsub("stddev", "", std_del), collapse = "|"), col_nms)]
    
    # Join together
    del_cols = c(manu_cols, pwr_cols, std_cols)
    del_cols = del_cols[!duplicated(del_cols)]
    
    # Remove whitelisted columns
    for (i in 1:length(white_cols)) {
      del_cols = del_cols[!(del_cols == white_cols[i])]
    }
    
    data_join = data_join %>% select(-one_of(del_cols))
  }
  
  # Reorder columns
  col_order = c(
    "chunk",
    "size",
    "time",
    "concentration",
    "frequency"
  )
  col_order = c(col_order, colnames(data_join)[!(colnames(data_join) %in% col_order)])
  
  data_join = data_join[, col_order]
  
  return(data_join)
  
  rm(data_t, i, nms, sz, tmstmp, chnk, time, conc, col_order)
}

read_test = function(file = NULL){
  # Open file
  if(is.null(file)){
    file = file.choose()
  }
  
  # Decode file name and folder structure, assumed to be: /Round/Analyte/(functionalization)_(ChipID).(test)
  round = basename(dirname(dirname(file)))
  analyte = basename(dirname(file))
  funct = strsplit(basename(file), split = "[_|.]")[[1]][1]
  chipID = strsplit(basename(file), split = "[_|.]")[[1]][2]
  test = strsplit(basename(file), split = "[_|.]")[[1]][3]
  
  data_raw = read.csv(file, sep = ";")
  
  # Cleanup
  data_join = clean_data(data_raw, deep_clean = T)
  data_join = extract_freq(data_join, target_phase = -45)
  
  data_join$Analyte = analyte
  data_join$Round = round
  data_join$Functionalization = funct
  data_join$ChipID = chipID
  data_join$Test = test
  
  data_join = data_join %>% select(Functionalization, Analyte, Round, ChipID, Test, everything())
  
  data_join$Analyte = as.factor(data_join$Analyte)
  data_join$Round = as.factor(data_join$Round)
  data_join$Functionalization = as.factor(data_join$Functionalization)
  data_join$ChipID = as.factor(data_join$ChipID)
  data_join$Test = as.factor(data_join$Test)
  
  
  return(data_join)
}

read_functionalization = function(folder = NULL, opts_multi){
  if(is.null(folder)){
    folder = choose.dir()
  }
  data_filenames = list.files(folder, pattern = "*.csv", full.names = T, recursive = T)
  
  if(length(data_filenames)==0){
    stop("No data files found in directory. Check that data files have '_data_' on the name, or pick a different folder")
  }
  
  data_mlist = list()
  for(i in 1:length(data_filenames)){
    data_join = read_test(data_filenames[i])
    
    data_mlist[[i]] = data_join
    names(data_mlist)[i] = paste(data_join$Analyte[1], data_join$ChipID[1], data_join$Test[1], collapse = "_")
    names(data_mlist) = make.names(names(data_mlist), unique = T)
  }
  return(data_mlist)
}

# Cleaning data
extract_freq = function(dataframe, target_phase = -45) {
  phase_deg = dataframe$phasez[dataframe$chunk == 0] * 180 / pi
  phase_close = order(abs(phase_deg - target_phase))
  
  freqs = dataframe$frequency[dataframe$chunk == 0]
  target_freq = freqs[phase_close[1]]
  
  data_target = dataframe[dataframe$frequency == target_freq, ]
  
  return(data_target)
}

frequency_sweep = function(df, phases, pca_opts){
  for (j in 1:length(phases)){
    data_temp = extract_freq(df, phases[j])
    data_temp = data_temp[,c("Analyte", "Test", "chunk", "time", "concentration", pca_opts$pca_fields)]
    vars = colnames(data_temp)[colnames(data_temp) %in% pca_opts$pca_fields]
    colnames(data_temp)[colnames(data_temp) %in% pca_opts$pca_fields] = paste0(vars, rep(phases[j], length(vars)))
    
    # If NAs from recording are carried over, delete them here
    data_temp = data_temp[!is.na(data_temp$chunk),]
    
    if(j == 1){
      data_mfreq = data_temp
    }else{
      data_mfreq = merge(data_mfreq, data_temp)
    }
  }
  data_freq = data_mfreq[order(data_mfreq$chunk),]
  
  return(data_freq)
}

# Data analysis
make_PCA_df = function(df, pca_opts){
  # Remove NAs
  for (i in 1:length(pca_opts$pca_fields)){
    data_pca = df[!is.na(df[,pca_opts$pca_fields[i]]),]
  }
  
  # Z-score normalization just in case
  if (pca_opts$scaling == "scaled") {
    for (i in length(pca_opts$pca_fields)) {
      data_pca[, pca_opts$pca_fields[i]] = (data_pca[, pca_opts$pca_fields[i]] - mean(data_pca[, pca_opts$pca_fields[i]])) / sd(data_pca[, pca_opts$pca_fields[i]])
    }
  }
  
  # Take in extraneous fields left behind
  # ignored_fields = c("chunk", "size", "time", "concentration", "frequency")
  if (pca_opts$fields == "reduced") {
    pca = prcomp(data_pca[, pca_opts$pca_fields])
  } else{
    pca = prcomp(data_pca[, -c(1:8)])
  }
  
  return(list(dataframe = data_pca, PCA = pca))
}

# Plot generation
cutie_layer = function() {
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = 'bottom',
    legend.key.width = unit(3, 'lines')
  )
}

customggsave = function(plot, upscale = 1.5, save_path = '', name = NULL) {
  save_path = paste0('./Plots', save_path)
  if (is.null(name)) {
    name = deparse(substitute(plot))
  }
  ggsave(
    paste0(name, ".png"),
    plot = plot,
    device = 'png',
    width = round(1920 * upscale),
    height = round(1080 * upscale),
    units = 'px',
    path = save_path
  )
}

# Data processing ----
# With functions
data_clean = read_test("./Data/JoeRound/APTES/Acetaldehyde/APTES_Joe.1.csv")

data_mlist = read_functionalization("./Data/JoeRound/APTES/")
data_multi = data_mlist[[1]]
for (i in 2:length(data_mlist)){
  data_multi = rbind(data_multi, data_mlist[[i]])
}

pca_list = make_PCA_df(data_multi, pca_opts) # Change data_multi for data_clean for a single analyte

# Make PCA plots
pca_loadings = signif(-1 * pca_list$PCA$rotation, digits = 3)
pca_obj = pca_list$PCA$x
plottype = c(paste(c(
  "Target:", "Fields:", "Scaling:", "Color:", "Label:"
), pca_opts[c("targeting", "fields", "scaling", "color", "text")]))

pca_plot = ggplot(pca_obj, aes(x = PC1, y = PC2, color = pca_list$dataframe$Analyte)) +
  geom_point(size = 2, alpha = .4, position = "jitter") +
  labs(x = paste0('PC1', ' (', (summary(pca_list$PCA)$importance[2, 1]), ')'),
       y = paste0('PC2', ' (', (summary(pca_list$PCA)$importance[2, 2]), ')'),
       color = "Analyte") +
  ggtitle(paste0("Multianalyte PCA: ", pca_list$dataframe$Functionalization[1], " functionalization")) +
  theme_bw() + cutie_layer()

pca_plot_tab = plot_grid(pca_plot, tableGrob(pca_loadings[, c("PC1", "PC2")]), rel_widths = c(3, 1)) +
  theme_bw()

# customggsave(pca_plot_tab,
#                upscale = 2,
#                save_path = paste0("/",functionalization),
#                name = paste(c("PCAmulti", pca_opts[c("targeting", "scaling", "color", "text")]), collapse = "_"))