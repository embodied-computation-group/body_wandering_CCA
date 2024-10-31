

# chordplot function (saves png)
chordplot <- function(cca_filename, threshold_percentile, chord2plot, savepath, region_labels, network_labels) {
  
  # install packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(readr,igraph,circlize,R.matlab,data.table, reshape2, tidyverse, ComplexHeatmap)
  
  ##############################################################################
  # Format cca results (thresholded adjacency list with regions & networks)
  ##############################################################################  
  # read in cca data  
  cca_results <- readMat(cca_filename) 
  
  Aload = cca_results$weightX 
  
  # rewrap function in R
  rewrap <- function(V) {
    # check for pacman package and install if not found
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(signal)
    nE <- length(V)
    Nn <- max(Re(roots(c(1, -1, -2*nE))))
    idx <- upper.tri(matrix(data=1,nrow=Nn,ncol=Nn), diag = FALSE)
    tmp <- matrix(data=0,nrow=Nn,ncol=Nn)
    tmp[idx] <- V
    W <- tmp + t(tmp)
    return(W)
  }
  
  atlas_size = 216
  Aload_matrix <- array(numeric(),c(atlas_size,atlas_size)) 
  
  Aload_matrix[ , ] <- rewrap(Aload[ ]);
  
  # threshold matrix 
  threshold = quantile(abs(Aload_matrix), threshold_percentile)
  idx = Aload_matrix[ , ] < threshold & Aload_matrix[ , ] > -threshold 
  Aload_matrix[ , ][idx] = NA; 
  
  # select which connectome results to plot (& remove lower triangle)
  connectivity = Aload_matrix 
  connectivity[lower.tri(connectivity)] <- NA  
  
  # add labels to connectivity matrix
  rownames(connectivity) = make.names(region_labels, unique=TRUE) # deal with duplicates atlas_labels by adding .1 .2 etc
  colnames(connectivity) = make.names(region_labels, unique=TRUE) # deal with duplicates atlas_labels by adding .1 .2 etc
  
  # matrix with removed NA's (values below threshold)
  connectivity_small <- connectivity[rowSums(is.na(connectivity))<ncol(connectivity),colSums(is.na(connectivity))<nrow(connectivity)] 
  
  # Make adjacency list with regions & networks
  my_cor_df <- melt(connectivity_small)
  # create adjacency list without NA
  my_adj_list <- my_cor_df %>% drop_na(value)
  names(my_adj_list) <- c('reg_A', 'reg_B', 'weight')
  
  # Network labels
  colnames(connectivity) = network_labels
  rownames(connectivity) = network_labels
  connectivity_small <- connectivity[rowSums(is.na(connectivity))<ncol(connectivity),colSums(is.na(connectivity))<nrow(connectivity)] # removerows/cols with all NAs
  cor_df_networks <- melt(connectivity_small)
  # create adjacency list without NA
  adj_list_networks <- cor_df_networks %>% drop_na(value)
  
  my_adj_list <- cbind(my_adj_list, adj_list_networks[ ,1:2])
  colnames(my_adj_list)[4] <- "mod_A"
  colnames(my_adj_list)[5] <- "mod_B"
  
  
  ##############################################################################
  # set chordplot details & plot!
  ##############################################################################  
  circos.clear
  
  # define network groups
  modul = c(structure(as.character(my_adj_list$mod_A), names=as.character(my_adj_list$reg_A)), structure(as.character(my_adj_list$mod_B), names=as.character(my_adj_list$reg_B)))
  modul = modul[!duplicated(names(modul))]
  modul = modul[order(modul, names(modul))]
  
  # define colour of networks
  mod_colors=c('maroon', 'darkslateblue', 'steelblue1', 'mediumpurple1', 'lightpink', 'slategray1', 'darkslategray1', 'hotpink1')
  modul_color = structure(mod_colors, names = unique(modul)) #mod_colors[1:8]
  
  # ensure network colours are the same in each chordplot
  for (idx in 1:length(modul_color)) {
    if (isTRUE(names(modul_color[idx]) == "SalVentAttn")) {
      modul_color[idx] = 'mediumpurple1' 
    } else if  (isTRUE(names(modul_color[idx]) == "Cont")) {
      modul_color[idx] = 'darkslateblue'
    } else if  (isTRUE(names(modul_color[idx]) == "DorsAttn")) {
      modul_color[idx] = 'darkslategray1'
    } else if  (isTRUE(names(modul_color[idx]) == "Subcortical")) {
      modul_color[idx] = 'lightpink'
    } else if  (isTRUE(names(modul_color[idx]) == "Vis")) {
      modul_color[idx] = 'maroon'
    } else if  (isTRUE(names(modul_color[idx]) == "Limbic")) {
      modul_color[idx] = 'slategray1'
    } else if  (isTRUE(names(modul_color[idx]) == "Default")) {
      modul_color[idx] = 'steelblue1'
    } else if  (isTRUE(names(modul_color[idx]) == "SomMot")) {
      modul_color[idx] = 'hotpink1' 
    }
  }
  #mod_colors = unique(modul_color)
  reg_color = structure(modul_color[as.factor(modul)], names = names(modul))
  
  # colours of links (positive & negative)
  cols = c('dodgerblue', 'brown1')
  
  # define gap between networks
  gap.after = do.call("c", lapply(table(modul), function(i) c(rep(2, i-1), 10)))
  circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0), canvas.ylim=c(-1.5,1.5)) # adjust gaps between regions & size of canvas (so text can fit)
  
  # SAVE           
  png(sprintf("%schordplot_CCAloadings_%.3fthresh_%.2fpercentile.png", savepath, round(threshold, 3), round(threshold_percentile, 3)), res = 150, height = 17, width = 17, units = "cm")
  
  # text size
  par(cex = 0.5) 
  
  # plot chordDiagram - with positive & negative links
  chordDiagram(my_adj_list, annotationTrack = "grid", grid.col = reg_color, order = names(modul), directional = 1, preAllocateTracks =  list(list(track.height = 0.02)), col=cols[(sign(my_adj_list$weight)+3)/2])
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  
  
  # add legend
  lgd_networks = Legend(labels = c(unique(modul)), legend_gp = gpar(fill = unique(reg_color)), title = "Network")
  lgd_list_vertical = packLegend(lgd_networks)
  draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
  
  dev.off()
  
  return(list(my_adj_list=my_adj_list, threshold=threshold))
}

