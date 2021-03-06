#' check for packages
if(!require("ggplot2")) stop("Please install ggplot2 packages")
suppressMessages(library(ggplot2))

#' process command line args
args = commandArgs(TRUE)
if( length(args) != 1 ) stop("please provide one argument (prefix)")

prefix <- args[1]

#' file names
seg_name <- paste0(prefix, ".roh")
segments_name <- paste0(prefix, ".rohseg.gz")
all_seg_name <- paste0(prefix, "allsegs.txt")

if( !(file.exists(seg_name) & file.exists(segments_name) & file.exists(all_seg_name)) ) stop("Missing RoH files")

# segments considered
all_seg <- read.table(all_seg_name, header = TRUE)
all_seg <- subset(all_seg, select = c(Chr, StartMB, StopMB))


# roh segments gz file
segments <- read.table(segments_name, header = TRUE)
segments <- subset(segments, select = c(FID, ID, Chr, StartMB, StopMB)) 

# roh file
roh <- read.table(seg_name, header = TRUE)
# 3 rd degree. Only draw the roh plots for samples with a high F_ROH (1/2^4.5)
roh_info <- roh[roh$F_ROH > 1/(2^4.5), c("FID","ID","F_ROH")] 
roh_info$FID <- as.character(roh_info$FID)
roh_info$ID <- as.character(roh_info$ID)

# generate plots
postscript(paste(prefix, "roh", "rplots.ps", sep = "_"), paper="letter", horizontal = T)

#### This loop is too slow, 0.5 second per iteration/sample. 
#### We need to consider parallel library (if all(required_packages %in% installed.packages()[,'Package'])) 
for (i in 1:nrow(roh_info)){
#### This loop makes roh_run_1 function, with g as return  
  #' get data for 1 individual & plot
  k <- subset(segments, FID == roh_info[i,1] & ID == roh_info[i,2])
  if(nrow(k) > 0) {
    theme_set(theme_bw(base_size = 18))
    f_roh <- roh_info[i,"F_ROH"]
    fid <- as.character(k[1,1])
    id <- as.character(k[1,2])
    g <- ggplot() +
      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), fill = 'white', color = "black", size = 0.85) + 
      geom_rect(data = k, aes(xmin = StartMB, xmax = StopMB, ymin = 0, ymax = 0.9), fill = "red") + 
      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), color = "black", alpha = 0, size = 0.85) +
      facet_grid(Chr ~ .) +
      scale_x_continuous(expand  = c(0, 0), limits = c(0, NA)) + 
      labs(x = "Position (MB)", y = "", title = bquote(paste('Run of Homozygosity for ', .(id), ' from FAM ', .(fid), ' in ', .(prefix), ' (F'['ROH']*' = ', .(f_roh), ')'))) +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title=element_text(size = 18) 
        # rel() is related to the parent. I try size = 18 and size = rel(1), the font sizes are the same.
      )
    print(g)
  }
}
#### Here is the idea of code for parallel computing:
#### if(all(c("parallel") %in% installed.packages()[,'Package'])){
####   cl <- makeCluster(detectCores()/2)
####   clusterExport(cl = cl)
####   glist <- parLapply(cl, 1:nrow(roh_info), roh_run_1)
####   stopCluster(cl)
#### }else{
####   glist <- lapply(cl, 1:nrow(roh_info), roh_run_1)
#### }
#### for(i in 1:nrow(roh_info)) g <- ggplot + glist[i]
#### print(g)
dev.off()
