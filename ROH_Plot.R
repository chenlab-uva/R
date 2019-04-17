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
segments <- subset(segments, select = c(ID, Chr, StartMB, StopMB))

# roh file
roh <- read.table(seg_name, header = TRUE)
# 3 rd degree. Only draw the roh plots for samples with a high F_ROH (1/2^4.5)
roh_id <- roh$ID[roh$F_ROH > 1/(2^4.5) ]

# generate plots
postscript(paste(prefix, "roh", "rplots.ps", sep = "_"), paper="letter", horizontal = T)
for (id in roh_id){
  #' get data for 1 individual & plot
  k <- subset(segments, ID == id)
  if(nrow(k) > 0) {
    theme_set(theme_bw(base_size = 18))
    f_roh <- roh[roh$ID==id,"F_ROH"]
#### Is there any way to start with 0 MB, i.e., getting rid of the gray margin on the left?    
    g <- ggplot() +
      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), fill = 'white', color = "black", size = 0.85) + 
      geom_rect(data = k, aes(xmin = StartMB, xmax = StopMB, ymin = 0, ymax = 0.9), fill = "red") + 
      geom_rect(data = all_seg, aes(xmin = StartMB, xmax = StopMB, ymin = 0, max = 0.9), color = "black", alpha = 0, size = 0.85) +
      facet_grid(Chr ~ .) +
      labs(x = "Position (MB)", y = "", title = bquote(paste('Run of Homozygosity for ', .(id), ' in ', .(prefix), ' (F'['ROH']*' = ', .(f_roh), ')'))) +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
    print(g)
  }
}
dev.off()
