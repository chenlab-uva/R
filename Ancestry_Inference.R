library(optparse)
library(e1071)
library(parallel)
library(iterators)
library(foreach)
library(doParallel)
library(ggplot2)
# input setting
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="PC file", metavar="character"),
  make_option(c("-p", "--phefile"), type="character", default=NULL,
              help="population dataset file name", metavar="character"),
  make_option(c("-n","--numberpc"), type="integer", default=10,
              help="number of PCs to be selected [default %default]", metavar="integer"),
  make_option(c("-k", "--kernel"), type="character", default="linear",
              help = "Kernel Option [default \"%default\"]", metavar="character"),
  make_option(c("--prefix"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$file) | is.null(opt$phe) | is.null(opt$prefix)){
  print_help(opt_parser)
  stop("PC file, population file and prefix  must be supplied", call.=FALSE)
}


pc <- read.table(opt$file, header = TRUE)
phe <- read.table(opt$phefile, header = TRUE)
number <- opt$numberpc


print(paste("Prepare the PC file and the reference file", sep=" "))
print(paste("Apply",number,"PCs to the analysis", sep=" "))
print(paste("Choose",opt$kernel,sep = " "))

pop <- phe[, c("IID", "population")]
train.data <- pc[pc$AFF == 1, c(2, 7:(number + 6))]
train.phe <- merge(train.data, pop, by = "IID")
test.data <- pc[pc$AFF == 2, c(1,2, 7:(number + 6))]

train.x <- train.phe[, !colnames(train.phe) %in% c("population", "IID")]
train.y <- train.phe$population

numCores <- detectCores()
registerDoParallel(cores = round((numCores/2)))

print(paste("Assign ",round((numCores/2))," cores for the grid search", sep = ""))
print("Grid search with a wide range")

if (opt$kernel == "linear") {
  # For Linear
  tuneresult <- foreach(cost = 2^(seq(-10, 10, by = 0.5)), .combine = rbind) %dopar% {
    set.seed(123)
    mod = tune(svm, train.x, as.factor(train.y), kernel = opt$kernel, 
               cost = cost, probability = TRUE)
    performance <- mod$performances[, c("error", "dispersion")]
    data.frame(cost = cost, performance = performance)
  }
  best.setting <- tuneresult[tuneresult$performance.error == min(tuneresult$performance.error), ]
  best.cost <- min(best.setting$cost)
  
  print(paste("The best cost is ", round(best.cost,6), " after the wide grid search", sep=""))
  print("Grid search with a small range")
  
  more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)
  tune.more <- foreach(cost = more.cost, .combine = rbind) %dopar% {
    set.seed(123)
    mod = tune(svm, train.x, as.factor(train.y), kernel = opt$kernel, 
               cost = cost, probability = TRUE)
    performance <- mod$performances[, c("error", "dispersion")]
    data.frame(cost = cost, performance = performance)
  }
  best.setting <- tune.more[tune.more$performance.error == min(tune.more$performance.error), ]
  best.cost <- min(best.setting$cost)
  print(paste("The best cost is ", round(best.cost,6), " after the small grid search", sep=""))
  # my final model
  set.seed(123)
  mymod <- svm(train.x, as.factor(train.y), cost = best.cost, kernel = opt$kernel,probability=TRUE)
} else {
   # For radial
  tuneresult <- foreach(gamma = 2^(seq(-15, 3, by = 0.5)), .combine = rbind) %dopar% {
    registerDoParallel(cores = 1)
    foreach(cost = 2^(seq(-10, 10, by = 0.5)), .combine = rbind) %dopar% 
    {
      set.seed(123)
      mod = tune(svm, train.x, as.factor(train.y), 
                 kernel = opt$kernel, ranges = list(cost = cost, gamma = gamma), 
                 probability = TRUE)
      performance <- mod$performances[, c("error", "dispersion")]
      data.frame(cost = cost, gamma = gamma, performance = performance)
    }
  }
  # write.table(tuneresult,"radial_check_tuneresult_1.txt",quote=FALSE, row.names=FALSE)
  best.setting <- tuneresult[tuneresult$performance.error == min(tuneresult$performance.error), ]
  best.cost <- best.setting$cost[1]
  best.gamma <- best.setting$gamma[1]
  
  print(paste("The best cost is ", round(best.cost,6), " after the wide grid search", sep=""))
  print(paste("The best gamma is ", round(best.gamma,6), " after the wide grid search", sep=""))
  print("Grid search with a small range")
  more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)
  more.gamma <- 2^seq(log2(best.gamma) - 0.5, log2(best.gamma) + 0.5, 
                      by = 0.05)
 registerDoParallel(cores = round((numCores/2)))  
 tune.more <- foreach(gamma= more.gamma, .combine = rbind) %dopar% {
   registerDoParallel(cores = 1)
   foreach(cost = more.cost, .combine = rbind) %dopar% {
      set.seed(123)
      mod = tune(svm, train.x, as.factor(train.y), kernel = opt$kernel, 
                 ranges = list(cost = cost, gamma = gamma), probability = TRUE)
      performance <- mod$performances[, c("error", "dispersion")]
      data.frame(cost = cost, gamma = gamma, performance = performance)
    }
  }
  
  best.setting <- tune.more[tune.more$performance.error == min(tune.more$performance.error), ]
  best.cost <- best.setting$cost[1]
  best.gamma <- best.setting$gamma[1]
  print(paste("The best cost is ", round(best.cost,6), " after the wide grid search", sep=""))
  print(paste("The best gamma is ", round(best.gamma,6), " after the wide grid search", sep=""))
  set.seed(123)
  mymod <- svm(train.x, as.factor(train.y), cost = best.cost, gamma = best.gamma, 
               kernel = opt$kernel,probability=TRUE)
}


print("Predict the ancestry information")

pred.pop <- predict(mymod, test.data[, !colnames(test.data) %in%c("FID","IID")],probability=TRUE)
test.data$PRED <- pred.pop
class.prob <- attr(pred.pop, "probabilities")


print("Prepare the summary file")

maxn <- function(n) function(x) order(x,decreasing = TRUE)[n]
fst_class.index <- apply(class.prob,1,maxn(1))
test.data$FST <- colnames(class.prob)[fst_class.index]
test.data$FST_PROB <- round(apply(class.prob,1,function(x) x[maxn(1)(x)]),3)

sec_class.index <- apply(class.prob,1,maxn(2))
test.data$SEC <- colnames(class.prob)[sec_class.index]
test.data$SEC_PROB <- round(apply(class.prob,1,function(x) x[maxn(2)(x)]),3)

third_class.index <- apply(class.prob,1,maxn(3))
test.data$THIRD <- colnames(class.prob)[third_class.index]
test.data$THIRD_PROB <- round(apply(class.prob,1,function(x) x[maxn(3)(x)]),3)

fourth_class.index <- apply(class.prob,1,maxn(4))
test.data$FOURTH <- colnames(class.prob)[fourth_class.index]
test.data$FOURTH_PROB <- round(apply(class.prob,1,function(x) x[maxn(4)(x)]),3)


pre.number.part <- test.data[,c("FST_PROB","SEC_PROB","THIRD_PROB","FOURTH_PROB")]
pre.pop.part <- test.data[,c("FST","SEC","THIRD","FOURTH")]
# convert it to character
pre.pop.part <- data.frame(lapply(pre.pop.part, as.character), stringsAsFactors=FALSE)
# For each sample, find the first several classes with prob > 0.65
check.number <- t(apply(pre.number.part,1,cumsum))
temp <- apply(check.number, 1, function(x) which(x > 0.65)[1])

#
for (i in 1:length(temp)){
  if (temp[i]==1) {
    test.data$PRED_CLASS[i] <- pre.pop.part[i,1]
    test.data$PRED_PROB[i] <- pre.number.part[i,1]
  }
  else {
    test.data$PRED_CLASS[i] <- paste(pre.pop.part[i,1:as.numeric(temp[i])],collapse=";")
    test.data$PRED_PROB[i] <- paste(pre.number.part[i,1:as.numeric(temp[i])],collapse=";")
  }
}


pred.out <- test.data[,c("FID","IID","PC1","PC2","FST","FST_PROB","SEC","SEC_PROB","PRED_CLASS","PRED_PROB")]
# 
write.table(pred.out, paste(opt$prefix,"_",opt$kernel,"_",opt$numberpc,"PCs.txt", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Results are saved to",paste(opt$prefix,"_",opt$kernel,"_",opt$numberpc,"PCs.txt", sep=""), sep=" "))

print("Generate ggplots")

uncertain.index <- (pred.out$FST_PROB <= 0.65)
pred.out$PRED_CLASS <- as.character(pred.out$PRED_CLASS)
pred.out$PRED_CLASS[uncertain.index] <- "UNCERTAIN"

Palette <- c("#0072B2", "#E69F00", "#CC79A7", "#009E73", "#F0E442", "#D55E00", 
    "#56B4E9", "#556b2f", "#0000FF", "#999999")
groups <- length(unique(pred.out$PRED_CLASS))
cbPalette <- Palette[c(1:(groups - 1), 10)]


pdf(paste(opt$prefix, "_", opt$kernel, "_", opt$numberpc, "PCs.pdf", sep = ""))
p <- ggplot(pred.out, aes(x = PC1, y = PC2, color = PRED_CLASS, label = IID))
p <- p + geom_point() + scale_colour_manual(values = cbPalette)
print(p)
dev.off()

# train.phe
pdf("Reference.pdf")
p <- ggplot(train.phe, aes(x = PC1, y = PC2, color = population, label = IID))
p <- p + geom_point() + scale_colour_manual(values = cbPalette)
print(p)
dev.off()
print(paste("Reference.pdf and ", opt$prefix, "_", opt$kernel, "_", opt$numberpc, "PCs.pdf are generated." , sep =""))
print("Done!")
