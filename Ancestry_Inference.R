print(date())

# check packages
required_packages <- c("e1071")
error_msg <- paste("Please install the following packages:", paste(required_packages, collapse = ", "))
if( !all(required_packages %in% installed.packages()[,'Package']) ) stop(error_msg)
suppressMessages(library(e1071))

# input setting
args = commandArgs(TRUE)
if( length(args) != 1 ) stop("please provide one arguments (prefix)")

prefix <- args[1]

pc <- read.table(paste0(prefix,"pc.txt"), header = TRUE)
phe <- read.table(paste0(prefix,"_popref.txt"), header = TRUE)

print(paste0("Prepare the PC file and the reference file, starts at ",date()))
print("Apply 10 PCs to the analysis")
print("Choose linear")

pop <- phe[, c("IID", "Population")]
train.data <- pc[pc$AFF == 1, c(2, 7:16)]
train.phe <- merge(train.data, pop, by = "IID")
test.data <- pc[pc$AFF == 2, c(1, 2, 7:16)]

train.x <- train.phe[, !colnames(train.phe) %in% c("Population", "IID")]
train.y <- train.phe[, "Population"]


# Linear Only
if (require("doParallel")) {
  numCores <- detectCores()
  registerDoParallel(cores = round((numCores/2)))
  tuneresults <- function(cost){
    tuneresult <- foreach(cost = cost, .combine = c) %dopar% 
    {
      set.seed(123)
      mod = tune(svm, train.x, as.factor(train.y), kernel = "linear", cost = cost, probability = TRUE, tunecontrol = tune.control(cross = 5))
      mod$performances[, c("error")]
    }
    best.cost <- cost[which.min(tuneresult)]
    return(best.cost)
  }
} else{
  numCores <- 2
  single.tune <- function(cost) {
    set.seed(123)
    mod = tune(svm, train.x, as.factor(train.y), kernel = "linear", 
               cost = cost, probability = TRUE, tunecontrol = tune.control(cross = 5))
    return(mod$performances[, c("error")])
  }
  tuneresults <- function(cost){
    return(cost[which.min(sapply(cost, single.tune))])
  }
}

print(paste0("Assign ", round((numCores/2)), " cores for the grid search."))
print(paste0("Grid search with a wide range, starts at ", date()))
best.cost <- tuneresults(2^(seq(-10, 10, by = 0.5)))
print(paste0("Grid search with a wide range, ends at ", date()))
print(paste0("The best cost is ", round(best.cost, 6), " after the wide grid search"))
print(paste0("Grid search with a small range, starts at ", date()))
more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)
best.cost <- tuneresults(more.cost)
print(paste0("Grid search with a small range, ends at ", date()))
print(paste0("The best cost is ", round(best.cost, 6), " after the small grid search"))


# final model
set.seed(123)
mymod <- svm(train.x, as.factor(train.y), cost = best.cost, kernel = "linear", probability=TRUE)

print(paste0("Predict ancestry information, starts at ", date()))
pred.pop <- predict(mymod, test.data[, !colnames(test.data) %in%c("FID","IID")], probability=TRUE)
test.data$PRED <- pred.pop
class.prob <- attr(pred.pop, "probabilities")

print(paste0("Prepare the summary file, starts at ", date()))
orders <- t(apply(class.prob, 1, function(x) order(x,decreasing =T)))
orders.class <- t(apply(orders, 1, function(x) colnames(class.prob)[x]))
orders.probs <- t(sapply(1:nrow(class.prob), function(x) class.prob[x, orders[x,]]))

# > 0.65
check.cumsum <- t(apply(orders.probs, 1, cumsum))
temp <- apply(check.cumsum, 1, function(x) which(x > 0.65)[1])

PRED_CLASS <- sapply(1:length(temp), function (x) paste(orders.class[x, 1:as.numeric(temp[x])], collapse = ";"))
PRED_PROB <- sapply(1:length(temp), function (x) paste(round(orders.probs[x, 1:as.numeric(temp[x])], 3), collapse = ";"))

pred.out <- cbind(test.data[, c("FID", "IID", "PC1", "PC2")], PRED_CLASS, PRED_PROB, orders.class[, 1:2], round(orders.probs[, 1:2], 3))
colnames(pred.out)[7:10] <- c("FST", "SEC", "FST_PROB", "SEC_PROB")

print(paste("summary file is ready ", date()))
write.table(pred.out, paste0(prefix, "_InferredAncestry.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Results are saved to", paste0(prefix, "_InferredAncestry.txt")))

print("Generate plots")
pred.out$PRED_CLASS <- as.character(pred.out$PRED_CLASS)
pred.out$PRED_CLASS[pred.out$FST_PROB <= 0.65] <- "UNCERTAIN"
Palette <- c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99","#999999")
groups <- length(unique(pred.out$PRED_CLASS))
cPalette <- Palette[c(1:(groups - 1), 13)]

if(!require("ggplot2")) {
  postscript(paste0(prefix, "_ancestryplot.ps"), paper = "letter", horizontal = T)
  # Predicted
  pred.out$PRED_CLASS <- as.factor(pred.out$PRED_CLASS)
  plot(pred.out$PC1, pred.out$PC2, col = cPalette[pred.out$PRED_CLASS], xlab = "PC1", ylab = "PC2",
       main = paste0("Inferred Populations as Ancestry in ", prefix), pch = 16)
  legend("topright", legend = unique(pred.out$PRED_CLASS), col = unique(cPalette[pred.out$PRED_CLASS]), pch = 16, cex = 1)
  
  # Reported
  plot(train.phe$PC1, train.phe$PC2, col = cPalette[train.phe$Population], xlab = "PC1", ylab = "PC2", 
       main = "Populations in Reference", pch = 16)
  legend("topright", legend = unique(train.phe$Population), 
         col = unique(cPalette[train.phe$Population]), pch = 16, cex = 1)
  dev.off()
  print(paste0(prefix, "_ancestryplot.ps is generated."))
  print("Done")
} else {
  postscript(paste0(prefix, "_ancestryplot.ps"), paper = "letter", horizontal = T)
  p <- ggplot(pred.out, aes(x = PC1, y = PC2, color = PRED_CLASS, label = IID)) 
  p <- p + geom_point() + labs(color = "") +  scale_colour_manual(values = cPalette) + ggtitle(paste0("Inferred Populations as Ancestry in ", prefix))
  print(p)
  p <- ggplot(train.phe, aes(x = PC1, y = PC2, color = Population, label = IID)) 
  p <- p + geom_point() + labs(color = "") + scale_colour_manual(values = cPalette) + ggtitle("Populations in Reference") + labs(fill = "")
  print(p)
  dev.off()
  print(paste0(prefix, "_ancestryplot.ps is generated."))
  print("Done")
}
print(date())