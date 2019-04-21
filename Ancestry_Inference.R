print(date())

# check packages
required_packages <- c("e1071", "doParallel")
error_msg <- paste("Please install the following packages:", paste(required_packages, collapse = ", "))
if( !all(required_packages %in% installed.packages()[,'Package']) ) stop(error_msg)
suppressMessages(library(e1071))
suppressMessages(library(doParallel))

# input setting
args = commandArgs(TRUE)
if( length(args) != 3 ) stop("please provide three arguments (pc file, population file and prefix)")

pc_name <- args[1]
phe_name <- args[2]
prefix <- args[3]

pc <- read.table(pc_name, header = TRUE)
phe <- read.table(phe_name, header = TRUE)

print(paste("Prepare the PC file and the reference file, starts at ",date(), sep = " "))
print(paste("Apply 10 PCs to the analysis", sep = " "))
print(paste("Choose linear", sep = " "))

pop <- phe[, c("IID", "Population")]
train.data <- pc[pc$AFF == 1, c(2, 7:16)]
train.phe <- merge(train.data, pop, by = "IID")
test.data <- pc[pc$AFF == 2, c(1, 2, 7:16)]

train.x <- train.phe[, !colnames(train.phe) %in% c("Population", "IID")]
train.y <- train.phe[,"Population"]

numCores <- detectCores()
registerDoParallel(cores = round((numCores/2)))

print(paste("Assign ", round((numCores/2)), " cores for the grid search, starts at ", date(), sep = ""))
print("Grid search with a wide range")

# Linear Only
tuneresult <- foreach(cost = 2^(seq(-10, 10, by = 0.5)), .combine = rbind) %dopar% {
  set.seed(123)
  mod = tune(svm, train.x, as.factor(train.y), kernel = "linear", cost = cost, probability = TRUE, 
             tunecontrol = tune.control(cross=5))
  performance <- mod$performances[, c("error", "dispersion")]
  data.frame(cost = cost, performance = performance)
}
best.cost <- tuneresult[which.min(tuneresult$performance.error), "cost"]

print(paste("The best cost is ", round(best.cost, 6), " after the wide grid search", sep=""))
print(paste("Grid search with a small range, starts at", date()))

more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)
tune.more <- foreach(cost = more.cost, .combine = rbind) %dopar% {
  set.seed(123)
  mod = tune(svm, train.x, as.factor(train.y), kernel = "linear", 
             cost = cost, probability = TRUE, tunecontrol = tune.control(cross=5))
  performance <- mod$performances[, c("error", "dispersion")]
  data.frame(cost = cost, performance = performance)
}
best.cost <- tune.more[which.min(tune.more$performance.error), "cost"]
print(paste("Grid search with a small range, ends at", date()))
print(paste("The best cost is ", round(best.cost, 6), " after the small grid search", sep=""))
# final model
set.seed(123)
mymod <- svm(train.x, as.factor(train.y), cost = best.cost, kernel = "linear", probability=TRUE)

print("Predict ancestry information")
print(date())

pred.pop <- predict(mymod, test.data[, !colnames(test.data) %in%c("FID","IID")], probability=TRUE)
test.data$PRED <- pred.pop
class.prob <- attr(pred.pop, "probabilities")


print(paste("Prepare the summary file, starts at", date()))

orders <- t(apply(class.prob, 1, function(x) order(x,decreasing =T)))
orders.class <- t(apply(orders, 1, function(x) colnames(class.prob)[x]))
orders.probs <- t(sapply(1:nrow(class.prob), function(x) class.prob[x, orders[x,]]))

# > 0.65
check.cumsum <- t(apply(orders.probs, 1, cumsum))
temp <- apply(check.cumsum, 1, function(x) which(x > 0.65)[1])

PRED_CLASS <- sapply(1:length(temp), function (x) paste(orders.class[x, 1:as.numeric(temp[x])], collapse = ";"))
PRED_PROB <- sapply(1:length(temp), function (x) paste(round(orders.probs[x, 1:as.numeric(temp[x])], 3), collapse = ";"))


pred.out <- cbind(test.data[, c("FID", "IID", "PC1", "PC2")], PRED_CLASS, PRED_PROB, orders.class[, 1:2], orders.probs[, 1:2])
colnames(pred.out)[7:10] <- c("FST", "SEC", "FST_PROB", "SEC_PROB")

print(paste("summary file is ready ", date()))
write.table(pred.out, paste(prefix, "_linear_", 10, "PCs.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("Results are saved to", paste(prefix, "_", "linear", "_", 10, "PCs.txt", sep = ""), sep = " "))

print("Generate ggplots")
pred.out$PRED_CLASS <- as.character(pred.out$PRED_CLASS)
pred.out$PRED_CLASS[pred.out$FST_PROB <= 0.65] <- "UNCERTAIN"
Palette <- c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928","#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F","#CAB2D6","#FFFF99","999999")
groups <- length(unique(pred.out$PRED_CLASS))
cPalette <- Palette[c(1:(groups - 1), 10)]

if(!require("ggplot2")) {
  pdf(paste(prefix, "_", "linear", "_", 10, "PCs.pdf", sep = ""))
  pred.out$PRED_CLASS <- as.factor(pred.out$PRED_CLASS)
  plot(pred.out$PC1, pred.out$PC2, col = cPalette[pred.out$PRED_CLASS], xlab = "PC1", ylab = "PC2", pch = 16)
  legend("topright", title = "Predicted", legend = unique(pred.out$PRED_CLASS), col = unique(cPalette[pred.out$PRED_CLASS]), pch = 16, cex = 1)
  dev.off()
  pdf("Reference.pdf")
  plot(train.phe$PC1, train.phe$PC2, col = cPalette[train.phe$Population], xlab = "PC1", ylab = "PC2", pch = 16)
  legend("topright", title = "Reported", legend = unique(train.phe$Population), 
         col = unique(cPalette[train.phe$Population]), pch = 16, cex = 1)
  print(paste("Reference.pdf and ", prefix, "_", "linear", "_", 10, "PCs.pdf are generated.", sep = ""))
  dev.off()
  print("Done")
} else {
  pdf(paste(prefix, "_", "linear", "_", 10, "PCs.pdf", sep = ""))
  p <- ggplot(pred.out, aes(x = PC1, y = PC2, color = PRED_CLASS, label = IID))
  p <- p + geom_point() + scale_colour_manual(values = cPalette)
  print(p)
  dev.off()
  
  pdf("Reference.pdf")
  p <- ggplot(train.phe, aes(x = PC1, y = PC2, color = Population, label = IID))
  p <- p + geom_point() + scale_colour_manual(values = cPalette)
  print(p)
  dev.off()
  print(paste("Reference.pdf and ", prefix, "_", "linear", "_", 10, "PCs.pdf are generated.", sep = ""))
  print("Done")
}
print(date())
