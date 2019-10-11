## king_ancestryplot.R for KING Ancestry plot, by Zhennan Zhu and Wei-Min Chen
library(e1071)
prefix = "king"
pc <- read.table(paste0(prefix, "pc.txt"), header = TRUE)
phe <- read.table(paste0(prefix, "_popref.txt"), header = TRUE)
print(paste("Prepare the PC file and the reference file, starts at ", date()))
pop <- phe[, c("IID", "Population")]
train.data <- pc[pc$AFF == 1, grep("IID|PC", colnames(pc))]
train.phe <- merge(train.data, pop, by = "IID")
test.data <- pc[pc$AFF == 2, grep("FID|IID|PC", colnames(pc))]
train.x <- train.phe[, !colnames(train.phe) %in% c("Population", "IID")]
train.y <- train.phe[, "Population"]
if (require("doParallel", quietly = TRUE)) {
  numCores <- detectCores()
  registerDoParallel(cores = min(round(numCores/2), 41))
  tuneresults <- function(cost) {
    tuneresult <- foreach(cost = cost, .combine = c) %dopar% {
      set.seed(123)
      mod = tune(svm, train.x, as.factor(train.y), kernel = "linear", cost = cost, 
                 probability = TRUE)
      mod$performances[, c("error")]
    }
    best.cost <- cost[which.min(tuneresult)]
    return(best.cost)
  }
} else {
  numCores <- 2
  tuneresults <- function(cost){
    set.seed(123) 
    tune.mod <- tune(svm, train.x, as.factor(train.y), kernel = "linear", ranges=(list(cost=cost)), 
                     probability = TRUE)
    return(tune.mod$best.parameters[1,1])
  }
}
print(paste0("Assign ", min(round(numCores/2), 41), " cores for the grid search."))
print(paste("Grid search with a wide range, starts at", date()))
best.cost <- tuneresults(2^(seq(-10, 10, by = 0.5)))
print(paste("Grid search with a wide range, ends at", date()))
print(paste0("The best cost is ", round(best.cost, 6), " after the wide grid search"))
print(paste("Grid search with a small range, starts at", date()))
more.cost <- 2^seq(log2(best.cost) - 0.5, log2(best.cost) + 0.5, by = 0.05)
best.cost <- tuneresults(more.cost)
print(paste("Grid search with a small range, ends at", date()))
print(paste0("The best cost is ", round(best.cost, 6), " after the small grid search"))
set.seed(123)
mymod <- svm(train.x, as.factor(train.y), cost = best.cost, kernel = "linear", probability = TRUE)
print(paste("Predict ancestry information, start at", date()))
pred.pop <- predict(mymod, test.data[, !colnames(test.data) %in% c("FID", "IID")], probability = TRUE)
test.data$PRED <- pred.pop
class.prob <- attr(pred.pop, "probabilities")
print(paste("Prepare the summary file, starts at", date()))
orders <- t(apply(class.prob, 1, function(x) order(x, decreasing = T)))
orders.class <- t(apply(orders, 1, function(x) colnames(class.prob)[x]))
orders.probs <- t(sapply(1:nrow(class.prob), function(x) class.prob[x, orders[x, ]]))
check.cumsum <- t(apply(orders.probs, 1, cumsum))
temp <- apply(check.cumsum, 1, function(x) which(x > 0.65)[1])
pred.class <- sapply(1:length(temp), function(x) paste(orders.class[x, 1:as.numeric(temp[x])], collapse = ";"))
pred.prob <- sapply(1:length(temp), function(x) paste(round(orders.probs[x, 1:as.numeric(temp[x])], 3), collapse = ";"))
pred.out <- cbind(test.data[, c("FID", "IID", "PC1", "PC2")], pred.class, pred.prob, 
                  orders.class[, 1], orders.class[, 2], round(orders.probs[, 1], 3), round(orders.probs[, 2],3))
colnames(pred.out)[5:10] <- c("Ancestry", "Pr_Anc", "Anc_1st", "Anc_2nd", "Pr_1st", "Pr_2nd")
print(paste("summary file is ready ", date()))
write.table(pred.out, paste0(prefix, "_InferredAncestry.txt"), quote = FALSE, row.names = FALSE)
print(paste("Results are saved to", paste0(prefix, "_InferredAncestry.txt")))
print("Generate plots")
pred.out$Ancestry <- as.character(pred.out$Ancestry)
pred.out$Ancestry[pred.out$Pr_1st <= 0.65] <- ">1 Pop"
Palette <- c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3", 
             "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#999999")

train.groups <- unique(train.phe$Population)
pred.colors <- rep(Palette[13], nrow(pred.out))

for (i in 1:length(train.groups)) {
  pred.colors[pred.out$Ancestry == train.groups[i]] <- Palette[i]
}

train.colors <- rep(0, nrow(train.phe))
for (i in 1:length(train.groups)) {
  train.colors[train.phe$Population == train.groups[i]] <- Palette[i]
}

# set xlim and ylim
x.adjust <- (max(train.phe$PC1, pred.out$PC1) - min(train.phe$PC1, pred.out$PC1))/10
x.low <- min(train.phe$PC1, pred.out$PC1) - x.adjust
x.high <- max(train.phe$PC1, pred.out$PC1) + x.adjust
y.adjust <- (max(train.phe$PC2, pred.out$PC2) - min(train.phe$PC2, pred.out$PC2))/10
y.low <- min(train.phe$PC2, pred.out$PC2) - y.adjust
y.high <- max(train.phe$PC2, pred.out$PC2) + y.adjust

# Generate plots
postscript(paste0(prefix, "_ancestryplot.ps"), paper = "letter", horizontal = T)
ncols <- min(3, ceiling(length(unique(pred.out$Ancestry))/2))
if (!require(ggplot2, quietly = TRUE)) {
  plot(pred.out$PC1, pred.out$PC2, col = pred.colors, xlab = "PC1", ylab = "PC2", xlim = c(x.low, x.high), 
       ylim = c(y.low, y.high), main = paste("Inferred Populations as Ancestry in", prefix), pch = 16)
  legend("topright", legend = sort(unique(pred.out$Ancestry)), col = unique(pred.colors)[order(unique(pred.out$Ancestry))], pch = 16, cex = 1)
  par(mfrow = c(2, ncols))
  for (i in sort(unique(pred.out$Ancestry))) {
    subdata <- subset(pred.out, Ancestry == i)
    plot(subdata$PC1, subdata$PC2, col = unique(pred.colors)[unique(pred.out$Ancestry) == i], 
         xlim = c(x.low, x.high), ylim = c(y.low, y.high), xlab = "PC1", ylab = "PC2", 
         main = paste0(i, " (N=", nrow(subdata), ")"))
  }
  par(mfrow = c(1, 1))
  plot(train.phe$PC1, train.phe$PC2, col = train.colors, xlim = c(x.low, x.high), 
       ylim = c(y.low, y.high), xlab = "PC1", ylab = "PC2", main = "Populations in Reference", pch = 16)
  legend("topright", legend = sort(unique(train.phe$Population)), 
         col = unique(train.colors)[order(unique(train.phe$Population))], pch = 16, cex = 1)
} else {
  p <- ggplot(pred.out, aes(x = PC1, y = PC2))
  p <- p + geom_point(aes(colour = factor(Ancestry, levels = sort(unique(Ancestry))))) + 
    xlim(x.low, x.high) + ylim(y.low, y.high) + labs(color = "") +
    scale_colour_manual(values = unique(pred.colors)[order(unique(pred.out$Ancestry))]) + 
    ggtitle(paste("Inferred Populations as Ancestry in", prefix))
  print(p)
  labels <- sapply(sort(unique(pred.out$Ancestry)), function(x) paste0(x, " (N=", sum(pred.out$Ancestry == x), ")"))
  p <- ggplot(pred.out, aes(x = PC1, y = PC2, colour = factor(Ancestry, levels = unique(Ancestry)))) + 
    scale_color_manual(values = unique(pred.colors)) + theme(legend.position = "none")
  p <- p + geom_point() + xlim(x.low, x.high) + ylim(y.low, y.high) + 
    facet_wrap(~factor(Ancestry, levels = sort(unique(Ancestry)), labels = labels), ncol = min(3, ncols))
  print(p)
  p <- ggplot(train.phe, aes(x = PC1, y = PC2))
  p <- p + geom_point(aes(colour = factor(Population, levels = sort(unique(Population))))) + 
    xlim(x.low, x.high) + ylim(y.low, y.high)
  p <- p + labs(color = "") + scale_colour_manual(values = unique(train.colors)[order(unique(train.phe$Population))]) + 
    ggtitle("Populations in Reference")
  print(p)
}
dev.off()
print(paste0(prefix, "_ancestryplot.ps is generated at ", date()))
