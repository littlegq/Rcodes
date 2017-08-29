# setwd("C:/Users/qgong/OneDrive - City of Hope National Medical Center/B.TME")
load("B.TME.R.RData")
if(FALSE){
library(reshape2)
dc <- read.csv("Steidl.2010.CHL.CIBERSORT.Output.csv")
cells <- colnames(dc)[2:23]
dc <- melt(dc[,1:23],id.vars = c('Input_Sample'))
colnames(dc) <- c('Samples', 'Celltypes', 'Content')

dg <- read.table("GDS4222.mt.txt", header=T)

# filter out genes lower 30% overall expression
rowSums <- apply(dg[, -1], 1, sum)
q30 <- quantile(rowSums,probs = 0.3)
filter <- apply(dg[,-1], 1, function(x) sum(x) > q30)
dg <- dg[filter, ]

# filter out genes with lower 70% variance
rowVar <- apply(dg[, -1], 1, var)
q70 <- quantile(rowVar, probs = 0.7)
filter <- apply(dg[, -1], 1, function(x) var(x) > q70)
dg <- dg[filter, ]

gcor <- c()
gcorp <- c()

# identify correlation
for (i in 1:length(cells)){
  for (j in 1:nrow(dg)){
    cell_ct <- dc[which(dc$Celltypes == cells[i]),]
    gev <- data.frame(Samples = colnames(dg[-1]), log2var = as.numeric(dg[j, -1]))
    #    colnames(gev)[2] <- as.character(dg[j,1])
    combine <- merge(cell_ct, gev, by = "Samples")
    out_name <- paste(dg[j, 1], cells[i], sep = ".")
    gcor[out_name] <- cor(combine$Content, combine$log2var)
    gcorp[out_name] <- cor.test(combine$Content, combine$log2var)$p.value
  }
}
dcor <- data.frame(names(gcor), gcor, gcorp)
for(i in 1:nrow(dcor)){
  gene_cell <- strsplit(as.character(dcor$names.gcor.[i]), ".", fixed = T)[[1]]
  dcor$gene[i] <- gene_cell[1]
  dcor$cell[i] <- gene_cell[2]
}
dcor <- dcor[, 2:5]
dcor <- dcor[order(dcor$gcorp),]
write.csv(dcor, "Gene_TME_correlation.csv", row.names = F)
save.image("B.TME.R.RData")
}

library(ggplot2)
library(Rmisc)
colnames(dc) <- c('Samples', 'Celltypes', 'Content')
pdf("CellularComposition.pdf", height = 8.5, width = 11)
## Define the color scale
colors <- c('brown', 'darkred', 'orange',
            'yellow', 'yellowgreen', 'greenyellow', 'lawngreen',
            'darkgreen', 'green', 'lightgreen',
            'cyan', 'darkcyan',
            'skyblue',
            'royalblue', 'blue', 'darkblue',
            'blueviolet', 'purple',
            'darkmagenta', 'magenta',
            'violet', 'pink')
names(colors) <- cells
p <- ggplot(data = dc, aes(x = Samples, y = Content, fill=Celltypes))
p1 <- p + geom_bar(stat="identity") + scale_fill_manual(values = colors) + theme_bw()
p1
dev.off()

pdf("CD68_MacrophagesM1.pdf", height = 8.5, width = 11)
cell <- "Macrophages_M1"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell),]
ev <- apply(dg[dg$GeneSymbol == gene, -1], 2, max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct, gev, by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = Macrophages_M1, y = CD68))
cor(combine$Macrophages_M1, combine$CD68)
p2 <- p + geom_point() + annotate("text", x = 0.2, y = 3.5, label = "r = 0.68") + theme_bw()
p2
dev.off()

pdf("CD68_MacrophagesM0.pdf", height = 8.5, width = 11)
cell <- "Macrophages_M0"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell),]
ev <- apply(dg[dg$GeneSymbol == gene, -1], 2, max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct, gev, by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = Macrophages_M0, y = CD68))
cor(combine$Macrophages_M0, combine$CD68)
p + geom_point() + annotate("text", x = 0.2, y = 3.5, label = "r = 0.45") + theme_bw()
dev.off()

pdf("CD68_CD4TMemoryResting.pdf", height = 8.5, width = 11)
cell <- "T_cells_CD4_memory_resting"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell),]
ev <- apply(dg[dg$GeneSymbol == gene, -1], 2, max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct,gev,by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = T_cells_CD4_memory_resting, y = CD68))
cor(combine$T_cells_CD4_memory_resting, combine$CD68)
p3 <- p + geom_point() + annotate("text", x = 0.2, y = 6.5, label = "r = -0.58") + theme_bw()
p3
dev.off()

pdf("CD68_CD8T.pdf", height = 8.5, width = 11)
cell <- "T_cells_CD8"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell), ]
ev <- apply(dg[dg$GeneSymbol == gene, -1], 2, max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct, gev, by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = T_cells_CD8, y = CD68))
cor(combine$T_cells_CD8, combine$CD68)
p4 <- p + geom_point() + annotate("text", x = 0.2, y = 3.5, label = "r = 0.56") + theme_bw()
p4
dev.off()


pdf("CD68_naiveB.pdf", height = 8.5, width = 11)
cell <- "B_cells_naive"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell),]
ev <- apply(dg[dg$GeneSymbol == gene, -1], 2, max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct, gev, by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = B_cells_naive, y = CD68))
cor(combine$B_cells_naive, combine$CD68)
p + geom_point() + annotate("text", x = 0.2, y = 6.5, label = "r = -0.41") + theme_bw()
dev.off()

pdf("CD68_Tfh.pdf", height = 8.5, width = 11)
cell <- "T_cells_follicular_helper"
gene <- "CD68,SNORA67"
cell_ct <- dc[which(dc$Celltypes == cell),]
ev <- apply(dg[dg$GeneSymbol == gene,-1],2,max)
gev <- data.frame(Samples = names(ev), log2var = as.numeric(ev))
combine <- merge(cell_ct, gev, by = "Samples")
colnames(combine)[3:4] <- c(cell, "CD68")
p <- ggplot(combine, aes(x = T_cells_follicular_helper, y = CD68))
cor(combine$T_cells_follicular_helper, combine$CD68)
p + geom_point() + annotate("text", x = 0.2, y = 6.5, label = "r = -0.41") + theme_bw()
dev.off()

library(Rmisc)
p <- list(p1, p2, p3, p4, p4)
pdf("FigureCombined.pdf", height = 8, width = 20)
lo <- matrix(
  c(1, 1, 1, 2, 3,
    1, 1, 1, 4, 5), 
  nrow = 2, byrow = T)
multiplot(plotlist = p, layout = lo)
dev.off()


# survivial analysis
library(survival)
library(rms)
library(ggplot2)
source("ggsurv.R")
tdc <- read.csv("Steidl.2010.CHL.CIBERSORT.Output.csv")
CDat <- read.table("Steidl.2010.CHL.CaseInfo.txt", head = T)
colnames(CDat)[19:21] <- c('OS', 'DSS', 'PFS')
tdc <- merge(tdc,CDat,by.x = "Input.Sample", by.y = "NcbiID")
pdf("PFS.plot.pdf")
CellTypes <- colnames(tdc)[2:23]
for (cell in CellTypes){
  MedCont <- median(tdc[, cell])
  pfs <- tdc[,c('Input.Sample', cell, 'CODENAPFS','PFS')]
  pfs$High <- pfs[, 2] > MedCont
  survp <- survdiff(Surv(PFS, CODENAPFS) ~ High, data = pfs)
  survp <- 1 - pchisq(survp$chisq, length(survp$n) - 1)
  if(survp < 0.05){
    High.surv <- survfit(Surv(PFS, CODENAPFS) ~ High, data = pfs)
    FigTitle <- paste0(cell,">",MedCont,"; p=",survp)
    print(ggsurv(High.surv) + ggtitle(FigTitle))
  }
}
dev.off()

pdf("OS.plot.pdf")
CellTypes <- colnames(tdc)[2:23]
for (cell in CellTypes){
  MedCont <- median(tdc[, cell])
  os <- tdc[, c('Input.Sample', cell, 'CODENAOS', 'OS')]
  os$High <- os[, 2] > MedCont
  survp <- survdiff(Surv(OS, CODENAOS) ~ High, data = os)
  survp <- 1 - pchisq(survp$chisq, length(survp$n) - 1)
  if(survp < 0.05){
    High.surv <- survfit(Surv(OS, CODENAOS) ~ High, data = os)
    FigTitle <- paste0(cell, ">", MedCont, "; p=", survp)
    print(ggsurv(High.surv) + ggtitle(FigTitle))
  }
}
dev.off()

