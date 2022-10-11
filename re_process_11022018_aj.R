# 11/02/2018
# asked to reprocess the data regressing out variables?
#
library(Seurat)
library(dplyr)
library(reshape2)
all<-readRDS(file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/source_seur_obj.RDS")
head(all@meta.data)

# Update the Seurat object
## check that Seurat package is at least v3.0
utils::packageVersion('Seurat') < 3 
## check version of Seurat object 
all@version < 3

## UpdateSeuratObject: Update old Seurat object to accomodate new features
all_updated <- UpdateSeuratObject(all)
class(all_updated)
head(all_updated@meta.data)
all <- all_updated

## What variable is the identity class of this seurat object
unique(Idents(object = all)) # orig.ident is identifier for each cell/library
## Assign another variable as identifier 
#all<-SetAllIdent(all, id = "sample_ID") # Deprecated
Idents(object = all) <- "sample_ID"
unique(Idents(object = all)) # Total 32 Samples 

#unique(all@ident) # 32 samples, remove 5 samples from batch 1, 2 samples from batch 3 # Deprecated

## Based on QC Remove 5 samples from batch 1, 2 samples from batch 3
all<-subset(all, idents=c("1443_sgACC_1","1443_DLPFC_1","2543_sgACC_1", "2543_DLPFC_1", "1443_sgACC_2", "2587_sgACC_1", "2587_sgACC_2"), invert=TRUE) # Depecated "SubsetData"
unique(Idents(object = all)) # Total 25 Samples 

#unique(all@ident) # 25 samples
#all<-SetAllIdent(all, id = "orig.ident") # Deprecated

# Visualize technical variables
pdata<-all@meta.data
pdf("pdf1.pdf")
# frac.mito
ggplot(pdata, aes(x = sample_ID, y = frac.mito )) +
  geom_violin() + geom_jitter(size=0.01) +
  geom_hline(yintercept=c(0.2), linetype="dashed", color = "red") + 
  theme(axis.text.x = element_text(angle=90))
# frac.ribo
ggplot(pdata, aes(x = sample_ID, y = frac.ribo )) +
  geom_violin() + geom_jitter(size=0.01) +
  theme(axis.text.x = element_text(angle=90))
dev.off()
# remove cells with frac.mito > 0.2; 
all_hiMitoFrac<-subset(all, subset= frac.mito > 0.2) # this removes about 800-900 cells
all<-subset(all, subset= frac.mito < 0.2) # Deprecated "FilterCells"




# remove mitochondrial genes
#expr.mat<-as(all@meta.data, "dgTMatrix") # Deprecated
expr.mat <- GetAssayData(all, slot = "counts")
#expr.mat <- as.matrix(expr.mat)
# remove 24 mitochondrially encoded genes
mgenes<-grep(pattern="^MT-|^MTRNR",x=rownames(expr.mat),value=T)
remove_rows<-which(rownames(expr.mat) %in% mgenes)
expr.mat.no.mito <-expr.mat[-remove_rows,]


# Decision: remove genes not expressed in at least 3 cells (5475 genes)
# decided to go with at least 3 cells
num.cells.expressing <-apply(expr.mat.no.mito, 1, function(x) sum(x>0))
pdf("pdf2.pdf")
# how many cells express each gene?
hist(num.cells.expressing, main="How many cells express each gene?", xlab="# of cells", 
     ylab="# of genes", col="darkmagenta", xlim=c(0,20), breaks=100000, ylim=c(0,5000))
hist(num.cells.expressing, main="How many cells express each gene?", xlab="# of cells", 
     ylab="# of genes", col="darkmagenta", breaks=100, ylim=c(0,5000))
plot(density(num.cells.expressing), main="How many cells express each gene?")
noise_genes<-which(num.cells.expressing<=3)
genes.to.be.removed<-rownames(expr.mat.no.mito)[noise_genes]
write.table(genes.to.be.removed, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/genes.removed",sep="\t")
expr.mat.no.mito.min<-expr.mat.no.mito[-noise_genes,]

# return matrix in dgCMatrix
expr.mat3<-as(expr.mat.no.mito.min, "dgCMatrix")
all@raw.data<-expr.mat3
# recalculate nGene, nUMI
cell.names<-rownames(all@meta.data)
dim(expr.mat)
expr.mat2<-expr.mat.no.mito.min[,cell.names]
nUMI<-colSums(expr.mat2)
nUMI<-as.data.frame(nUMI)
nGene<-apply(expr.mat2, 2, function(x) sum(x!=0))
nGene<-as.data.frame(nGene)
identical(rownames(nUMI),rownames(nGene))
nUMI$nGene<-nGene$nGene

all@meta.data$nUMI <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "nUMI"]
all@meta.data$nGene <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "nGene"]

# recalculate frac.ribo (177 ribo protein genes)
rpgenes<-grep(pattern="^MRPL|^MRPS|^RPS|^RPL", x=rownames(expr.mat.no.mito.min), value=T)
rp_rows<-which(rownames(expr.mat.no.mito.min) %in% rpgenes)
rp.expr<-expr.mat.no.mito.min[rp_rows, cell.names]
rp.tot.counts<-colSums(rp.expr)
rp.tot.counts<-as.data.frame(rp.tot.counts)
identical(rownames(nUMI), rownames(rp.tot.counts))
nUMI$ribo.tot.counts<-rp.tot.counts$rp.tot.counts
nUMI$frac.ribo.2<-(nUMI$ribo.tot.counts/nUMI$nUMI)
head(nUMI)
class(nUMI$frac.ribo.2)
all@meta.data$ribo.tot.counts <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "ribo.tot.counts"]
all@meta.data$frac.ribo.2 <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "frac.ribo.2"]


plot(density(all@meta.data$frac.ribo))
lines(density(all@meta.data$frac.ribo.2))

# recalculate frac.mito (24 mito encoded genes)
# calculate mito counts for each cell
mito_rows<-remove_rows
mito.expr<-expr.mat[mito_rows,cell.names]
mito.tot.counts<-colSums(mito.expr)
mito.tot.counts<-as.data.frame(mito.tot.counts)

identical(rownames(nUMI),rownames(mito.tot.counts))
nUMI$mito.tot.counts<-mito.tot.counts$mito.tot.counts
nUMI$frac.mito.2<-(nUMI$mito.tot.counts/nUMI$nUMI)
head(nUMI)
class(nUMI$frac.mito.2)
all@meta.data$mito.tot.counts <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "mito.tot.counts"]
all@meta.data$frac.mito.2 <- nUMI[match(rownames(nUMI), rownames(all@meta.data)), "frac.mito.2"]
max(all@meta.data$frac.mito.2)

plot(density(all@meta.data$frac.mito))
plot(lines(all@meta.data$frac.mito.2))

pdata<-all@meta.data
# nGene
ggplot(pdata, aes(x = sample_ID, y = nGene )) +
  geom_violin() + 
  geom_hline(yintercept=c(200,5000), linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle=90))
# nUMI
ggplot(pdata, aes(x = sample_ID, y = nUMI)) +
  geom_violin() +
  geom_hline(yintercept=c(200,15000), linetype="dashed", color = "red") +
  theme(axis.text.x = element_text(angle=90))
######## Filter cells
all<-FilterCells(all, subset.names=c("nGene","nUMI"), 
                  low.thresholds = c(200,200),
                  high.thresholds = c(5000,15000))
saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.RDS")
all<-readRDS(file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.RDS")
################################################# Normalize and scale data based on highly variable genes ##########################################################3
# like CPM, allow you to compare different library sizes
all<-NormalizeData(all,normalization.method="LogNormalize",scale.factor=10000)
all<-FindVariableGenes(all,mean.function=ExpMean,dispersion.function = LogVMR,x.low.cutoff = 0.015,x.high.cutoff = 4,y.cutoff = 0.65)
length(all@var.genes)
# X variable genes
# Scale data, if using negative binomial, regresses on the UMI count data
all<-ScaleData(all, vars.to.regress=c("nUMI","batch","frac.mito.2"),
               model.use = "negbinom", genes.use=all@var.genes, do.par=TRUE)
saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.RDS")
# get error/warning message from ScaleData that some genes could not be regressed to neg binom
# "the following genes failed with glm.nb and fell back to scale(log(y+1))"
# "LINC01090, ZNF385D-AS2, FAM19A1, KCCAT211, KCNIP4, AC011288.2, NPSR1, RP11-30J20.1, RP11-314P15.2"
hvgs<-all@var.genes
failed.genes<-c("LINC01090", "ZNF385D-AS2", "FAM19A1", "KCCAT211", "KCNIP4", "AC011288.2", "NPSR1", "RP11-30J20.1", "RP11-314P15.2")
hvgs<-hvgs[!hvgs %in% failed.genes]
# 4443 - 9 genes failed = 4434 hvgs
length(grep(pattern="^MRPL|^MRPS|^RPS|^RPL", x=hvgs, value=T))
# 102 rp genes

# Perform PCA - linear dimensional reduction
all<-RunPCA(all, pc.genes=hvgs, do.print=TRUE, pcs.print=1:5, genes.print=5, pcs.compute=100)
PCElbowPlot(all,num.pc=100)
# plot heatmap of PCs
PCHeatmap(all, pc.use=c(1:9), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(10:18), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(19:26), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(27:35), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(36:44), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(45:53), cells.use=500, num.genes=30, do.balanced = TRUE)
PCHeatmap(all, pc.use=c(54:62), cells.use=500, num.genes=30, do.balanced = TRUE)

# Plot PCA by variables
PCAPlot(object=all,dim.1=1,dim.2=2, group.by="orig.ident")
PCAPlot(object=all,dim.1=1,dim.2=2, group.by="br")
PCAPlot(object=all,dim.1=3,dim.2=4, group.by="br")
PCAPlot(object=all,dim.1=3,dim.2=4, group.by="batch")
#

#######
saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_11032018.RDS")
all<-readRDS(file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_11032018.RDS")
#######
########################################## Find Clusters and Remove Mixes/Unclassifiable ##########################################
all<-FindClusters(all,reduction.type="pca",dims.use=1:50, resolution=3, save.SNN=TRUE)
saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.RDS")
all<-readRDS(file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.RDS")
plot_vln <- function(t, my.genes3) {
  d <- as.matrix(t@data[intersect(my.genes3, rownames(t@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$tree.ident <- t@ident[dd$cell]
  str(dd$tree.ident)
  dd$gene <- factor(dd$gene, levels = intersect(my.genes3, rownames(t@data)))
  ggplot(dd, aes(tree.ident, value, fill = tree.ident)) + geom_violin(scale = "width", trim = T, alpha = 0.8, adjust = 1) + facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") + theme(strip.background = element_blank(), strip.placement="outside", axis.text.y=element_blank(), axis.title.y=element_blank(), strip.text.y = element_text(colour="red", angle=360, size=10), legend.position="none", panel.grid=element_blank(), panel.border = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.9)), legend.position = "none") + xlab("")
}

# mk gene panels
#
plot_vln(all,c("KCNK1","KIAA1324","LNX1","NELL1","COBL","SLITRK1","DPYSL5","C14orf37",
                 "DLX1","DLX6","GLRA2","DLX2","DLX5","SLC10A4","EGFR","SST","PNOC",
                 "NXPH1","BCL11A","DCN","TMEM130","CNTN4","CDO1","NFASC","LRRTM3",
                 "GRIA3","RELN","SNAP25","SYP","RBFOX3","SNHG11",
                 "SLC17A7","SATB2","GAD1","GAD2","SLC32A1","SLC6A1"))

plot_vln(all, c("DAAM2","ASPA","MAL","SEC14L5","MAP6D1","DPYD","PPP1R14A",
                  "GJB1", "FA2H", "MAG", "CDK18", "LGI3", "SHC4", "UGT8", "KLK6",
                  "KCNH8", "GPR37", "MOBP", "LPAR1", "ADAMTS4", "ERMN", "OPALIN",
                  "CLDN11", "PLEKHB1", "GSN", "GRM3", "CNP", "MBP", "PLP1",
                  "MBP", "MOG", "ENPP6", "MOG1"))
plot_vln(sgACC, c("SLC17A7","SATB2", "GAD1","GAD2","SLC1A2","GJA1","AQP4","ATP1A2","MOG","PLP1", "PDGFRA","GPR17","CD74","CX3CR1","FLT1","PECAM1"))

plot_vln(all, c("SLC14A1", "GLIS3", "GLI3", "PPP1R3C", "CHRDL1", "CYBRD1", "CTH", "SORCS2",
                  "ITGB4", "RNF43", "NWD1", "PAQR6", "C16orf89", "ALDH1L1", "TRIM66", "HGF",
                  "CBS", "ITGA7", "SLC30A10", "FGFR", "BMPR1B", "ATP13A4", "AQP4",
                  "ATP1A2", "GJA1", "SLC1A2"))
#microglia
plot_vln(all, c("GPR183","CCL4","CD83","LAPTM5","CSF1R","HLA-DRA","BCL2A1","CD14","CCL2","CTSS","ITGAM","PTPRC","CD74","CX3CR1","SPI1"))
#endo
plot_vln(all, c("APOLD1","TM4SF1","FLT1","A2M","PECAM1","TEK","MYL9","PDGFRB"))
#opc
plot_vln(all, c("PDGFRA", "LHFPL3", "MEGF11", "PCDH15", "CSPG4", "GPR17"))

ggplot(sgACC@meta.data, aes(x=interim_id, y=nGene)) + geom_violin() + geom_jitter(size=0.01) + 
  theme(axis.text.x = element_text(angle=90))


# NSC
plot_vln(all, c("SPAG17","WDR96","C9orf117", "CAPS", "SLC47A2", "DTHD1", "EGR1",
                  "CCDC42B", "LRRIQ1", "CYR61", "WDR63", "ANXA1", "TTC18", "TOB11",
                  "MAP3K19", "WDR52", "RSPH1", "FHAD1", "CCDC146", "ARMC3", "FOS1",
                  "FOXJ1", "C12orf55", "TCTEX1D1", "CP", "SPAG8",
                  "ID4","VIM","SOX2","AQP4"))



plot_vln(all, c("SYP","RBFOX3","CNTN4",
                 "KCNK1","SLC17A7","SATB2",
                 "GAD1","GAD2","SLC32A1","SLC6A1",
                 "SLC14A1","GLI3","ALDH1L1","BMPR1B","ATP13A4","AQP4","SLC1A2",
                 "MOG","PLP1","ASPA","OPALIN","CLDN11","CNP","PLP1",
                 "PDGFRA","MEGF11","GPR17",
                 "GPR183","CCL4","LAPTM5","CSF1R","CD14","CD74","PTPRC",
                 "APOLD1","TM4SF1","FLT1","PECAM1","MYL9"))
plot_vln(all, c("SLC17A7","SATB2", "GAD1","GAD2","SLC1A2","GJA1","AQP4","ATP1A2","MOG","PLP1","PDGFRA","GPR17","CD74","CX3CR1","FLT1","PECAM1"))




############################################## ID Clusters #######################################
IDs<-data.frame(cluster.mem=c(0,2,3,5,8,9,10,11,15,18,23,28,30,32,34,35,42,43,46,47,51,12,16,20,26,
                              1,7,13,14,17,19,21,36,37,38,45,49,50,
                              4,24,25,33,41,
                              22,31,39,
                              6,44,
                              29,
                              27,40,56,
                              48,52,53,54,55,57),
                cluster.name=c(paste("Ex",1:25,sep=""),paste("In",1:13,sep=""),
                               paste("Asc",1:5,sep=""),paste("Odc",1:3,sep=""),
                               paste("Opc",1:2,sep=""),paste("Mg",1,sep=""),
                               paste("End",1:3,sep=""),paste("Mix",1:6,sep="")),
                neuron.name=c(rep("Neuron",38,sep=""),
                              rep("Non-neuron",14,sep=""),
                              rep("Mix/Unknown",6,sep="")),
                broad.class.name=c(rep("Excitatory Neurons",25,sep=""),
                                   rep("Inhibitory Neurons",13,sep=""),
                                   rep("Astrocytes",5,sep=""),
                                   rep("Oligodendrocytes",3,sep=""),
                                   rep("Oligodendrocyte Precursor Cells",2,sep=""),
                                   rep("Microglia",1,sep=""),
                                   rep("Endothelial Cells",3,sep=""),
                                   rep("Mixtures",6,sep="")),
                broad.class.color=c(rep("orange",25,sep=""),
                                    rep("yellow",13,sep=""),
                                    rep("dodgerblue1",5,sep=""),
                                    rep("green2",3,sep=""),
                                    rep("green4",2,sep=""),
                                    rep("gray",1,sep=""),
                                    rep("dodgerblue4",3,sep=""),
                                    rep("darkviolet",6,sep="")))
neuron.name<-c("Neuron", "Non-neuron", "Mix/Unknown")
neuron.name.colors<-c("red","blue","gray")
broad.class<-c("Excitatory Neurons", "Inhibitory Neurons", "Astrocytes", "Oligodendrocytes", 
               "Oligodendrocyte Precursor Cells","Microglia", "Endothelial Cells", "Mixtures")
broad.class.colors<-c("orange","yellow","dodgerblue1","green2","green4","gray","dodgerblue4","darkviolet")

id_clusters<-function(seur_obj, ID_df){
  temp<-seur_obj@meta.data
  temp$interim_id<-plyr::mapvalues(temp$res.3, ID_df$cluster.mem, as.character(ID_df$cluster.name))
  temp$neuron<-plyr::mapvalues(temp$res.3, ID_df$cluster.mem, as.character(ID_df$neuron.name))
  temp$broad.class<-plyr::mapvalues(temp$res.3, ID_df$cluster.mem, as.character(ID_df$broad.class.name))
  head(temp)
  seur_obj@meta.data<-temp
  seur_obj<-SetAllIdent(seur_obj, id = "interim_id")
  seur_obj@ident<-factor(seur_obj@ident, levels=as.character(ID_df$cluster.name))
  return(seur_obj)
}

all<-id_clusters(all, IDs)


pdata<-all@meta.data
pdata$interim_id<-factor(pdata$interim_id, levels=IDs$cluster.name)
pdata$neuron<-factor(pdata$neuron, levels=neuron.name)
pdata$broad.class<-factor(pdata$broad.class, levels=broad.class)

ggplot(pdata, aes(x=interim_id,y=nUMI, fill=broad.class.colors)) + 
  geom_violin() + geom_jitter(size=0.05) + theme(axis.text.x=element_text(angle=90))
ggplot(pdata, aes(x=interim_id,y=nUMI)) + geom_violin() + theme(axis.text.x=element_text(angle=90))

# interim id
ggplot(pdata, aes(x=interim_id,y=frac.ribo, fill=broad.class)) + geom_violin() + 
  theme(axis.text.x=element_text(angle=90)) + 
  scale_fill_manual("Broad Cell Class", 
                    values=broad.class.colors)
# broad cell class
ggplot(pdata, aes(x=broad.class,y=frac.ribo, fill=broad.class)) + geom_violin() + 
  theme(axis.text.x=element_text(angle=90)) + 
  scale_fill_manual("Broad Cell Class", 
                    values=broad.class.colors)
# neuron or non-neuron
ggplot(pdata, aes(x=interim_id,y=nUMI, fill=neuron)) + geom_violin() + 
  theme(axis.text.x=element_text(angle=90)) + 
  scale_fill_manual("Neuronal or Non-neuronal", 
                    values=neuron.name.colors)


################################################## tSNE ################################################################ 
all<-RunTSNE(all, reduction.use="pca", dims.use=1:50, do.fast=TRUE)

# tSNE visualization
TSNEPlot(all, do.label=TRUE, pt.size=0.1)
all@meta.data$broad.class<-factor(all@meta.data$broad.class, levels=broad.class)
TSNEPlot(all, do.label=TRUE, group.by="broad.class", pt.size=0.01, colors.use=broad.class.colors)

TSNEPlot(object = all, group.by="batch", pt.size=0.01)
TSNEPlot(object = all, group.by="orig.ident", pt.size = 0.01)
TSNEPlot(object = all, group.by="br", pt.size=0.01)

FeaturePlot(all, features.plot="nUMI")
FeaturePlot(all, features.plot="nGene")
FeaturePlot(all, features.plot="frac.mito.2")
FeaturePlot(all, features.plot="frac.ribo.2" )

tsne.cell.embeddings<-all@dr$tsne@cell.embeddings
pdata$tSNE_1<-tsne.cell.embeddings[match(rownames(tsne.cell.embeddings), row.names(pdata)), "tSNE_1"]
pdata$tSNE_2<-tsne.cell.embeddings[match(rownames(tsne.cell.embeddings), row.names(pdata)), "tSNE_2"]
ggplot(pdata, aes(x=tSNE_1, y=tSNE_2, color = broad.class)) + geom_point(size=0.01) + 
  scale_color_manual("Broad Cell Class", 
                     values=broad.class.colors)


####################################################### Compute AVG EXPRESSION for cor & Hclust ###############################################

measurements.Exp1<-as(all@data, "dgTMatrix")
measurements.Exp1<-as.matrix(measurements.Exp1)
measurements.Exp1<-exp(measurements.Exp1)-1

membership.Exp1<-all@meta.data[,"interim_id"]
manual.Experiment1.order<-unique(membership.Exp1)

means.Exp1 = sapply(manual.Experiment1.order,function(group) rowMeans(measurements.Exp1[,membership.Exp1 == group]))
avg.all<-means.Exp1
rm(measurements.Exp1, membership.Exp1, manual.Experiment1.order, means.Exp1)


length(all@var.genes) #3907
avg.all2<-avg.all[all@var.genes,]

######################################################### CORRELATION HEATMAP #############################################
get_cor_df<-function(avg.data.1, avg.data.2)
{
  genelist<-intersect(unique(rownames(avg.data.1)), rownames(avg.data.2)) #only consider co-expressed markers
  print(paste("Number of highly variable genes compared: ",length(genelist),sep=""))
  df1<-avg.data.1[genelist,]
  df2<-avg.data.2[genelist,]
  
  df1<-log1p(df1)
  df2<-log1p(df2)
  num<-as.numeric(dim(df1)[2]*dim(df2)[2])
  df.cor<-data.frame(matrix(ncol=3,nrow=num))
  
  k=0
  for(i in colnames(df1)){
    for(j in colnames(df2)){
      k=k+1
      df.cor[k,]<-c(i,j,cor(df1[,i],df2[,j],use="complete.obs"))
      #print(c(i,j,cor(df2[[i]],df1[[j]])))
      #print(k)
    }
  }
  df.cor$X3<-as.numeric(df.cor$X3)
  return(df.cor)
}
comp.n<-get_cor_df(avg.all2,avg.all2)
head(comp.n)
# for res 2
comp.n$X1<-factor(comp.n$X1, levels=IDs$cluster.name)
comp.n$X2<-factor(comp.n$X2, levels=IDs$cluster.name)
ggplot(data=comp.n,aes(x=X2,y=X1,fill=X3))+geom_tile() + 
  scale_fill_gradientn(colors=rev(c(colorRampPalette(c("red", "white"))(5)[-5],colorRampPalette(c("white", "blue"))(6)))) +
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), axis.ticks.x=element_blank(),axis.title.x=element_blank())+ xlab("")+
  ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + ggtitle("Compare to Self; Res 2")
################################### hierarchical clustering analysis - dendrogram ############################
make_hc<-function(avg.exp.mat){
  sc.avg<-t(scale(t(avg.exp.mat),scale=T))
  dim(sc.avg)
  d<-dist(as.matrix(t(sc.avg)))
  hcc<-hclust(d,method="ward.D2")
  return(hcc)
}
hc_all<-make_hc(avg.all2)
plot(hc_all,hang=-1,main="Cluster Dendrogram")

library(grid)
library(ggdendro)
hc_all<-as.dendrogram(hc_all)
data <- dendro_data(hc_all)
dendro.plot <- ggplot() + geom_segment(data = segment(data), 
                                       aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
  geom_text(data = label(data), 
            aes_string(x = "x", y = "y", label = "label"), hjust = 1, angle = 0) + 
  scale_y_continuous(expand=c(0.2, 0)) + coord_flip() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.line.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
#dendro.plot<-ggdendrogram(hc_all, rotate=TRUE)+ theme(axis.text.x=element_text(size=8))
# re order dendrogram
h.order<-order.dendrogram(hc_all)
comp.n<-get_cor_df(avg.all2,avg.all2)# for res 2
comp.n$X1<-factor(comp.n$X1, levels=colnames(avg.all2)[h.order], ordered=TRUE)
comp.n$X2<-factor(comp.n$X2, levels=colnames(avg.all2)[h.order], ordered=TRUE)
mid<-mean(comp.n$X3)
heatmap.plot<-ggplot(data=comp.n, aes(x=X1, y=X2, fill=X3)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white",high="red", midpoint=mid, space="Lab") + 
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), 
        axis.title.x=element_blank()) + 
  xlab("") + ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + 
  ggtitle("Compare to self; Res 2")
# width=1000 save default w:h ratio
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.31, y = 0.5, width = 0.65, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.80, y = 0.55, width = 0.4, height = 0.91))


#### Remove suspected technical clusters to see effects on hierarchical clustering
##
################################ remove mixes and unknowns   #################################
# remove Mix1-6
remove.mix.unk<-avg.all2[, -which(colnames(avg.all2) %in% c(paste("Mix",1:6,sep="")))]
hc_remove_mixes<-make_hc(remove.mix.unk)
plot(hc_remove_mixes,hang=-1,main="Cluster Dendrogram, Remove Mixes and Unknowns")

comp.n<-get_cor_df(remove.mix.unk, remove.mix.unk)
head(comp.n)
temp<-c(paste("Ex",1:25,sep=""),paste("Inh",1:13,sep=""),
        paste("Asc",1:5,sep=""),paste("Odc",1:3,sep=""),
        paste("Opc",1:2,sep=""),paste("Mg",1,sep=""),
        paste("End",1:3,sep=""))

comp.n$X1<-factor(comp.n$X1, levels=temp)
comp.n$X2<-factor(comp.n$X2, levels=temp)

library(grid)
library(ggdendro)
hc_remove_mixes<-as.dendrogram(hc_remove_mixes)
data <- dendro_data(hc_remove_mixes)
dendro.plot <- ggplot() + geom_segment(data = segment(data), 
                                       aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
  geom_text(data = label(data), 
            aes_string(x = "x", y = "y", label = "label"), hjust = 1, angle = 0) + 
  scale_y_continuous(expand=c(0.2, 0)) + coord_flip() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.line.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
#dendro.plot<-ggdendrogram(hc_remove_mixes, rotate=TRUE)+ theme(axis.text.x=element_text(size=8))
# re order dendrogram
h.order<-order.dendrogram(hc_remove_mixes)
comp.n<-get_cor_df(remove.mix.unk,remove.mix.unk)# for res 2
comp.n$X1<-factor(comp.n$X1, levels=colnames(remove.mix.unk)[h.order], ordered=TRUE)
comp.n$X2<-factor(comp.n$X2, levels=colnames(remove.mix.unk)[h.order], ordered=TRUE)
mid<-mean(comp.n$X3)
heatmap.plot<-ggplot(data=comp.n, aes(x=X1, y=X2, fill=X3)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white",high="red", midpoint=mid, space="Lab") + 
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), 
        axis.title.x=element_blank()) + 
  xlab("") + ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + 
  ggtitle("Compare to self; Res 2")
# width=1000 save default w:h ratio
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.31, y = 0.5, width = 0.65, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.80, y = 0.55, width = 0.4, height = 0.91))




############################### remove Ex4, Ex23:25 #################################
remove.ex<-remove.mix.unk[, -which(colnames(remove.mix.unk) %in% c("Ex4","Ex23","Ex24","Ex25"))]

hc_remove_ex<-make_hc(remove.ex)
plot(hc_remove_ex,hang=-1,main="Cluster Dendrogram, Remove Mixes and Unknowns")

comp.n<-get_cor_df(remove.ex, remove.ex)
head(comp.n)
temp<-c(paste("Ex",1:3,sep=""),
        paste("Ex",5:22,sep=""),
        paste("Inh",1:13,sep=""),
        paste("Asc",1:5,sep=""),paste("Odc",1:3,sep=""),
        paste("Opc",1:2,sep=""),paste("Mg",1,sep=""),
        paste("End",1:3,sep=""))

comp.n$X1<-factor(comp.n$X1, levels=temp)
comp.n$X2<-factor(comp.n$X2, levels=temp)

library(grid)
library(ggdendro)
hc_remove_ex<-as.dendrogram(hc_remove_ex)
data <- dendro_data(hc_remove_ex)
dendro.plot <- ggplot() + geom_segment(data = segment(data), 
                                       aes_string(x = "x", y = "y", xend = "xend", yend = "yend")) + 
  geom_text(data = label(data), 
            aes_string(x = "x", y = "y", label = "label"), hjust = 1, angle = 0) + 
  scale_y_continuous(expand=c(0.2, 0)) + coord_flip() + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.line.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
#dendro.plot<-ggdendrogram(hc_remove_ex, rotate=TRUE)+ theme(axis.text.x=element_text(size=8))
# re order dendrogram
h.order<-order.dendrogram(hc_remove_ex)
comp.n<-get_cor_df(remove.ex,remove.ex)# for res 2
comp.n$X1<-factor(comp.n$X1, levels=colnames(remove.ex)[h.order], ordered=TRUE)
comp.n$X2<-factor(comp.n$X2, levels=colnames(remove.ex)[h.order], ordered=TRUE)
mid<-mean(comp.n$X3)
heatmap.plot<-ggplot(data=comp.n, aes(x=X1, y=X2, fill=X3)) + 
  geom_tile() + 
  scale_fill_gradient2(low="blue",mid="white",high="red", midpoint=mid, space="Lab") + 
  theme(axis.text.x=element_text(angle = 45,hjust=1,vjust=1), 
        axis.title.x=element_blank()) + 
  xlab("") + ylab("") + theme(legend.position="bottom",legend.title=element_blank()) + 
  ggtitle("Compare to self; Res 2")
# width=1000 save default w:h ratio
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.31, y = 0.5, width = 0.65, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.80, y = 0.55, width = 0.4, height = 0.91))

####################################  Individual contribution to each cluster #########################
# Contribution of each subject to a cluster
#for id in s_IDs$cluster.name
make_fraction<-function(seur_obj, id){
  temp<-seur_obj@meta.data
  temp2<-temp %>% 
    group_by(interim_id, orig.ident) %>% 
    dplyr::summarise(nuclei.total=n())
  temp2<-as.data.frame(temp2)
  colnames(temp2)<-c("cluster", "br_number", "nuclei")
  temp3<-temp2 %>% group_by(br_number) %>% mutate(fraction=nuclei/sum(nuclei))
  temp3<-as.data.frame(temp3)
  temp3$cluster<-factor(temp3$cluster, levels=id$cluster.name)
  return(temp3)
}
fr.all<-make_fraction(all, IDs)
head(fr.all)

ggplot(fr.all, aes(x=cluster, y=fraction, fill=br_number)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# br contribution to each cluster
make_fraction<-function(seur_obj, id){
  temp<-seur_obj@meta.data
  temp2<-temp %>% 
    group_by(interim_id, br) %>% 
    dplyr::summarise(nuclei.total=n())
  temp2<-as.data.frame(temp2)
  colnames(temp2)<-c("cluster", "br_region", "nuclei")
  temp3<-temp2 %>% group_by(br_region) %>% mutate(fraction=nuclei/sum(nuclei))
  temp3<-as.data.frame(temp3)
  #temp3$cluster<-factor(temp3$cluster, levels=id$cluster.name)
  return(temp3)
}
fr.all<-make_fraction(all, IDs)
head(fr.all)
fr.all$cluster<-factor(fr.all$cluster, levels=IDs$cluster.name)
ggplot(fr.all, aes(x=cluster, y=fraction, fill=br_region)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))



####################################  Individual contribution to each cluster #########################
# Contribution of each subject to a cluster
#for id in s_IDs$cluster.name
make_fraction<-function(seur_obj, id){
  temp<-seur_obj@meta.data
  temp2<-temp %>% 
    group_by(interim_id, orig.ident) %>% 
    dplyr::summarise(nuclei.total=n())
  temp2<-as.data.frame(temp2)
  colnames(temp2)<-c("cluster", "br_number", "nuclei")
  temp3<-temp2 %>% group_by(br_number) %>% mutate(fraction=nuclei/sum(nuclei))
  temp3<-as.data.frame(temp3)
  temp3$cluster<-factor(temp3$cluster, levels=id$cluster.name)
  return(temp3)
}
fr.all<-make_fraction(all, IDs)
head(fr.all)

ggplot(fr.all, aes(x=cluster, y=fraction, fill=br_number)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# batch contribution to each cluster
make_fraction<-function(seur_obj, id){
  temp<-seur_obj@meta.data
  temp2<-temp %>% 
    group_by(interim_id, batch) %>% 
    dplyr::summarise(nuclei.total=n())
  temp2<-as.data.frame(temp2)
  colnames(temp2)<-c("cluster", "batch", "nuclei")
  temp3<-temp2 %>% group_by(batch) %>% mutate(fraction=nuclei/sum(nuclei))
  temp3<-as.data.frame(temp3)
  #temp3$cluster<-factor(temp3$cluster, levels=id$cluster.name)
  return(temp3)
}
fr.all<-make_fraction(all, IDs)
head(fr.all)
fr.all$cluster<-factor(fr.all$cluster, levels=IDs$cluster.name)
ggplot(fr.all, aes(x=cluster, y=fraction, fill=batch)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))



saveRDS(all, file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.tSNE.RDS")
all<-readRDS(file="/data/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_postfilt.scale.clust.tSNE.RDS")