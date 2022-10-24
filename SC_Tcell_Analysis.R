# Collaboration with Maya and Donna Farber (190828_DONNA_MEI_6_HUMAN_PLATE)
# Isolated total T cells from tissue of spleen, BM, lung and blood from one donor
# 21 plates
# Stimulated with peptide pools for CMV, CD4, CMV CD8, flu CD8, flu HA, flu internal proteins (others)
# =24hr stimulations (differs from first experiment)
# Sorted by production of IFgamma (representing stimulated cells)
# Goal: to see response to viral antigen flu vs. cmv in tissue sites, across sites in each
library(Seurat)
library(ggfortify)
library(stringr)
library(dplyr)
setwd("~/Dropbox/Maya_Tcell/")

#filter for protein coding
proteins <- read.table("protein_coding.csv", header=T, sep=",")
umi.counts <- read.table("Tcell_matrix.csv",header = T,row.names = 1,sep=",")
umi.counts <- subset(umi.counts, rownames(umi.counts) %in% proteins$Gene)
meta <- read.table("Tcell_meta.csv",header = T,sep=",",row.names=1)
meta <- cbind(meta,data.frame(str_split_fixed(meta$Identity, "_", 3)))
meta$ID <- rownames(meta)
colnames(meta) <- c(colnames(meta)[1:4],"Tissue","Virus","Cell","ID")
meta <- subset(meta, meta$Identity != "EMPTY")
umi.counts <- subset(umi.counts, select=colnames(umi.counts) %in% meta$ID)

# Seurat object creation and general filtering
tcell <- CreateSeuratObject(counts = umi.counts, project = "All_Plates", names.delim = "", meta.data = meta)
tcell@meta.data[["orig.ident"]] <- tcell@meta.data[["Identity"]]
tcell <- subset(tcell, subset = Identity !=  "EMPTY")
tcell[["percent.mt"]] <- PercentageFeatureSet(tcell, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(tcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(tcell, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tcell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Filter
tcell <- subset(tcell, subset = nFeature_RNA > 200 & percent.mt < 10)
tcell <- SetIdent(object = tcell, value = tcell@meta.data[["Identity"]])
VlnPlot(tcell, features = c("nFeature_RNA", "percent.mt"), ncol = 2)

#just CD8
tcell <- subset(tcell, subset = Cell ==  "CD8")
#CD8 outliers from PCA
tcell <- subset(tcell, subset = ID !="DM037_CGTGTCTA" & ID !="DM037_CTAGCTGT"
                & ID != "DM035_CTGTGTAC")
#just CD4
tcell <- subset(tcell, subset = Cell ==  "CD4")

#Normalize by log
tcell <- NormalizeData(tcell, normalization.method = "LogNormalize", scale.factor = 10000)                
#scale data for PCA etc
all.genes <- rownames(tcell)
tcell <- ScaleData(tcell, features = all.genes) #,vars.to.regress = "percent.mt")
tcell <- FindVariableFeatures(tcell, selection.method = "vst")
tcell <- RunPCA(tcell, features = VariableFeatures(object = tcell))
DimPlot(tcell, reduction = "pca", group.by = "Identity", pt.size=2.5)

#Find HVG
tcell <- FindVariableFeatures(tcell, selection.method = "vst")
#PCA loadings
VizDimLoadings(tcell, dims = 1:2, reduction = "pca")
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tcell), 10)
# plot variable features with and without labels
LabelPoints(plot = VariableFeaturePlot(tcell), points = top10, repel = TRUE)

#Visualizing key genes
FeaturePlot(tcell, features = c("CCL4","IFNG","GZMB","CCL3","FABP5","DNAJB1"), pt.size=2)
FeaturePlot(tcell, features = c("CCL4","IFNG","GZMB","FABP5","SCGB1A1","IL2"), pt.size=2)
#CD8 HVG
FeaturePlot(tcell, features = c("CCL4","XCL1","CCL4L2","XCL1","XCL2","GZMB","IFNG","DNAJB1","GNLY"), pt.size=2)
#CD4 HVG
FeaturePlot(tcell, features = c("CCL4","IFNG","CCL4L2","IL2","FABP5","CCL3","IL22","XCL1","CXCL10"), pt.size=2)
#Tmem signature
FeaturePlot(tcell, features = c("CXCR6","ITGA1","ITGAE","CX3CR1","PDCD1","CRTAM"), pt.size=2)
#PCA variance 
pca = tcell@reductions[["pca"]]
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
varExplained
mat <- Seurat::GetAssayData(tcell, assay = "RNA", slot = "scale.data")
pca <- tcell[["pca"]]
# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
varExplained

#determining the best PC for UMAP
PCAcoord <- tcell[['pca']]@cell.embeddings[,1:2]
ElbowPlot(tcell, ndims=50)
DimHeatmap(tcell, dims = c(1:3, 10:18), cells = 500, balanced = TRUE)
tcell <- JackStraw(tcell, num.replicate = 100)
tcell <- ScoreJackStraw(tcell, dims = 1:20)
JackStrawPlot(tcell, dims = 1:20)

#creating UMAP
tcell <- RunUMAP(tcell, dims = 1:15)
DimPlot(tcell, reduction = "umap", group.by = "Tissue", cols = "Set2",pt.size=3)
DimPlot(tcell, reduction = "umap", group.by = "Identity", pt.size=3)
DimPlot(tcell, reduction = "umap", group.by = "Virus", pt.size=3, cols= c("blue","red"))
tcell <- FindNeighbors(tcell, dims = 1:15)
# Findings clusters in the UMAP
tclust <- FindClusters(tcell, resolution = 0.5)
DimPlot(tclust, reduction = "umap", pt.size=3, cols = "Set1")
#Finding markers per cluster
tclust.markers <- FindAllMarkers(tclust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
toptclust <- tclust.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
DoHeatmap(tclust, features = toptclust$gene, size = 3, label =F, group.colors = c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"), raster=F)
#write.table(56dtclust@meta.data, file='Cell_Cluster.csv', quote=FALSE, sep=',', col.names = T)
fullmeta <- tclust@meta.data
subset(fullmeta,fullmeta$seurat_clusters == "4") %>% group_by(Tissue) %>% tally()

#UMAP Cluster Markers
#Cluster0
FeaturePlot(tcell, features = c("CCL5", "KLRB1", "IFITM1", "ANXA1", "IL7R", "STAT1", "IRF1", "TPT1", "XAF1"), pt.size=0.25)
#Cluster1
FeaturePlot(tcell, features = c("XCL2", "XCL1", "CCL4", "EGR2", "BATF", "NFKB2", "TNFRSF9", "TNFRSF4", "IL2RG", "CRTAM", "TIGIT", "BCL2A1"), pt.size=0.25)
#Cluster2
FeaturePlot(tcell, features = c("NME1", "GZMB", "IL2RA", "IFNG", "DDX21", "ABCE1", "ILF3", "BANF1", "UBE2N", "HNRNPA2B1"), pt.size=0.25)
#Cluster3
FeaturePlot(tcell, features = c("IFITM2", "IFITM3", "IFITM1", "IL32", "LTB", "MIF", "CD226", "LGALS1"), pt.size=0.25)
#Cluster4 + GZMH
FeaturePlot(tcell, features = c("CTLA4", "CXCR4", "BNIP3L","GZMH"), pt.size=1)

FeaturePlot(tcell, features = c("NME1"), pt.size=2)
FeaturePlot(tcell, features = c("GZMB"), pt.size=2)
FeaturePlot(tcell, features = c("IL2RA"), pt.size=2)
FeaturePlot(tcell, features = c("IL7R"), pt.size=2)
FeaturePlot(tcell, features = c("CCL5"), pt.size=2)
FeaturePlot(tcell, features = c("IFITM1"), pt.size=2)
FeaturePlot(tcell, features = c("STAT1"), pt.size=2)
FeaturePlot(tcell, features = c("XCL1"), pt.size=2)
FeaturePlot(tcell, features = c("XCL2"), pt.size=2)
FeaturePlot(tcell, features = c("CCL4"), pt.size=2)
FeaturePlot(tcell, features = c("TNFRSF9"), pt.size=2)
FeaturePlot(tcell, features = c("TNFRSF4"), pt.size=2)
FeaturePlot(tcell, features = c("CRTAM"), pt.size=2)
FeaturePlot(tcell, features = c("LGALS1"), pt.size=2)
FeaturePlot(tcell, features = c("KLRB1"), pt.size=2)

# find all markers between BM groups
cluster1.markers <- FindMarkers(tclust,ident.1 = 0, ident.2 = 2, min.pct = 0.25)
top1 <- cluster1.markers %>% top_n(n = -50, wt = p_val_adj)
DoHeatmap(tclust, features = rownames(top1), cells = WhichCells(tclust,idents=c(0,2)),size = 3, label =F)
#write.table(tclust.markers, "Cluster_MarkersCD4.csv", quote=F, sep = ",")

# overall markers for CD8 or CD4
tcell.markers <- FindAllMarkers(tcell, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
tcell.markers %>% group_by(cluster) %>% top_n(n = -2, wt = p_val_adj)
top10 <- tcell.markers %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
DoHeatmap(tcell, features = top10$gene, size = 3, label =F,, raster=F)

filtmarkers <- tcell.markers %>% group_by(cluster) %>% top_n(n = -200, wt = p_val_adj)
filtmarkers <- subset(filtmarkers,filtmarkers$p_val_adj<0.05)
#write.table(filtmarkers, "MarkersCD4.csv", quote=F, sep = ",")

###################### LFC Ratio Analysis ###################################
tcell <- SetIdent(object = tcell, value = tcell@meta.data[["Tissue"]])
lungspleen.markers <- FindMarkers(tcell, ident.1 = c("lung"), 
                                  ident.2 = c("spleen"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(lungspleen.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("spleen","lung"))) + NoLegend()
sig_lungspleen <- subset(lungspleen.markers,p_val_adj < 0.1)
#write.table(lungspleen.markers, "LungSpleenDiff_CD4.csv", quote=F, sep = ",")

llnspleen.markers <- FindMarkers(tcell, ident.1 = c("LLN"), 
                                 ident.2 = c("spleen"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(llnspleen.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("spleen","LLN"))) + NoLegend()
sig_llnspleen <- subset(llnspleen.markers,p_val_adj < 0.1)
#write.table(llnspleen.markers, "LLNSpleenDiff_CD4.csv", quote=F, sep = ",")

BMspleen.markers <- FindMarkers(tcell, ident.1 = c("BM"), 
                                ident.2 = c("spleen"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(BMspleen.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("spleen","BM"))) + NoLegend()
sig_BMspleen <- subset(BMspleen.markers,p_val_adj < 0.1)
#write.table(BMspleen.markers, "BMSpleenDiff_CD4.csv", quote=F, sep = ",")

#### union of rownames from Lung, Spleen, LLN, BM
lfcspleen <- merge(lungspleen.markers, llnspleen.markers, by="row.names")
rownames(lfcspleen) <- lfcspleen$Row.names
lfcspleen <- merge(lfcspleen, BMspleen.markers, by="row.names")
rownames(lfcspleen) <- lfcspleen$Row.names
lfcspleen <- lfcspleen[,c(4,9,14)]
colnames(lfcspleen) <- c("lung_logFC","lln_logFC","BM_logFC")
lfcspleen$keyLLN <- ifelse(rownames(lfcspleen) %in% rownames(sig_lungspleen), "Sig Diff. in Lung",
                           ifelse(rownames(lfcspleen) %in% rownames(sig_llnspleen) , "Sig Diff. in LLN","NS"))
lfcspleen[which(rownames(lfcspleen) %in% rownames(sig_lungspleen) &
                  rownames(lfcspleen) %in% rownames(sig_llnspleen), arr.ind=TRUE), 4] <- "Sig Diff. in Both"
lfcspleen$keyBM <- ifelse(rownames(lfcspleen) %in% rownames(sig_lungspleen), "Sig Diff. in Lung",
                          ifelse(rownames(lfcspleen) %in% rownames(sig_BMspleen) , "Sig Diff. in BM","NS"))
lfcspleen[which(rownames(lfcspleen) %in% rownames(sig_lungspleen) &
                  rownames(lfcspleen) %in% rownames(sig_BMspleen), arr.ind=TRUE), 5] <- "Sig Diff. in Both"
lfcspleen$keyLLNBM <- ifelse(rownames(lfcspleen) %in% rownames(sig_llnspleen), "Sig Diff. in LLN",
                             ifelse(rownames(lfcspleen) %in% rownames(sig_BMspleen) , "Sig Diff. in BM","NS"))
lfcspleen[which(rownames(lfcspleen) %in% rownames(sig_llnspleen) &
                  rownames(lfcspleen) %in% rownames(sig_BMspleen), arr.ind=TRUE), 6] <- "Sig Diff. in Both"
#write.table(lfcspleen, "DE_SpleenLFC_CD8.csv", quote=F, sep = ",")

####### lung LLN
ggplot(lfcspleen %>% arrange(match(keyLLN, c("NS","Sig Diff. in Lung","Sig Diff. in LLN", "Sig Diff. in Both"))), aes(x=lung_logFC, y=lln_logFC)) +
  geom_hline(yintercept=0) +geom_vline(xintercept = 0)+geom_point(aes(color=keyLLN, shape=keyLLN), size=2) +
  labs(x="LFC(Lung/Spleen)", y="LFC(LLN/Spleen)")+
  scale_color_manual(breaks=c("Sig Diff. in Both","Sig Diff. in Lung","Sig Diff. in LLN", "NS"),
                     values = c("NS" = "gray","Sig Diff. in Lung" = "blue","Sig Diff. in LLN" = "red", "Sig Diff. in Both" = "green"))+
  scale_shape_manual(breaks=c("Sig Diff. in Both","Sig Diff. in Lung","Sig Diff. in LLN", "NS"),
                     values = c(16,17,17,16)) + theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + xlim(-2.7,2.7)+ylim(-2.6,2.6) #xlim(-2.7,2.7)+ylim(-2,2)# + xlim(-2.5,2.5)+ylim(-2.6,2.6)

####### BM lung
ggplot(lfcspleen %>% arrange(match(keyBM, c("NS","Sig Diff. in Lung","Sig Diff. in BM", "Sig Diff. in Both"))),
       aes(x=lung_logFC, y=BM_logFC)) +geom_hline(yintercept=0) +geom_vline(xintercept = 0)+ geom_point(aes(color=keyBM, shape=keyBM), size=2) +
  labs(x="LFC(Lung/Spleen)", y="LFC(BM/Spleen)")+
  scale_color_manual(breaks=c("Sig Diff. in Both","Sig Diff. in Lung","Sig Diff. in BM", "NS"),
                     values = c("NS" = "gray","Sig Diff. in Lung" = "blue","Sig Diff. in BM" = "red", "Sig Diff. in Both" = "green"))+
  scale_shape_manual(breaks=c("Sig Diff. in Both","Sig Diff. in Lung","Sig Diff. in BM", "NS"),
                     values = c(16,17,17,16)) + theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + xlim(-2.7,2.7)+ylim(-2.6,2.6)# + xlim(-2.5,2.5)+ylim(-2.6,2.6)

####### BM LLN
ggplot(lfcspleen %>% arrange(match(keyLLNBM, c("NS","Sig Diff. in LLN","Sig Diff. in BM", "Sig Diff. in Both"))),
       aes(x=lln_logFC, y=BM_logFC)) +geom_hline(yintercept=0) +geom_vline(xintercept = 0)+ geom_point(aes(color=keyLLNBM, shape=keyLLNBM), size=2) +
  labs(x="LFC(LLN/Spleen)", y="LFC(BM/Spleen)")+
  scale_color_manual(breaks=c( "Sig Diff. in Both","Sig Diff. in LLN","Sig Diff. in BM","NS"),
                     values = c("NS" = "gray","Sig Diff. in LLN" = "blue","Sig Diff. in BM" = "red", "Sig Diff. in Both" = "green"))+
  scale_shape_manual(breaks=c( "Sig Diff. in Both","Sig Diff. in LLN","Sig Diff. in BM","NS"),
                     values = c(16,17,17,16)) + theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) + xlim(-2.7,2.7)+ylim(-2.6,2.6) #+ xlim(-2.7,2.7)+ylim(-2,2)


###################### CD8 DE Analyses ######################################

spleenCD8.markers <- FindMarkers(tcell, ident.1 = c("spleen_flu_CD8"), 
                                 ident.2 = c("spleen_CMV_CD8"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(spleenCD8.markers)[1:50], size = 3,angle=0,raster=FALSE,
          cells = WhichCells(tcell,idents=c("spleen_flu_CD8","spleen_CMV_CD8"))) + NoLegend()
sig_spleenCD8 <- subset(spleenCD8.markers,spleenCD8.markers$p_val_adj < 0.1)
#write.table(spleenCD8.markers, "ALLspleenDiffCD8.csv", quote=F, sep = ",")

bloodCD8.markers <- FindMarkers(tcell, ident.1 = c("blood_flu_CD8"), 
                                ident.2 = c("blood_CMV_CD8"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(bloodCD8.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("blood_flu_CD8","blood_CMV_CD8"))) + NoLegend()
sig_bloodCD8 <- subset(bloodCD8.markers,p_val_adj < 0.1)
#write.table(bloodCD8.markers, "ALLbloodDiffCD8.csv", quote=F, sep = ",")

lungCD8.markers <- FindMarkers(tcell, ident.1 = c("lung_flu_CD8"), 
                               ident.2 = c("lung_CMV_CD8"),min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(lungCD8.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("lung_flu_CD8","lung_CMV_CD8"))) + NoLegend()
sig_lungCD8 <- subset(lungCD8.markers,p_val_adj < 0.1)
#write.table(lungCD8.markers, "ALLlungDiffCD8.csv", quote=F, sep = ",")

BMCD8.markers <- FindMarkers(tcell, ident.1 = c("BM_flu_CD8"), 
                             ident.2 = c("BM_CMV_CD8"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(BMCD8.markers)[1:15], size = 3,angle=0, raster=F,
          cells = WhichCells(tcell,idents=c("BM_flu_CD8","BM_CMV_CD8"))) + NoLegend()
sig_BMCD8 <- subset(BMCD8.markers,p_val_adj < 0.1)
#write.table(BMCD8.markers, "ALLBMDiffCD8.csv", quote=F, sep = ",")

LLNCD8.markers <- FindMarkers(tcell, ident.1 = c("LLN_flu_CD8"), 
                              ident.2 = c("LLN_CMV_CD8"), min.pct = 0,logfc.threshold = 0, min.cells.feature = 1,min.cells.group = 3)
DoHeatmap(tcell, features = rownames(LLNCD8.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("LLN_flu_CD8","LLN_CMV_CD8"))) + NoLegend()
sig_LLNCD8 <- subset(LLNCD8.markers,p_val_adj < 0.1)
#write.table(LLNCD8.markers, "ALLLLNDiffCD8.csv", quote=F, sep = ",")

allsig <- list(rownames(sig_spleenCD8),rownames(sig_lungCD8),rownames(sig_BMCD8),rownames(sig_bloodCD8),rownames(sig_LLNCD8))
combn(allsig, 2, function(x) intersect(x[[1]], x[[2]]), simplify = F)
tb <- table(unlist(allsig))
tb[tb>1]

###################### CD4 DE Analyses ######################################

spleenCD4.markers <- FindMarkers(tcell, ident.1 = c("spleen_flu_CD4"), 
                                 ident.2 = c("spleen_CMV_CD4"), min.pct = 0.25)
DoHeatmap(tcell, features = rownames(spleenCD4.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("spleen_flu_CD4","spleen_CMV_CD4"))) + NoLegend()

sig_spleenCD4 <- subset(spleenCD4.markers,spleenCD4.markers$p_val_adj < 0.1)
#write.table(spleenCD4.markers, "SpleenDiffCD4.csv", quote=F, sep = ",")

bloodCD4.markers <- FindMarkers(tcell, ident.1 = c("blood_flu_CD4"), 
                                ident.2 = c("blood_CMV_CD4"), min.pct = 0.25)
DoHeatmap(tcell, features = rownames(bloodCD4.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("blood_flu_CD4","blood_CMV_CD4"))) + NoLegend()
sig_bloodCD4 <- subset(bloodCD4.markers,p_val_adj < 0.1)
#write.table(bloodCD4.markers, "BloodDiffCD4.csv", quote=F, sep = ",")

lungCD4.markers <- FindMarkers(tcell, ident.1 = c("lung_flu_CD4"), 
                               ident.2 = c("lung_CMV_CD4"), min.pct = 0.25)
DoHeatmap(tcell, features = rownames(lungCD4.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("lung_flu_CD4","lung_CMV_CD4"))) + NoLegend()
sig_lungCD4 <- subset(lungCD4.markers,p_val_adj < 0.1)
#write.table(lungCD4.markers, "LungDiffCD4.csv", quote=F, sep = ",")

BMCD4.markers <- FindMarkers(tcell, ident.1 = c("BM_flu_CD4"), 
                             ident.2 = c("BM_CMV_CD4"), min.pct = 0.25)
DoHeatmap(tcell, features = rownames(BMCD4.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("BM_flu_CD4","BM_CMV_CD4"))) + NoLegend()
sig_BMCD4 <- subset(BMCD4.markers,p_val_adj < 0.1)
#write.table(BMCD4.markers, "BMDiffCD4.csv", quote=F, sep = ",")

LLNCD4.markers <- FindMarkers(tcell, ident.1 = c("LLN_flu_CD4"), 
                              ident.2 = c("LLN_CMV_CD4"), min.pct = 0.25)
DoHeatmap(tcell, features = rownames(LLNCD4.markers)[1:50], size = 3,angle=0,
          cells = WhichCells(tcell,idents=c("LLN_flu_CD4","LLN_CMV_CD4"))) + NoLegend()
sig_LLNCD4 <- subset(LLNCD4.markers,p_val_adj < 0.1)
#write.table(LLNCD4.markers, "LLNDiffCD4.csv", quote=F, sep = ",")

allsig <- list(rownames(sig_spleenCD4),rownames(sig_lungCD4),rownames(sig_BMCD4),rownames(sig_bloodCD4),rownames(sig_LLNCD4))
combn(allsig, 2, function(x) intersect(x[[1]], x[[2]]), simplify = F)
tb <- table(unlist(allsig))
tb[tb>1]

###################### CD8 Cross-Tissue Flu Analyses ######################################

#just flu CD8
tdiff <- subset(tcell, subset = Cell ==  "CD8" & Virus == "flu")
#Normalize by log
tdiff <- NormalizeData(tdiff, normalization.method = "LogNormalize", scale.factor = 10000)                
#scale data for PCA etc
all.genes <- rownames(tdiff)
tdiff <- ScaleData(tdiff, features = all.genes) #,vars.to.regress = "percent.mt")

bloodCD8.markers <- FindMarkers(tdiff, ident.1 = c("blood_flu_CD8"), 
                                ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(bloodCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_bloodCD8 <- subset(bloodCD8.markers,p_val_adj < 0.1)
#write.table(bloodCD8.markers, "BloodFluCD8.csv", quote=F, sep = ",")

BMCD8.markers <- FindMarkers(tdiff, ident.1 = c("BM_flu_CD8"), 
                             ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(BMCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_BMCD8 <- subset(BMCD8.markers,p_val_adj < 0.1)
#write.table(BMCD8.markers, "BMFluCD8.csv", quote=F, sep = ",")

LLNCD8.markers <- FindMarkers(tdiff, ident.1 = c("LLN_flu_CD8"), 
                              ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(LLNCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_LLNCD8 <- subset(LLNCD8.markers,p_val_adj < 0.1)
#write.table(LLNCD8.markers, "LLNFluCD8.csv", quote=F, sep = ",")

lungCD8.markers <- FindMarkers(tdiff, ident.1 = c("lung_flu_CD8"), 
                               ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(lungCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_lungCD8 <- subset(lungCD8.markers,p_val_adj < 0.1)
#write.table(lungCD8.markers, "lungFluCD8.csv", quote=F, sep = ",")

spleenCD8.markers <- FindMarkers(tdiff, ident.1 = c("spleen_flu_CD8"), 
                                 ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(spleenCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_spleenCD8 <- subset(spleenCD8.markers,p_val_adj < 0.1)
#write.table(spleenCD8.markers, "spleenFluCD8.csv", quote=F, sep = ",")

###################################### CD8 CMV ########################################

#just CMV CD8
tdiff <- subset(tcell, subset = Cell ==  "CD8" & Virus == "CMV")
#Normalize by log
tdiff <- NormalizeData(tdiff, normalization.method = "LogNormalize", scale.factor = 10000)                
#scale data for PCA etc
all.genes <- rownames(tdiff)
tdiff <- ScaleData(tdiff, features = all.genes) #,vars.to.regress = "percent.mt")

bloodCD8.markers <- FindMarkers(tdiff, ident.1 = c("blood_CMV_CD8"), 
                                ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(bloodCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_bloodCD8 <- subset(bloodCD8.markers,p_val_adj < 0.1)
#write.table(bloodCD8.markers, "BloodCMVCD8.csv", quote=F, sep = ",")

BMCD8.markers <- FindMarkers(tdiff, ident.1 = c("BM_CMV_CD8"), 
                             ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(BMCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_BMCD8 <- subset(BMCD8.markers,p_val_adj < 0.1)
#write.table(BMCD8.markers, "BMCMVCD8.csv", quote=F, sep = ",")

LLNCD8.markers <- FindMarkers(tdiff, ident.1 = c("LLN_CMV_CD8"), 
                              ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(LLNCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_LLNCD8 <- subset(LLNCD8.markers,p_val_adj < 0.1)
#write.table(LLNCD8.markers, "LLNCMVCD8.csv", quote=F, sep = ",")

lungCD8.markers <- FindMarkers(tdiff, ident.1 = c("lung_CMV_CD8"), 
                               ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(lungCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_lungCD8 <- subset(lungCD8.markers,p_val_adj < 0.1)
#write.table(lungCD8.markers, "lungCMVCD8.csv", quote=F, sep = ",")

spleenCD8.markers <- FindMarkers(tdiff, ident.1 = c("spleen_CMV_CD8"), 
                                 ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(spleenCD8.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_spleenCD8 <- subset(spleenCD8.markers,p_val_adj < 0.1)
#write.table(spleenCD8.markers, "spleenCMVCD8.csv", quote=F, sep = ",")

###################### CD4 Cross-Tissue Flu Analyses ######################################

#just flu CD4
tdiff <- subset(tcell, subset = Cell ==  "CD4" & Virus == "flu")
#Normalize by log
tdiff <- NormalizeData(tdiff, normalization.method = "LogNormalize", scale.factor = 10000)                
#scale data for PCA etc
all.genes <- rownames(tdiff)
tdiff <- ScaleData(tdiff, features = all.genes) #,vars.to.regress = "percent.mt")

bloodCD4.markers <- FindMarkers(tdiff, ident.1 = c("blood_flu_CD4"), 
                                ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(bloodCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_bloodCD4 <- subset(bloodCD4.markers,p_val_adj < 0.1)
#write.table(bloodCD4.markers, "BloodFluCD4.csv", quote=F, sep = ",")

BMCD4.markers <- FindMarkers(tdiff, ident.1 = c("BM_flu_CD4"), 
                             ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(BMCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_BMCD4 <- subset(BMCD4.markers,p_val_adj < 0.1)
#write.table(BMCD4.markers, "BMFluCD4.csv", quote=F, sep = ",")

LLNCD4.markers <- FindMarkers(tdiff, ident.1 = c("LLN_flu_CD4"), 
                              ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(LLNCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_LLNCD4 <- subset(LLNCD4.markers,p_val_adj < 0.1)
#write.table(LLNCD4.markers, "LLNFluCD4.csv", quote=F, sep = ",")

lungCD4.markers <- FindMarkers(tdiff, ident.1 = c("lung_flu_CD4"), 
                               ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(lungCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_lungCD4 <- subset(lungCD4.markers,p_val_adj < 0.1)
#write.table(lungCD4.markers, "lungFluCD4.csv", quote=F, sep = ",")

spleenCD4.markers <- FindMarkers(tdiff, ident.1 = c("spleen_flu_CD4"), 
                                 ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(spleenCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_spleenCD4 <- subset(spleenCD4.markers,p_val_adj < 0.1)
#write.table(spleenCD4.markers, "spleenFluCD4.csv", quote=F, sep = ",")

###################################### CD4 CMV ########################################

#just CMV CD4
tdiff <- subset(tcell, subset = Cell ==  "CD4" & Virus == "CMV")
#Normalize by log
tdiff <- NormalizeData(tdiff, normalization.method = "LogNormalize", scale.factor = 10000)                
#scale data for PCA etc
all.genes <- rownames(tdiff)
tdiff <- ScaleData(tdiff, features = all.genes) #,vars.to.regress = "percent.mt")

bloodCD4.markers <- FindMarkers(tdiff, ident.1 = c("blood_CMV_CD4"), 
                                ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(bloodCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_bloodCD4 <- subset(bloodCD4.markers,p_val_adj < 0.1)
#write.table(bloodCD4.markers, "BloodCMVCD4.csv", quote=F, sep = ",")

BMCD4.markers <- FindMarkers(tdiff, ident.1 = c("BM_CMV_CD4"), 
                             ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(BMCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_BMCD4 <- subset(BMCD4.markers,p_val_adj < 0.1)
#write.table(BMCD4.markers, "BMCMVCD4.csv", quote=F, sep = ",")

LLNCD4.markers <- FindMarkers(tdiff, ident.1 = c("LLN_CMV_CD4"), 
                              ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(LLNCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_LLNCD4 <- subset(LLNCD4.markers,p_val_adj < 0.1)
#write.table(LLNCD4.markers, "LLNCMVCD4.csv", quote=F, sep = ",")

lungCD4.markers <- FindMarkers(tdiff, ident.1 = c("lung_CMV_CD4"), 
                               ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(lungCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_lungCD4 <- subset(lungCD4.markers,p_val_adj < 0.1)
#write.table(lungCD4.markers, "lungCMVCD4.csv", quote=F, sep = ",")

spleenCD4.markers <- FindMarkers(tdiff, ident.1 = c("spleen_CMV_CD4"), 
                                 ident.2 = NULL, min.pct = 0.25)
DoHeatmap(tdiff, features = rownames(spleenCD4.markers)[1:50], size = 3,angle=0) + NoLegend()
sig_spleenCD4 <- subset(spleenCD4.markers,p_val_adj < 0.1)
#write.table(spleenCD4.markers, "spleenCMVCD4.csv", quote=F, sep = ",")