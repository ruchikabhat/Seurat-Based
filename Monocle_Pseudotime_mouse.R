library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(patchwork)
library('reticulate')
library(Matrix)
library(cowplot)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(patchwork)
library(BiocManager)
#library(multtest)
#library(metap)
#library(SingleR)
library(EnhancedVolcano)
memory.limit(size=7800000000000)
library(RColorBrewer)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
a<-DiscretePalette(n, palette = NULL)
library(ggplot2)
#BiocManager::install("monocle")
library(monocle)

######load the seurat object first ######
setwd ("D:/HPC_DATA_MOUSE_NEW COVID_April2022/Rrun")
#mouse<-readRDS("Allmouse_updatedcelltypes_14June.rds") #just one version before newest
mouse<-readRDS("mouseAnnotated_16thMarch2023.rds")     ##NEWEST ##
FibroMono<-readRDS("Pseudotime_onlyFRCs.rds")
Idents(mouse)<-"seurat_clusters"

FRC<-subset(Frcd,idents= c('FRCs'), invert=FALSE)
Idents(FRC)<-"seurat_clusters"
cluster_letters <- Idents(object = mouse)
names(cluster_letters) <- colnames(x = mouse)
FRC <- AddMetaData(
  object = mouse,
  metadata = cluster_letters,
  col.name = 'Clusters'
)
View(FRC@meta.data)
Idents(FRC)<-'Clusters'
FRC <- RenameIdents(FRC, '0' = 'Ccl19low Type1', '1' = 'Cd34+ Type2', '3'='Inmt+', '5'='Nr4a1+','7'='Ccl19high','4'='Ccl19low Type2','17'='Cd34+ Type1','6'='Enthothelial Cells','18'='NerveCells','9'='SmoothMuscle','12'='Pericytes','2'='CapEC1','8'='fLECs','10'='Art','11'='cLECs','13'='TrEC+HEV','14'='CapEC2','15'='Macro/Ptx3 LECs','16'='Cd45+Bcells','19'='Cd45+Tcells','20'='CapIfn','21'='Cd45+Bcells')
cluster_letters <- Idents(object = FRC)
names(cluster_letters) <- colnames(x = FRC)
FRC <- AddMetaData(
  object = FRC,
  metadata = cluster_letters,
  col.name = 'Clusters'
)
Idents(Frcd)<-'Condition'
Ad<-subset(Frcd,idents= c('Adult'), invert=FALSE)
Od<-subset(Frcd,idents= c('Old'), invert=FALSE)
saveRDS(FRC, "mouseAnnotated_16thMarch2023.rds")
View(FRC@meta.data)
#seuratobject<-readRDS(YOURPATH/SEURATRDSfile)
FRC<-Od
Idents(FRC)<-"Clusters"
data <- as(as.matrix(FRC@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = FRC@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))# save list of gene names
fd <- new('AnnotatedDataFrame', data = fData)

#install.packages("vctrs")

#cds <- new_cell_data_set(data,                         cell_metadata = pd                         gene_metadata = fd)

#FibroMono <- newCellDataSet(data,
         #                   phenoData  = pd,
          #                  featureData = fd,expressionFamily = negbinomial.size())
FibroMono <- newCellDataSet(data,
                            phenoData  = pd,
                            featureData = fd,expressionFamily = uninormal())

pData(FibroMono)
fData(FibroMono)


#Run ordering algorithm
var_genes <- FRC[["RNA"]]@var.features
ordering_genes <- var_genes

FibroMono <- setOrderingFilter(FibroMono, ordering_genes)
print(dim(exprs(FibroMono)))


## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
FibroMono <- reduceDimension(FibroMono,norm_method="none", 
                             reduction_method="DDRTree",
                             max_components=4,
                             scaling=TRUE,
                             verbose=TRUE,
                             pseudo_expr=0)

# First decide what you want to color your cells by
print(head(pData(FibroMono)))

## order cells change colors and theta to match your plot
FibroMono <- orderCells(FibroMono)




plot_cell_trajectory(FibroMono, 
                     color_by = "updated_subtypes",
                     theta = -15,
                     show_branch_points = TRUE,show_backbone = TRUE,
                     backbone_color = "black",
                     show_tree = TRUE,use_color_gradient = FALSE,
                     cell_size = 1.5, values=a) + scale_color_manual(values = a, name = "Only FRCs Pseudotime")+ theme(legend.position = "top")

saveRDS(FibroMono, "Pseudotime.rds")
# now onwards to pseudotemporal plots
genex <- c("Inmt")

sig_gene_names <- (genex)
head(sig_gene_names)
pseudotemporalplot <- plot_pseudotime_heatmap(FibroMono[my_genes],
                                              num_clusters = 9, 
                                              cores = 4,
                                              hmcols = NULL,
                                              show_rownames = T)

pseudotemporalplot 
print(pData(cds_subset)$)

###################################### OLd one ##################3

####ESTIMATE DISPERSIONS AND SIZE FACTORS############
FibroMono <- estimateSizeFactors(FibroMono)
FibroMono <- estimateDispersions(FibroMono)


#######FILTER LOW QUALITY GENES/CELLS##########
FibroMono <- detectGenes(FibroMono, min_expr = 0.1)
print(head(fData(FibroMono)))

expressed_genes <- row.names(subset(fData(FibroMono),
                                    num_cells_expressed >= 10))
######DEG TEST##########


diff_test_res <- differentialGeneTest(FibroMono[1:100,],
                                      fullModelFormulaStr = "~updated_subtypes")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

########### NEW PART FROM ANOTHER SCRIPT ##############
var.genes <- VariableFeatures(mouse)
FibroMono <- setOrderingFilter(FibroMono, var.genes)
p1 <- plot_ordering_genes(FibroMono)

disp_table <- dispersionTable(FibroMono)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
FibroMono <- setOrderingFilter(FibroMono, disp.genes)
p2 <- plot_ordering_genes(FibroMono)

p1
p2
#######################################
FibroMono<- setOrderingFilter(FibroMono, ordering_genes)
FibroMono <- reduceDimension(FibroMono, max_components = 2,method = 'DDRTree', auto_param_selection = F)
FibroMono <- orderCells(FibroMono)
png(filename="Trajectory.png")
plot_ordering_genes(FibroMono)
dev.off()
png(filename="TrajectorybyGroup.png")
plot_cell_trajectory(FibroMono, color_by = "group")
dev.off()


###SIGNIFICANT GENES CHANGING WITH PSEUDOTIME##########
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
diff_test_res <- differentialGeneTest(FibroMono[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
FibroMono <- orderCells(FibroMono, root_state = GM_state(FibroMono))

###PLOT####
png(filename="TrajectorybyPSEUDOTIME.png")
plot_cell_trajectory(FibroMono, color_by = "Conditon")
dev.off()
png(filename="PSEUDOTIMEHEATMAP.png")
plot_pseudotime_heatmap(FibroMono[my_genes,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
dev.off()


data(FibroMono)
my_genes <- row.names(subset(fData(FibroMono), gene_short_name %in% c("Inmt", "Nr4a1", "Ccl19")))
my_genes
cds_subset <- FibroMono[my_genes,]
cds_subset
plot_genes_in_pseudotime(cds_subset, color_by="updated")
