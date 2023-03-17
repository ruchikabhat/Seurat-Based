#devtools::install_github("sqjin/CellChat")
#install.packages('NMF')
#devtools::install_github("jokergoo/circlize")
#devtools::install_github("jokergoo/ComplexHeatmap")
######################################         Script starts          ##############################
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
load("D:/HPC_DATA_MOUSE_NEW COVID_April2022/Rrun/data_humanSkin_CellChat.rda")
data.input = data_humanSkin$data # normalized data matrix
meta = data_humanSkin$meta # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
View(data_humanSkin$meta)
View(data.input@Dimnames)
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
View(meta)
unique(meta$labels)


setwd ("D:/HPC_DATA_MOUSE_NEW COVID_April2022/Rrun")
mouse<-readRDS("Allmouse_updatedcelltypes_14June.rds") 
cell.use = rownames(mouse)[mouse$Condition == "IL7_treated"]
meta = data.frame(labels = mouse$updated_subtypes[cell.use], row.names = colnames(mouse)) # manually create a dataframe consisting of the cell labels

Idents(mouse)<-"Condition"
data.input<-subset(mouse,idents= c('IL7_treated'), invert=FALSE)
View(data.input@meta.data)
cellchat <- createCellChat(object = data.input,group.by = "updated_subtypes")
cellchat <- setIdent(cellchat, ident.use = "updated_subtypes")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
View(groupSize)
CellChatDB<-CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# Extra line:
CellChatDB.use <- CellChatDB
## extraline ends
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM_Receptor")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat@data.signaling
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

saveRDS(cellchat, file = "cellchat_IL7treated_ECMR.rds")

adult<-readRDS("cellchat_Adult.rds")
old<-readRDS("cellchat_Old.rds")
Il7<-readRDS("cellchat_IL7treated.rds")
View(adult@meta)
levels(cellchat@meta$updated_subtypes)
object.list <- list(Adult = adult, Old = old, IL7=Il7)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

group.cellType <- c(rep("FRCs", 4), rep("LECs", 4), rep("BECs", 4))
group.cellType <- factor(group.cellType, levels = c("FRCs", "LECs", "BECs"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
gg
df$labels = factor(df$labels,levels = levels(object@idents))
df$labels <- factor(df$labels, levels = names(incoming.cells))
cellchatMO <- computeCommunProbPathway(cellchatMO)
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

cellchatMO <- computeNetSimilarityPairwise(cellchatMO, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchatMO <- netEmbedding(cellchatMO, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchatMO <- netClustering(cellchatMO, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchatMO, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


#####################33333333




























gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}




gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,  comparison = c(2,1,3))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,  comparison = c(2,1,3))
gg1 + gg2

netVisual_bubble(cellchat, signaling = "MK",  comparison = c(2,1,3), angle.x = 45)
#> Comparing communications on a merged object




par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]],signaling = "VPTN", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("VPTN ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]],signaling = "PTN", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("PTN ", names(object.list)[1]))
netVisual_chord_gene(object.list[[3]],signaling = "PTN", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("PTN ", names(object.list)[3]))


netVisual_chord_gene(object.list[[2]],signaling = "PSAP", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("PSAP ", names(object.list)[1]))

netVisual_chord_gene(object.list[[1]],signaling = "PDGF", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("PDGF ", names(object.list)[1]))


netVisual_chord_gene(object.list[[3]],signaling = "MK", slot.name = 'net', lab.cex = 0.8, small.gap = 2, title.name = paste0("MK ", names(object.list)[3]))



netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



#########################################################################################################3

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
#> 
#> # Chord diagram
###group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
##names(group.cellType) <- levels(cellchat@idents)
###netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> 
 netAnalysis_contribution(cellchat, signaling = pathways.show)
 
 # Access all the signaling pathways showing significant communications
 pathways.show.all <- cellchat@netP$pathways
 # check the order of cell identity to set suitable vertex.receiver
 levels(cellchat@idents)
 vertex.receiver = seq(1,4)
 for (i in 1:length(pathways.show.all)) {
   # Visualize communication network associated with both signaling pathway and individual L-R pairs
   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
 }
 
 # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
 netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
 #> Comparing communications on a single object
 # show all the significant interactions (L-R pairs) associated with certain signaling pathways
 netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
 #> Comparing communications on a single object
 #> 
 # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
 # show all the interactions sending from Inflam.FIB
 netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
 netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
 plotGeneExpression(cellchat, signaling = "CXCL")
 cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
 netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
 # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
 gg1 <- netAnalysis_signalingRole_scatter(cellchat)
 #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
 # Signaling role analysis on the cell-cell communication networks of interest
 gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
 #> Signaling role analysis on the cell-cell communication network from user's input
 gg1 + gg2
 ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
 ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
 ht1 + ht2
 nPatterns = 3
 cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
 # river plot
 library(ggalluvial)
 netAnalysis_river(cellchat, pattern = "outgoing")
 #> Please make sure you have load `library(ggalluvial)` when running this function 
 netAnalysis_river(cellchat, pattern = "incoming")
 nPatterns = 4
 cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
 netAnalysis_dot(cellchat, pattern = "incoming") 
 