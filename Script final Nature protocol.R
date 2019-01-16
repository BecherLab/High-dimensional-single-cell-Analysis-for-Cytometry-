#Script adapted from Nowicka et al.

#Load the content of the metadata.xlsx file into R (see Nowicka et al for description of "metadata")
metadata_filename <- "metadata.xlsx"
library(readxl)
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$condition <- factor(md$condition, levels = c("HD","Blood", "Tumor"))
head(data.frame(md))

# Define colors for conditions 
color_conditions <- c( "#278210","#fc0008", "#020202") 
  names(color_conditions) <- levels(md$condition)

# Load the content of the .fcs files into R
library(flowCore)
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

# Load the content of the channel.names.xlsx file into R
channel.names <- read_excel("channel.names.xlsx")
new.names <- channel.names$Antigen
colnames(fcs_raw) <-new.names
panel_filename <- "channel.names.xlsx"
panel <- read_excel(panel_filename)
head(data.frame(panel))

# Replace problematic characters (to eliminate potential error source)
panel$Antigen <- gsub("-", "_", panel$Antigen)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
panel_fcs$desc <- gsub("-", "_", panel_fcs$name)

# Define Lineage markers (or other classification)
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Define Functional markers (or other classification)
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Define All markers in your panel (for global analysis)
(all_markers <- panel$Antigen[panel$All == 1])

# Spot checks for markers
all(lineage_markers %in% panel_fcs$name)
all(functional_markers %in% panel_fcs$name)
all(all_markers %in% panel_fcs$name)

# Column subsetting
fcs <- fsApply((fcs_raw[, c(all_markers)]), exprs)
head(fcs)

#---------------------------------------------------------------------

# Optional, if you are not using CyTOBANK/MATLAB it's possible to do arcsinh transformation of the channels included for the clustering in R
# Please, check if the files were correctly exported in FlowJo, sometines files have a scaling bug
# PLease compared the transformed data to the original .fcs file to check if correct transformation factor was choosen


# Cofactor 100 
asinh_scale_100 <- c("CCR8")
cofactor_1 <- 100
fcs[,asinh_scale_100] <- asinh(fcs[,asinh_scale_100] / cofactor_1)

# Cofactor 200 
asinh_scale_200 <- c("GNLY")
cofactor_2 <- 200
fcs[,asinh_scale_200] <- asinh(fcs[,asinh_scale_200] / cofactor_2)

# Cofactor 300 
asinh_scale_300 <- c("PD1")
cofactor_3 <- 300
fcs[,asinh_scale_300] <- asinh(fcs[,asinh_scale_300] / cofactor_3)

# Cofactor 400 
asinh_scale_400 <- c("CD39", "CD244")
cofactor_4 <- 400
fcs[,asinh_scale_400] <- asinh(fcs[,asinh_scale_400] / cofactor_4)

# Cofactor 500
asinh_scale_500 <- c("CD27", "CCR7", "CD25", "CD4", "CD95", "BTLA")
cofactor_5 <- 500
fcs[,asinh_scale_500] <- asinh(fcs[,asinh_scale_500] / cofactor_5)

# Cofactor 800
asinh_scale_800 <- c("CD127", "CD45RA", "CD57", "CCR4")
cofactor_6 <- 800
fcs[,asinh_scale_800] <- asinh(fcs[,asinh_scale_800] / cofactor_6)

# Cofactor 1000
asinh_scale_1000 <- c("CXCR5", "CD38", "CD71")
cofactor_7 <- 1000
fcs[,asinh_scale_1000] <- asinh(fcs[,asinh_scale_1000] / cofactor_7)

# Cofactor 1100
asinh_scale_1100 <- c("ICOS")
cofactor_8 <- 1100
fcs[,asinh_scale_1100] <- asinh(fcs[,asinh_scale_1100] / cofactor_8)

# Cofactor 1500
asinh_scale_1500 <- c("HLADR", "Tbet")
cofactor_9 <- 1500
fcs[,asinh_scale_1500] <- asinh(fcs[,asinh_scale_1500] / cofactor_9)

# Cofactor 1800
asinh_scale_1800 <- c("Ki67", "CD8")
cofactor_10 <- 1800
fcs[,asinh_scale_1800] <- asinh(fcs[,asinh_scale_1800] / cofactor_10)

# Cofactor 3000
asinh_scale_3000 <- c("CD28")
cofactor_11 <- 3000
fcs[,asinh_scale_3000] <- asinh(fcs[,asinh_scale_3000] / cofactor_11)

#---------------------------------------------------------------------

# Transform the markers to values between 0-1 (Normalization for FlowSOM algorithm)
library(matrixStats)
rng <- colQuantiles(fcs, probs = c(0.01, 0.99))
expr01 <- t((t(fcs) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

# Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))
library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, fcs)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

theme_bw2 <- function (base_size = 14, base_family = "", face="bold") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.background = element_rect(fill = "white", 
                                          colour = NA), panel.border = element_rect(fill = NA, 
                                                                                    colour = "NA"), panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "NA", 
                                          size = 0.25), strip.background = element_rect(fill = "NA", 
                                                                                        colour = "NA"), legend.key = element_rect(fill = "white", 
                                                                                                                                  colour = NA), complete = TRUE)
}

plot1 <- ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") + theme_bw2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions)
plot1

ggsave(filename = "expression_matrix.png", plot = plot1, 
       scale = 1, width = 10, height = 7, units = c("in"))

# The next three steps are a data normalization for running a tSNE analysis
# Adjust data to start from (approximately) 0
q.vector <- apply(fcs[,all_markers], 2, function(x) quantile(x, 0.01, names = F))
q.vector
data.shift <- fcs
data.shift[,all_markers] <- sweep(data.shift[,all_markers], 2, q.vector)
head(data.shift)

# Normalize to have everything from 0 to 1
per.vector <- apply(data.shift[,all_markers], 2, function(x) quantile(x, 0.99999, names = F))
per.vector
data.shift[,all_markers] <- t(t(data.shift[,all_markers]) / as.numeric(per.vector))

# Check whether you adjusted the range approximately from 0 to 1
apply(data.shift, 2, min)
apply(data.shift, 2, max)

# Check successful transformation for all markers
# Example "CD4" vs "CD8"
# Show biaxial plot
ggplot(data = data.frame(data.shift[1:90000,]), aes(x =CD4, y = CD8)) + geom_point(alpha=1/2)


# Find and skip duplicates
dups <- which(!duplicated(data.shift[, lineage_markers]))

# Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

# Optional for real big data sets: How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 10000)

# Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

expr <- as.data.frame(data.shift)
expr$cell.id <- 1:nrow(expr)

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]

## Run t-SNE

if(!require(devtools)) install.packages("devtools") 

# In this step, you can eventually run Rtsne multicore If not already installed, use devtools::install_github("RGLab/Rtsne.multicore") <-- not availale for macOS.
# In the script you should replate Rtsne.multicore for Rtsne

library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, dims = 2, perplexity = 100, theta = 0.5, 
                            max_iter = 2000, verbose = T, pca = F, check_duplicates=F)

# Plot t-SNE colored by single expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                 expr[tsne_inds, lineage_markers])


# Prepare the expression data
expr$cell.id  <- as.factor(expr$cell.id)
data.ix.df <- data.frame(expr[tsne_inds,c(lineage_markers, "cell.id")])
library(reshape2)
data.melt <- melt(data.ix.df, variable.name = "antigen", value.name = "expression")
dr$cell.id <-expr[tsne_inds,c("cell.id")]
dim(dr)

# Create a vector called dr2 with tSNE1, tSNE2, "cell.id" values
dr2 <-dr[,c(1:2, 27)]
joined.expr <- merge(data.melt, dr2, by = "cell.id")

# Plot tSNE plot with with blue color palette
Blue_palette <- c('#2F2C62', '#42399B', '#4A52A7', '#A7DA64',
             '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')

plot2 <- ggplot(joined.expr, aes(x = tSNE1, y = tSNE2, color = expression)) +
  geom_point(size = 0.04) +
  theme_bw2() +
  scale_color_gradientn(colours = Blue_palette,
                        limits = c(-0, 1)) +
  facet_wrap(~ antigen, ncol = 4, scales = "free") 
plot2

ggsave(filename = "tSNE_allmarkers.png", plot = plot2, 
       scale = 1, width = 10, height = 11, units = c("in"))


## Flowsome

library(FlowSOM)
head(fcs)

fsom <- ReadInput(flowFrame(exprs = fcs, desc = list(FIL = 1)), transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

# Metaclustering with ConsensusClusterPlus
library(ConsensusClusterPlus)
codes <- som$map$codes
plot_outdir <- "consensus_plots"

# Metaclustering K30
nmc <- 30
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", distance = "euclidean", seed = 1234)

# Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

color_clusters <- c("#56ff0d", "#1965B0", "#7BAFDE", "#DC050C", "#882E72",
                    "#FF7F00", "#B17BA6", "#E7298A","#FDB462",  "#E78AC3",
                    "#0a0d05", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#ffff00", "#33A02C", 
                    "#8B3800", "#4E8500", "#33A02C" )  

plot_clustering_heatmap_wrapper <- function(fcs, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL) {
  
# Calculate the median expression
  expr_median <- data.frame(fcs, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  
# Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
# This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(fcs)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
# Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering)) 
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }  
  
# Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 0.001,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)
}

plot_clustering_heatmap_wrapper(fcs = fcs[, lineage_markers],
                                expr01 = expr01[, lineage_markers],
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)



# Plot Flowsome 
dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

# Plot t-SNE colored by clusters
plot3 <-ggplot(dr, aes(x = tSNE1, y = tSNE2,color = cell_clustering1)) +
  geom_density_2d(data = dr[,c(1,2)], aes(x = tSNE1, y = tSNE2), colour = "lightgrey", size =0.5, bins=30) +
  geom_point(size = 1) +
  theme_bw2() + coord_fixed(ratio = 1)+  
  ylim(-45,45) + xlim(-45,45) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
plot3

# Facet per condition
plot3 + facet_wrap(~ condition)

# and you count k30 cluster
counts_table30 <- table(cell_clustering1, sample_ids)
props_table30 <- t(t(counts_table30) / colSums(counts_table30)) * 100
counts30 <- as.data.frame.matrix(counts_table30)
props30 <- as.data.frame.matrix(props_table30)
props30$cluster <- c(1:30)

ggdf <- melt(props30, id.vars = "cluster", value.name = "Frequency", variable.name = "sample_id")


# Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$cluster <- factor(ggdf$cluster, levels = c(1:length(ggdf$cluster)), ordered = T)

plot4 <-ggplot(ggdf, aes(x = sample_id, y = Frequency, fill = cluster)) + geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = color_clusters)
plot4 

## Manual metaclustering
#Load the content of the cluster_mergings.xlsx file into R (see Nowicka et al for description of "cluster_mergings")

cluster_merging_filename <- "cluster_mergings.xlsx" 
cluster_merging1 <- read_excel(cluster_merging_filename) 
data.frame(cluster_merging1)
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

# You update the t-SNE plot with the new annotated cell populations.
dr$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])

color_clusters2 <- c("#DC050C", "#33A02C", "#1965B0", "#7BAFDE", "#882E72",
                     "#FF7F00", "#B17BA6", "#FDB462", "#E7298A", "#E78AC3",
                     "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                     "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                     "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                     "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
                     "#8B3800", "#0a0d05", "#4E8500")


plot5 <-ggplot(dr, aes(x = tSNE1, y = tSNE2,color = cell_clustering1m)) +
  geom_density_2d(data = dr[,c(1,2)], aes(x = tSNE1, y = tSNE2), colour = "lightgrey", size =0.5, bins=30) +
  geom_point(size = 0.8) +
  theme_bw2() + coord_fixed(ratio = 1)+  
  ylim(-45,45) + xlim(-45,45) +
  scale_color_manual(values = color_clusters2) +
  guides(color = guide_legend(override.aes = list(size = 4)))
plot5

# Facet per condition
plot5 + facet_wrap(~ condition)


# Heatmap Metaclastering
dr_sub <-dr[,c(lineage_markers, "cell_clustering1m")]
head(dr_sub)

median_subset <- dr_sub %>%
  group_by(cell_clustering1m) %>%
  summarize_all(funs(median(., na.rm=TRUE)))

class(median_subset)
median_subset <- as.data.frame(median_subset)
median_subset_2 <- as.matrix(sapply(median_subset[, -1], as.numeric))
rownames(median_subset_2) <- median_subset$cell_clustering1m
class(median_subset_2)

# Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(3)

pheatmap(median_subset_2, color = color,
         number_color = "black", fontsize_number = 5, clustering_method = "average",
         cluster_cols = FALSE,
         border_color = "black",fontsize = 10,
         cellwidth = 20, cellheight = 15,
         display_numbers = matrix(ifelse(median_subset_2 > 5, "*", ""), nrow(median_subset_2)))
dev.off()


# and you count cluster mergin
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "Frequency", variable.name = "sample_id")

# Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])

ggplot(ggdf, aes(x = sample_id, y = Frequency, fill = cluster)) + geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = color_clusters2)

ggdf$patient_id <- factor(md$patient_id[mm])

# example sub-sample "CD8+ CM T cells"
ggdf_CM <- ggdf[ggdf$cluster %in% c("CD8+ CM T cells"), ]

ggplot(ggdf_CM) +
  geom_boxplot(aes(x = condition, y = Frequency, color = condition,
                   fill = condition), position = position_dodge(), alpha = 0.5,
               outlier.color = NA) +
  geom_point(aes(x = condition, y = Frequency, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) + facet_wrap(~ cluster, scales = "free", nrow = 2) +
  theme_bw2() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 10, face="bold")) + scale_color_manual(values = color_conditions) + scale_fill_manual(values = color_conditions) + 
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33))


## Differential analysis of marker expression stratified by cell population

# Get median marker expression per sample and cluster (normalized to 0-1)
expr_median_sample_cluster_tbl <- data.frame(expr[, functional_markers],
                                             sample_id = sample_ids, cluster = cell_clustering1m) %>%
  group_by(sample_id, cluster) %>%
  summarize_each(funs(median))
expr_median_sample_cluster_melt <- melt(expr_median_sample_cluster_tbl,
                                        id.vars = c("sample_id", "cluster"), value.name = "median_expression",
                                        variable.name = "antigen")

# Rearange so the rows represent clusters and markers
expr_median_sample_cluster <- dcast(expr_median_sample_cluster_melt,
                                    cluster + antigen ~ sample_id, value.var = "median_expression") 
rownames(expr_median_sample_cluster) <- paste0(expr_median_sample_cluster$cluster,
                                               "_", expr_median_sample_cluster$antigen)

# Optionl: Eliminate clusters with low frequency (e.g. outliers)
clusters_keep <- names(which((rowSums(counts < 0) == 0)))
keepLF <- expr_median_sample_cluster$cluster %in% clusters_keep
expr_median_sample_cluster <- expr_median_sample_cluster[keepLF, ]

# Eliminate cases with zero expression in all samples
keep0 <- rowSums(expr_median_sample_cluster[, md$sample_id]) > 0
expr_median_sample_cluster <- expr_median_sample_cluster[keep0, ]
ggdf <- expr_median_sample_cluster_melt[expr_median_sample_cluster_melt$cluster 
                                        %in% clusters_keep, ]

# Add info about samples
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$condition[mm])
ggdf$patient_id <- factor(md$patient_id[mm])

ggdf_CD8CM <- ggdf[ggdf$cluster %in% c("CD8+ CM T cells"), ]
ggdf_PD1CM <- ggdf_CD8CM[ggdf_CD8CM$antigen %in% c("PD1"), ]
ggplot(ggdf_PD1CM) +
  geom_boxplot(aes(x = antigen, y = median_expression,
                   color = condition, fill = condition),
               position = position_dodge(), alpha = 0.5, outlier.color = NA) +
  geom_point(aes(x = antigen, y = median_expression, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge(),
             size = 1.1) +
  facet_wrap(~ cluster, scales = "free_y", ncol=5) +
  theme_bw2() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face="bold"), strip.text = element_text(size = 10, face="bold")) +  scale_color_manual(values = color_conditions) + 
  scale_fill_manual(values = color_conditions) + scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39))


save.image(file = "natprot.RData", compress = T)
