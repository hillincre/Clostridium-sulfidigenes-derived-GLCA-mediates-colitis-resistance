setwd("path")

library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(patchwork)

out_dir <- "./pseudotime_ILC3"
dir.create(out_dir, showWarnings = FALSE)

group_cols <- c("Con" = "#7b95c6", "GLCA" = "#f0486f")

obj <- readRDS("./anchored_data_processed.rds")
DefaultAssay(obj) <- "RNA"
ly <- Layers(obj, assay = "RNA")
if (!"counts" %in% ly && any(grepl("^counts\\.", ly))) obj <- JoinLayers(obj, assay = "RNA")

print(table(obj$group, obj$celltype))

TRAJ_CELLTYPES <- c("ILC3")

sub <- subset(obj, celltype %in% TRAJ_CELLTYPES)
sub$celltype <- factor(sub$celltype)
Idents(sub) <- sub$celltype
print(table(sub$group))

sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
hvg <- VariableFeatures(sub)

ap_program_genes <- c(
  "Cd74","H2-Aa","H2-Ab","Cd86",
  "Cd69","Cd44","Icos","Il2ra","Cd38","Cd2","Thy1",
  "Mki67",
  "Il23r","Il1r1","Il18r1","Tnfrsf25"
)

cytokine_genes <- c("Il17a","Il17f","Il22","Csf2","Tnf","Ifng","Lta","Ltb","Il13")

tf_genes <- c("Rorc","Ahr","Id2","Tcf7","Tox","Gata3","Ikzf1","Ikzf3",
              "Notch1","Hif1a","Myc")

receptor_genes <- c("Il7r","Kit","Il23r","Il1r1","Il18r1","Il2ra","Il2rb",
                    "Ccr6","Cxcr5","Ncr1","Icos","Cd2","Thy1",
                    "H2-Aa","H2-Ab","Tnfrsf25","Tnfsf11","Tnfsf4","Klrd1")

activation_markers <- c("Cd69","Cd44","Mki67","Cd38","Cd63","Cd74","Cd86",
                        "Klrg1","Cd200")

resting_genes <- c("Sell","Ccr7","Lef1","Klf2","S1pr1","Tcf7","Il7r")

activation_score_genes <- ap_program_genes
cytokine_score_genes  <- cytokine_genes
resting_score_genes   <- resting_genes

all_panel <- unique(c(ap_program_genes, cytokine_genes, tf_genes,
                      receptor_genes, activation_markers, resting_genes))
all_panel <- intersect(all_panel, rownames(sub))

detect_genes <- c("Cd74","H2-Aa","H2-Ab","Cd86","Cd69","Cd44","Icos","Il2ra",
                  "Mki67","Il23r","Il17a","Il17f","Il22","Csf2","Tnf","Ifng",
                  "Il7r","Sell")
detect_genes <- intersect(detect_genes, rownames(sub))
m_counts <- GetAssayData(sub, layer = "counts")
detect_pct <- sapply(detect_genes, function(g) round(100 * mean(m_counts[g, ] > 0), 1))
print(detect_pct)

genes_keep <- unique(c(hvg, all_panel))

sub <- AddModuleScore(sub, features = list(intersect(activation_score_genes, rownames(sub))),
                      name = "Activation_", nbin = 24, ctrl = 100)
sub <- AddModuleScore(sub, features = list(intersect(resting_score_genes, rownames(sub))),
                      name = "Resting_", nbin = 24, ctrl = 100)
sub <- AddModuleScore(sub, features = list(intersect(cytokine_score_genes, rownames(sub))),
                      name = "Cytokine_", nbin = 24, ctrl = 100)
colnames(sub@meta.data)[colnames(sub@meta.data) == "Activation_1"] <- "activation_score"
colnames(sub@meta.data)[colnames(sub@meta.data) == "Resting_1"]    <- "resting_score"
colnames(sub@meta.data)[colnames(sub@meta.data) == "Cytokine_1"]   <- "cytokine_score"

counts <- as.matrix(GetAssayData(sub, layer = "counts")[genes_keep, ])
pd <- new("AnnotatedDataFrame", data = sub@meta.data)
fd <- new("AnnotatedDataFrame",
          data = data.frame(gene_short_name = genes_keep, row.names = genes_keep))

cds <- newCellDataSet(counts, phenoData = pd, featureData = fd,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, ordering_genes = hvg[hvg %in% genes_keep])

cds <- reduceDimension(cds, max_components = 2, method = "DDRTree",
                       norm_method = "log", pseudo_expr = 1, verbose = FALSE)
cds <- orderCells(cds)

root_cell  <- rownames(pData(cds))[which.max(pData(cds)$resting_score)]
root_state <- pData(cds)[root_cell, "State"]
cds <- orderCells(cds, root_state = root_state)

cor_act <- cor(pData(cds)$Pseudotime, pData(cds)$activation_score, use = "complete.obs")
cor_rest <- cor(pData(cds)$Pseudotime, pData(cds)$resting_score, use = "complete.obs")
cor_cyto <- cor(pData(cds)$Pseudotime, pData(cds)$cytokine_score, use = "complete.obs")

state_pt <- tapply(pData(cds)$Pseudotime, pData(cds)$State, median, na.rm = TRUE)
state_pt <- sort(state_pt)
print(state_pt)

grp <- cut(seq_along(state_pt), breaks = 3,
           labels = c("1", "2", "3"))
map <- setNames(as.character(grp), names(state_pt))
print(map)

pData(cds)$State <- factor(
  map[as.character(pData(cds)$State)],
  levels = c("1", "2", "3")
)
print(table(pData(cds)$State))

p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = TRUE,
                           cell_size = 1.5) + theme_bw(base_size = 8) + labs(title = "Pseudotime")
ggsave(file.path(out_dir, "1_trajectory_pseudotime.pdf"), p1,
       device = cairo_pdf, width = 7, height = 5, units = "cm")

p2 <- plot_cell_trajectory(cds, color_by = "group", show_branch_points = TRUE,
                           cell_size = 1.5) + scale_color_manual(values = group_cols) +
  theme_bw(base_size = 8) + labs(title = "Group")
ggsave(file.path(out_dir, "2_trajectory_group.pdf"), p2,
       device = cairo_pdf, width = 7, height = 5, units = "cm")

p3 <- plot_cell_trajectory(cds, color_by = "State", show_branch_points = TRUE,
                           cell_size = 1.5) + theme_bw(base_size = 8) + labs(title = "State")
ggsave(file.path(out_dir, "3_trajectory_state.pdf"), p3,
       device = cairo_pdf, width = 7, height = 5, units = "cm")

p1/p2/p3
ggsave(file.path(out_dir, "1+2+3_trajectory.pdf"), p1/p2/p3,
       device = cairo_pdf, width = 6, height = 12, units = "cm")


sig_genes <- rownames(deg_pt)[1:min(50, nrow(deg_pt))]
pdf(file.path(out_dir, "4_pseudotime_heatmap.pdf"), width = 10/2.53, height = 13/2.53)
plot_pseudotime_heatmap(cds[sig_genes, ], num_clusters = 3, cores = 1, show_rownames = TRUE)
dev.off()

saveRDS(cds, file.path(out_dir, "ILC3_pseudotime_cds.rds"))
write.csv(d[, c("celltype","group","orig.ident","Pseudotime","State",
                "activation_score","resting_score","cytokine_score")],
          file.path(out_dir, "ILC3_pseudotime_cell_scores.csv"))