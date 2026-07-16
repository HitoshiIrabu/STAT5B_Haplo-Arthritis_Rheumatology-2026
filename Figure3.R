#!/usr/bin/env Rscript

INPUT_RDS          <- "seurat_object.rds"
CELL_META_CSV      <- "cell_metadata_with_tls_dist.csv"
JIA_EXAMPLE_SAMPLE <- "JIA5"   # representative JIA section for the spatial overlays (Panels H, K)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(RANN)
  library(patchwork)
  library(scales)
  library(ggrastr)
  library(ggrepel)
  library(Matrix)
  library(grid)
})

run_nn2 <- function(data_xy, query_xy, k) {
  data_xy  <- as.matrix(data_xy)
  query_xy <- as.matrix(query_xy)
  RANN::nn2(data = data_xy, query = query_xy, k = min(as.integer(k), nrow(data_xy)))
}

normalize_cd4_label <- function(x) {
  x2 <- as.character(x)
  ifelse(grepl("Th17|TH17", x2, ignore.case = TRUE), "TH17",
  ifelse(grepl("Treg", x2, ignore.case = TRUE), "Treg",
  ifelse(grepl("TPH|TFH", x2, ignore.case = TRUE), "TPH_TFH",
  ifelse(grepl("na[iï]ve|naive", x2, ignore.case = TRUE), "CD4+ naïve",
  ifelse(grepl("CD4", x2, ignore.case = TRUE), "CD4+ effector", NA_character_)))))
}

# ------------------------------------------------------------------------------
# Load and harmonize main object
# ------------------------------------------------------------------------------
ste0 <- readRDS(INPUT_RDS)

keep_mask <- ste0@meta.data[["disease"]] %in% c("JIA", "STAT5b LOF")
ste <- subset(ste0, cells = rownames(ste0@meta.data)[keep_mask])

mix_sub1   <- grepl("MIX", ste@meta.data[["Annotation_sub1"]], ignore.case = TRUE)
mix_detail <- grepl("MIX", ste@meta.data[["Annotation_detail"]], ignore.case = TRUE)
ste <- subset(ste, cells = colnames(ste)[!(mix_sub1 | mix_detail)])

meta <- ste@meta.data
coords <- do.call(rbind, lapply(Images(ste), function(im) GetTissueCoordinates(ste, image = im)))
rownames(coords) <- coords$cell
meta$x <- coords[rownames(meta), "x"]
meta$y <- coords[rownames(meta), "y"]
meta$disease_group <- ifelse(grepl("STAT5b|stat5b", meta[["disease"]], ignore.case = TRUE), "STAT5b", "JIA")
meta$disease_plot  <- ifelse(meta$disease_group == "STAT5b", "Patient", "JIA")

meta$tls_bin <- ifelse(is.na(meta[["TLS_pos"]]), "TLS-", "TLS+")
meta$is_th17 <- meta[["Annotation_detail"]] == "Th17"

ste@meta.data <- meta

meta_dt <- as.data.table(meta)
meta_dt[, bc := rownames(meta)]
meta_dt_tls <- meta_dt[section_id %in% unique(meta_dt[tls_bin == "TLS+"]$section_id)]

# Representative JIA section and the STAT5b LOF patient, used in the spatial overlays.
patient_sample  <- unique(meta$sample_id[meta$disease == "STAT5b LOF"])
overlay_samples <- c(JIA_EXAMPLE_SAMPLE, patient_sample)

# ------------------------------------------------------------------------------
# Panel A: global and disease-specific UMAPs
# ------------------------------------------------------------------------------
umap_emb <- as.data.frame(Embeddings(ste, reduction = "umap"))
colnames(umap_emb) <- c("UMAP_1", "UMAP_2")
umap_df <- cbind(
  umap_emb,
  meta[, c("Annotation_detail", "Annotation_sub1", "disease", "disease_plot"), drop = FALSE],
  cell = rownames(meta)
)
cell_counts <- table(umap_df[["Annotation_detail"]])
umap_df$draw_order <- as.numeric(cell_counts[as.character(umap_df[["Annotation_detail"]])])
umap_df <- umap_df[sample(nrow(umap_df)), ]
umap_df <- umap_df[order(umap_df$draw_order), ]

base_umap_theme <- theme_void() +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.35, "cm"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

pA_total <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = .data[["Annotation_detail"]])) +
  geom_point(size = 0.03, alpha = 0.6) +
  ggtitle(sprintf("total (%s cells)", format(nrow(umap_df), big.mark = ","))) +
  base_umap_theme

make_disease_umap <- function(dis_name) {
  idx <- umap_df$disease_plot == dis_name
  ggplot() +
    geom_point(data = umap_df, aes(UMAP_1, UMAP_2), color = "grey88", size = 0.02, alpha = 0.5) +
    geom_point(data = umap_df[idx, , drop = FALSE], aes(UMAP_1, UMAP_2, color = .data[["Annotation_detail"]]), size = 0.04, alpha = 0.7) +
    ggtitle(sprintf("%s (%s cells)", dis_name, format(sum(idx), big.mark = ","))) +
    coord_cartesian(xlim = range(umap_df$UMAP_1), ylim = range(umap_df$UMAP_2)) +
    base_umap_theme
}
pA_jia <- make_disease_umap("JIA")
pA_patient <- make_disease_umap("Patient")
pA_combined <- (pA_total | pA_jia | pA_patient) + plot_layout(guides = "collect")

# ------------------------------------------------------------------------------
# Panel B / C / D: CD4-specific pipeline
# ------------------------------------------------------------------------------
cd4_mask <- grepl("CD4|Th17|Treg|TPH|TFH", meta[["Annotation_detail"]], ignore.case = TRUE) |
            grepl("CD4", meta[["Annotation_sub2"]], ignore.case = TRUE)
cd4_cells <- rownames(meta)[cd4_mask]
cd4_obj <- subset(ste, cells = cd4_cells)

cd4_obj$cd4_display <- normalize_cd4_label(cd4_obj[["Annotation_detail"]][, 1])
cd4_obj <- subset(cd4_obj, cells = colnames(cd4_obj)[!is.na(cd4_obj$cd4_display)])
cd4_obj$cd4_display <- factor(cd4_obj$cd4_display,
                              levels = c("CD4+ effector", "CD4+ naïve", "Treg", "TPH_TFH", "TH17"))

DefaultAssay(cd4_obj) <- "SCT"

cd4_pca_genes <- c(
  "RORC","CCR6","IL23R",
  "FOXP3","IL2RA","CTLA4","IKZF2",
  "PDCD1","CXCL13","CXCR5","BCL6",
  "CCR7","SELL","TCF7","LEF1",
  "CD3E","CD3D","CD4","CD2","CD5","CD28","ICOS",
  "IFNG","TNF","IL2","GZMB","GZMK","GZMA","PRF1",
  "BATF","IRF4","MAF","TBX21","EOMES",
  "CXCR3","CCR4","ITGA4","ITGAL",
  "TIGIT","ENTPD1","TNFRSF18","TNFRSF4",
  "ZAP70","LCK","FYN","STAT3","STAT5A","STAT5B",
  "TOX","TOX2","ZNF683","HOPX","PRDM1","LAG3","HAVCR2",
  "IL7R","CD44","CD69","KLRB1","KLRG1"
)
cd4_use_genes <- intersect(cd4_pca_genes, rownames(cd4_obj[["SCT"]]))

cd4_obj <- ScaleData(cd4_obj, assay = "SCT", features = cd4_use_genes, verbose = FALSE)
cd4_obj <- RunPCA(cd4_obj, assay = "SCT", features = cd4_use_genes, npcs = 20, verbose = FALSE)
cd4_obj <- FindNeighbors(cd4_obj, reduction = "pca", dims = 1:10, verbose = FALSE)
cd4_obj <- RunUMAP(cd4_obj, reduction = "pca", dims = 1:10, n.neighbors = 10, min.dist = 0.1, verbose = FALSE)

cd4_umap <- as.data.frame(Embeddings(cd4_obj, "umap"))
colnames(cd4_umap) <- c("UMAP_1", "UMAP_2")
cd4_umap$cd4_display <- cd4_obj$cd4_display
cd4_umap$disease_plot <- ifelse(grepl("STAT5b|stat5b", cd4_obj[["disease"]][, 1], ignore.case = TRUE), "Patient", "JIA")
rownames(cd4_umap) <- colnames(cd4_obj)

cd4_umap <- cd4_umap[sample(nrow(cd4_umap)), ]
df_effector <- cd4_umap %>% filter(cd4_display == "CD4+ effector")
df_other    <- cd4_umap %>% filter(cd4_display != "CD4+ effector")

pB <- ggplot() +
  ggrastr::geom_point_rast(data = df_effector, aes(UMAP_1, UMAP_2, color = cd4_display), size = 1.0, alpha = 0.35, raster.dpi = 300) +
  ggrastr::geom_point_rast(data = df_other,    aes(UMAP_1, UMAP_2, color = cd4_display), size = 1.4, alpha = 0.85, raster.dpi = 300) +
  theme_void(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggtitle("CD4")

make_cd4_split <- function(dis) {
  df_dis <- cd4_umap %>% filter(disease_plot == dis)
  df_eff <- df_dis %>% filter(cd4_display == "CD4+ effector")
  df_oth <- df_dis %>% filter(cd4_display != "CD4+ effector")
  ggplot() +
    ggrastr::geom_point_rast(data = cd4_umap, aes(UMAP_1, UMAP_2), color = "#E0E0E0", size = 0.5, alpha = 0.3, raster.dpi = 300) +
    ggrastr::geom_point_rast(data = df_eff, aes(UMAP_1, UMAP_2, color = cd4_display), size = 1.0, alpha = 0.35, raster.dpi = 300) +
    ggrastr::geom_point_rast(data = df_oth, aes(UMAP_1, UMAP_2, color = cd4_display), size = 1.4, alpha = 0.85, raster.dpi = 300) +
    coord_cartesian(xlim = range(cd4_umap$UMAP_1), ylim = range(cd4_umap$UMAP_2)) +
    theme_void(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
    ggtitle(dis)
}
pB_jia     <- make_cd4_split("JIA")
pB_patient <- make_cd4_split("Patient")

# Panel C: stacked bar
cd4_meta <- cd4_obj@meta.data
cd4_meta$cd4_display <- cd4_obj$cd4_display
cd4_meta$disease_plot <- ifelse(grepl("STAT5b|stat5b", cd4_meta[["disease"]], ignore.case = TRUE), "Patient", "JIA")
panelC_df <- cd4_meta %>%
  dplyr::count(disease_plot, cd4_display, name = "n") %>%
  dplyr::group_by(disease_plot) %>%
  dplyr::mutate(pct = n / sum(n) * 100) %>%
  dplyr::ungroup()
panelC_df$cd4_display  <- factor(panelC_df$cd4_display,  levels = rev(levels(cd4_obj$cd4_display)))
panelC_df$disease_plot <- factor(panelC_df$disease_plot, levels = c("JIA", "Patient"))

pC <- ggplot(panelC_df, aes(disease_plot, pct, fill = cd4_display)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_y_continuous(labels = percent_format(scale = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "", y = "Proportion (%)", title = "CD4 subtype") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 12, face = "bold"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Panel D: CD4 FeaturePlots with fixed per-gene expression ranges
panelD_params <- list(
  CCR7  = list(breaks = c(0, 2,   5), lim = c(0, 5)),
  SELL  = list(breaks = c(0, 3.5, 5), lim = c(0, 5)),
  TCF7  = list(breaks = c(0, 3.5, 5), lim = c(0, 5)),
  CCR6  = list(breaks = c(0, 3,   4), lim = c(0, 4)),
  RORC  = list(breaks = c(0, 3,   4), lim = c(0, 4)),
  IL23R = list(breaks = c(0, 3,   4), lim = c(0, 4)),
  PRDM1 = list(breaks = c(0, 2,   5), lim = c(0, 5)),
  LAG3  = list(breaks = c(0, 2.5, 5), lim = c(0, 5)),
  CXCL13= list(breaks = c(0, 3.5, 5), lim = c(0, 5)),
  PDCD1 = list(breaks = c(0, 2.5, 5), lim = c(0, 5)),
  FOXP3 = list(breaks = c(0, 3,   5), lim = c(0, 5)),
  IL2RA = list(breaks = c(0, 3.5, 5), lim = c(0, 5))
)

make_fp_d <- function(gene, params) {
  p <- params[[gene]]
  br  <- p$breaks
  lim <- p$lim

  plot_df <- cd4_umap
  plot_df$expr <- as.numeric(GetAssayData(cd4_obj, assay = "SCT", layer = "data")[gene, rownames(plot_df)])
  plot_df$expr_plot <- pmin(pmax(plot_df$expr, lim[1]), lim[2])
  plot_df <- plot_df[order(plot_df$expr_plot), , drop = FALSE]

  ggplot(plot_df, aes(UMAP_1, UMAP_2, color = expr_plot)) +
    ggrastr::geom_point_rast(size = 1.0, alpha = 0.90, raster.dpi = 300) +
    scale_color_gradientn(
      colors = c("grey90", "grey90", "#08519C"),
      values = scales::rescale(br, from = lim),
      limits = lim,
      breaks = lim,
      labels = as.character(lim),
      oob    = scales::squish,
      name   = "Expression"
    ) +
    coord_cartesian(xlim = range(cd4_umap$UMAP_1), ylim = range(cd4_umap$UMAP_2)) +
    ggtitle(gene) +
    theme_void(base_size = 9) +
    theme(
      plot.title        = element_text(size = 11, face = "bold.italic", hjust = 0.5),
      legend.position   = "right",
      legend.key.width  = unit(0.25, "cm"),
      legend.key.height = unit(0.6,  "cm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 7)
    )
}

panelD_genes <- intersect(names(panelD_params), rownames(cd4_obj[["SCT"]]))
fp_list <- lapply(panelD_genes, make_fp_d, params = panelD_params)
pD <- wrap_plots(fp_list, ncol = 6)

# ------------------------------------------------------------------------------
# Panel E / K / L: common preprocessing (Th17 / SL environment)
# ------------------------------------------------------------------------------
IL17_CANDIDATES <- c("CCL20","CXCL2","IL6","TRAF3IP2","ZC3H12A","IL17RA")
TH17_LABEL_HM  <- "Th17"
SL_SUB2_LABEL  <- "SL"

hm_md <- as.data.table(ste@meta.data)
hm_md[, cell_id := rownames(ste@meta.data)]

prev_hm <- fread(CELL_META_CSV, select = c("cell_id", "th17_tls_label"))
hm_md   <- merge(hm_md, prev_hm, by = "cell_id", all.x = TRUE)

hm_md[, disease_label := fcase(
  disease == "STAT5b LOF", "STAT5b_LOF",
  disease == "JIA",        "JIA",
  default = NA_character_
)]
hm_md <- hm_md[!is.na(disease_label)]

n_jia_s    <- length(unique(hm_md[disease_label == "JIA"]$sample_id))
il17_genes <- intersect(IL17_CANDIDATES, rownames(ste[["SCT"]]))

th17_md <- hm_md[Annotation_detail == TH17_LABEL_HM & !is.na(th17_tls_label)]
sl_md   <- hm_md[Annotation_sub2 == SL_SUB2_LABEL]

th17_md[, win_label_coarse := fcase(
  th17_tls_label == "TLS-associated", "TLS-assoc",
  default = "non-TLS-assoc"
)]

# Mean log-normalized IL-17 pathway score per SL fibroblast
sl_norm_hm <- GetAssayData(ste, assay = "SCT", layer = "data")[il17_genes, sl_md$cell_id, drop = FALSE]
sl_norm_dt_hm <- as.data.table(t(as.matrix(sl_norm_hm)))
sl_norm_dt_hm[, cell_id := sl_md$cell_id]
sl_norm_dt_hm[, IL17_norm_score := rowMeans(.SD, na.rm = TRUE), .SDcols = il17_genes]
sl_md <- merge(sl_md, sl_norm_dt_hm[, .(cell_id, IL17_norm_score)], by = "cell_id", all.x = TRUE)

# ------------------------------------------------------------------------------
# Panel E: Th17 pseudobulk volcano (STAT5b LOF vs JIA)
# ------------------------------------------------------------------------------
VOLCANO_GENES_HM <- c("AHR","CCR6","IL26","IL23R","BATF","RORC","TYK2","TBX21","BACH2")

th17_pb_md_hm <- hm_md[Annotation_detail == TH17_LABEL_HM]

all_expr_hm <- GetAssayData(ste, assay = "Spatial", layer = "counts")[, th17_pb_md_hm$cell_id, drop = FALSE]
samples_pb_hm <- unique(th17_pb_md_hm$sample_id)
pb_mat_hm <- sapply(samples_pb_hm, function(s) {
  cells <- th17_pb_md_hm[sample_id == s, cell_id]
  if (length(cells) == 1) as.numeric(all_expr_hm[, cells])
  else rowSums(all_expr_hm[, cells, drop = FALSE])
})
rownames(pb_mat_hm) <- rownames(all_expr_hm)

samp_dis_hm  <- th17_pb_md_hm[, .(disease_label = disease_label[1]), by = "sample_id"]
jia_s_hm     <- samp_dis_hm[disease_label == "JIA"       ][["sample_id"]]
stat5b_s_hm  <- samp_dis_hm[disease_label == "STAT5b_LOF"][["sample_id"]]
lib_size_hm  <- colSums(pb_mat_hm)
cpm_mat_hm   <- sweep(pb_mat_hm, 2, lib_size_hm, "/") * 1e6
jia_cpm_hm   <- cpm_mat_hm[, intersect(jia_s_hm,   colnames(cpm_mat_hm)), drop = FALSE]
stat5_cpm_hm <- cpm_mat_hm[, intersect(stat5b_s_hm, colnames(cpm_mat_hm)), drop = FALSE]

eps_hm        <- 0.1
jia_mean_hm   <- rowMeans(jia_cpm_hm,  na.rm = TRUE)
jia_sd_hm     <- apply(jia_cpm_hm, 1, sd, na.rm = TRUE)
stat5_mean_hm <- rowMeans(stat5_cpm_hm, na.rm = TRUE)
log2fc_hm     <- log2((stat5_mean_hm + eps_hm) / (jia_mean_hm + eps_hm))
zscore_hm     <- ifelse(jia_sd_hm > 0, (stat5_mean_hm - jia_mean_hm) / jia_sd_hm, NA_real_)

pb_dt_hm <- data.table(gene   = rownames(cpm_mat_hm),
                        log2FC = round(log2fc_hm, 4),
                        zscore = round(zscore_hm, 4))
pb_dt_hm <- pb_dt_hm[order(-abs(zscore))]

pb_vol_hm <- pb_dt_hm[!is.na(zscore) & !is.na(log2FC) & is.finite(zscore) & is.finite(log2FC)]
pb_vol_hm[, abs_z    := abs(zscore)]
pb_vol_hm[, is_label := toupper(gene) %in% toupper(VOLCANO_GENES_HM)]
pb_vol_hm[, label    := ifelse(is_label, gene, "")]
pb_vol_hm[, color_group := fcase(
  is_label & log2FC > 0, "STAT5b high",
  is_label & log2FC < 0, "JIA high",
  default = "neutral"
)]

bg_dt_hm    <- pb_vol_hm[color_group == "neutral"]
label_dt_hm <- pb_vol_hm[is_label == TRUE]

pE <- ggplot() +
  geom_point(data = bg_dt_hm,
             aes(x = log2FC, y = abs_z),
             color = "grey80", size = 0.8, alpha = 0.5) +
  geom_point(data = label_dt_hm,
             aes(x = log2FC, y = abs_z, color = color_group),
             size = 2.8, alpha = 0.95) +
  ggrepel::geom_text_repel(data = label_dt_hm,
                aes(x = log2FC, y = abs_z, label = label, color = color_group),
                size = 3.2, fontface = "italic",
                max.overlaps = Inf, box.padding = 0.5, force = 18,
                segment.size = 0.3, segment.color = "grey50", show.legend = FALSE) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(
    title    = "Th17 pseudobulk: STAT5b LOF vs JIA",
    subtitle = sprintf(
      "X: log2FC (STAT5b/JIA)  |  Y: |Z-score|  |  JIA: %d samples, STAT5b LOF: 1  |  CPM-normalized",
      length(jia_s_hm)),
    x = "log2 Fold Change (STAT5b LOF / JIA)", y = "|Z-score|"
  )

# ------------------------------------------------------------------------------
# Panel F / G: Th17 distance-band neighbor analysis
# ------------------------------------------------------------------------------
panel_G_genes <- c("CXCR4", "CCL19", "CCR6", "CCR7", "CXCL13", "CXCL12", "CXCR5", "CCL20")
expr_genes_for_G <- intersect(panel_G_genes, rownames(ste[["SCT"]]))
expr_mat_G <- GetAssayData(ste, assay = "SCT", layer = "data")[expr_genes_for_G, , drop = FALSE]
for (g in expr_genes_for_G) {
  meta_dt[[paste0("expr_", g)]] <- as.numeric(expr_mat_G[g, meta_dt$bc])
}

neighbor_bands <- c(20, 40, 60, 80, 100, 150, 200)
immune_pattern <- "CD4|CD8|TRM|Th17|Treg|TPH_TFH|NK cells|gdT|Bcell|B cell|conventional B|plasmablast|Plasma cell"

compute_neighbor_band_summary <- function(dt_in, label_col, with_expr = FALSE) {
  res <- list()
  for (sec_id in unique(dt_in$section_id)) {
    ds <- copy(dt_in[section_id == sec_id])
    idx_th17 <- which(ds$is_th17 == TRUE)
    if (length(idx_th17) < 2) next

    if (with_expr) {
      idx_other <- which(!ds$is_th17 & grepl(immune_pattern, ds[["Annotation_detail"]], ignore.case = TRUE))
    } else {
      idx_other <- which(!ds$is_th17)
    }
    if (length(idx_other) < 10) next

    nn_res <- run_nn2(ds[idx_other, .(x, y)], ds[idx_th17, .(x, y)], k = min(200, length(idx_other)))

    for (b in seq_along(neighbor_bands)) {
      r_max <- neighbor_bands[b]
      r_min <- if (b == 1) 0 else neighbor_bands[b - 1]
      for (j in seq_along(idx_th17)) {
        dists <- nn_res$nn.dists[j, ]
        idxs  <- nn_res$nn.idx[j, ]
        in_band <- which(dists > r_min & dists <= r_max)
        if (length(in_band) == 0) next
        use_rows <- idx_other[idxs[in_band]]
        if (with_expr) {
          for (g in expr_genes_for_G) {
            vals <- ds[[paste0("expr_", g)]][use_rows]
            if (all(is.na(vals))) next
            res[[length(res) + 1]] <- data.table(
              section_id = sec_id, band = paste0(r_min, "-", r_max),
              r_min = r_min, r_max = r_max, label = g,
              pos_frac = mean(vals > 0, na.rm = TRUE), n_neighbors = length(use_rows)
            )
          }
        } else {
          lbl_use <- as.character(ds[[label_col]][use_rows])
          ct_tab <- table(lbl_use)
          for (ct in names(ct_tab)) {
            res[[length(res) + 1]] <- data.table(
              section_id = sec_id, band = paste0(r_min, "-", r_max),
              r_min = r_min, r_max = r_max, label = ct, count = as.integer(ct_tab[[ct]])
            )
          }
        }
      }
    }
  }
  if (length(res) == 0) return(NULL)
  rbindlist(res)
}

F_matched <- c("Bcell", "Plasma cell", "plasmablast", "CD4 memory", "TPH_TFH")

F_label_map <- c(
  "Bcell"       = "B cell",
  "Plasma cell" = "Plasma cell",
  "plasmablast" = "Plasmablast",
  "CD4 memory"  = "CD4+ effector",
  "TPH_TFH"     = "TPH/TFH"
)
F_display_levels <- c("B cell", "Plasma cell", "Plasmablast", "CD4+ effector", "TPH/TFH")

# Neighbor denominator = all non-Th17 cells in each band;
# selected labels are extracted only after the band denominator is computed.
F_raw_all <- compute_neighbor_band_summary(meta_dt, label_col = "Annotation_detail", with_expr = FALSE)

F_all_sum <- F_raw_all[, .(total = sum(count)), by = .(band, r_min, r_max, label)]
F_band_total <- F_all_sum[, .(band_total = sum(total)), by = .(band, r_min, r_max)]
F_sum <- merge(F_all_sum, F_band_total, by = c("band", "r_min", "r_max"), all.x = TRUE)
F_sum[, observed_pct := ifelse(band_total > 0, total / band_total * 100, NA_real_)]
F_sum$band <- factor(F_sum$band, levels = unique(F_sum[order(r_min)]$band))

F_plot_sum <- F_sum[label %in% F_matched]

# Baseline = 0-20um band
F_base <- F_plot_sum[r_min == 0, .(label, baseline = observed_pct)]
F_plot_sum <- merge(F_plot_sum, F_base, by = "label", all.x = TRUE)
F_plot_sum <- F_plot_sum[!is.na(baseline)]

F_plot_sum[, rel_change_pct := ifelse(
  is.na(baseline) | baseline <= 0, NA_real_, (observed_pct / baseline - 1) * 100)]

F_plot_sum[, label_plot := factor(F_label_map[as.character(label)], levels = F_display_levels)]

plot_F_data <- F_plot_sum[!is.na(rel_change_pct) & !is.na(label_plot)]
pF <- ggplot(plot_F_data, aes(band, rel_change_pct, color = label_plot, group = label_plot)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(title = "Selected TLS-related cell types by distance from TH17",
       x = "Distance band (µm)", y = "Relative change vs 0-20µm (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

G_raw <- compute_neighbor_band_summary(meta_dt, label_col = "Annotation_detail", with_expr = TRUE)
G_sum <- G_raw[, .(mean_pos = mean(pos_frac, na.rm = TRUE)), by = .(band, r_min, r_max, label)]
G_base <- G_sum[r_min == 0, .(label, baseline = mean_pos)]
G_sum <- merge(G_sum, G_base, by = "label", all.x = TRUE)
G_sum[, rel_change_pct := ifelse(is.na(baseline) | baseline <= 0, NA_real_, (mean_pos / baseline - 1) * 100)]
G_sum$band  <- factor(G_sum$band, levels = unique(G_sum[order(r_min)]$band))
G_sum$label <- factor(G_sum$label, levels = expr_genes_for_G)
pG <- ggplot(G_sum[!is.na(rel_change_pct)], aes(band, rel_change_pct, color = label, group = label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(title = "Selected chemokines by distance from TH17",
       x = "Distance band (µm)", y = "Relative change vs 0-20µm (%)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Panel H: TH17/TLS overlay (representative JIA and STAT5b LOF patient)
# ------------------------------------------------------------------------------
MAX_OVERLAY_CELLS <- 50000L

overlay_list <- list()
for (sid in overlay_samples) {
  ds <- meta_dt[sample_id == sid]
  if (nrow(ds) > MAX_OVERLAY_CELLS) ds <- ds[sample(.N, MAX_OVERLAY_CELLS)]

  p <- ggplot(ds, aes(x, y)) +
    geom_point(data = ds[tls_bin == "TLS-" & !is_th17], color = "grey78", size = 0.01, alpha = 0.35) +
    geom_point(data = ds[tls_bin == "TLS+"],             color = "#2166AC", size = 0.04, alpha = 0.40) +
    geom_point(data = ds[is_th17 == TRUE],               color = "#E41A1C", size = 0.70, alpha = 0.90) +
    coord_equal() +
    theme_void(base_size = 6) +
    labs(
      title    = paste0(sid, " (", unique(ds$disease_plot), ")"),
      subtitle = paste0("Red = TH17 (", sum(ds$is_th17), ") | Blue = TLS (", sum(ds$tls_bin == "TLS+"), ")")
    ) +
    theme(plot.title = element_text(face = "bold", size = 8), plot.subtitle = element_text(size = 6))

  overlay_list[[sid]] <- p
}
pH <- wrap_plots(overlay_list, ncol = 2) +
  plot_annotation(title = "TH17/TLS overlay (Red = TH17, Blue = TLS)")

# ------------------------------------------------------------------------------
# Panels I / J: TH17 distance-to-TLS analysis
# ------------------------------------------------------------------------------
DIST_CAP       <- 1500
DIST_PLOT_YMAX <- 1500

th17_edge_rows <- list()
for (sec_id in unique(meta_dt_tls$section_id)) {
  ds <- meta_dt_tls[section_id == sec_id]
  tls_sec  <- ds[tls_bin == "TLS+"]
  th17_out <- ds[is_th17 == TRUE & tls_bin == "TLS-"]
  if (nrow(tls_sec) == 0 || nrow(th17_out) == 0) next
  nn_res <- run_nn2(tls_sec[, .(x, y)], th17_out[, .(x, y)], k = 1)
  th17_edge_rows[[length(th17_edge_rows) + 1]] <- data.table(
    bc = th17_out$bc, section_id = sec_id,
    sample_id = th17_out[["sample_id"]],
    disease = th17_out$disease_group,
    disease_plot = th17_out$disease_plot,
    edge_dist = as.numeric(nn_res$nn.dists[, 1])
  )
}
th17_edge_dt <- rbindlist(th17_edge_rows, fill = TRUE)
th17_edge_dt <- th17_edge_dt[!is.na(edge_dist) & edge_dist <= DIST_CAP]

ecdf_sample <- th17_edge_dt[, .(
  n_th17_outside = .N,
  mean_dist   = mean(edge_dist, na.rm = TRUE),
  median_dist = median(edge_dist, na.rm = TRUE)
), by = .(sample_id, disease, disease_plot)]

box_df <- copy(th17_edge_dt)
sample_order <- box_df[, .(med = median(edge_dist, na.rm = TRUE)), by = sample_id][order(med)]$sample_id
box_df$sample_id <- factor(box_df$sample_id, levels = sample_order)

pI_left <- ggplot(box_df, aes(sample_id, edge_dist, fill = disease_plot)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.75) +
  coord_cartesian(ylim = c(0, DIST_PLOT_YMAX)) +
  labs(title = "Distance of TLS-external TH17 to nearest TLS boundary", x = "", y = "Distance to TLS boundary (µm)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")

pI_right <- ggplot(ecdf_sample, aes(disease_plot, median_dist, color = disease_plot, shape = disease_plot)) +
  geom_jitter(aes(size = n_th17_outside), width = 0.12, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "black", linewidth = 0.6) +
  scale_shape_manual(values = c("JIA" = 16, "Patient" = 17)) +
  scale_size_continuous(range = c(2, 6), name = "Outside TH17") +
  coord_cartesian(ylim = c(0, DIST_PLOT_YMAX)) +
  labs(title = "Median distance", x = "", y = "Median distance (µm)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

pI <- pI_left + pI_right + plot_layout(widths = c(3.3, 1.0))

x_cap <- DIST_CAP
pJ <- ggplot(th17_edge_dt[edge_dist <= x_cap], aes(x = edge_dist, color = disease_plot, group = sample_id)) +
  stat_ecdf(data = th17_edge_dt[disease_plot == "JIA"     & edge_dist <= x_cap], aes(group = sample_id), alpha = 0.35, linewidth = 0.7) +
  stat_ecdf(data = th17_edge_dt[disease_plot == "Patient" & edge_dist <= x_cap], aes(group = sample_id), linewidth = 1.4) +
  coord_cartesian(xlim = c(0, x_cap), ylim = c(0, 1)) +
  labs(title = "Cumulative fraction of TLS-external TH17",
       x = "Distance to nearest TLS boundary (µm)", y = "Cumulative fraction of TLS-external TH17") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# Panel K: IL-17 hotspot overlay (representative JIA and STAT5b LOF patient)
# ------------------------------------------------------------------------------
IL17_GENES_OVERLAY <- c("CCL20","CXCL2","IL6","TRAF3IP2","ZC3H12A","IL17RA")

norm_mat_ov   <- GetAssayData(ste, assay = "SCT", layer = "data")
use_genes_ov  <- intersect(IL17_GENES_OVERLAY, rownames(norm_mat_ov))
ste$IL17_score_ov <- as.numeric(Matrix::colMeans(norm_mat_ov[use_genes_ov, , drop = FALSE]))
meta_ov <- ste@meta.data

# TLS outlines as the convex hull of the TLS+ cells of each TLS
make_tls_layer_ov <- function(df_tls) {
  ids <- unique(df_tls$TLS_id)
  ids <- ids[!is.na(ids)]
  hull_df <- do.call(rbind, lapply(ids, function(tid) {
    sub <- df_tls[df_tls$TLS_id == tid, ]
    idx <- chull(sub$x, sub$y)
    data.frame(x = sub$x[c(idx, idx[1])], y = sub$y[c(idx, idx[1])],
               tls_id = as.character(tid))
  }))
  geom_polygon(data = hull_df, aes(x = x, y = y, group = tls_id),
               fill = "#1565C0", color = "#0D47A1", linewidth = 0.5, alpha = 0.4)
}

pK_list <- list()

for (smp in overlay_samples) {
  cells_smp_ov <- rownames(meta_ov)[meta_ov[["sample_id"]] == smp]

  coords_ov <- GetTissueCoordinates(ste, image = smp)
  rownames(coords_ov) <- coords_ov$cell
  cells_common_ov <- intersect(coords_ov$cell, cells_smp_ov)
  coords_sub_ov   <- coords_ov[cells_common_ov, ]

  meta_img_ov <- meta_ov[cells_common_ov,
                         c("Annotation_sub2", "Annotation_detail", "TLS_pos", "TLS_id"),
                         drop = FALSE]
  meta_img_ov$IL17_score <- ste$IL17_score_ov[cells_common_ov]
  meta_img_ov$x          <- coords_sub_ov[["x"]]
  meta_img_ov$y          <- coords_sub_ov[["y"]]

  df_all_ov  <- meta_img_ov
  df_sl_ov   <- meta_img_ov[meta_img_ov$Annotation_sub2 == "SL", ]
  df_th17_ov <- meta_img_ov[meta_img_ov$Annotation_detail == "Th17", ]
  df_tls_ov  <- meta_img_ov[!is.na(meta_img_ov$TLS_pos), ]

  df_sl_pos_ov <- df_sl_ov[df_sl_ov$IL17_score > 0, ]

  df_sl_pos_ov$z <- Nebulosa:::calculate_density(
    df_sl_pos_ov$IL17_score,
    df_sl_pos_ov[, c("x","y")],
    method = "wkde"
  )
  p_hot_ov <- ggplot() +
    geom_point(data = df_all_ov, aes(x = x, y = y),
               color = "grey88", size = 0.12, alpha = 0.35) +
    geom_point(data = df_sl_pos_ov, aes(x = x, y = y, color = z),
               size = 0.7, alpha = 0.85) +
    scale_color_gradientn(
      colours = c("grey90","#C7E9C0","#74C476","#238B45","#00441B"),
      values  = scales::rescale(c(
        quantile(df_sl_pos_ov$z, 0.25, na.rm = TRUE),
        quantile(df_sl_pos_ov$z, 0.55, na.rm = TRUE),
        quantile(df_sl_pos_ov$z, 0.75, na.rm = TRUE),
        quantile(df_sl_pos_ov$z, 0.90, na.rm = TRUE),
        quantile(df_sl_pos_ov$z, 0.99, na.rm = TRUE)
      )),
      name  = "IL-17\nwKDE",
      guide = guide_colorbar(barwidth = 0.8, barheight = 6,
                             title.position = "top")
    ) +
    geom_point(data = df_th17_ov, aes(x = x, y = y),
               color = "#D62728", size = 0.7, alpha = 0.85) +
    make_tls_layer_ov(df_tls_ov) +
    labs(title    = "IL-17 Hotspot (wKDE) - SL Fibroblast",
         subtitle = paste0("Sample: ", smp),
         x = "X (µm)", y = "Y (µm)",
         caption  = paste0("Genes: ", paste(use_genes_ov, collapse = ", "))) +
    coord_equal() +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          plot.caption  = element_text(size = 7, color = "grey55"),
          panel.grid    = element_blank(),
          legend.position = "right")

  pK_list[[smp]] <- p_hot_ov
}

# ------------------------------------------------------------------------------
# Panel L: SL fibroblast IL-17 score vs distance to nearest Th17 (LOESS)
# ------------------------------------------------------------------------------
DIS2_BIN_WS <- c(20, 10)
DIS2_SPAN   <- 0.6
DIS2_DISP   <- 200
DIS2_MIN_N  <- 5L
DIS2_ORDER  <- c("JIA", "STAT5b_LOF")

build_dis2_master_v3 <- function(th17_subset, label) {
  pair_list_d2 <- list()
  for (sec in unique(th17_subset$section_id)) {
    th17_sec_d2 <- th17_subset[section_id == sec]
    sl_sec_d2   <- sl_md[section_id == sec]
    if (nrow(th17_sec_d2) == 0 || nrow(sl_sec_d2) == 0) next
    th17_xy_d2 <- as.matrix(th17_sec_d2[, .(x, y)])
    sl_xy_d2   <- as.matrix(sl_sec_d2[,   .(x, y)])
    k_max_d2   <- min(nrow(th17_sec_d2), 500L)
    nn_res_d2  <- RANN::nn2(th17_xy_d2, sl_xy_d2, k = k_max_d2,
                             searchtype = "radius", radius = DIS2_DISP)
    for (i in seq_len(nrow(sl_sec_d2))) {
      hits  <- nn_res_d2$nn.idx[i, ]
      dsts  <- nn_res_d2$nn.dists[i, ]
      valid <- hits > 0 & dsts <= DIS2_DISP
      if (!any(valid)) next
      pair_list_d2[[paste(sec, i, sep = "_")]] <- data.table(
        sl_cell_id           = sl_sec_d2$cell_id[i],
        dist_to_nearest_th17 = min(dsts[valid]),
        section_id           = sec,
        sample_id            = sl_sec_d2$sample_id[i],
        disease_label        = sl_sec_d2$disease_label[i]
      )
    }
  }
  pairs_d2 <- rbindlist(pair_list_d2)
  pairs_d2 <- merge(pairs_d2, sl_md[, .(cell_id, IL17_norm_score)],
                    by.x = "sl_cell_id", by.y = "cell_id", all.x = TRUE)
  pairs_d2[, analysis_label := label]
  pairs_d2
}

dis2_bin_agg_v3 <- function(dt, bw, min_n = DIS2_MIN_N) {
  dt2 <- copy(dt)
  dt2[, bin_lo     := floor(dist_to_nearest_th17 / bw) * bw]
  dt2[, bin_center := bin_lo + bw / 2]
  res <- dt2[!is.na(IL17_norm_score) & is.finite(IL17_norm_score), {
    x  <- IL17_norm_score
    n  <- length(x)
    mn <- mean(x)
    se <- if (n > 1) sd(x) / sqrt(n) else 0
    .(bin_center = bin_center[1], wt_mean = mn, wt_se = se, n_cells = n)
  }, by = .(disease_label, sample_id, bin_lo)]
  res[n_cells >= min_n]
}

dis2_grp_agg_v3 <- function(samp_dt) {
  samp_dt[, .(
    grp_mean  = mean(wt_mean, na.rm = TRUE),
    grp_se    = sd(wt_mean,   na.rm = TRUE) / sqrt(sum(!is.na(wt_mean))),
    n_samples = sum(!is.na(wt_mean))
  ), by = .(disease_label, bin_lo, bin_center)]
}

make_dis2_loess_v3 <- function(grp_dt, bw, span, y_max, title, subtitle) {
  dt <- grp_dt[is.finite(grp_mean)]
  dt[, disease_label := factor(disease_label, levels = DIS2_ORDER)]
  ggplot(dt, aes(x = bin_center, y = grp_mean,
                 color = disease_label, fill = disease_label,
                 group = disease_label)) +
    geom_smooth(method = "loess", formula = y ~ x,
                span = span, se = TRUE, level = 0.95,
                linewidth = 1.1, alpha = 0.15) +
    scale_x_continuous(breaks = seq(0, DIS2_DISP, max(bw * 2, 40)),
                       limits = c(0, DIS2_DISP), expand = c(0.01, 0)) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme_bw(base_size = 12) +
    theme(legend.position  = "bottom",
          legend.key.width = unit(1.4, "cm"),
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill = "white", color = NA)) +
    labs(title = title, subtitle = subtitle,
         x = "Distance to nearest Th17 (um)",
         y = "Mean log-normalized IL-17 score")
}

dis2_all_dt <- build_dis2_master_v3(th17_md, "All-Th17")
dis2_non_dt <- build_dis2_master_v3(
  th17_md[win_label_coarse == "non-TLS-assoc"], "non-TLS-assoc-Th17-only")

y_vals_d2 <- c()
for (dt_obj in list(dis2_all_dt, dis2_non_dt)) {
  for (bw in DIS2_BIN_WS) {
    g <- dis2_grp_agg_v3(dis2_bin_agg_v3(dt_obj, bw))
    y_vals_d2 <- c(y_vals_d2, g$grp_mean + 1.96 * g$grp_se)
  }
}
dis2_y_max <- max(y_vals_d2, na.rm = TRUE) * 1.05

grp_m <- dis2_grp_agg_v3(dis2_bin_agg_v3(dis2_all_dt, 10))
pM <- make_dis2_loess_v3(
  grp_m, 10, DIS2_SPAN, dis2_y_max,
  title    = "SL fibroblast IL-17 score: JIA vs STAT5b  [All Th17 pooled]",
  subtitle = sprintf("LOESS span=%.2f | bin=10um | 1 SL = 1 cell | JIA n=%d | STAT5b n=1",
                     DIS2_SPAN, n_jia_s)
)

grp_n <- dis2_grp_agg_v3(dis2_bin_agg_v3(dis2_non_dt, 10))
pN <- make_dis2_loess_v3(
  grp_n, 10, DIS2_SPAN, dis2_y_max,
  title    = "SL fibroblast IL-17 score: JIA vs STAT5b  [non-TLS-assoc Th17 only]",
  subtitle = sprintf("LOESS span=%.2f | bin=10um | 1 SL = 1 cell | non-TLS-assoc Th17 | JIA n=%d | STAT5b n=1",
                     DIS2_SPAN, n_jia_s)
)
