#!/usr/bin/env Rscript

library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(grid)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(scCustomize)
library(magrittr)
library(dplyr)
library(stringr)
library(hdf5r)
library(harmony)
library(ggrastr)
library(Nebulosa)
library(scDblFinder)
library(MOFA2)
library(MOFAdata)
library(biomaRt)
library(StabMap)
library(ggalluvial)
library(pheatmap)
library(Matrix)
library(ggraph)
library(circlize)
library(scales)
library(CellChat)
library(ComplexHeatmap)
library(parallel)
library(RANN)
library(gridExtra)
library(scp)
library(randomcoloR)
library(fgsea)
library(msigdbr)
library(glmGamPoi)
library(data.table)
library(lisi)
library(SoupX)
library(stats)

# Load the Xenium data
ReadXenium <- function (data.dir, outs = c("matrix", "microns"), type = "centroids",
                        mols.qv.threshold = 20)
{
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                    several.ok = TRUE)
  outs <- c(outs, type)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir,
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      cell_info <- as.data.frame(data.table::fread(file.path(data.dir,
                                                             "cells.csv.gz")))
      cell_centroid_df <- data.frame(x = cell_info$x_centroid,
                                     y = cell_info$y_centroid, cell = cell_info$cell_id,
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir,
                                                                      "cell_boundaries.csv.gz")))
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location,
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0) {
  my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
  nn.ranked <- my.knn$nn.idx

  snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
  rownames(snn_res) <- row.names(data.use)
  colnames(snn_res) <- row.names(data.use)
  return(snn_res)
}
environment(BuildSNNSeurat) <- asNamespace("Seurat")

# Moran's I for a row-normalized spatial weight matrix W and a vector x
compute_morans_i <- function(W, x) {
  N <- length(x)
  x_centered <- x - mean(x)
  num <- as.numeric(t(x_centered) %*% W %*% x_centered)
  den <- sum(x_centered^2)
  I <- (N / sum(W)) * (num / den)
  return(I)
}

compute_gearys_c <- function(W, x) {
  N <- length(x)
  x_mean <- mean(x)
  var_x <- sum((x - x_mean)^2) / (N - 1)
  num <- 0.5 * sum(W * outer(x, x, FUN = function(a, b) (a - b)^2))
  C <- (N - 1) * num / sum(W) / (N * var_x)
  return(C)
}

# Permutation p-value: one-sided for Moran's I, deviation from 1 for Geary's C
permutation_test <- function(obs_stat, W, x, mode = "moran", n_perms = 999, seed = 42) {
  if (!is.null(seed)) set.seed(seed)

  perm_stats <- replicate(n_perms, {
    x_perm <- sample(x)
    if (mode == "moran") {
      compute_morans_i(W, x_perm)
    } else {
      compute_gearys_c(W, x_perm)
    }
  })

  if (mode == "moran") {
    pval <- mean(perm_stats >= obs_stat)
  } else {
    pval <- mean(abs(perm_stats - 1) >= abs(obs_stat - 1))
  }

  return(list(
    perm_stats = perm_stats,
    pval_perm = pval
  ))
}

# Global spatial autocorrelation (Moran's I or Geary's C) per gene.
# Spatial neighbors are built per image and combined block-diagonally.
spatial_autocorr <- function(
    seurat_obj,
    neighbors.k = 30,
    connectivity_key = "nn",
    genes = NULL,
    mode = c("moran", "geary"),
    transformation = TRUE,
    n_perms = 500,
    corr_method = "BH",
    assay = "SCT",
    attr = "data",
    seed = 1938493,
    copy = FALSE
) {

  all_nn <- list()
  all_snn <- list()
  cell_id <- vector()
  for (name in names(seurat_obj@images)) {
    coords <- seurat_obj[[name]]$centroids@coords %>%
      as.data.frame() %>%
      dplyr::mutate(cell = Cells(seurat_obj[[name]]))
    cells <- coords$cell
    rownames(coords) <- cells
    coords <- as.matrix(coords[, c("x", "y")])

    neighbors <- FindNeighbors(coords, k.param = neighbors.k, verbose = FALSE)

    all_nn[[name]] <- neighbors$nn
    all_snn[[name]] <- neighbors$snn
    cell_id <- c(cell_id, cells)
  }

  if (connectivity_key == "nn") {
    W <- bdiag(all_nn)
  } else {
    W <- bdiag(all_snn)
  }
  rownames(W) <- colnames(W) <- cell_id

  # Row-normalize the weights by the degree of each cell
  if (transformation) {
    degrees <- Matrix::colSums(W) + 1
    W <- W / degrees
  }

  mat <- Seurat::GetAssayData(seurat_obj, layer = attr, assay = assay)
  all_genes <- rownames(mat)

  if (is.null(genes)) {
    hvf <- VariableFeatures(seurat_obj)
    genes <- if (length(hvf) > 0) hvf else all_genes
  }
  genes <- intersect(genes, all_genes)

  vals <- as.matrix(mat[genes, , drop = FALSE])
  n_cells <- ncol(vals)

  mode <- match.arg(mode)

  results_list <- lapply(seq_along(genes), function(i) {
    x_vec <- vals[i, ]

    if (mode == "moran") {
      stat_value <- compute_morans_i(W, x_vec)
      expected <- -1.0 / (n_cells - 1)
    } else {
      stat_value <- compute_gearys_c(W, x_vec)
      expected <- 1.0
    }

    perm_res <- permutation_test(
      obs_stat = stat_value,
      W = W,
      x = x_vec,
      mode = mode,
      n_perms = n_perms,
      seed = seed
    )

    data.frame(
      feature  = genes[i],
      score    = stat_value,
      expected = expected,
      pval_sim = perm_res$pval_perm
    )
  })

  df <- do.call(rbind, results_list)
  rownames(df) <- df$feature
  df$pval_sim_adj <- p.adjust(df$pval_sim, method = corr_method)

  df <- df[order(df$score, decreasing = (mode == "moran")), ]

  if (copy) {
    return(df)
  } else {
    if (mode == "moran") {
      seurat_obj@misc$moranI <- df
    } else {
      seurat_obj@misc$gearyC <- df
    }
    return(seurat_obj)
  }
}
