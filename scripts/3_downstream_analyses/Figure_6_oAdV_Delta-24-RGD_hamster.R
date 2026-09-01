# ==============================================================================
# Figures 6 and S5 analyses using hamster Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts/3_downstream_analyses
#
# Required input files: 
# ../../hamster/analysis_output_hamster_early_28500/Global_with_final_harmonized_labels_28500.rds
# ../../hamster/analysis_output_hamster_terminal_28488/Global_with_final_harmonized_labels_28488.rds
#
# Output directories:
# ../../plots
# ../../output

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

suppressPackageStartupMessages({
  
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(RANN)
  library(lme4)
  library(lmerTest)
  library(spatstat.geom)
  library(spatstat.random)
  library(patchwork)
  library(ggpubr)
  library(ggrastr)
  
})

set.seed(42)
dir.create("../../plots", showWarnings = FALSE)
dir.create("../../output", showWarnings = FALSE)

# ============================================================
# Figure settings
# ============================================================

mm_to_in <- function(mm) mm / 25.4

P_Width <- 58
P_Height <- 50

pt_size <- 5.8

theme_fig <- function() {
  theme_classic(base_size = 6, base_family = "Arial") +
    theme(
      plot.title   = element_blank(),
      axis.title   = element_text(size = pt_size),
      axis.text    = element_text(size = pt_size),
      legend.text  = element_text(size = pt_size),
      legend.title = element_text(size = pt_size),
      strip.text   = element_text(size = pt_size)
    )
}

save_pdf <- function(filename, plot_obj, w_mm, h_mm) {
  
  ggsave(
    filename,
    plot_obj,
    width  = mm_to_in(w_mm),
    height = mm_to_in(h_mm),
    device = cairo_pdf,
    limitsize = FALSE
  )
  
}

format_p_plotmath <- function(p) {
  
  p_str <- base::format.pval(p, digits = 2)
  
  ifelse(
    grepl("<", p_str),
    paste0("italic(P) < '", sub("<\\s*", "", p_str), "'"),
    paste0("italic(P) == '", p_str, "'")
  )
  
}

pal_best <- c(
  "Myeloid"                   = "#D55E00",
  "Oligodendrocytes/OPC-like" = "#E69F00",
  "Astrocytes/AC-like"        = "#F0E442",
  "MES-like/Hypoxic Tumor"    = "#7F7F7F",
  "Vasculature/Stroma"        = "#0072B2",
  "T/NK Cells"                = "#009E73",
  "Neurons"                   = "#7cec37",
  "B/Plasma Cells"            = "#CC79A7",
  "Other"                     = "#000000"
)

# ============================================================
# Parameters
# ============================================================

MY_R <- "Myeloid"
B_R  <- "B/Plasma"
T_R  <- "T/NK"
VASC_R <- "Vasculature/Stroma"

r_tri <- 30

NEAR_MAX <- 60
FAR_MIN  <- 300

ASSAY_USE <- "Xenium"

CELL_TYPE_COL <- "final_harmonized_label"

VIRAL_GENES <- c(
  "E1A","E1B","E2A","E2B","E3","E4","IVa2",
  "IX","L1","L2","L3","L4","L5"
)

TUMOR_LABELS <- c(
  "Oligodendrocytes/OPC-like",
  "Astrocytes/AC-like",
  "MES-like/Hypoxic Tumor"
)

# ============================================================
# Data loading
# ============================================================

EARLY_RDS <- paste0(
  "../../hamster/analysis_output_hamster_early_28500/",
  "Global_with_final_harmonized_labels_28500.rds"
)

TERMINAL_RDS <- paste0(
  "../../hamster/analysis_output_hamster_terminal_28488/",
  "Global_with_final_harmonized_labels_28488.rds"
)

xe_early <- readRDS(EARLY_RDS)
xe_terminal <- readRDS(TERMINAL_RDS)

xe_early <- UpdateSeuratObject(xe_early)
xe_terminal <- UpdateSeuratObject(xe_terminal)

DefaultAssay(xe_early)    <- ASSAY_USE
DefaultAssay(xe_terminal) <- ASSAY_USE

xe_early$Timepoint <- "Early"
xe_terminal$Timepoint <- "Terminal"

# ============================================================
# Merge objects
# ============================================================

xe_hamster <- merge(
  xe_early,
  y = xe_terminal,
  add.cell.ids = c(
    "Early",
    "Terminal"
  )
)

# ============================================================
# Build metadata with coordinates
# ============================================================

raw_coords <- do.call(
  
  rbind,
  
  lapply(
    Images(xe_hamster),
    function(fov) {
      
      GetTissueCoordinates(
        xe_hamster[[fov]]
      ) |>
        dplyr::mutate(
          FOV = fov
        )
      
    }
  )
  
)

meta <- xe_hamster@meta.data |>
  tibble::rownames_to_column(
    "cell"
  ) |>
  dplyr::left_join(
    raw_coords,
    by = "cell"
  )

# ============================================================
# Metadata
# ============================================================

meta$cell_type <- as.character(
  meta[[CELL_TYPE_COL]]
)

meta$Timepoint <- factor(
  meta$Timepoint,
  levels = c(
    "Early",
    "Terminal"
  )
)

meta$spot_id <- as.character(
  meta$spot_id
)

meta$HamsterID <- as.character(
  meta$HamsterID
)

meta$FOV <- as.character(
  meta$FOV
)

meta$Treatment <- ifelse(
  meta$Condition == "RGD",
  "Virus",
  as.character(meta$Condition)
)

meta <- meta %>%
  dplyr::filter(
    Treatment %in% c(
      "PBS",
      "Virus"
    )
  )

meta$Treatment <- factor(
  meta$Treatment,
  levels = c(
    "PBS",
    "Virus"
  )
)

# ============================================================
# Unified biological sample identifier
# ============================================================

meta$Sample <- NA_character_

# Terminal samples
meta$Sample[
  meta$Timepoint == "Terminal"
] <- meta$HamsterID[
  meta$Timepoint == "Terminal"
]

# Early samples
early_idx <- meta$Timepoint == "Early"

meta$Sample[early_idx] <-
  gsub(
    ".*?(PBS|RGD)\\s+([0-9]+).*",
    "\\1 \\2",
    meta$spot_id[early_idx]
  )

# ============================================================
# Diagnostics
# ============================================================

message(
  "\nCell type counts:\n"
)

print(
  sort(
    table(meta$cell_type),
    decreasing = TRUE
  )
)

message(
  "\nSample summary:\n"
)

print(
  meta %>%
    dplyr::distinct(
      Sample,
      Timepoint,
      Treatment
    ) %>%
    dplyr::arrange(
      Timepoint,
      Treatment,
      Sample
    )
)

message(
  "\nSamples by condition:\n"
)

print(
  meta %>%
    dplyr::distinct(
      Sample,
      Timepoint,
      Treatment
    ) %>%
    dplyr::count(
      Timepoint,
      Treatment
    )
)

meta %>%
  filter(Timepoint == "Early") %>%
  distinct(
    spot_id,
    Sample
  ) %>%
  arrange(Sample)

# ============================================================
# Virus-positive calls
# ============================================================

viral_expr <- FetchData(
  xe_hamster,
  vars = intersect(
    VIRAL_GENES,
    rownames(xe_hamster)
  ),
  cells = meta$cell,
  layer = "counts"
)

meta$Virus_Pos <- rowSums(
  viral_expr > 0
) >= 2


# ============================================================
# Helper functions
# ============================================================

# KDE function
gauss_kde_at <- function(
    query_xy,
    anchor_xy,
    bw = 100,
    k = 100
) {
  
  if (
    nrow(query_xy) == 0 ||
    nrow(anchor_xy) == 0
  ) {
    
    return(
      rep(0, nrow(query_xy))
    )
    
  }
  
  k_use <- min(
    k,
    nrow(anchor_xy)
  )
  
  nn <- RANN::nn2(
    as.matrix(anchor_xy),
    as.matrix(query_xy),
    k = k_use
  )
  
  rowSums(
    exp(
      -(nn$nn.dists^2) /
        (2 * bw^2)
    )
  )
  
}

# Nearest-neighbor distance
nn_dists <- function(
    query_xy,
    ref_xy
) {
  
  if (
    nrow(query_xy) == 0 ||
    nrow(ref_xy) == 0
  ) {
    
    return(
      rep(Inf, nrow(query_xy))
    )
    
  }
  
  RANN::nn2(
    as.matrix(ref_xy),
    as.matrix(query_xy),
    k = 1
  )$nn.dists[, 1]
  
}

# Distance to nearest reference set
dist_to_set <- function(
    query_df,
    ref_df
) {
  
  if (
    nrow(query_df) == 0 ||
    nrow(ref_df) == 0
  ) {
    
    return(
      rep(Inf, nrow(query_df))
    )
    
  }
  
  RANN::nn2(
    data  = as.matrix(
      ref_df[, c("x", "y")]
    ),
    query = as.matrix(
      query_df[, c("x", "y")]
    ),
    k = 1
  )$nn.dists[, 1]
  
}

# Define myeloid triad anchors
triad_flags_myeloid <- function(df) {
  
  A <- subset(
    df,
    grepl(
      MY_R,
      cell_type,
      ignore.case = TRUE
    )
  )
  
  B <- subset(
    df,
    grepl(
      B_R,
      cell_type,
      ignore.case = TRUE
    )
  )
  
  T <- subset(
    df,
    grepl(
      T_R,
      cell_type,
      ignore.case = TRUE
    )
  )
  
  if (
    nrow(A) == 0 ||
    nrow(B) == 0 ||
    nrow(T) == 0
  ) {
    
    return(
      data.frame(
        cell = character(),
        Triad = logical()
      )
    )
    
  }
  
  dB <- RANN::nn2(
    B[, c("x", "y")],
    A[, c("x", "y")],
    k = 1
  )$nn.dists[, 1]
  
  dT <- RANN::nn2(
    T[, c("x", "y")],
    A[, c("x", "y")],
    k = 1
  )$nn.dists[, 1]
  
  data.frame(
    cell = A$cell,
    Triad = (dB <= r_tri) &
      (dT <= r_tri)
  )
  
}

# Return nearby cells for each anchor
cells_within_radius <- function(
    anchors,
    targets,
    radius = r_tri,
    k_max = 256
) {
  
  if (
    nrow(anchors) == 0 ||
    nrow(targets) == 0
  ) {
    
    return(
      vector(
        "list",
        nrow(anchors)
      )
    )
    
  }
  
  k_use <- min(
    k_max,
    nrow(targets)
  )
  
  nn <- RANN::nn2(
    data = as.matrix(
      targets[, c("x", "y")]
    ),
    query = as.matrix(
      anchors[, c("x", "y")]
    ),
    k = k_use
  )
  
  out <- vector(
    "list",
    nrow(anchors)
  )
  
  for (i in seq_len(nrow(anchors))) {
    
    keep <- which(
      nn$nn.dists[i, ] <= radius
    )
    
    out[[i]] <- if (
      length(keep)
    ) {
      
      targets$cell[
        nn$nn.idx[i, keep]
      ]
      
    } else {
      
      character(0)
      
    }
    
  }
  
  out
  
}

# Mean module score across genes
module_score_vec <- function(
    object,
    cells,
    genes,
    layer = "data"
) {
  
  genes <- intersect(
    genes,
    rownames(object)
  )
  
  if (
    length(genes) == 0 ||
    length(cells) == 0
  ) {
    
    return(
      setNames(
        rep(
          NA_real_,
          length(cells)
        ),
        cells
      )
    )
    
  }
  
  expr <- FetchData(
    object,
    vars  = genes,
    cells = cells,
    layer = layer
  )
  
  setNames(
    rowMeans(
      expr,
      na.rm = TRUE
    ),
    cells
  )
  
}

# Bootstrap CI for Near vs Far comparisons
bootstrap_delta_ci <- function(
    x,
    group,
    B = 2000
) {
  
  keep <- is.finite(x) &
    !is.na(group)
  
  x <- x[keep]
  
  group <- as.character(
    group[keep]
  )
  
  if (
    length(unique(group)) < 2 ||
    sum(group == "Near") < 5 ||
    sum(group == "Far") < 5
  ) {
    
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
    
  }
  
  near_vals <- x[
    group == "Near"
  ]
  
  far_vals <- x[
    group == "Far"
  ]
  
  deltas <- replicate(
    B,
    median(
      sample(
        near_vals,
        replace = TRUE
      )
    ) -
      median(
        sample(
          far_vals,
          replace = TRUE
        )
      )
  )
  
  quantile(
    deltas,
    c(0.025, 0.975),
    na.rm = TRUE
  )
  
}

# Cell-type frequency helper
celltype_freq <- function(df) {
  
  df %>%
    count(
      Timepoint,
      cell_type
    ) %>%
    group_by(Timepoint) %>%
    mutate(
      freq = n / sum(n)
    ) %>%
    ungroup()
  
}

# ============================================================
# Identify triad anchors
# ============================================================
# A myeloid cell is a triad anchor if:
# nearest B/Plasma cell <= r_tri
# AND
# nearest T/NK cell <= r_tri

triads <- meta |>
  dplyr::group_by(
    Sample,
    FOV
  ) |>
  dplyr::group_modify(
    ~ triad_flags_myeloid(.x)
  ) |>
  dplyr::ungroup()

# ============================================================
# Store triad calls
# ============================================================

meta <- meta |>
  dplyr::select(
    -dplyr::any_of("Triad")
  ) |>
  dplyr::left_join(
    dplyr::select(
      triads,
      cell,
      Triad
    ),
    by = "cell"
  )

meta$Triad[
  is.na(meta$Triad)
] <- FALSE

# Triad anchor cell IDs
triad_anchors <- triads |>
  dplyr::filter(Triad) |>
  dplyr::pull(cell)

N_TRIADS <- length(
  triad_anchors
)

message(
  "\nTotal triad anchors: ",
  N_TRIADS
)

triad_xy <- meta %>%
  dplyr::filter(
    cell %in% triad_anchors
  ) %>%
  dplyr::select(
    cell,
    Sample,
    FOV,
    x,
    y
  )

# ============================================================
# Label triad niche
# ============================================================
# A target cell is considered "In triad" if its nearest
# triad-anchor myeloid cell is within r_tri.

label_specific_niche <- function(
    df,
    anchors,
    target_regex
){
  
  anc <- subset(
    df,
    cell %in% anchors
  )
  
  tgt <- subset(
    df,
    grepl(
      target_regex,
      cell_type,
      ignore.case = TRUE
    )
  )
  
  # no target cells
  if (nrow(tgt) == 0) {
    
    return(
      data.frame(
        cell = character(),
        Specific_Niche = character()
      )
    )
    
  }
  
  # no triad anchors in this Sample/FOV
  # -> every target cell is outside a triad
  if (nrow(anc) == 0) {
    
    return(
      data.frame(
        cell = tgt$cell,
        Specific_Niche = "Outside triad"
      )
    )
    
  }
  
  d <- nn_dists(
    tgt[, c("x", "y")],
    anc[, c("x", "y")]
  )
  
  data.frame(
    cell = tgt$cell,
    Specific_Niche = ifelse(
      d <= r_tri,
      "In triad",
      "Outside triad"
    )
  )
  
}

# Label Myeloid niche
myeloid_niche <- meta |>
  dplyr::group_by(Sample, FOV) |>
  dplyr::group_modify(
    ~label_specific_niche(
      .x,
      anchors = triad_anchors,
      target_regex = MY_R
    )
  ) |>
  dplyr::ungroup()

# Label T/NK niche
tcell_niche <- meta |>
  dplyr::group_by(Sample, FOV) |>
  dplyr::group_modify(
    ~label_specific_niche(
      .x,
      anchors = triad_anchors,
      target_regex = T_R
    )
  ) |>
  dplyr::ungroup()

# Label B/Plasma niche
bcell_niche <- meta |>
  dplyr::group_by(Sample, FOV) |>
  dplyr::group_modify(
    ~label_specific_niche(
      .x,
      anchors = triad_anchors,
      target_regex = B_R
    )
  ) |>
  dplyr::ungroup()

# ============================================================
# Store niche labels
# ============================================================

meta <- meta |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "Myeloid_Niche",
        "Tcell_Niche",
        "Bcell_Niche"
      )
    )
  ) |>
  dplyr::left_join(
    myeloid_niche |>
      dplyr::select(
        cell,
        Specific_Niche
      ) |>
      dplyr::rename(
        Myeloid_Niche = Specific_Niche
      ),
    by = "cell"
  ) |>
  dplyr::left_join(
    tcell_niche |>
      dplyr::select(
        cell,
        Specific_Niche
      ) |>
      dplyr::rename(
        Tcell_Niche = Specific_Niche
      ),
    by = "cell"
  ) |>
  dplyr::left_join(
    bcell_niche |>
      dplyr::select(
        cell,
        Specific_Niche
      ) |>
      dplyr::rename(
        Bcell_Niche = Specific_Niche
      ),
    by = "cell"
  )

# ============================================================
# UMAP
# ============================================================

# Early timepoint UMAP
red_name <- intersect(
  c("full.umap", "umap"),
  Reductions(xe_early)
)[1]

stopifnot(!is.na(red_name))

emb <- Embeddings(
  xe_early,
  red_name
)

umap_df_early <- data.frame(
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  final_harmonized_label =
    xe_early$final_harmonized_label
)

p_umap_early <- ggplot(
  umap_df_early,
  aes(
    UMAP_1,
    UMAP_2,
    color = final_harmonized_label
  )
) +
  geom_point_rast(
    size = 0.32,
    raster.dpi = 600
  ) +
  scale_color_manual(
    values = pal_best
  ) +
  theme_fig() +
  labs(
    title = "Early",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "none"
  )

# Terminal timepoint UMAP
red_name <- intersect(
  c("full.umap", "umap"),
  Reductions(xe_terminal)
)[1]

stopifnot(!is.na(red_name))

emb <- Embeddings(
  xe_terminal,
  red_name
)

umap_df_terminal <- data.frame(
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  final_harmonized_label =
    xe_terminal$final_harmonized_label
)

p_umap_terminal <- ggplot(
  umap_df_terminal,
  aes(
    UMAP_1,
    UMAP_2,
    color = final_harmonized_label
  )
) +
  geom_point_rast(
    size = 0.32,
    raster.dpi = 600
  ) +
  scale_color_manual(
    values = pal_best
  ) +
  theme_fig() +
  labs(
    title = "Terminal",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "none"
  )

# Combined
p_umaps <-
  p_umap_early +
  p_umap_terminal +
  plot_layout(ncol = 2)

# Legend
cell_order <- c(
  "Oligodendrocytes/OPC-like",
  "Astrocytes/AC-like",
  "MES-like/Hypoxic Tumor",
  "Neurons",
  "Vasculature/Stroma",
  "B/Plasma Cells",
  "T/NK Cells",
  "Myeloid",
  "Other"
)

cell_df <- data.frame(
  Label = factor(
    cell_order,
    levels = cell_order
  ),
  x = seq_along(cell_order),
  y = 1
)

p_cell <- ggplot(
  cell_df,
  aes(
    x,
    y,
    fill = Label
  )
) +
  
  geom_point(
    shape = 22,
    size = 3.5,
    colour = "black",
    stroke = 0.3
  ) +
  
  scale_fill_manual(
    values = pal_best,
    breaks = cell_order,
    name = NULL
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "right",
    legend.key.height = unit(0.18,"cm"),
    legend.key.width  = unit(0.18,"cm")
  )

cell_legend <- cowplot::get_legend(
  p_cell
)

legend_plot <- wrap_elements(
  full = cowplot::ggdraw(cell_legend)
)

final_panel <-
  p_umap_early +
  p_umap_terminal +
  legend_plot +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 0.6)
  )

save_pdf(
  "../../plots/FigS5_UMAPs.pdf",
  final_panel,
  2.5 * P_Width,
  P_Height
)

# ============================================================
# Cell-type composition
# ============================================================

cell_order <- c(
  "Other",
  "Myeloid",
  "T/NK Cells",
  "B/Plasma Cells",
  "Vasculature/Stroma",
  "Neurons",
  "MES-like/Hypoxic Tumor",
  "Astrocytes/AC-like",
  "Oligodendrocytes/OPC-like"
)

comp_bio <- meta %>%
  dplyr::filter(
    Treatment %in% c("PBS", "Virus")
  ) %>%
  dplyr::group_by(
    Sample,
    Timepoint,
    Treatment,
    cell_type
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    Sample,
    Timepoint,
    Treatment
  ) %>%
  dplyr::mutate(
    Freq = n / sum(n)
  ) %>%
  dplyr::ungroup()

comp_bio$Timepoint <- factor(
  comp_bio$Timepoint,
  levels = c(
    "Early",
    "Terminal"
  )
)

comp_bio$Treatment <- factor(
  comp_bio$Treatment,
  levels = c(
    "PBS",
    "Virus"
  )
)

comp_bio$cell_type <- factor(
  comp_bio$cell_type,
  levels = rev(cell_order)
)

comp_early <- comp_bio %>%
  filter(Timepoint == "Early")

comp_terminal <- comp_bio %>%
  filter(Timepoint == "Terminal")

p_comp_early <- ggplot(
  comp_early,
  aes(
    x = Sample,
    y = Freq,
    fill = cell_type
  )
) +
  geom_col(
    width = 0.85,
    position = "fill"
  ) +
  facet_grid(
    ~ Treatment,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = pal_best,
    drop = FALSE
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 4
    ),
    legend.position = "none"
  ) +
  labs(
    title = "Early",
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1,
      margin = margin(
        b = 5
      )
    )
  )

p_comp_terminal <- ggplot(
  comp_terminal,
  aes(
    x = Sample,
    y = Freq,
    fill = cell_type
  )
) +
  geom_col(
    width = 0.85,
    position = "fill"
  ) +
  facet_grid(
    ~ Treatment,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = pal_best,
    drop = FALSE
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 4
    ),
    legend.position = "none"
  ) +
  labs(
    title = "Terminal",
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1,
      margin = margin(
        b = 5
      )
    )
  )

p_comp_pair <-
  p_comp_early +
  p_comp_terminal +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

save_pdf(
  "../../plots/Fig6_Composition.pdf",
  p_comp_pair,
  1.5 * P_Width,
  P_Height
)

# ============================================================
# Triad fraction
# ============================================================

triad_summary <- meta %>%
  dplyr::filter(
    grepl(
      MY_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::group_by(
    Sample,
    Timepoint,
    Treatment
  ) %>%
  dplyr::summarise(
    n_myeloid  = dplyr::n(),
    n_triad    = sum(Triad),
    triad_frac = n_triad / n_myeloid,
    .groups = "drop"
  )

print(triad_summary)

triad_summary %>%
  dplyr::count(
    Timepoint,
    Treatment
  )

triad_summary$Group <- factor(
  paste(
    triad_summary$Timepoint,
    triad_summary$Treatment
  ),
  levels = c(
    "Early PBS",
    "Early Virus",
    "Terminal PBS",
    "Terminal Virus"
  )
)

# Overall Kruskal-Wallis test
p_kw <- kruskal.test(
  triad_frac ~ Group,
  data = triad_summary
)$p.value

# Pairwise Wilcoxon test
pw <- pairwise.wilcox.test(
  triad_summary$triad_frac,
  triad_summary$Group,
  p.adjust.method = "BH"
)

pw$p.value

sig_df <- data.frame(
  group1 = c(
    "Early PBS",
    "Terminal PBS",
    "Early Virus"
  ),
  group2 = c(
    "Early Virus",
    "Terminal Virus",
    "Terminal Virus"
  ),
  p = c(
    pw$p.value["Early Virus", "Early PBS"],
    pw$p.value["Terminal Virus", "Terminal PBS"],
    pw$p.value["Terminal Virus", "Early Virus"]
  ),
  y.position = c(
    max(triad_summary$triad_frac) * 1.05,
    max(triad_summary$triad_frac) * 1.15,
    max(triad_summary$triad_frac) * 1.25
  )
)

sig_df$p.adj.signif <- dplyr::case_when(
  sig_df$p < 0.001 ~ "***",
  sig_df$p < 0.01  ~ "**",
  sig_df$p < 0.05  ~ "*",
  TRUE ~ ""
)

sig_df <- sig_df %>%
  dplyr::filter(
    p.adj.signif != ""
  )

p_tri <- ggplot(
  triad_summary,
  aes(
    x = Group,
    y = triad_frac,
    fill = Treatment
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.6
  ) +
  scale_fill_manual(
    values = c(
      PBS = "#BDBDBD",
      Virus = "#08306B"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Early PBS"      = "Early\nPBS",
      "Early Virus"    = "Early\nVirus",
      "Terminal PBS"   = "Terminal\nPBS",
      "Terminal Virus" = "Terminal\nVirus"
    )
  )

p_tri_final <- p_tri +
  stat_pvalue_manual(
    sig_df,
    label = "p.adj.signif",
    tip.length = 0.01
  ) +
  theme_fig() +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Triad fraction among myeloid cells"
  )

save_pdf(
  "../../plots/Fig6_Triad_Fraction.pdf",
  p_tri_final,
  P_Width,
  P_Height
)

# ============================================================
# Triad-associated immune programs
# ============================================================

MODULES <- list(
  
  # Myeloid
  "Inflammatory chemokine" = c(
    "Cxcl9","Cxcl10","Cxcl11","Ccl5"
  ),
  
  "Antigen presentation" = c(
    "Cd74","Ctss","Cd80","Cd86"
  ),
  
  # T/NK
  "Stemness" = c(
    "Tcf7","Il7r"
  ),
  
  "Cytotoxicity" = c(
    "Prf1","Nkg7"
  ),
  
  "Exhaustion" = c(
    "Lag3","Tox"
  ),
  
  "Regulatory T" = c(
    "Foxp3","Il2ra","Ctla4","Tigit"
  ),
  
  # B cells
  "Lymphoid organization" = c(
    "Cxcl13","Ccr7"
  ),
  
  "Regulatory B" = c(
    "Il10","Tgfb1"
  )
  
)

MODULES <- lapply(
  MODULES,
  intersect,
  rownames(xe_hamster)
)

MODULES <- MODULES[
  sapply(MODULES, length) >= 2
]

module_info <- tibble::tribble(
  
  ~Module,                  ~Regex, ~Group,
  
  "Inflammatory chemokine", MY_R,   "Myeloid",
  "Antigen presentation",   MY_R,   "Myeloid",
  
  "Stemness",               T_R,    "T/NK Cells",
  "Cytotoxicity",           T_R,    "T/NK Cells",
  "Exhaustion",             T_R,    "T/NK Cells",
  "Regulatory T",           T_R,    "T/NK Cells",
  
  "Lymphoid organization",  B_R,    "B/Plasma Cells",
  "Regulatory B",           B_R,    "B/Plasma Cells"
  
)

plot_data_list <- list()
ext_lmm_list <- list()
sample_effects_list <- list()

for(tp in c("Early","Terminal")) {
  
  message("Processing: ", tp)
  
  meta_tp <- meta %>%
    dplyr::filter(
      Treatment == "Virus",
      Timepoint == tp
    )
  
  for(mod_name in names(MODULES)) {
    
    info <- module_info %>%
      dplyr::filter(
        Module == mod_name
      )
    
    genes <- MODULES[[mod_name]]
    
    if(length(genes) < 2)
      next
    
    niche_df <- switch(
      
      info$Group,
      
      "Myeloid" =
        meta_tp %>%
        dplyr::select(
          cell,
          Specific_Niche = Myeloid_Niche
        ),
      
      "T/NK Cells" =
        meta_tp %>%
        dplyr::select(
          cell,
          Specific_Niche = Tcell_Niche
        ),
      
      "B/Plasma Cells" =
        meta_tp %>%
        dplyr::select(
          cell,
          Specific_Niche = Bcell_Niche
        )
      
    )
    
    df_lmm <- meta_tp %>%
      dplyr::filter(
        grepl(
          info$Regex,
          cell_type,
          ignore.case = TRUE
        )
      ) %>%
      dplyr::select(
        cell,
        Sample,
        FOV
      ) %>%
      dplyr::inner_join(
        niche_df,
        by = "cell"
      ) %>%
      dplyr::filter(
        Specific_Niche %in%
          c(
            "Outside triad",
            "In triad"
          )
      )
    
    if(nrow(df_lmm) < 100)
      next
    
    expr_mat <- FetchData(
      xe_hamster,
      vars = genes,
      cells = df_lmm$cell,
      layer = "data"
    )
    
    module_score <- rowMeans(
      expr_mat,
      na.rm = TRUE
    )
    
    df_lmm$Module_Score <- module_score
    
    df_lmm$Score_Z <- as.numeric(
      scale(module_score)
    )
    
    df_lmm$Timepoint <- tp
    df_lmm$Module <- mod_name
    
    plot_data_list[[paste(tp, mod_name, sep = "_")]] <- df_lmm
    
    df_lmm$Specific_Niche <- factor(
      df_lmm$Specific_Niche,
      levels = c(
        "Outside triad",
        "In triad"
      )
    )
    
    # Percent change (raw score)
    mean_outside <- mean(
      df_lmm$Module_Score[
        df_lmm$Specific_Niche=="Outside triad"
      ],
      na.rm = TRUE
    )
    
    mean_inside <- mean(
      df_lmm$Module_Score[
        df_lmm$Specific_Niche=="In triad"
      ],
      na.rm = TRUE
    )
    
    pct_change <- ifelse(
      abs(mean_outside) > 1e-8,
      100 * (mean_inside - mean_outside) /
        abs(mean_outside),
      NA_real_
    )
    
    fit <- tryCatch(
      
      lmer(
        Score_Z ~
          Specific_Niche +
          (1 | Sample/FOV),
        data = df_lmm
      ),
      
      error = function(e) NULL
      
    )
    
    if(is.null(fit))
      next
    
    coef_tab <- summary(fit)$coefficients
    
    rn <- rownames(coef_tab)[
      grepl(
        "Specific_Niche",
        rownames(coef_tab)
      )
    ][1]
    
    beta <- coef_tab[rn,"Estimate"]
    se   <- coef_tab[rn,"Std. Error"]
    pval <- coef_tab[rn,"Pr(>|t|)"]
    
    ext_lmm_list[[paste(tp, mod_name, sep="_")]] <-
      data.frame(
        
        Timepoint = tp,
        Module = mod_name,
        Group = info$Group,
        
        Mean_Outside = mean_outside,
        Mean_InTriad = mean_inside,
        
        Percent_Change = pct_change,
        
        Beta = beta,
        P = pval,
        
        CI_low  = beta - 1.96 * se,
        CI_high = beta + 1.96 * se,
        
        stringsAsFactors = FALSE
      )
    
    # Sample-specific estimates
    for(sm in unique(df_lmm$Sample)) {
      
      tmp <- df_lmm %>%
        dplyr::filter(
          Sample == sm
        )
      
      if(nrow(tmp) < 100)
        next
      
      if(length(unique(tmp$Specific_Niche)) < 2)
        next
      
      mean_out_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche=="Outside triad"
        ],
        na.rm = TRUE
      )
      
      mean_in_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche=="In triad"
        ],
        na.rm = TRUE
      )
      
      pct_sm <- ifelse(
        abs(mean_out_sm) > 1e-8,
        100 * (mean_in_sm - mean_out_sm) /
          abs(mean_out_sm),
        NA_real_
      )
      
      sample_effects_list[[paste(tp,mod_name,sm,sep="_")]] <- 
        data.frame(
        
        Timepoint = tp,
        Module = mod_name,
        Sample = sm,
        
        Percent_Change = pct_sm,
        
        stringsAsFactors = FALSE
      )
    }
    
  }
  
}

# Combine results
ext_stats <- bind_rows(
  ext_lmm_list
)

sample_effects <- bind_rows(
  sample_effects_list
)

# Consistency metrics
consistency_summary <- sample_effects %>%
  group_by(
    Timepoint,
    Module
  ) %>%
  summarise(
    
    N_Pos = sum(
      Percent_Change > 0,
      na.rm = TRUE
    ),
    
    N_Neg = sum(
      Percent_Change < 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

ext_stats <- ext_stats %>%
  dplyr::left_join(
    consistency_summary,
    by = c(
      "Timepoint",
      "Module"
    )
  )

ext_stats$Consistency <- ifelse(
  
  ext_stats$N_Pos >= ext_stats$N_Neg,
  
  paste0(
    ext_stats$N_Pos,
    "/",
    ext_stats$N_Pos + ext_stats$N_Neg,
    "+"
  ),
  
  paste0(
    ext_stats$N_Neg,
    "/",
    ext_stats$N_Pos + ext_stats$N_Neg,
    "-"
  )
)

# Multiple testing correction
ext_stats <- ext_stats %>%
  dplyr::group_by(
    Timepoint
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(
      P,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

ext_stats$Sig <- dplyr::case_when(
  ext_stats$FDR < 0.001 ~ "***",
  ext_stats$FDR < 0.01  ~ "**",
  ext_stats$FDR < 0.05  ~ "*",
  TRUE ~ ""
)

# Sample-level validation
specimen_validation <- list()

for(nm in names(plot_data_list)) {
  
  df_plot <- plot_data_list[[nm]]
  
  tp <- unique(df_plot$Timepoint)
  
  mod <- strsplit(nm,"_")[[1]]
  mod <- paste(mod[-1], collapse = "_")
  
  sample_summary <- df_plot %>%
    group_by(
      Sample,
      Specific_Niche
    ) %>%
    summarise(
      Mean_Score = mean(
        Module_Score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Specific_Niche,
      values_from = Mean_Score
    ) %>%
    mutate(
      Delta = `In triad` -
        `Outside triad`
    )
  
  wilcox_p <- tryCatch(
    
    wilcox.test(
      sample_summary$Delta,
      mu = 0,
      exact = FALSE
    )$p.value,
    
    error = function(e) NA_real_
  )
  
  specimen_validation[[nm]] <- data.frame(
    
    Timepoint = tp,
    
    Module = mod,
    
    N_Samples = nrow(sample_summary),
    
    Positive = sum(
      sample_summary$Delta > 0,
      na.rm = TRUE
    ),
    
    Negative = sum(
      sample_summary$Delta < 0,
      na.rm = TRUE
    ),
    
    Wilcox_P = wilcox_p,
    
    stringsAsFactors = FALSE
  )
}

specimen_validation_df <-
  bind_rows(specimen_validation)

specimen_validation_df$Wilcox_FDR <-
  p.adjust(
    specimen_validation_df$Wilcox_P,
    method = "BH"
  )

# Save output
final_results <- ext_stats %>%
  left_join(
    specimen_validation_df,
    by = c(
      "Timepoint",
      "Module"
    )
  )

write.csv(
  final_results,
  "../../output/Fig6_Triad_Programs.csv",
  row.names = FALSE
)

# Ordering
module_order <- c(
  
  "Inflammatory chemokine",
  "Antigen presentation",
  
  "Stemness",
  "Cytotoxicity",
  "Exhaustion",
  "Regulatory T",
  
  "Lymphoid organization",
  "Regulatory B"
  
)

ext_stats$Module <- factor(
  ext_stats$Module,
  levels = rev(module_order)
)

sample_effects$Module <- factor(
  sample_effects$Module,
  levels = rev(module_order)
)

ext_stats$Group <- factor(
  ext_stats$Group,
  levels = c(
    "Myeloid",
    "T/NK Cells",
    "B/Plasma Cells"
  )
)

# Module summary dot plot
module_dot <- ext_stats %>%
  dplyr::select(
    Module,
    Group,
    Beta,
    FDR
  )

pct_df <- lapply(
  names(plot_data_list),
  function(nm) {
    
    df_plot <- plot_data_list[[nm]]
    
    data.frame(
      
      Timepoint = unique(df_plot$Timepoint),
      
      Module = unique(df_plot$Module),
      
      Percent_Positive =
        mean(
          df_plot$Module_Score > 0,
          na.rm = TRUE
        ) * 100
    )
  }
) %>%
  bind_rows()

module_dot <- ext_stats %>%
  dplyr::select(
    Timepoint,
    Module,
    Group,
    Beta,
    FDR,
    Sig,
    Consistency
  ) %>%
  left_join(
    pct_df,
    by = c(
      "Timepoint",
      "Module"
    )
  )

module_dot$Group <- factor(
  module_dot$Group,
  levels = c(
    "Myeloid",
    "T/NK Cells",
    "B/Plasma Cells"
  )
)

module_dot$Module <- factor(
  module_dot$Module,
  levels = rev(c(
    "Inflammatory chemokine",
    "Antigen presentation",
    
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T",
    
    "Lymphoid organization",
    "Regulatory B"
  ))
)

beta_lim <- max(
  abs(module_dot$Beta),
  na.rm = TRUE
)

p_module_summary <- ggplot(
  module_dot,
  aes(
    x = 1,
    y = Module
  )
) +
  
  geom_point(
    aes(
      fill = Beta
    ),
    size = 5,
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +
  
  geom_text(
    aes(
      x = 1.2,
      label = Sig
    ),
    size = pt_size / 2
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -beta_lim,
      beta_lim
    ),
    name = "β"
  ) +
  
  scale_x_continuous(
    limits = c(
      0.85,
      1.25
    ),
    breaks = NULL
  ) +
  
  facet_grid(
    Group ~ Timepoint,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(1, "cm"),
      barwidth = unit(0.3, "cm")
    )
  ) +
  
  theme_fig() +
  
  theme(
    
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),
    
    axis.text.x = element_blank(),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      size = pt_size * 0.7
    ),
    
    legend.text = element_text(
      size = pt_size * 0.7
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

save_pdf(
  "../../plots/Fig6_Triad_Immune_Modules.pdf",
  p_module_summary,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Tumor cells
# ============================================================

tumor_meta <- meta %>%
  dplyr::filter(
    cell_type %in% TUMOR_LABELS
  )

tumor_meta$Tumor_State <- factor(
  tumor_meta$cell_type,
  levels = TUMOR_LABELS
)

tumor_meta$cell <- as.character(
  tumor_meta$cell
)

AP_GENE_POOL <- c(
  "Ctss",
  "Ctsb",
  "Cd74"
)

IFN_GENE_POOL <- c(
  "Cxcl9",
  "Cxcl10",
  "Isg15",
  "Ifit2",
  "Ifit3",
  "Mx1"
)

HYPOXIA_GENE_POOL <- c(
  "Vegfa",
  "Flt1",
  "Mmp9"
)

TUMOR_MODULES <- list(
  
  "Antigen processing" =
    AP_GENE_POOL,
  
  "Interferon response" =
    IFN_GENE_POOL,
  
  "Hypoxia/angiogenesis" =
    HYPOXIA_GENE_POOL
  
)

TUMOR_MODULES <- lapply(
  TUMOR_MODULES,
  intersect,
  rownames(xe_hamster)
)

TUMOR_MODULES <- TUMOR_MODULES[
  sapply(TUMOR_MODULES, length) >= 2
]

print(
  sapply(
    TUMOR_MODULES,
    length
  )
)

module_info <- tibble::tribble(
  
  ~Module,                  ~Group,
  
  "Antigen processing",     "Tumor",
  "Interferon response",    "Tumor",
  "Hypoxia/angiogenesis",   "Tumor"
  
)

# Compute tumor triad density
tumor_meta$triad_density <- 0

for(sid in unique(tumor_meta$Sample)) {
  
  for(fov in unique(
    tumor_meta$FOV[
      tumor_meta$Sample == sid
    ]
  )) {
    
    idx <- which(
      tumor_meta$Sample == sid &
        tumor_meta$FOV == fov
    )
    
    tri_idx <- which(
      triad_xy$Sample == sid &
        triad_xy$FOV == fov
    )
    
    if(
      length(idx) == 0 ||
      length(tri_idx) == 0
    ) {
      next
    }
    
    tumor_meta$triad_density[idx] <-
      gauss_kde_at(
        
        query_xy =
          tumor_meta[
            idx,
            c("x","y")
          ],
        
        anchor_xy =
          triad_xy[
            tri_idx,
            c("x","y")
          ]
        
      )
    
  }
  
}

tumor_meta$TriadDensity_Z <-
  as.numeric(
    scale(
      log1p(
        tumor_meta$triad_density
      )
    )
  )

# Mixed model analysis
tumor_plot_data <- list()

tumor_lmm_list <- list()

tumor_sample_effects <- list()

for(tp in c(
  "Early",
  "Terminal"
)) {
  
  message(
    "Processing: ",
    tp
  )
  
  tumor_tp <- tumor_meta %>%
    dplyr::filter(
      Timepoint == tp,
      Treatment == "Virus"
    )
  
  for(mod_name in names(
    TUMOR_MODULES
  )) {
    
    genes <- TUMOR_MODULES[[mod_name]]
    
    expr_mat <- FetchData(
      xe_hamster,
      vars = genes,
      cells = tumor_tp$cell,
      layer = "data"
    )
    
    score <- rowMeans(
      expr_mat,
      na.rm = TRUE
    )
    
    if(sd(score, na.rm = TRUE) == 0)
      next
    
    df_lmm <- tumor_tp
    
    df_lmm$Module_Score <- score
    
    df_lmm$Score_Z <- as.numeric(
      scale(score)
    )
    
    df_lmm$Timepoint <- tp
    df_lmm$Module <- mod_name
    
    tumor_plot_data[[paste(tp, mod_name, sep = "_")]] <- df_lmm
    
    fit <- tryCatch(
      
      lmer(
        Score_Z ~
          TriadDensity_Z +
          (1 | Sample),
        data = df_lmm
      ),
      
      error = function(e)
        NULL
      
    )
    
    if(
      is.null(fit)
    ) next
    
    coef_tab <-
      summary(fit)$coefficients
    
    beta <- coef_tab[
      "TriadDensity_Z",
      "Estimate"
    ]
    
    se <- coef_tab[
      "TriadDensity_Z",
      "Std. Error"
    ]
    
    pval <- coef_tab[
      "TriadDensity_Z",
      "Pr(>|t|)"
    ]
    
    tumor_lmm_list[[
      paste(
        tp,
        mod_name,
        sep = "_"
      )
    ]] <- data.frame(
      
      Timepoint = tp,
      
      Module = mod_name,
      
      Group = "Tumor Cells",
      
      N_cells = nrow(
        df_lmm
      ),
      
      Beta = beta,
      SE = se,
      P = pval,
      
      CI_low =
        beta - 1.96*se,
      
      CI_high =
        beta + 1.96*se,
      
      stringsAsFactors = FALSE
      
    )
    
    for(sm in unique(
      df_lmm$Sample
    )) {
      
      tmp <- df_lmm %>%
        dplyr::filter(
          Sample == sm
        )
      
      if(
        nrow(tmp) < 30
      ) next
      
      fit_sm <- tryCatch(
        
        lm(
          Module_Score ~ TriadDensity_Z,
          data = tmp
        ),
        
        error = function(e) NULL
      )
      
      if(
        is.null(fit_sm)
      ) next
      
      cf <- summary(
        fit_sm
      )$coefficients
      
      if(
        !"TriadDensity_Z" %in%
        rownames(cf)
      ) {
        next
      }
      
      tumor_sample_effects[[
        paste(
          tp,
          mod_name,
          sm,
          sep = "_"
        )
      ]] <- data.frame(
        
        Timepoint = tp,
        Module = mod_name,
        Sample = sm,
        
        Effect =
          cf[
            "TriadDensity_Z",
            "Estimate"
          ],
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
  }
  
}

# Combine
tumor_stats <- bind_rows(
  tumor_lmm_list
)

sample_effects <- bind_rows(
  tumor_sample_effects
)

# Consistency metrics
consistency_summary <-
  sample_effects %>%
  group_by(
    Timepoint,
    Module
  ) %>%
  summarise(
    
    N_Pos = sum(
      Effect > 0,
      na.rm = TRUE
    ),
    
    N_Neg = sum(
      Effect < 0,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

tumor_stats <- tumor_stats %>%
  left_join(
    consistency_summary,
    by = c(
      "Timepoint",
      "Module"
    )
  )

# Sample-level consistency
tumor_stats$Consistency <- ifelse(
  
  tumor_stats$N_Pos >= tumor_stats$N_Neg,
  
  paste0(
    tumor_stats$N_Pos,
    "/",
    tumor_stats$N_Pos +
      tumor_stats$N_Neg,
    "+"
  ),
  
  paste0(
    tumor_stats$N_Neg,
    "/",
    tumor_stats$N_Pos +
      tumor_stats$N_Neg,
    "-"
  )
)

tumor_stats <- tumor_stats %>%
  dplyr::group_by(
    Timepoint
  ) %>%
  dplyr::mutate(
    FDR = p.adjust(
      P,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

tumor_stats$Sig <- case_when(
  
  tumor_stats$FDR < 0.001 ~ "***",
  tumor_stats$FDR < 0.01  ~ "**",
  tumor_stats$FDR < 0.05  ~ "*",
  
  TRUE ~ ""
  
)

# Percentage positive cell
pct_df_tumor <- lapply(
  names(tumor_plot_data),
  function(nm) {
    
    df_plot <- tumor_plot_data[[nm]]
    
    data.frame(
      
      Timepoint = unique(df_plot$Timepoint),
      
      Module = unique(df_plot$Module),
      
      Percent_Positive =
        mean(
          df_plot$Module_Score > 0,
          na.rm = TRUE
        ) * 100,
      
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

# Plot order
module_order <- c(
  
  "Antigen processing",
  
  "Interferon response",
  
  "Hypoxia/angiogenesis"
  
)

tumor_stats$Module <- factor(
  tumor_stats$Module,
  levels = rev(
    module_order
  )
)

tumor_stats$Group <-
  "Tumor Cells"


# Annotation
tumor_stats$Annotation <- ifelse(
  tumor_stats$Sig == "",
  tumor_stats$Consistency,
  paste0(
    tumor_stats$Sig,
    "  ",
    tumor_stats$Consistency
  )
)

# Sample-level validation
tumor_validation <- list()

for(tp in c("Early","Terminal")) {
  
  for(mod_name in names(TUMOR_MODULES)) {
    
    tmp <- sample_effects %>%
      filter(
        Timepoint == tp,
        Module == mod_name
      )
  
  if (nrow(tmp) == 0)
    next
  
  wilcox_p <- tryCatch(
    
    wilcox.test(
      tmp$Effect,
      mu = 0,
      exact = FALSE
    )$p.value,
    
    error = function(e) NA_real_
    
  )
  
  tumor_validation[[paste(
    tp,
    mod_name,
    sep = "_"
  )]] <- data.frame(
    
    Timepoint = tp,
    
    Module = mod_name,
    
    N_Samples = nrow(tmp),
    
    Positive = sum(
      tmp$Effect > 0,
      na.rm = TRUE
    ),
    
    Negative = sum(
      tmp$Effect < 0,
      na.rm = TRUE
    ),
    
    Wilcox_P = wilcox_p,
    
    stringsAsFactors = FALSE
  )
  }
}

tumor_validation_df <- bind_rows(
  tumor_validation
)

tumor_validation_df$Wilcox_FDR <- p.adjust(
  tumor_validation_df$Wilcox_P,
  method = "BH"
)

# Save output
tumor_final <- tumor_stats %>%
  
  left_join(
    tumor_validation_df,
    by = c(
      "Timepoint",
      "Module"
    )
  ) %>%
  
  select(
    Timepoint,
    Module,
    
    Beta,
    CI_low,
    CI_high,
    
    P,
    FDR,
    Sig,
    
    Consistency,
    
    N_Samples,
    Positive,
    Negative,
    
    Wilcox_P,
    Wilcox_FDR
  )

write.csv(
  tumor_final,
  "../../output/Fig6_Tumor_Program_Statistics.csv",
  row.names = FALSE
)

# Module summary dot plot
tumor_dot <- tumor_final %>%
  select(
    Timepoint,
    Module,
    Beta,
    FDR,
    Sig
  ) %>%
  left_join(
    pct_df_tumor,
    by = c(
      "Timepoint",
      "Module"
    )
  )

tumor_dot$Group <- "Tumor Cells"

module_order_tumor <- c(
  
  "Antigen processing",
  "Interferon response",
  "Hypoxia/angiogenesis"
)

tumor_dot$Module <- factor(
  tumor_dot$Module,
  levels = rev(module_order_tumor)
)

beta_lim_tumor <- max(
  abs(tumor_dot$Beta),
  na.rm = TRUE
)

p_tumor <- ggplot(
  tumor_dot,
  aes(
    x = 1,
    y = Module
  )
) +
  
  geom_point(
    
    aes(
      fill = Beta
    ),
    size = 5,
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +
  
  geom_text(
    
    aes(
      x = 1.12,
      label = Sig
    ),
    
    size = pt_size / 2
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -beta_lim_tumor,
      beta_lim_tumor
    ),
    name = "β"
  ) +
  
  scale_x_continuous(
    
    limits = c(
      0.85,
      1.15
    ),
    
    breaks = NULL
  ) +
  
  facet_grid(
    Group ~ Timepoint,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(1.5, "cm"),
      barwidth = unit(0.35, "cm"),
      order = 2
    )
  ) +
  
  theme_fig() +
  
  theme(
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),
    
    axis.text.x = element_blank(),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      size = pt_size * 0.7
    ),
    
    legend.text = element_text(
      size = pt_size * 0.7
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

save_pdf(
  "../../plots/Fig6_Triad_Tumor_Modules.pdf",
  p_tumor,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Cell-type marker dot plot
# ============================================================

markers <- list(
  "Other" = c(
    "Prom1"
  ),
  "Myeloid" = c(
    "Cd68", "Itgam", "Csf1r", "Trem2"
  ),
  "T/NK Cells" = c(
    "Cd3e", "Cd8b", "Nkg7", "Prf1"
  ),
  "B/Plasma Cells" = c(
    "Ms4a1", "Cd79a", "Mzb1", "Jchain"
  ),
  "Vasculature/Stroma" = c(
    "Flt1", "Rgs5"
  ),
  "Neurons" = c(
    "Gria1", "Gria2", "Eno2"
  ),
  "MES-like/Hypoxic Tumor" = c(
    "Cd44", "Vegfa", "Ndrg1"
  ),
  "Astrocytes/AC-like" = c(
    "Aqp4", "Slc1a3"
  ),
  "Oligodendrocytes/OPC-like" = c(
    "Pdgfra", "Olig1", "Olig2"
  )
)

cell_order <- c(
  "Other",
  "Myeloid",
  "T/NK Cells",
  "B/Plasma Cells",
  "Vasculature/Stroma",
  "Neurons",
  "MES-like/Hypoxic Tumor",
  "Astrocytes/AC-like",
  "Oligodendrocytes/OPC-like"
)

gene_order <- unique(unlist(markers))

# Extract expression
cells_use <- rownames(xe_hamster@meta.data)[
  xe_hamster$final_harmonized_label %in% cell_order
]

expr_df <- FetchData(
  xe_hamster,
  vars = gene_order,
  cells = cells_use
)

expr_df$final_harmonized_label <- xe_hamster$final_harmonized_label[
  rownames(expr_df)
]

# Per-cell-type summaries
dp <- bind_rows(
  
  lapply(
    cell_order,
    function(ct) {
      
      df_ct <- expr_df[
        expr_df$final_harmonized_label == ct,
        ,
        drop = FALSE
      ]
      
      bind_rows(
        
        lapply(
          gene_order,
          function(g) {
            
            data.frame(
              final_harmonized_label = ct,
              gene = g,
              avg = mean(df_ct[[g]], na.rm = TRUE),
              pct = mean(df_ct[[g]] > 0, na.rm = TRUE) * 100
            )
            
          }
        )
        
      )
      
    }
  )
  
)

# Scale average expression within marker
dp <- dp |>
  group_by(gene) |>
  mutate(
    avg_scaled = as.numeric(scale(avg))
  ) |>
  ungroup()

dp$gene <- factor(
  dp$gene,
  levels = gene_order
)

dp$final_harmonized_label <- factor(
  dp$final_harmonized_label,
  levels = cell_order
)

# Plot
p_dot <- ggplot(
  dp,
  aes(gene, final_harmonized_label)
) +
  geom_point(
    aes(
      size = pct,
      color = avg_scaled
    )
  ) +
  scale_size(
    range  = c(0.4, 2.2),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name   = "% cells"
  ) +
  scale_color_gradientn(
    colors = c(
      "grey80",
      "#4393C3",
      "#08306B"
    ),
    limits = c(-1.5, 2.5),
    oob = scales::squish,
    name = "Scaled\nExpr"
  ) +
  guides(
    size = guide_legend(
      title.position = "top"
    ),
    colour = guide_colorbar(
      title.position = "top",
      barheight = unit(1.2, "cm"),
      barwidth  = unit(0.25, "cm")
    )
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text  = element_text(size = 4.5),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm")
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_pdf(
  "../../plots/FigS5_DotPlot_CellTypeMarkers.pdf",
  p_dot,
  2 * P_Width,
  P_Height
)

# ============================================================
# Representative virus-near triads and virus-distant triads
# ============================================================

ham <- meta %>%
  filter(
    Timepoint == "Early",
    Treatment == "Virus",
    Sample == "RGD 606"
  )

# Helper
dist_to_set <- function(query_df, ref_df){
  
  if(nrow(query_df) == 0 || nrow(ref_df) == 0){
    return(rep(Inf, nrow(query_df)))
  }
  
  RANN::nn2(
    data  = as.matrix(ref_df[, c("x","y")]),
    query = as.matrix(query_df[, c("x","y")]),
    k = 1
  )$nn.dists[,1]
  
}

# Virus coordinates
virus_xy <- ham %>%
  filter(Virus_Pos) %>%
  select(x, y)

# Triad coordinates
tri_xy <- ham %>%
  filter(cell %in% triad_anchors) %>%
  select(cell, x, y) %>%
  rename(anchor = cell)

# Distance to nearest virus+ cell
tri_xy$DistTriadToVirus <- dist_to_set(
  tri_xy,
  virus_xy
)

# Near / Mid / Far
tri_xy <- tri_xy %>%
  mutate(
    Zone = case_when(
      DistTriadToVirus <= NEAR_MAX ~ "Virus-near",
      DistTriadToVirus >= FAR_MIN  ~ "Virus-distant",
      TRUE                         ~ "Mid"
    )
  )

tri_xy$Zone <- factor(
  tri_xy$Zone,
  levels = c("Virus-near","Mid","Virus-distant")
)

print(table(tri_xy$Zone))

# Manual region selection
cx <- 950   # manual x coordinate
cy <- 14000   # manual y coordinate

WINDOW <- 500

xmin <- cx - WINDOW
xmax <- cx + WINDOW

ymin <- cy - WINDOW
ymax <- cy + WINDOW


# Data in plotting window
plot_cells <- ham %>%
  filter(
    x >= xmin,
    x <= xmax,
    y >= ymin,
    y <= ymax
  )

plot_virus <- plot_cells %>%
  filter(Virus_Pos)

plot_triads <- tri_xy %>%
  filter(
    x >= xmin,
    x <= xmax,
    y >= ymin,
    y <= ymax,
    Zone %in% c("Virus-near","Virus-distant")
  )

# Scale bar
SCALE_LEN <- 200

x0 <- xmax - SCALE_LEN - 50
x1 <- x0 + SCALE_LEN

y0 <- ymin + 50

# Combined dataframe for plotting + legend
plot_highlight <- bind_rows(
  
  plot_virus %>%
    mutate(Type = "Virus+"),
  
  plot_triads %>%
    filter(Zone == "Virus-near") %>%
    mutate(Type = "Near triad"),
  
  plot_triads %>%
    filter(Zone == "Virus-distant") %>%
    mutate(Type = "Distant triad")
  
)

plot_highlight$Type <- factor(
  plot_highlight$Type,
  levels = c(
    "Near triad",
    "Distant triad",
    "Virus+"
  )
)

# Plot
p_triad_zoom <- ggplot() +
  
  geom_point_rast(
    data = plot_cells,
    aes(x, y),
    color = "grey88",
    size = 0.15,
    alpha = 0.8,
    raster.dpi = 600
  ) +
  
  geom_point(
    data = plot_highlight,
    aes(
      x = x,
      y = y,
      color = Type,
      shape = Type
    ),
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "Near triad"    = "blue",
      "Distant triad" = "black",
      "Virus+"        = "red3"
    ),
    name = NULL
  ) +
  
  scale_shape_manual(
    values = c(
      "Near triad"    = 17,
      "Distant triad" = 17,
      "Virus+"        = 16
    ),
    name = NULL
  ) +
  
  ggplot2::annotate(
    "segment",
    x = x0,
    xend = x1,
    y = y0,
    yend = y0,
    linewidth = 0.8
  ) +
  
  coord_fixed(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  
  labs(
    title = "Early: RGD 606"
  ) + 
  
  theme_fig() +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.margin = margin(
      t = -1,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.box.margin = margin(
      t = -1,
      r = 0,
      b = 0,
      l = 0
    ),
    
    legend.title = element_blank(),
    
    legend.key.width  = unit(0.35, "cm"),
    legend.key.height = unit(0.25, "cm"),
    
    plot.margin = margin(0, 0, 0, 0)
  ) +
  
  guides(
    shape = guide_legend(
      nrow = 1,
      override.aes = list(
        size = 2,
        alpha = 1
      )
    )
  )

save_pdf(
  "../../plots/Fig6_Hamster_Spatial_Virus_Near_Distant_Triads.pdf",
  p_triad_zoom,
  P_Width,
  P_Height
)

# ============================================================
# Virus-near vs virus-distant triad programs
# ============================================================

module_info <- data.frame(
  Module = c(
    "Inflammatory chemokine",
    "Antigen presentation",
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T",
    "Lymphoid organization",
    "Regulatory B"
  ),
  Group = c(
    "Myeloid",
    "Myeloid",
    "T/NK Cells",
    "T/NK Cells",
    "T/NK Cells",
    "T/NK Cells",
    "B/Plasma Cells",
    "B/Plasma Cells"
  ),
  stringsAsFactors = FALSE
)

# Helper
mean_over <- function(ids, vec){
  
  ids <- intersect(
    ids,
    names(vec)
  )
  
  if(!length(ids))
    return(NA_real_)
  
  mean(
    vec[ids],
    na.rm = TRUE
  )
  
}

# Virus-positive cells
virus_xy <- meta %>%
  dplyr::filter(
    Virus_Pos
  ) %>%
  dplyr::select(
    FOV,
    x,
    y
  )

# Triad anchors
tri_nf <- triad_xy %>%
  dplyr::rename(
    anchor = cell
  ) %>%
  dplyr::select(
    -any_of(c("Sample","Timepoint"))
  ) %>%
  dplyr::left_join(
    meta %>%
      dplyr::select(
        cell,
        Sample,
        Timepoint
      ) %>%
      dplyr::distinct(),
    by = c("anchor" = "cell")
  )

colnames(df_lmm)

"Sample" %in% colnames(df_lmm)
"FOV" %in% colnames(df_lmm)

tri_nf$DistTriadToVirus <- Inf

for(fv in unique(tri_nf$FOV)) {
  
  tri_idx <- which(
    tri_nf$FOV == fv
  )
  
  virus_sub <- virus_xy %>%
    dplyr::filter(
      FOV == fv
    )
  
  tri_nf$DistTriadToVirus[
    tri_idx
  ] <- dist_to_set(
    
    tri_nf[
      tri_idx,
      c("x","y")
    ],
    
    virus_sub
    
  )
  
}

tri_nf <- tri_nf %>%
  dplyr::mutate(
    Zone = dplyr::case_when(
      DistTriadToVirus <= NEAR_MAX ~ "Virus-near",
      DistTriadToVirus >= FAR_MIN  ~ "Virus-distant",
      TRUE                         ~ "Mid"
    )
  ) %>%
  dplyr::filter(
    Timepoint == "Early",
    Zone %in% c("Virus-near", "Virus-distant")
  )

tri_nf$Zone <- factor(
  tri_nf$Zone,
  levels = c(
    "Virus-distant",
    "Virus-near"
  )
)

message(
  "\n=== Early triad counts ==="
)

print(
  table(
    tri_nf$Zone
  )
)

# T/B members within each triad
T_all <- meta %>%
  dplyr::filter(
    grepl(
      T_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    cell,
    x,
    y
  )

B_all <- meta %>%
  dplyr::filter(
    grepl(
      B_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::select(
    cell,
    x,
    y
  )

tri_nf$T_members <- cells_within_radius(
  tri_nf,
  T_all,
  radius = r_tri
)

tri_nf$B_members <- cells_within_radius(
  tri_nf,
  B_all,
  radius = r_tri
)

# Module analysis
virus_list <- list()
sample_effects_list <- list()

for(i in seq_len(
  nrow(module_info)
)) {
  
  mod_name <- module_info$Module[i]
  group_name <- module_info$Group[i]
  
  genes <- MODULES[[mod_name]]
  
  if(length(genes) < 2)
    next
  
  # Myeloid
  if(group_name == "Myeloid") {
    
    score_vec <- module_score_vec(
      xe_hamster,
      tri_nf$anchor,
      genes
    )
    
    score <- score_vec[
      tri_nf$anchor
    ]
    
    # T/NK
  } else if(group_name == "T/NK Cells") {
    
    ids <- unique(
      unlist(
        tri_nf$T_members
      )
    )
    
    if(length(ids) < 20)
      next
    
    score_vec <- module_score_vec(
      xe_hamster,
      ids,
      genes
    )
    
    score <- vapply(
      tri_nf$T_members,
      mean_over,
      numeric(1),
      vec = score_vec
    )
    
    # B/Plasma
  } else {
    
    ids <- unique(
      unlist(
        tri_nf$B_members
      )
    )
    
    if(length(ids) < 20)
      next
    
    score_vec <- module_score_vec(
      xe_hamster,
      ids,
      genes
    )
    
    score <- vapply(
      tri_nf$B_members,
      mean_over,
      numeric(1),
      vec = score_vec
    )
    
  }
  
  df_lmm <- tri_nf %>%
    dplyr::mutate(
      Score = score
    ) %>%
    dplyr::filter(
      is.finite(Score)
    )
  
  if(nrow(df_lmm) < 20)
    next
  
  df_lmm$Score_Z <- as.numeric(
    scale(
      df_lmm$Score
    )
  )
  
  fit <- tryCatch(
    
    lmer(
      Score_Z ~
        Zone +
        (1 | Sample/FOV),
      data = df_lmm
    ),
    
    error = function(e)
      NULL
    
  )
  
  if(is.null(fit))
    next
  
  coef_tab <- summary(
    fit
  )$coefficients
  
  rn <- grep(
    "^Zone",
    rownames(coef_tab),
    value = TRUE
  )[1]
  
  if(is.na(rn))
    next
  
  beta <- coef_tab[
    rn,
    "Estimate"
  ]
  
  se <- coef_tab[
    rn,
    "Std. Error"
  ]
  
  pval <- coef_tab[
    rn,
    "Pr(>|t|)"
  ]
  
  virus_list[[mod_name]] <- data.frame(
    
    Module = mod_name,
    Group = group_name,
    
    N_Triads = nrow(df_lmm),
    
    Beta = beta,
    SE = se,
    P = pval,
    
    CI_low =
      beta - 1.96 * se,
    
    CI_high =
      beta + 1.96 * se,
    
    stringsAsFactors = FALSE
    
  )
  
  # Sample consistency
    for(sm in unique(
    df_lmm$Sample
  )) {
    
    tmp <- df_lmm %>%
      dplyr::filter(
        Sample == sm
      )
    
    if(nrow(tmp) < 20)
      next
    
    if(length(unique(tmp$Zone)) < 2)
      next
    
    fit_sm <- lm(
      Score_Z ~ Zone,
      data = tmp
    )
    
    cf <- summary(
      fit_sm
    )$coefficients
    
    rn_sm <- grep(
      "^Zone",
      rownames(cf),
      value = TRUE
    )[1]
    
    sample_effects_list[[paste(
      mod_name,
      sm,
      sep = "_"
    )]] <- data.frame(
      
      Module = mod_name,
      Sample = sm,
      
      Beta =
        cf[
          rn_sm,
          "Estimate"
        ],
      
      stringsAsFactors = FALSE
      
    )
    
  }
  
}

# Combine results
virus_stats <- bind_rows(
  virus_list
)

sample_effects <- bind_rows(
  sample_effects_list
)

# Consistency metrics
consistency_summary <- sample_effects %>%
  dplyr::group_by(
    Module
  ) %>%
  dplyr::summarise(
    
    N_Pos =
      sum(Beta > 0),
    
    N_Neg =
      sum(Beta < 0),
    
    Prop_Pos =
      mean(Beta > 0),
    
    Median_Beta =
      median(
        Beta,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )

virus_stats <- virus_stats %>%
  dplyr::left_join(
    consistency_summary,
    by = "Module"
  )

# FDR
virus_stats$FDR <- p.adjust(
  virus_stats$P,
  method = "BH"
)

virus_stats$Stars <- dplyr::case_when(
  virus_stats$FDR < 0.001 ~ "***",
  virus_stats$FDR < 0.01  ~ "**",
  virus_stats$FDR < 0.05  ~ "*",
  TRUE ~ ""
)

# Ordering
module_order <- c(
  
  "Inflammatory chemokine",
  "Antigen presentation",
  
  "Stemness",
  "Cytotoxicity",
  "Exhaustion",
  "Regulatory T",
  "Chemokine receptor",
  
  "Lymphoid organization",
  "Regulatory B"
  
)

virus_stats$Module <- factor(
  virus_stats$Module,
  levels = rev(module_order)
)

sample_effects$Module <- factor(
  sample_effects$Module,
  levels = rev(module_order)
)

virus_stats$Group <- factor(
  virus_stats$Group,
  levels = c(
    "Myeloid",
    "T/NK Cells",
    "B/Plasma Cells"
  )
)

write.csv(
  virus_stats,
  "../../output/Fig6_Virus_Near_Triad_Programs_Early_Only.csv",
  row.names = FALSE
)

# Plot
module_dot <- virus_stats %>%
  
  dplyr::select(
    Module,
    Group,
    Beta,
    FDR,
    Stars
  )

module_dot$Group <- factor(
  module_dot$Group,
  levels = c(
    "Myeloid",
    "T/NK Cells",
    "B/Plasma Cells"
  )
)

module_dot$Module <- factor(
  module_dot$Module,
  levels = rev(module_order)
)


beta_lim <- max(
  abs(module_dot$Beta),
  na.rm = TRUE
)

p_triad_modules <- ggplot(
  module_dot,
  aes(
    x = 1,
    y = Module
  )
) +
  
  geom_point(
    aes(
      fill = Beta
    ),
    size = 5,
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +
  
  geom_text(
    aes(
      x = 1.15,
      label = Stars
    ),
    size = pt_size / 2
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -beta_lim,
      beta_lim
    ),
    name = expression(beta)
  ) +
  
  scale_x_continuous(
    limits = c(
      0.85,
      1.18
    ),
    breaks = NULL
  ) +
  
  facet_grid(
    Group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(1, "cm"),
      barwidth = unit(0.3, "cm")
    )
  ) +
  
  theme_fig() +
  
  theme(
    
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),

    axis.text.x = element_blank(),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    legend.box = "vertical",
    
    legend.title = element_text(
      size = pt_size * 0.7
    ),
    
    legend.text = element_text(
      size = pt_size * 0.7
    )
  )

save_pdf(
  "../../plots/Fig6_Hamster_Triad_Module_Virus_Distance.pdf",
  p_triad_modules,
  1.2 * P_Width,
  P_Height
)

# ============================================================
# Triad enrichment near virus
# (% myeloid cells that are triad anchors)
# ============================================================

# Early myeloid cells
myeloid_zone_df <- meta %>%
  
  dplyr::filter(
    Timepoint == "Early"
  ) %>%
  
  dplyr::filter(
    grepl(
      MY_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  
  dplyr::mutate(
    DistToVirus = dist_to_set(
      .,
      virus_xy
    ),
    Virus_Proximity = dplyr::case_when(
      DistToVirus <= NEAR_MAX
      ~ "Virus-near",
      DistToVirus >= FAR_MIN
      ~ "Virus-distant",
      TRUE
      ~ NA_character_
    ),
    
    Triad = cell %in% triad_anchors
  ) %>%
  dplyr::filter(
    !is.na(Virus_Proximity)
  )

# Mixed-effects logistic regression
myeloid_zone_df$Virus_Proximity <- factor(
  myeloid_zone_df$Virus_Proximity,
  levels = c(
    "Virus-distant",
    "Virus-near"
  )
)

fit_glmm <- glmer(
  Triad ~ Virus_Proximity +
    (1 | Sample/FOV),
  data = myeloid_zone_df,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa"
  )
)

coef_tab <- summary(fit_glmm)$coefficients

beta <- coef_tab[
  "Virus_ProximityVirus-near",
  "Estimate"
]

se <- coef_tab[
  "Virus_ProximityVirus-near",
  "Std. Error"
]

pval <- coef_tab[
  "Virus_ProximityVirus-near",
  "Pr(>|z|)"
]

or_result <- tibble(
  Comparison = "Virus-near",
  OR      = exp(beta),
  CI_Low  = exp(beta - 1.96 * se),
  CI_High = exp(beta + 1.96 * se),
  P       = pval
) %>%
  mutate(
    Sig = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

print(or_result)

write.csv(
  or_result,
  "../../output/Fig6_Early_Triad_Enrichment_OR.csv",
  row.names = FALSE
)

# Prevalence
triad_prevalence <- myeloid_zone_df %>%
  dplyr::group_by(
    Virus_Proximity
  ) %>%
  dplyr::summarise(
    Fraction_Triad =
      mean(Triad),
    Triad_Cells =
      sum(Triad),
    Total_Cells =
      n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Virus_Proximity = factor(
      Virus_Proximity,
      levels = c(
        "Virus-distant",
        "Virus-near"
      )
    )
  )

print(
  triad_prevalence
)

write.csv(
  triad_prevalence,
  "../../output/Fig6_Early_Triad_Prevalence.csv",
  row.names = FALSE
)

# Plot
p_triad_enrichment <- ggplot(
  triad_prevalence,
  aes(
    x = Virus_Proximity,
    y = Fraction_Triad,
    fill = Virus_Proximity
  )
) +
  
  geom_col(
    width = 0.65,
    colour = "black",
    linewidth = 0.3
  ) +
  
  geom_text(
    aes(
      label =
        scales::percent(
          Fraction_Triad,
          accuracy = 0.1
        )
    ),
    vjust = -0.4,
    size = pt_size / 2.2
  ) +
  
  scale_fill_manual(
    values = c(
      "Virus-distant" = "black",
      "Virus-near" = "blue"
    )
  ) +
  
  ggplot2::annotate(
    "text",
    x = 1.5,
    y =
      max(
        triad_prevalence$Fraction_Triad
      ) * 1.4,
    
    label = paste0(
      "OR = ",
      round(or_result$OR, 2),
      "\n95% CI ",
      round(or_result$CI_Low, 2),
      "–",
      round(or_result$CI_High, 2)
    ),
    
    size = pt_size / 2.5,
    lineheight = 0.9
  ) +
  
  ggplot2::annotate(
    "text",
    x = 1.5,
    y =
      max(
        triad_prevalence$Fraction_Triad
      ) * 1.1,
    label = or_result$Sig,
    size =
      pt_size / 1.4,
  ) +
  
  facet_wrap(
    ~ factor("Early"),
    nrow = 1
  ) +
  
  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    expand = expansion(
      mult = c(0, 0.35)
    )
  ) +
  
  labs(
    x = NULL,
    y = "% of myeloid cells participating in triads"
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "none",
    strip.text = element_text(
      size = pt_size + 1
    )
  )

save_pdf(
  "../../plots/Fig6_Early_Triad_Enrichment_Near_Virus.pdf",
  p_triad_enrichment,
  0.8 * P_Width,
  P_Height
)
