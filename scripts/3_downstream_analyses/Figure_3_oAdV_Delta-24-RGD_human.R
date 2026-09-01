# ==============================================================================
# Figures 3 and S3 analyses using Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts/3_downstream_analyses
#
# Required input files: 
# ../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds
# ../../NCT00805376_rGBM_oAdV/clinical_data.csv
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
  
})

set.seed(24)
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
  "B/Plasma Cells"            = "#CC79A7"
)

# ============================================================
# Parameters
# ============================================================

MY_R   <- "Myeloid"
B_R    <- "B/Plasma"
T_R    <- "T/NK"
VASC_R <- "Vasculature/Stroma"

r_tri  <- 30


TUMOR_LABELS <- c(
  "Oligodendrocytes/OPC-like",
  "Astrocytes/AC-like",
  "MES-like/Hypoxic Tumor"
)

# ============================================================
# Data loading
# ============================================================

xe_global <- readRDS(
  "../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds"
)

xe_global <- UpdateSeuratObject(xe_global)

DefaultAssay(xe_global) <- "Xenium"

# ============================================================
# Build metadata with coordinates
# ============================================================

raw_coords <- do.call(
  rbind,
  lapply(Images(xe_global), function(fov) {
    
    GetTissueCoordinates(xe_global[[fov]]) |>
      dplyr::mutate(FOV = fov)
    
  })
)

meta <- xe_global@meta.data |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(raw_coords, by = "cell")

# Patient/sample identifiers
meta$sample_id <- stringr::str_extract(
  meta$orig.ident,
  "P[0-9]+-(post|pre)"
)

meta$Patient_Match <- stringr::str_extract(
  meta$sample_id,
  "P[0-9]+"
)

# ============================================================
# Define triad anchors
# ============================================================
# A myeloid cell is a triad anchor if:
# nearest one (k=1) B/Plasma distance <= 30 um
# AND
# nearest one (k=1) T/NK distance     <= 30 um

triad_flags_myeloid <- function(df) {
  
  A <- subset(
    df,
    grepl(MY_R, cell_type, ignore.case = TRUE)
  )
  
  B <- subset(
    df,
    grepl(B_R, cell_type, ignore.case = TRUE)
  )
  
  T <- subset(
    df,
    grepl(T_R, cell_type, ignore.case = TRUE)
  )
  
  if (nrow(A) == 0 ||
      nrow(B) == 0 ||
      nrow(T) == 0) {
    
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
    Triad = (dB <= r_tri) & (dT <= r_tri)
  )
}

# ============================================================
# Identify triad anchors
# ============================================================

triads <- meta |>
  dplyr::group_by(sample_id, FOV) |>
  dplyr::group_modify(
    ~triad_flags_myeloid(.x)
  ) |>
  dplyr::ungroup()

# Store triad calls in metadata
meta <- meta |>
  dplyr::left_join(
    dplyr::select(triads, cell, Triad),
    by = "cell"
  )

meta$Triad[is.na(meta$Triad)] <- FALSE

# Triad anchor cell IDs
triad_anchors <- triads |>
  dplyr::filter(Triad) |>
  dplyr::pull(cell)


N_TRIADS <- length(triad_anchors)

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
    sample_id,
    FOV,
    x,
    y
  )

# ============================================================
# Label triad niche
# ============================================================
# A target cell is considered "In triad" if its nearest
# triad-anchor myeloid cell is within 30 µm.

label_specific_niche <- function(df,
                                 anchors,
                                 target_regex) {
  
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
  
  if (nrow(anc) == 0 ||
      nrow(tgt) == 0) {
    
    return(
      data.frame(
        cell = character(),
        Specific_Niche = character()
      )
    )
  }
  
  d <- RANN::nn2(
    anc[, c("x", "y")],
    tgt[, c("x", "y")],
    k = 1
  )$nn.dists[, 1]
  
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
  dplyr::group_by(sample_id, FOV) |>
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
  dplyr::group_by(sample_id, FOV) |>
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
  dplyr::group_by(sample_id, FOV) |>
  dplyr::group_modify(
    ~label_specific_niche(
      .x,
      anchors = triad_anchors,
      target_regex = B_R
    )
  ) |>
  dplyr::ungroup()

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

# Distance function
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

# ============================================================
# Distance of myeloid cells to vasculature/stroma
# ============================================================

df_my_dist <- meta |>
  dplyr::filter(
    grepl(MY_R, cell_type, ignore.case = TRUE)
  ) |>
  dplyr::group_by(sample_id, FOV) |>
  dplyr::group_modify(~{
    
    vas <- meta |>
      dplyr::filter(
        sample_id == .y$sample_id,
        FOV == .y$FOV,
        grepl(VASC_R, cell_type, ignore.case = TRUE)
      ) |>
      dplyr::select(x, y)
    
    if (nrow(vas) == 0)
      return(tibble())
    
    d <- RANN::nn2(
      vas,
      .x[, c("x","y")],
      k = 1
    )$nn.dists[,1]
    
    tibble(
      cell = .x$cell,
      Dist = d,
      Niche = ifelse(
        .x$Triad,
        "In triad",
        "Outside triad"
      )
    )
    
  }) |>
  dplyr::ungroup()

# Sample-level medians
stats_struct_long <- df_my_dist |>
  dplyr::group_by(sample_id, Niche) |>
  dplyr::summarise(
    Median_Dist = median(Dist),
    .groups = "drop"
  )

stats_struct <- stats_struct_long |>
  tidyr::pivot_wider(
    names_from = Niche,
    values_from = Median_Dist
  ) |>
  dplyr::filter(
    is.finite(`In triad`),
    is.finite(`Outside triad`)
  )

# Paired Wilcoxon test
p_paired <- wilcox.test(
  stats_struct$`In triad`,
  stats_struct$`Outside triad`,
  paired = TRUE
)$p.value

# Plot
p_struct <- ggplot(
  stats_struct_long,
  aes(x = Niche, y = Median_Dist)
) +
  geom_boxplot(
    aes(fill = Niche),
    outlier.shape = NA,
    width = 0.5,
    alpha = 0.85
  ) +
  geom_line(
    aes(group = sample_id),
    color = "grey50",
    alpha = 0.6
  ) +
  geom_point(
    aes(fill = Niche),
    shape = 21,
    size = 1.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Outside triad" = "grey70",
      "In triad" = "#D55E00"
    )
  ) +
  theme_fig() +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Median distance to vasculature/stroma (µm)"
  )

# Add p-value
y_annot <- max(
  stats_struct_long$Median_Dist,
  na.rm = TRUE
) * 1.05

p_struct <- p_struct +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = y_annot,
    label = format_p_plotmath(p_paired),
    parse = TRUE,
    size = pt_size / 2.2
  )

# Save plot
save_pdf(
  "../../plots/Fig3_Myeloid_Distance_to_Vasculature.pdf",
  p_struct,
  P_Width,
  P_Height
)

# ============================================================
# Triad anchor spatial plot (P15)
# ============================================================

anchors_all <- meta %>%
  dplyr::filter(Triad) %>%
  dplyr::select(
    cell,
    sample_id,
    FOV,
    x,
    y
    )

# Representative ROI (manually selected)
sid <- "P15-post"
fv  <- "fov"

cx <- 3300
cy <- 2300

ROI_WIN   <- 2000
SCALE_LEN <- 500

# Data
df_bg <- meta %>%
  dplyr::filter(
    sample_id == sid,
    FOV == fv
  )

df_v <- df_bg %>%
  dplyr::filter(
    grepl(
      VASC_R,
      cell_type,
      ignore.case = TRUE
    )
  )

df_a <- anchors_all %>%
  dplyr::filter(
    sample_id == sid,
    FOV == fv
  )

center <- tryCatch(
  pick_overlap_center(
    df_a,
    df_v[, c("x","y")]
  ),
  error = function(e)
    c(
      mean(df_a$x),
      mean(df_a$y)
    )
)

cx <- center[1]
cy <- center[2]

message(
  sprintf(
    "ROI center = (%.1f, %.1f)",
    cx,
    cy
  )
)

# Crop ROI
df_bg_roi <- df_bg %>%
  dplyr::filter(
    x >= cx - ROI_WIN,
    x <= cx + ROI_WIN,
    y >= cy - ROI_WIN,
    y <= cy + ROI_WIN
  )

df_v_roi <- df_v %>%
  dplyr::filter(
    x >= cx - ROI_WIN,
    x <= cx + ROI_WIN,
    y >= cy - ROI_WIN,
    y <= cy + ROI_WIN
  ) %>%
  dplyr::mutate(
    Layer = "Vasculature/Stroma"
  )

df_a_roi <- df_a %>%
  dplyr::filter(
    x >= cx - ROI_WIN,
    x <= cx + ROI_WIN,
    y >= cy - ROI_WIN,
    y <= cy + ROI_WIN
  ) %>%
  dplyr::mutate(
    Layer = "Triad anchors"
  )

df_v_roi$Layer <- factor(
  df_v_roi$Layer,
  levels = c(
    "Vasculature/Stroma",
    "Triad anchors"
  )
)

df_a_roi$Layer <- factor(
  df_a_roi$Layer,
  levels = c(
    "Vasculature/Stroma",
    "Triad anchors"
  )
)

# Scale bar
sb_x1 <- (cx + ROI_WIN) - 40 - SCALE_LEN
sb_x2 <- sb_x1 + SCALE_LEN
sb_y  <- (cy - ROI_WIN) + 40

# Plot
p_fig3_roi <- ggplot() +
  
  ggrastr::geom_point_rast(
    data = df_bg_roi,
    aes(x, y),
    colour = "grey93",
    size = 0.06,
    alpha = 0.30,
    raster.dpi = 600,
    show.legend = FALSE
  ) +
  
  ggrastr::geom_point_rast(
    data = df_v_roi,
    aes(
      x, y,
      colour = Layer,
      shape  = Layer
    ),
    size = 0.35,
    alpha = 0.60,
    stroke = 0,
    raster.dpi = 600
  ) +
  
  geom_point(
    data = df_a_roi,
    aes(
      x, y,
      colour = Layer,
      shape  = Layer
    ),
    size = 1.5,
    alpha = 0.95,
    stroke = 0.15
  ) +
  
  ggplot2::annotate(
    "segment",
    x = sb_x1,
    xend = sb_x2,
    y = sb_y,
    yend = sb_y,
    linewidth = 0.45,
    colour = "black"
  ) +
  
  coord_fixed(
    ratio = 1,
    xlim = c(
      cx - ROI_WIN,
      cx + ROI_WIN
    ),
    ylim = c(
      cy - ROI_WIN,
      cy + ROI_WIN
    ),
    clip = "on"
  ) +
  
  scale_shape_manual(
    values = c(
      "Vasculature/Stroma" = 16,
      "Triad anchors"      = 17
    ),
    name = NULL
  ) +
  
  scale_colour_manual(
    values = c(
      "Vasculature/Stroma" = "#0072B2",
      "Triad anchors"      = "#D55E00"
    ),
    name = NULL
  ) +
  
  labs(
    title = "P15"
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
  "../../plots/Fig3_Triad_Spatial_Map_P15.pdf",
  p_fig3_roi,
  P_Width,
  P_Height
)

# ============================================================
# Triad-anchor clustering versus CSR
# ============================================================

if (nrow(anchors_all) > 0) {
  
  B <- 200
  
  # Observed NN distances
  obs_list <- anchors_all %>%
    dplyr::group_by(sample_id, FOV) %>%
    dplyr::group_modify(~{
      
      if (nrow(.x) < 10) {
        return(
          tibble(
            NN   = numeric(0),
            kind = character(0)
          )
        )
      }
      
      tibble(
        NN = RANN::nn2(
          as.matrix(.x[, c("x", "y")]),
          k = 2
        )$nn.dists[, 2],
        kind = "Observed"
      )
      
    }) %>%
    dplyr::ungroup()
  
  # CSR simulations
  csr_list <- anchors_all %>%
    dplyr::group_by(sample_id, FOV) %>%
    dplyr::group_modify(~{
      
      n <- nrow(.x)
      
      if (n < 10) {
        return(
          tibble(
            NN   = numeric(0),
            kind = character(0)
          )
        )
      }
      
      cells_fov <- meta %>%
        dplyr::filter(
          sample_id == .y$sample_id[1],
          FOV == .y$FOV[1]
        )
      
      win <- tryCatch(
        spatstat.geom::convexhull.xy(
          x = cells_fov$x,
          y = cells_fov$y
        ),
        error = function(e)
          spatstat.geom::owin(
            xrange = range(cells_fov$x),
            yrange = range(cells_fov$y)
          )
      )
      
      sims_nn <- unlist(
        replicate(
          B,
          {
            sim_pts <- spatstat.random::runifpoint(n, win)
            
            RANN::nn2(
              cbind(sim_pts$x, sim_pts$y),
              k = 2
            )$nn.dists[, 2]
          },
          simplify = FALSE
        )
      )
      
      tibble(
        NN   = sims_nn,
        kind = "CSR"
      )
      
    }) %>%
    dplyr::ungroup()
  
  # Combine
  kde_df <- dplyr::bind_rows(
    obs_list,
    csr_list
  )
  
  if (nrow(kde_df) > 0) {
    
    X_MAX <- max(
      quantile(
        obs_list$NN,
        0.99,
        na.rm = TRUE
      ),
      quantile(
        csr_list$NN,
        0.99,
        na.rm = TRUE
      ),
      na.rm = TRUE
    )
    
    # Sample-level paired comparison
    med_obs <- obs_list %>%
      dplyr::group_by(sample_id, FOV) %>%
      dplyr::summarise(
        M = median(NN),
        .groups = "drop"
      )
    
    med_csr <- csr_list %>%
      dplyr::group_by(sample_id, FOV) %>%
      dplyr::summarise(
        M = median(NN),
        .groups = "drop"
      )
    
    med_samp <- dplyr::inner_join(
      med_obs,
      med_csr,
      by = c("sample_id", "FOV")
    ) %>%
      dplyr::group_by(sample_id) %>%
      dplyr::summarise(
        O = median(M.x),
        C = median(M.y),
        .groups = "drop"
      )
    
    p_wil <- tryCatch(
      wilcox.test(
        med_samp$O,
        med_samp$C,
        paired = TRUE
      )$p.value,
      error = function(e) NA_real_
    )
    
    # Global medians
    med_obs_global <- median(
      obs_list$NN,
      na.rm = TRUE
    )
    
    med_csr_global <- median(
      csr_list$NN,
      na.rm = TRUE
    )
    
    # Plot
    p_kde <- ggplot(
      kde_df,
      aes(
        x = NN,
        fill = kind,
        colour = kind
      )
    ) +
      
      geom_density(
        alpha = 0.35,
        linewidth = 0.30,
        adjust = 1
      ) +
      
      geom_vline(
        xintercept = med_obs_global,
        colour = pal_best["Myeloid"],
        linetype = 2,
        linewidth = 0.35
      ) +
      
      geom_vline(
        xintercept = med_csr_global,
        colour = "#6E6E6E",
        linetype = 2,
        linewidth = 0.35
      ) +
      
      coord_cartesian(
        xlim = c(0, X_MAX)
      ) +
      
      scale_fill_manual(
        values = c(
          "Observed" = unname(pal_best["Myeloid"]),
          "CSR"      = "#BDBDBD"
        ),
        breaks = c(
          "Observed",
          "CSR"
        ),
        name = NULL
      ) +
      
      scale_colour_manual(
        values = c(
          "Observed" = unname(pal_best["Myeloid"]),
          "CSR"      = "#6E6E6E"
        ),
        breaks = c(
          "Observed",
          "CSR"
        ),
        name = NULL
      ) +
      
      ggplot2::annotate(
        "text",
        x = X_MAX * 0.6,
        y = Inf,
        hjust = 1,
        vjust = 1.4,
        label = format_p_plotmath(p_wil),
        parse = TRUE,
        size = pt_size / 2.2
      ) +
      
      labs(
        x = "Nearest-neighbor distance (µm)",
        y = "Density"
      ) +
      
      theme_fig() +
      
      theme(
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank()
      )
    
    save_pdf(
      "../../plots/Fig3_TriadAnchor_NN_CSR.pdf",
      p_kde,
      P_Width,
      P_Height
    )
    
  }
}

# ============================================================
# Triad-associated immune programs
# ============================================================

# Niche labels
lab_my <- myeloid_niche
lab_t  <- tcell_niche
lab_b  <- bcell_niche

# Module definitions
MODULES_EXT <- list(
  
  # Myeloid
  "Inflammatory chemokine" = c("CXCL9","CXCL10","CXCL11","CCL5"),
  "Antigen presentation"   = c("CD74","CTSS","CD80","CD86"),
  "IFNγ response"          = c("STAT1","IRF1","IFNGR1"),
  
  # T/NK
  "Stemness"     = c("TCF7","IL7R","SELL"),
  "Cytotoxicity" = c("PRF1","NKG7","GNLY","GZMB"),
  "Exhaustion"   = c("PDCD1","LAG3","HAVCR2","TOX"),
  "Regulatory T" = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  
  # B/Plasma
  "Lymphoid organization" = c("CXCL13","LTB","CCR7","CCL19"),
  "Regulatory B"          = c("CD274","IL10","TGFB1")
)

# Keep only genes present in panel
MODULES_EXT <- lapply(
  MODULES_EXT,
  intersect,
  rownames(xe_global)
)

# Keep modules with >= 2 genes
MODULES_EXT <- MODULES_EXT[
  sapply(MODULES_EXT, length) >= 2
]

# Module metadata
module_info <- tibble::tribble(
  ~Module,                  ~Regex, ~Group,
  
  "Inflammatory chemokine", MY_R,   "Myeloid",
  "Antigen presentation",   MY_R,   "Myeloid",
  "IFNγ response",          MY_R,   "Myeloid",
  
  "Stemness",               T_R,    "T/NK Cells",
  "Cytotoxicity",           T_R,    "T/NK Cells",
  "Exhaustion",             T_R,    "T/NK Cells",
  "Regulatory T",           T_R,    "T/NK Cells",
  
  "Lymphoid organization",  B_R,    "B/Plasma Cells",
  "Regulatory B",           B_R,    "B/Plasma Cells"
)

# Mixed model analysis
plot_data_list <- list()
ext_lmm_list <- list()
sample_effects_list <- list()

for (mod_name in names(MODULES_EXT)) {
  
  message("Processing: ", mod_name)
  
  info <- module_info %>%
    dplyr::filter(Module == mod_name)
  
  genes <- MODULES_EXT[[mod_name]]
  
  if (length(genes) < 2)
    next
  
  labs_use <- switch(
    info$Group,
    
    "Myeloid"        = lab_my,
    "T/NK Cells"     = lab_t,
    "B/Plasma Cells" = lab_b
  )
  
  df_lmm <- meta %>%
    dplyr::filter(
      grepl(
        info$Regex,
        cell_type,
        ignore.case = TRUE
      )
    ) %>%
    dplyr::select(
      cell,
      sample_id,
      FOV
    ) %>%
    dplyr::inner_join(
      labs_use %>%
        dplyr::select(
          cell,
          Specific_Niche
        ),
      by = "cell"
    )
  
  if (nrow(df_lmm) < 100)
    next
  
  expr_mat <- FetchData(
    xe_global,
    vars = genes,
    cells = df_lmm$cell,
    layer = "data"
  )
  
  module_score <- rowMeans(
    expr_mat,
    na.rm = TRUE
  )
  
  if (sd(module_score, na.rm = TRUE) == 0)
    next
  
  df_lmm$Module_Score <- module_score
  
  df_lmm$Score_Z <- as.numeric(
    scale(module_score)
  )
  
  plot_data_list[[mod_name]] <- df_lmm
  
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
      df_lmm$Specific_Niche == "Outside triad"
    ],
    na.rm = TRUE
  )
  
  mean_inside <- mean(
    df_lmm$Module_Score[
      df_lmm$Specific_Niche == "In triad"
    ],
    na.rm = TRUE
  )
  
  pct_change <- ifelse(
    abs(mean_outside) > 1e-8,
    100 * (mean_inside - mean_outside) /
      abs(mean_outside),
    NA_real_
  )
  
  # LMM on standardized score
  fit <- tryCatch(
    
    lmer(
      Score_Z ~
        Specific_Niche +
        (1 | sample_id/FOV),
      data = df_lmm
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit))
    next
  
  coef_tab <- summary(fit)$coefficients
  
  coef_name <- rownames(coef_tab)[
    grepl(
      "Specific_Niche",
      rownames(coef_tab)
    )
  ][1]
  
  beta <- coef_tab[coef_name, "Estimate"]
  se   <- coef_tab[coef_name, "Std. Error"]
  pval <- coef_tab[coef_name, "Pr(>|t|)"]
  
  ci_low  <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  
  ext_lmm_list[[mod_name]] <- data.frame(
    
    Module = mod_name,
    Group = info$Group,
    
    Mean_Outside = mean_outside,
    Mean_InTriad = mean_inside,
    
    Percent_Change = pct_change,
    
    Beta = beta,
    P = pval,
    
    CI_low  = ci_low,
    CI_high = ci_high,
    
    stringsAsFactors = FALSE
  )
  
  # Sample-level direction consistency
  for (sm in unique(df_lmm$sample_id)) {
    
    tmp <- df_lmm %>%
      dplyr::filter(sample_id == sm)
    
    if (nrow(tmp) < 100)
      next
    
    if (length(unique(tmp$Specific_Niche)) < 2)
      next
    
    mean_out_sm <- mean(
      tmp$Module_Score[
        tmp$Specific_Niche == "Outside triad"
      ],
      na.rm = TRUE
    )
    
    mean_in_sm <- mean(
      tmp$Module_Score[
        tmp$Specific_Niche == "In triad"
      ],
      na.rm = TRUE
    )
    
    pct_sm <- ifelse(
      abs(mean_out_sm) > 1e-8,
      100 * (mean_in_sm - mean_out_sm) /
        abs(mean_out_sm),
      NA_real_
    )
    
    sample_effects_list[[paste0(mod_name, "_", sm)]] <-
      data.frame(
        Module = mod_name,
        sample_id = sm,
        Percent_Change = pct_sm,
        stringsAsFactors = FALSE
      )
  }
}

# Combine results
ext_stats <- bind_rows(
  ext_lmm_list
)

sample_effects <- bind_rows(
  sample_effects_list
)

# Patient-level consistency
consistency_summary <- sample_effects %>%
  group_by(Module) %>%
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
  left_join(
    consistency_summary,
    by = "Module"
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
ext_stats$FDR <- p.adjust(
  ext_stats$P,
  method = "BH"
)

ext_stats <- ext_stats %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Plot ordering
module_order <- c(
  
  "Inflammatory chemokine",
  "Antigen presentation",
  "IFNγ response",
  
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

ext_stats$Group <- factor(
  ext_stats$Group,
  levels = c(
    "Myeloid",
    "T/NK Cells",
    "B/Plasma Cells"
  )
)

# Annotation
ext_stats$Annotation <- ifelse(
  ext_stats$Sig == "",
  ext_stats$Consistency,
  paste0(
    ext_stats$Sig,
    "  ",
    ext_stats$Consistency
  )
)

# Sample-level validation
specimen_validation <- list()

for (mod_name in names(plot_data_list)) {
  
  df_plot <- plot_data_list[[mod_name]]
  
  sample_summary <- df_plot %>%
    group_by(
      sample_id,
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
      Delta = `In triad` - `Outside triad`
    )
  
  wilcox_p <- tryCatch(
    
    wilcox.test(
      sample_summary$Delta,
      mu = 0,
      exact = FALSE
    )$p.value,
    
    error = function(e) NA_real_
  )
  
  specimen_validation[[mod_name]] <- data.frame(
    
    Module = mod_name,
    
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

specimen_validation_df <- bind_rows(
  specimen_validation
)

specimen_validation_df$Wilcox_FDR <- p.adjust(
  specimen_validation_df$Wilcox_P,
  method = "BH"
)

# Save output
final_results <- ext_stats %>%
  left_join(
    specimen_validation_df,
    by = "Module"
  ) %>%
  dplyr::select(
    
    Module,
    Group,
    
    Mean_Outside,
    Mean_InTriad,
    Percent_Change,
    
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
  final_results,
  "../../output/Fig3_Triad_Module_Statistics.csv",
  row.names = FALSE
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
  function(mod_name) {
    
    df_plot <- plot_data_list[[mod_name]]
    
    data.frame(
      Module = mod_name,
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
    Module,
    Group,
    Beta,
    FDR,
    Sig,
    Consistency
  ) %>%
  left_join(
    pct_df,
    by = "Module"
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
    "IFNγ response",
    
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
      x = 1.15,
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
    ),
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 5,
      l = 5
    )
  )

save_pdf(
  "../../plots/Fig3_Triad_Immune_Modules.pdf",
  p_module_summary,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Triad module scores vs overall survival
# ============================================================

clinical_df <- readr::read_csv(
  "../../NCT00805376_rGBM_oAdV/clinical_data.csv",
  show_col_types = FALSE
)

os_df <- clinical_df |>
  transmute(
    Patient_Match = str_extract(Patient_ID, "P[0-9]+"),
    OS_months = as.numeric(OS_Months)
  ) |>
  filter(
    !is.na(Patient_Match),
    is.finite(OS_months)
  ) |>
  distinct(Patient_Match, .keep_all = TRUE)

modules_use <- c(
  "Inflammatory chemokine",
  "Stemness",
  "Lymphoid organization"
)

# Patient-level triad module scores
module_os_df <- bind_rows(
  
  lapply(
    modules_use,
    function(mod_name){
      
      df_plot <- plot_data_list[[mod_name]]
      
      df_plot %>%
        dplyr::filter(
          Specific_Niche == "In triad"
        ) %>%
        dplyr::group_by(sample_id) %>%
        dplyr::summarise(
          Module_Score = mean(
            Module_Score,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Module = mod_name,
          Patient_Match = stringr::str_extract(
            sample_id,
            "P[0-9]+"
          )
        )
    }
  )
  
)

module_os_df <- module_os_df %>%
  dplyr::left_join(
    os_df,
    by = "Patient_Match"
  ) %>%
  dplyr::filter(
    is.finite(OS_months)
  )

module_os_df %>%
  dplyr::select(
    sample_id,
    Patient_Match,
    Module,
    Module_Score,
    OS_months
  )

module_os_df <- module_os_df %>%
  dplyr::filter(
    Patient_Match != "P20"
  )

module_os_df <- module_os_df %>%
  filter(grepl("post", sample_id, ignore.case = TRUE))

os_results <- bind_rows(
  
  lapply(
    modules_use,
    function(mod_name){
      
      tmp <- module_os_df %>%
        dplyr::filter(
          Module == mod_name
        )
      
      ct <- cor.test(
        tmp$Module_Score,
        tmp$OS_months,
        method = "spearman",
        exact = FALSE
      )
      
      data.frame(
        Module = mod_name,
        N = nrow(tmp),
        Rho = unname(ct$estimate),
        P = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
  
)

os_results$FDR <- p.adjust(
  os_results$P,
  method = "BH"
)

print(os_results)

write.csv(
  os_results,
  "../../output/Fig3_Module_OS_Correlations.csv",
  row.names = FALSE
)

# ============================================================
# Myeloid subtype analysis for triads
# ============================================================

# Define myeloid subtype signatures
SIG_LIST <- list(
  
  Microglia = c(
    "P2RY12",
    "CX3CR1",
    "GPR34"
  ),
  
  Macrophage = c(
    "CD163",
    "MARCO",
    "VSIG4",
    "FCGR3A"
  ),
  
  Monocyte = c(
    "CCR2",
    "CD14",
    "VCAN",
    "FCGR1A",
    "S100A9"
  ),
  
  DC = c(
    "CD1C",
    "FCER1A",
    "CLEC10A",
    "CD83",
    "CCR7",
    "LAMP3"
  )
  
)

# Fetch expression
all_genes <- unique(unlist(SIG_LIST))

expr <- FetchData(
  xe_global,
  vars = all_genes,
  layer = "data"
)

expr$cell <- rownames(expr)

# Score signatures
score_mean <- function(df, genes){
  
  genes <- intersect(
    genes,
    colnames(df)
  )
  
  rowMeans(
    df[, genes, drop = FALSE]
  )
  
}

expr_sub <- expr

for(nm in names(SIG_LIST)){
  
  expr_sub[[paste0(nm, "_Score")]] <-
    score_mean(
      expr_sub,
      SIG_LIST[[nm]]
    )
  
}

score_cols <- c(
  "Microglia_Score",
  "Macrophage_Score",
  "Monocyte_Score",
  "DC_Score"
)

score_mat <- as.matrix(
  expr_sub[, score_cols]
)

expr_sub$MaxScore <-
  apply(score_mat, 1, max)

expr_sub$Myeloid_Subtype <-
  gsub(
    "_Score",
    "",
    score_cols[
      max.col(score_mat)
    ]
  )

expr_sub$Myeloid_Subtype[
  expr_sub$MaxScore < 0.05
] <- "Unassigned"

myeloid_df <- expr_sub %>%
  
  select(
    cell,
    Myeloid_Subtype,
    all_of(all_genes)
  ) %>%
  
  left_join(
    meta %>%
      select(
        cell,
        cell_type,
        Triad,
        sample_id,
        FOV
      ),
    by = "cell"
  ) %>%
  
  filter(
    grepl(
      "Myeloid",
      cell_type
    )
  )

subtypes <- c(
  "Microglia",
  "Macrophage",
  "Monocyte",
  "DC"
)

# Keep only assigned myeloid states
myeloid_df_lmm <- myeloid_df %>%
  filter(
    Myeloid_Subtype %in% subtypes
  )

or_tbl <- lapply(
  
  subtypes,
  
  function(st){
    
    message("Processing: ", st)
    
    df <- myeloid_df_lmm
    
    # Binary outcome:
    # Is the cell of this subtype?
    df$Target <- as.integer(
      df$Myeloid_Subtype == st
    )
    
    fit <- tryCatch(
      
      glmer(
        Target ~
          Triad +
          (1 | sample_id/FOV),
        data = df,
        family = binomial,
        control = glmerControl(
          optimizer = "bobyqa"
        )
      ),
      
      error = function(e){
        
        message(
          "Failed for ",
          st,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    if (is.null(fit))
      return(NULL)
    
    coef_tab <- summary(fit)$coefficients
    
    beta <- coef_tab[
      "TriadTRUE",
      "Estimate"
    ]
    
    se <- coef_tab[
      "TriadTRUE",
      "Std. Error"
    ]
    
    pval <- coef_tab[
      "TriadTRUE",
      "Pr(>|z|)"
    ]
    
    data.frame(
      
      Subtype = st,
      
      OR = exp(beta),
      
      CI_low = exp(
        beta - 1.96 * se
      ),
      
      CI_high = exp(
        beta + 1.96 * se
      ),
      
      Beta = beta,
      
      P = pval,
      
      stringsAsFactors = FALSE
    )
  }
  
) %>%
  
  bind_rows()

# Multiple-testing correction
or_tbl <- or_tbl %>%
  
  mutate(
    
    FDR = p.adjust(
      P,
      method = "BH"
    ),
    
    Sig = case_when(
      
      FDR < 0.001 ~ "***",
      
      FDR < 0.01 ~ "**",
      
      FDR < 0.05 ~ "*",
      
      TRUE ~ ""
    )
  )

# Sample-level consistency
sample_consistency <- lapply(
  
  subtypes,
  
  function(st){
    
    myeloid_df_lmm %>%
      
      group_by(
        sample_id,
        Triad
      ) %>%
      
      summarise(
        Fraction =
          mean(
            Myeloid_Subtype == st
          ),
        .groups = "drop"
      ) %>%
      
      tidyr::pivot_wider(
        names_from = Triad,
        values_from = Fraction
      ) %>%
      
      mutate(
        
        Delta =
          `TRUE` - `FALSE`,
        
        Subtype = st
      )
  }
  
) %>%
  
  bind_rows()

consistency_tbl <- sample_consistency %>%
  
  group_by(Subtype) %>%
  
  summarise(
    
    N_Pos = sum(
      Delta > 0,
      na.rm = TRUE
    ),
    
    N_Neg = sum(
      Delta < 0,
      na.rm = TRUE
    ),
    
    N = N_Pos + N_Neg,
    
    Consistency =
      ifelse(
        N_Pos >= N_Neg,
        paste0(
          N_Pos,
          "/",
          N,
          "+"
        ),
        paste0(
          N_Neg,
          "/",
          N,
          "-"
        )
      ),
    
    .groups = "drop"
  )

or_tbl <- or_tbl %>%
  
  left_join(
    consistency_tbl,
    by = "Subtype"
  )

write.csv(
  or_tbl,
  "../../output/FigS3_Myeloid_OR_LMM.csv",
  row.names = FALSE
)

print(or_tbl)

# Forest plot: myeloid subtype enrichment in triads
plot_or <- or_tbl %>%
  
  arrange(OR) %>%
  
  mutate(
    
    Subtype = factor(
      Subtype,
      levels = Subtype
    ),
    
    Label = paste0(
      sprintf("%.2f", OR),
      " ",
      Sig,
      "   "
    )
  )

p_myeloid_or <-
  
  ggplot(
    plot_or,
    aes(
      x = OR,
      y = Subtype,
      fill = Subtype
    )
  ) +
  
  # 95% CI
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.25,
    linewidth = 1
  ) +
  
  # OR estimate
  geom_point(
    shape = 21,
    size = 3,
    stroke = 0.6,
    color = "black"
  ) +
  
  # null effect
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  
  # labels
  geom_text(
    aes(
      x = CI_high + 0.75,
      label = Label
    ),
    hjust = 0,
    size = 2.8
  ) +
  
  scale_fill_manual(
    values = c(
      "Microglia"  = "#4DAF4A",
      "Macrophage" = "#E41A1C",
      "Monocyte"   = "#FF7F00",
      "DC"         = "#377EB8"
    )
  ) +
  
  labs(
    x = "Odds ratio (triad enrichment)",
    y = NULL
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      max(plot_or$CI_high) * 1.45
    ),
    clip = "off"
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "none",
    plot.margin = margin(
      5,
      35,
      5,
      5
    )
  )

save_pdf(
  "../../plots/FigS3_Myeloid_OR_Forest.pdf",
  p_myeloid_or,
  P_Width,
  P_Height
)

# Myeloid subtype marker validation
marker_order <- c(
  
  # Microglia
  "P2RY12",
  "CX3CR1",
  "GPR34",
  
  # Macrophage
  "CD163",
  "MARCO",
  "VSIG4",
  "FCGR3A",
  
  # Monocyte
  "CCR2",
  "CD14",
  "VCAN",
  "FCGR1A",
  "S100A9",
  
  # DC
  "CD1C",
  "FCER1A",
  "CLEC10A",
  "CD83",
  "CCR7",
  "LAMP3"
  
)

subtype_order <- c(
  "Microglia",
  "Macrophage",
  "Monocyte",
  "DC"
)

dot_df <- myeloid_df %>%
  
  filter(
    Myeloid_Subtype != "Unassigned"
  ) %>%
  
  pivot_longer(
    cols = all_of(marker_order),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  
  group_by(
    Myeloid_Subtype,
    Gene
  ) %>%
  
  summarise(
    avg = mean(Expression),
    pct = mean(Expression > 0) * 100,
    .groups = "drop"
  )

# Scale within gene (same as cell-type marker plot)

dot_df <- dot_df %>%
  
  group_by(Gene) %>%
  
  mutate(
    avg_scaled =
      as.numeric(
        scale(avg)
      )
  ) %>%
  
  ungroup()

dot_df$Gene <- factor(
  dot_df$Gene,
  levels = marker_order
)

dot_df$Myeloid_Subtype <- factor(
  dot_df$Myeloid_Subtype,
  levels = subtype_order
)

write.csv(
  dot_df,
  "../../output/FigS3_Myeloid_Marker_Dotplot_Data.csv",
  row.names = FALSE
)

# Plot
p_dot <- ggplot(
  dot_df,
  aes(
    Gene,
    Myeloid_Subtype
  )
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
    
    legend.title = element_text(
      size = 5
    ),
    
    legend.text = element_text(
      size = 4.5
    ),
    
    legend.key.height =
      unit(0.25, "cm"),
    
    legend.key.width =
      unit(0.25, "cm")
  ) +
  
  labs(
    x = NULL,
    y = NULL
  )

save_pdf(
  "../../plots/FigS3_Myeloid_Marker_DotPlot.pdf",
  p_dot,
  1.2 * P_Width,
  P_Height
)

# ============================================================
# Myeloid subtype triad-associated programs
# ============================================================

MODULES_MY <- MODULES_EXT[
  c(
    "Inflammatory chemokine",
    "Antigen presentation",
    "IFNγ response"
  )
]

myeloid_module_plot_data <- list()

myeloid_module_list <- list()

myeloid_module_samples_list <- list()

for (subtype in subtypes) {
  
  message("Processing: ", subtype)
  
  cells_use <- myeloid_df_lmm$cell[
    myeloid_df_lmm$Myeloid_Subtype == subtype
  ]
  
  for (mod_name in names(MODULES_MY)) {
    
    genes <- MODULES_MY[[mod_name]]
    
    if (length(genes) < 2)
      next
    
    df_lmm <- myeloid_niche %>%
      dplyr::filter(
        cell %in% cells_use
      )
    
    colnames(df_lmm)[
      colnames(df_lmm) == "sample_id.x"
    ] <- "sample_id"
    
    colnames(df_lmm)[
      colnames(df_lmm) == "FOV.x"
    ] <- "FOV"
    
    keep_cols <- c(
      "cell",
      "sample_id",
      "FOV",
      "Specific_Niche"
    )
    
    df_lmm <- df_lmm[
      ,
      intersect(
        keep_cols,
        colnames(df_lmm)
      )
    ]
    
    if (nrow(df_lmm) < 30)
      next
    
    if (length(unique(df_lmm$Specific_Niche)) < 2)
      next
    
    expr_mat <- FetchData(
      xe_global,
      vars = genes,
      cells = df_lmm$cell,
      layer = "data"
    )
    
    module_score <- rowMeans(
      expr_mat,
      na.rm = TRUE
    )
    
    if (sd(module_score, na.rm = TRUE) == 0)
      next
    
    # Raw score for biological interpretation
    df_lmm$Module_Score <- module_score
    
    # Z-score for statistical testing
    df_lmm$Score_Z <- as.numeric(
      scale(module_score)
    )
    
    # Save for plotting / validation
    myeloid_module_plot_data[[
      paste0(
        subtype,
        "_",
        mod_name
      )
    ]] <- df_lmm
    
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
        df_lmm$Specific_Niche == "Outside triad"
      ],
      na.rm = TRUE
    )
    
    mean_inside <- mean(
      df_lmm$Module_Score[
        df_lmm$Specific_Niche == "In triad"
      ],
      na.rm = TRUE
    )
    
    pct_change <- ifelse(
      abs(mean_outside) > 1e-8,
      100 * (mean_inside - mean_outside) /
        abs(mean_outside),
      NA_real_
    )
    
    # LMM on standardized score
    fit <- tryCatch(
      
      lmer(
        Score_Z ~
          Specific_Niche +
          (1 | sample_id/FOV),
        data = df_lmm
      ),
      
      error = function(e) NULL
    )
    
    if (is.null(fit))
      next
    
    coef_tab <- summary(fit)$coefficients
    
    coef_name <- rownames(coef_tab)[
      grepl(
        "Specific_Niche",
        rownames(coef_tab)
      )
    ][1]
    
    if (is.na(coef_name))
      next
    
    beta <- coef_tab[coef_name, "Estimate"]
    se   <- coef_tab[coef_name, "Std. Error"]
    pval <- coef_tab[coef_name, "Pr(>|t|)"]
    
    ci_low  <- beta - 1.96 * se
    ci_high <- beta + 1.96 * se
    
    myeloid_module_list[[paste0(
      subtype,
      "_",
      mod_name
    )]] <- data.frame(
      
      Subtype = subtype,
      Module = mod_name,
      
      N_Cells = nrow(df_lmm),
      
      Mean_Outside = mean_outside,
      Mean_InTriad = mean_inside,
      
      Percent_Change = pct_change,
      
      Beta = beta,
      P = pval,
      
      CI_low = ci_low,
      CI_high = ci_high,
      
      stringsAsFactors = FALSE
    )
    
    # Sample-level direction consistency
    for (sm in unique(df_lmm$sample_id)) {
      
      tmp <- df_lmm %>%
        dplyr::filter(
          sample_id == sm
        )
      
      if (nrow(tmp) < 10)
        next
      
      if (length(unique(tmp$Specific_Niche)) < 2)
        next
      
      mean_out_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche == "Outside triad"
        ],
        na.rm = TRUE
      )
      
      mean_in_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche == "In triad"
        ],
        na.rm = TRUE
      )
      
      pct_sm <- ifelse(
        abs(mean_out_sm) > 1e-8,
        100 * (mean_in_sm - mean_out_sm) /
          abs(mean_out_sm),
        NA_real_
      )
      
      myeloid_module_samples_list[[paste0(
        subtype,
        "_",
        mod_name,
        "_",
        sm
      )]] <- data.frame(
        
        Subtype = subtype,
        Module = mod_name,
        sample_id = sm,
        
        Percent_Change = pct_sm,
        
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine results
myeloid_module_stats <- bind_rows(
  myeloid_module_list
)

myeloid_module_samples <- bind_rows(
  myeloid_module_samples_list
)

# Patient-level consistency
consistency_summary_myeloid_module <- myeloid_module_samples %>%
  dplyr::group_by(
    Subtype,
    Module
  ) %>%
  dplyr::summarise(
    
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

myeloid_module_stats <- myeloid_module_stats %>%
  left_join(
    consistency_summary_myeloid_module,
    by = c(
      "Subtype",
      "Module"
    )
  )

myeloid_module_stats$Consistency <- ifelse(
  
  myeloid_module_stats$Beta > 0,
  
  paste0(
    myeloid_module_stats$N_Pos,
    "/",
    myeloid_module_stats$N_Pos +
      myeloid_module_stats$N_Neg,
    "+"
  ),
  
  paste0(
    myeloid_module_stats$N_Neg,
    "/",
    myeloid_module_stats$N_Pos +
      myeloid_module_stats$N_Neg,
    "-"
  )
)

# Multiple testing correction
myeloid_module_stats$FDR <- p.adjust(
  myeloid_module_stats$P,
  method = "BH"
)

myeloid_module_stats <- myeloid_module_stats %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Sample-level validation
myeloid_module_validation <- list()
myeloid_module_delta_list <- list()

for (nm in names(myeloid_module_plot_data)) {
  
  tmp <- myeloid_module_plot_data[[nm]]
  
  parts <- strsplit(
    nm,
    "_"
  )[[1]]
  
  subtype <- parts[1]
  
  module <- paste(
    parts[-1],
    collapse = "_"
  )
  
  sample_summary <- tmp %>%
    
    group_by(
      sample_id,
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
      Delta =
        `In triad` -
        `Outside triad`
    )
  
  sample_summary$Subtype <- subtype
  sample_summary$Module <- module
  
  myeloid_module_delta_list[[nm]] <- sample_summary
  
  wilcox_p <- tryCatch(
    
    wilcox.test(
      sample_summary$Delta,
      mu = 0,
      exact = FALSE
    )$p.value,
    
    error = function(e) NA_real_
  )
  
  myeloid_module_validation[[nm]] <- data.frame(
    
    Subtype = subtype,
    Module = module,
    
    N_Samples = nrow(
      sample_summary
    ),
    
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

myeloid_module_validation_df <- bind_rows(
  myeloid_module_validation
)

myeloid_module_validation_df$Wilcox_FDR <- p.adjust(
  myeloid_module_validation_df$Wilcox_P,
  method = "BH"
)

myeloid_module_deltas <- bind_rows(
  myeloid_module_delta_list
)

# Save output
myeloid_module_final <- myeloid_module_stats %>%
  
  left_join(
    myeloid_module_validation_df,
    by = c(
      "Subtype",
      "Module"
    )
  ) %>%
  
  select(
    
    Subtype,
    Module,
    
    Mean_Outside,
    Mean_InTriad,
    
    Percent_Change,
    
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
  myeloid_module_final,
  "../../output/Fig3_Myeloid_Subtype_Statistics.csv",
  row.names = FALSE
)

# Heatmap: subtype-specific triad programs
heat_df <- myeloid_module_final %>%
  
  dplyr::select(
    Subtype,
    Module,
    Beta,
    Sig
  )

heat_df$Subtype <- factor(
  heat_df$Subtype,
  levels = c(
    "DC",
    "Monocyte",
    "Macrophage",
    "Microglia"
  )
)

heat_df$Module <- factor(
  heat_df$Module,
  levels = c(
    "Inflammatory chemokine",
    "Antigen presentation",
    "IFNγ response"
  )
)

beta_lim <- max(
  abs(heat_df$Beta),
  na.rm = TRUE
)

p_myeloid_heat <- ggplot(
  heat_df,
  aes(
    x = Module,
    y = Subtype,
    fill = Beta
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
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
    name = expression(beta)
  ) +
  
  scale_x_discrete(
    labels = c(
      "Inflammatory chemokine" = "Inflammatory\nchemokine",
      "Antigen presentation"   = "Antigen\npresentation",
      "IFNγ response"          = "IFNγ\nresponse"
    )
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_fig() +
  
  theme(
    
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = pt_size * 0.8,
    ),
    
    legend.text = element_text(
      size = pt_size * 0.8
    )
    
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(1.5, "cm"),
      barwidth = unit(0.35, "cm")
    )
  )

save_pdf(
  "../../plots/Fig3_Myeloid_Module_Heatmap.pdf",
  p_myeloid_heat,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# CD4 vs CD8 triad-associated programs
# ============================================================

# Define CD4 and CD8 cells
expr_cd <- FetchData(
  xe_global,
  vars = c("CD4", "CD8A", "CD8B"),
  layer = "data"
)

expr_cd$cell <- rownames(expr_cd)

t_meta <- meta %>%
  filter(
    grepl(
      T_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  left_join(
    expr_cd,
    by = "cell"
  ) %>%
  mutate(
    
    CD4_only =
      CD4 > 0 &
      CD8A == 0 &
      CD8B == 0,
    
    CD8_only =
      CD4 == 0 &
      (CD8A > 0 | CD8B > 0)
    
  )

message(
  "CD4-only cells: ",
  sum(t_meta$CD4_only)
)

message(
  "CD8-only cells: ",
  sum(t_meta$CD8_only)
)

# T cell modules
MODULES_CD <- MODULES_EXT[
  c(
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T"
  )
]

# Mixed-model analysis
cd4_cd8_plot_data <- list()

cd4_cd8_list <- list()
cd4_cd8_samples_list <- list()

for (lineage in c("CD4", "CD8")) {
  
  message("Processing: ", lineage)
  
  cells_use <- switch(
    
    lineage,
    
    "CD4" =
      t_meta$cell[
        t_meta$CD4_only
      ],
    
    "CD8" =
      t_meta$cell[
        t_meta$CD8_only
      ]
  )
  
  for (mod_name in names(MODULES_CD)) {
    
    genes <- MODULES_CD[[mod_name]]
    
    if (length(genes) < 2)
      next
    
    df_lmm <- tcell_niche %>%
      dplyr::filter(
        cell %in% cells_use
      )
    
    colnames(df_lmm)[
      colnames(df_lmm) == "sample_id.x"
    ] <- "sample_id"
    
    colnames(df_lmm)[
      colnames(df_lmm) == "FOV.x"
    ] <- "FOV"
    
    keep_cols <- c(
      "cell",
      "sample_id",
      "FOV",
      "Specific_Niche"
    )
    
    df_lmm <- df_lmm[
      ,
      intersect(
        keep_cols,
        colnames(df_lmm)
      )
    ]
    
    if (nrow(df_lmm) < 100)
      next
    
    if (length(unique(df_lmm$Specific_Niche)) < 2)
      next
    
    expr_mat <- FetchData(
      xe_global,
      vars = genes,
      cells = df_lmm$cell,
      layer = "data"
    )
    
    module_score <- rowMeans(
      expr_mat,
      na.rm = TRUE
    )
    
    if (sd(module_score, na.rm = TRUE) == 0)
      next
    
    # Raw score for biological interpretation
    df_lmm$Module_Score <- module_score
    
    # Z-score for statistical testing
    df_lmm$Score_Z <- as.numeric(
      scale(module_score)
    )
    
    # Save for plotting / validation
    cd4_cd8_plot_data[[
      paste0(
        lineage,
        "_",
        mod_name
      )
    ]] <- df_lmm
    
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
        df_lmm$Specific_Niche == "Outside triad"
      ],
      na.rm = TRUE
    )
    
    mean_inside <- mean(
      df_lmm$Module_Score[
        df_lmm$Specific_Niche == "In triad"
      ],
      na.rm = TRUE
    )
    
    pct_change <- ifelse(
      abs(mean_outside) > 1e-8,
      100 * (mean_inside - mean_outside) /
        abs(mean_outside),
      NA_real_
    )
    
    # LMM on standardized score
    fit <- tryCatch(
      
      lmer(
        Score_Z ~
          Specific_Niche +
          (1 | sample_id/FOV),
        data = df_lmm
      ),
      
      error = function(e) NULL
    )
    
    if (is.null(fit))
      next
    
    coef_tab <- summary(fit)$coefficients
    
    coef_name <- rownames(coef_tab)[
      grepl(
        "Specific_Niche",
        rownames(coef_tab)
      )
    ][1]
    
    if (is.na(coef_name))
      next
    
    beta <- coef_tab[coef_name, "Estimate"]
    se   <- coef_tab[coef_name, "Std. Error"]
    pval <- coef_tab[coef_name, "Pr(>|t|)"]
    
    ci_low  <- beta - 1.96 * se
    ci_high <- beta + 1.96 * se
    
    cd4_cd8_list[[paste0(
      lineage,
      "_",
      mod_name
    )]] <- data.frame(
      
      Lineage = lineage,
      Module = mod_name,
      
      Mean_Outside = mean_outside,
      Mean_InTriad = mean_inside,
      
      Percent_Change = pct_change,
      
      Beta = beta,
      P = pval,
      
      CI_low  = ci_low,
      CI_high = ci_high,
      
      stringsAsFactors = FALSE
    )
    
    # Sample-level direction consistency
    for (sm in unique(df_lmm$sample_id)) {
      
      tmp <- df_lmm %>%
        dplyr::filter(
          sample_id == sm
        )
      
      if (nrow(tmp) < 20)
        next
      
      if (length(unique(tmp$Specific_Niche)) < 2)
        next
      
      mean_out_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche == "Outside triad"
        ],
        na.rm = TRUE
      )
      
      mean_in_sm <- mean(
        tmp$Module_Score[
          tmp$Specific_Niche == "In triad"
        ],
        na.rm = TRUE
      )
      
      pct_sm <- ifelse(
        abs(mean_out_sm) > 1e-8,
        100 * (mean_in_sm - mean_out_sm) /
          abs(mean_out_sm),
        NA_real_
      )
      
      cd4_cd8_samples_list[[paste0(
        lineage,
        "_",
        mod_name,
        "_",
        sm
      )]] <- data.frame(
        
        Lineage = lineage,
        Module = mod_name,
        sample_id = sm,
        
        Percent_Change = pct_sm,
        
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine results
cd4_cd8_stats <- bind_rows(
  cd4_cd8_list
)

cd4_cd8_samples <- bind_rows(
  cd4_cd8_samples_list
)

# Patient-level consistency
consistency_summary_cd4cd8 <- cd4_cd8_samples %>%
  dplyr::group_by(
    Lineage,
    Module
  ) %>%
  dplyr::summarise(
    
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

cd4_cd8_stats <- cd4_cd8_stats %>%
  left_join(
    consistency_summary_cd4cd8,
    by = c(
      "Lineage",
      "Module"
    )
  )

cd4_cd8_stats$Consistency <- ifelse(
  
  cd4_cd8_stats$Beta > 0,
  
  paste0(
    cd4_cd8_stats$N_Pos,
    "/",
    cd4_cd8_stats$N_Pos +
      cd4_cd8_stats$N_Neg,
    "+"
  ),
  
  paste0(
    cd4_cd8_stats$N_Neg,
    "/",
    cd4_cd8_stats$N_Pos +
      cd4_cd8_stats$N_Neg,
    "-"
  )
)

# Multiple testing correction
cd4_cd8_stats$FDR <- p.adjust(
  cd4_cd8_stats$P,
  method = "BH"
)

cd4_cd8_stats <- cd4_cd8_stats %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Sample-level validation
cd4_cd8_validation <- list()
cd4_cd8_delta_list <- list()

for (nm in names(cd4_cd8_plot_data)) {
  
  tmp <- cd4_cd8_plot_data[[nm]]
  
  parts <- strsplit(
    nm,
    "_"
  )[[1]]
  
  lineage <- parts[1]
  
  module <- paste(
    parts[-1],
    collapse = "_"
  )
  
  sample_summary <- tmp %>%
    
    group_by(
      sample_id,
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
      Delta =
        `In triad` -
        `Outside triad`
    )
  
  sample_summary$Lineage <- lineage
  sample_summary$Module <- module
  
  cd4_cd8_delta_list[[nm]] <- sample_summary
  
  wilcox_p <- tryCatch(
    
    wilcox.test(
      sample_summary$Delta,
      mu = 0,
      exact = FALSE
    )$p.value,
    
    error = function(e) NA_real_
  )
  
  cd4_cd8_validation[[nm]] <- data.frame(
    
    Lineage = lineage,
    Module = module,
    
    N_Samples = nrow(
      sample_summary
    ),
    
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

cd4_cd8_validation_df <- bind_rows(
  cd4_cd8_validation
)

cd4_cd8_validation_df$Wilcox_FDR <- p.adjust(
  cd4_cd8_validation_df$Wilcox_P,
  method = "BH"
)

cd4_cd8_deltas <- bind_rows(
  cd4_cd8_delta_list
)

# Save output
cd4_cd8_final <- cd4_cd8_stats %>%
  
  left_join(
    cd4_cd8_validation_df,
    by = c(
      "Lineage",
      "Module"
    )
  ) %>%
  
  select(
    
    Lineage,
    Module,
    
    Mean_Outside,
    Mean_InTriad,
    
    Percent_Change,
    
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
  cd4_cd8_final,
  "../../output/Fig3_CD4_CD8_Statistics.csv",
  row.names = FALSE
)

# Heatmap
heat_df_cd4_cd8 <- cd4_cd8_final %>%
  dplyr::select(
    Lineage,
    Module,
    Beta,
    Sig
  )

heat_df_cd4_cd8$Lineage <- factor(
  heat_df_cd4_cd8$Lineage,
  levels = c(
    "CD8",
    "CD4"
  )
)

heat_df_cd4_cd8$Module <- factor(
  heat_df_cd4_cd8$Module,
  levels = c(
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T"
  )
)

beta_lim <- max(
  abs(heat_df_cd4_cd8$Beta),
  na.rm = TRUE
)

p_cd4_cd8_heat <- ggplot(
  heat_df_cd4_cd8,
  aes(
    x = Module,
    y = Lineage,
    fill = Beta
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
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
    name = expression(beta)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_fig() +
  
  theme(
    
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = pt_size * 0.8,
    ),
    
    legend.text = element_text(
      size = pt_size * 0.8
    )
    
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(1.5, "cm"),
      barwidth = unit(0.35, "cm")
    )
  )

save_pdf(
  "../../plots/Fig3_CD4_CD8_Module_Heatmap.pdf",
  p_cd4_cd8_heat,
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
  "HLA-B", "PSMB10", "CTSS", "CTSL", "CTSB", "CD74"
)

IFN_GENE_POOL <- c(
  "STAT1", "STAT2", "IRF1", "CXCL9", "CXCL10", "ISG15", 
  "IFIT2", "IFIT3", "IFITM3", "MX1"
)

HYPOXIA_GENE_POOL <- c(
  "VEGFA", "HIF1A", "FLT1", "FGF2", "MMP9"
)

TUMOR_MODULES <- list(
  "Antigen processing" = AP_GENE_POOL,
  "Interferon response" = IFN_GENE_POOL,
  "Hypoxia/angiogenesis" = HYPOXIA_GENE_POOL
)

TUMOR_MODULES <- lapply(
  TUMOR_MODULES,
  intersect,
  rownames(xe_global)
)

TUMOR_MODULES <- TUMOR_MODULES[
  sapply(TUMOR_MODULES, length) >= 2
]

message("\n=== Tumor modules ===")

print(
  sapply(
    TUMOR_MODULES,
    length
  )
)

module_info <- tibble::tribble(
  
  ~Module,                  ~Group,
  
  "Antigen processing",   "Tumor",
  "Interferon response",    "Tumor",
  "Hypoxia/angiogenesis", "Tumor"
  
)

# Compute tumor triad density
tumor_meta$triad_density <- 0
for (sid in unique(tumor_meta$sample_id)) {
  
  for (fov in unique(
    tumor_meta$FOV[
      tumor_meta$sample_id == sid
    ]
  )) {
    
    idx <- which(
      tumor_meta$sample_id == sid &
        tumor_meta$FOV == fov
    )
    
    tri_idx <- which(
      triad_xy$sample_id == sid &
        triad_xy$FOV == fov
    )
    
    if (
      length(idx) == 0 ||
      length(tri_idx) == 0
    ) {
      next
    }
    
    tumor_meta$triad_density[idx] <-
      gauss_kde_at(
        query_xy =
          tumor_meta[idx, c("x", "y")],
        
        anchor_xy =
          triad_xy[tri_idx, c("x", "y")]
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

message(
  "\n=== Tumor triad density summary ==="
)

print(
  summary(
    tumor_meta$triad_density
  )
)

print(
  table(
    tumor_meta$triad_density > 0
  )
)

density_summary <- tumor_meta %>%
  summarise(
    N_cells = n(),
    Median = median(triad_density),
    Mean = mean(triad_density),
    SD = sd(triad_density)
  )

write.csv(
  density_summary,
  "../../output/Fig3_TriadDensity_Summary.csv",
  row.names = FALSE
)

# Triad-associated tumor programs
tumor_plot_data <- list()

tumor_lmm_list <- list()

tumor_sample_effects <- list()

for (mod_name in names(TUMOR_MODULES)) {
  
  message("Processing: ", mod_name)
  
  genes <- TUMOR_MODULES[[mod_name]]
  
  if (length(genes) < 2)
    next
  
  expr_mat <- FetchData(
    xe_global,
    vars = genes,
    cells = tumor_meta$cell,
    layer = "data"
  )
  
  module_score <- rowMeans(
    expr_mat,
    na.rm = TRUE
  )
  
  if (sd(module_score, na.rm = TRUE) == 0)
    next
  
  df_lmm <- tumor_meta
  
  df_lmm$Module_Score <- module_score
  
  df_lmm$Score_Z <- as.numeric(
    scale(module_score)
  )
  
  tumor_plot_data[[mod_name]] <- df_lmm
  
  fit <- tryCatch(
    
    lmer(
      Score_Z ~
        TriadDensity_Z +
        (1 | sample_id/FOV),
      data = df_lmm
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(fit))
    next
  
  coef_tab <- summary(fit)$coefficients
  
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
  
  ci_low  <- beta - 1.96 * se
  ci_high <- beta + 1.96 * se
  
  tumor_lmm_list[[mod_name]] <- data.frame(
    
    Module = mod_name,
    
    Beta = beta,
    P = pval,
    
    CI_low  = ci_low,
    CI_high = ci_high,
    
    stringsAsFactors = FALSE
  )
  
  # Sample-level direction consistency
  for (sm in unique(df_lmm$sample_id)) {
    
    tmp <- df_lmm %>%
      dplyr::filter(
        sample_id == sm
      )
    
    if (nrow(tmp) < 100)
      next
    
    fit_sm <- tryCatch(
      
      lm(
        Module_Score ~ TriadDensity_Z,
        data = tmp
      ),
      
      error = function(e) NULL
    )
    
    if (is.null(fit_sm))
      next
    
    cf <- summary(fit_sm)$coefficients
    
    if (!"TriadDensity_Z" %in% rownames(cf))
      next
    
    tumor_sample_effects[[paste0(mod_name, "_", sm)]] <-
      data.frame(
        
        Module = mod_name,
        sample_id = sm,
        
        Effect = cf[
          "TriadDensity_Z",
          "Estimate"
        ],
        
        stringsAsFactors = FALSE
      )
  }
}

# Combine results
tumor_stats <- bind_rows(
  tumor_lmm_list
)

tumor_sample_effects <- bind_rows(
  tumor_sample_effects
)

tumor_stats$Group <- "Tumor Cells"

# Patient-level consistency
consistency_summary <- tumor_sample_effects %>%
  group_by(Module) %>%
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
    by = "Module"
  )

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

# Multiple testing correction
tumor_stats$FDR <- p.adjust(
  tumor_stats$P,
  method = "BH"
)

tumor_stats <- tumor_stats %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Percentage positive cell
pct_df_tumor <- lapply(
  names(tumor_plot_data),
  function(mod_name) {
    
    df_plot <- tumor_plot_data[[mod_name]]
    
    data.frame(
      
      Module = mod_name,
      
      Percent_Positive =
        mean(
          df_plot$Module_Score > 0,
          na.rm = TRUE
        ) * 100
    )
  }
) %>%
  bind_rows()

# Plot ordering
module_order <- c(
  
  "Antigen processing",
  "Interferon response",
  "Hypoxia/angiogenesis"
  
)

tumor_stats$Module <- factor(
  tumor_stats$Module,
  levels = rev(module_order)
)

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

for (mod_name in names(TUMOR_MODULES)) {
  
  tmp <- tumor_sample_effects %>%
    dplyr::filter(
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
  
  tumor_validation[[mod_name]] <- data.frame(
    
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
    by = "Module"
  ) %>%
  
  select(
    
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
  "../../output/Fig3_Tumor_Program_Statistics.csv",
  row.names = FALSE
)

# Module summary dot plot
tumor_dot <- tumor_final %>%
  dplyr::select(
    Module,
    Beta,
    FDR,
    Sig
  ) %>%
  left_join(
    pct_df_tumor,
    by = "Module"
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
    size = 7,
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
  "../../plots/Fig3_Triad_Tumor_Modules.pdf",
  p_tumor,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Triad density effect with and without vessel-distance adjustment
# ============================================================

hypoxia_genes <- TUMOR_MODULES[["Hypoxia/angiogenesis"]]

expr_hyp <- FetchData(
  xe_global,
  vars = hypoxia_genes,
  cells = tumor_meta$cell,
  layer = "data"
)

tumor_meta$Hypoxia_Z <- as.numeric(
  scale(
    rowMeans(
      expr_hyp,
      na.rm = TRUE
    )
  )
)

# Vessel coordinates
vessel_xy <- meta %>%
  dplyr::filter(
    cell_type == "Vasculature/Stroma"
  ) %>%
  dplyr::select(
    sample_id,
    FOV,
    x,
    y
  )

# Distance to nearest vessel
tumor_meta$Dist_Vessel <- Inf

for (sid in unique(tumor_meta$sample_id)) {
  
  for (fov in unique(
    tumor_meta$FOV[
      tumor_meta$sample_id == sid
    ]
  )) {
    
    idx <- which(
      tumor_meta$sample_id == sid &
        tumor_meta$FOV == fov
    )
    
    v_idx <- which(
      vessel_xy$sample_id == sid &
        vessel_xy$FOV == fov
    )
    
    if (
      length(idx) == 0 ||
      length(v_idx) == 0
    ) next
    
    tumor_meta$Dist_Vessel[idx] <-
      nn_dists(
        tumor_meta[idx, c("x", "y")],
        vessel_xy[v_idx, c("x", "y")]
      )
  }
}

tumor_meta$Dist_Vessel_Z <- as.numeric(
  scale(
    log1p(tumor_meta$Dist_Vessel)
  )
)

# Pooled models
fit_vessel <- lmer(
  Hypoxia_Z ~
    Dist_Vessel_Z +
    (1 | sample_id/FOV),
  data = tumor_meta
)

fit_unadj <- lmer(
  Hypoxia_Z ~
    TriadDensity_Z +
    (1 | sample_id/FOV),
  data = tumor_meta
)

fit_adj <- lmer(
  Hypoxia_Z ~
    TriadDensity_Z +
    Dist_Vessel_Z +
    (1 | sample_id/FOV),
  data = tumor_meta
)

# Extract coefficients
extract_coef <- function(
    fit,
    term,
    model_name
) {
  
  cf <- summary(fit)$coefficients
  
  beta <- cf[
    term,
    "Estimate"
  ]
  
  se <- cf[
    term,
    "Std. Error"
  ]
  
  p <- cf[
    term,
    "Pr(>|t|)"
  ]
  
  data.frame(
    
    Model = model_name,
    
    Beta = beta,
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    
    P = p,
    
    stringsAsFactors = FALSE
  )
}

# Summary table
hypoxia_control <- bind_rows(
  
  extract_coef(
    fit_vessel,
    "Dist_Vessel_Z",
    "Vessel distance only"
  ) %>%
    mutate(Term = "Dist_Vessel_Z"),
  
  extract_coef(
    fit_unadj,
    "TriadDensity_Z",
    "Triad density only"
  ) %>%
    mutate(Term = "TriadDensity_Z"),
  
  extract_coef(
    fit_adj,
    "TriadDensity_Z",
    "Triad density +\nvessel distance"
  ) %>%
    mutate(Term = "TriadDensity_Z")
)

write.csv(
  hypoxia_control,
  "../../output/FigS3_Hypoxia_Vessel_Adjustment.csv",
  row.names = FALSE
)

cor.test(
  tumor_meta$TriadDensity_Z,
  tumor_meta$Dist_Vessel_Z,
  method = "spearman"
)

triad_unadj <- hypoxia_control$Beta[
  hypoxia_control$Model ==
    "Triad density only"
]

triad_adj <- hypoxia_control$Beta[
  hypoxia_control$Model ==
    "Triad density +\nvessel distance"
]

attenuation_pct <-
  100 *
  (abs(triad_unadj) - abs(triad_adj)) /
  abs(triad_unadj)

attenuation_summary <- data.frame(
  Unadjusted_Beta = triad_unadj,
  Adjusted_Beta = triad_adj,
  Percent_Attenuation = attenuation_pct
)

write.csv(
  attenuation_summary,
  "../../output/FigS3_Hypoxia_Vessel_Adjustment_Attenuation.csv",
  row.names = FALSE
)





