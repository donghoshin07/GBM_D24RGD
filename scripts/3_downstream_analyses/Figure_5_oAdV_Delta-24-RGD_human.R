# ==============================================================================
# Figure 5 analysis using Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts/3_downstream_analyses
#
# Required input files: 
# ../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds
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
  library(ggrastr)
  library(readr)
  library(ggpubr)
  library(mgcv)
  library(broom)
  
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

ASSAY_USE <- "Xenium"

NEAR_MAX    <- 60
FAR_MIN     <- 300
CONTACT_MAX <- 10

SCALE_LEN <- 1000

VIRAL_GENES <- c(
  "HAdVC-E1A",
  "HAdVC-E1B-55k",
  "HAdVC-Fiber",
  "HAdVC-Hexon"
)

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

meta$Timepoint <- ifelse(
  grepl("post", meta$sample_id),
  "post",
  "pre"
)

# Virus+ rule: >= 2 distinct viral genes 
viral_expr <- FetchData(
  xe_global,
  vars = intersect(
    VIRAL_GENES,
    rownames(xe_global)
  ),
  cells = meta$cell,
  layer = "counts"
)

meta$Virus_Pos <- rowSums(
  viral_expr > 0
) >= 2

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

message(
  "\nTotal triad anchors: ",
  length(triad_anchors)
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
# Figure 5 case selection
# ============================================================

P17_FOV <- "fov.9"

p17 <- meta |>
  dplyr::filter(
    Patient_Match == "P17",
    Timepoint == "post",
    FOV == P17_FOV
  )

stopifnot(
  nrow(p17) > 0
)

p17 %>%
  count(Timepoint)

p17 %>%
  count(cell_type, sort = TRUE)

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

# Average module score
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
    rowMeans(expr, na.rm = TRUE),
    cells
  )
  
}

# Return cells within radius of each anchor
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
    data  = as.matrix(
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

# Bootstrap CI for Near vs Far median difference
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

# Distance to nearest set (virus+ cell)
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

# ============================================================
# Spatial distribution of virus-positive cells
# ============================================================

x0 <- max(p17$x) - (SCALE_LEN + 40)
x1 <- x0 + SCALE_LEN
y0 <- min(p17$y) + 30

virus_label <- sprintf(
  "Virus-positive cells (n = %s)",
  sum(p17$Virus_Pos)
)

p17$Virus_Label <- virus_label

p_spatial_virus <- ggplot() +
  
  geom_point_rast(
    data = p17,
    aes(x, y),
    color = "grey90",
    size = 0.05,
    alpha = 0.35,
    raster.dpi = 600
  ) +
  
  geom_point_rast(
    data = filter(p17, Virus_Pos),
    aes(
      x,
      y,
      color = Virus_Label
    ),
    size = 0.1,
    alpha = 0.95,
    raster.dpi = 600
  ) +
  
  scale_color_manual(
    values = c("red3"),
    name = NULL
  ) +
  
  ggplot2::annotate(
    "segment",
    x = x0,
    xend = x1,
    y = y0,
    yend = y0,
    linewidth = 0.33
  ) +
  
  coord_fixed() +
  
  labs(
    title = "P17"
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
  "../../plots/Fig5_Virus_Positive_Spatial.pdf",
  p_spatial_virus,
  P_Width,
  P_Height
)

# ============================================================
# Composition of virus-positive cells
# ============================================================

viral_df <- p17 %>%
  filter(Virus_Pos)

comp_tbl <- viral_df %>%
  count(
    cell_type,
    name = "N"
  ) %>%
  mutate(
    Fraction = N / sum(N)
  ) %>%
  arrange(Fraction)

write_csv(
  comp_tbl,
  "../../output/Fig5_Virus_Composition.csv"
)

p_comp <- ggplot(
  comp_tbl,
  aes(
    x = Fraction,
    y = reorder(cell_type, Fraction),
    fill = cell_type
  )
) +
  
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = pal_best
  ) +
  
  scale_x_continuous(
    labels = scales::percent_format()
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "none"
  ) +
  
  labs(
    x = "Fraction of virus+ cells",
    y = NULL
  )

save_pdf(
  "../../plots/Fig5_Virus_Composition.pdf",
  p_comp,
  P_Width,
  P_Height
)

# ============================================================
# Virus+ vs Virus- module heatmap in P17-post tumor states 
# ============================================================

MIN_POS_PER_TYPE_TEST    <- 5
MIN_POS_PER_TYPE_DISPLAY <- 3

MODULES <- list(
  
  "Interferon response" = c(
    "STAT1","STAT2","IRF1",
    "ISG15","IFIT2","IFIT3",
    "IFITM3","MX1","LY6E"
  ),
  
  "Inflammatory chemokine" = c(
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CCL5"
  ),
  
  "Antigen processing" = c(
    "HLA-B",
    "PSMB10",
    "ICAM1"
  ),
  
  "Stress response" = c(
    "JUN",
    "JUNB",
    "FOS",
    "HIF1A",
    "XBP1",
    "SOCS1",
    "SOCS3"
  ),
  
  "CXCL16" = c(
    "CXCL16"
  )
)

# Function to compute module deltas
compute_module_deltas <- function(df_case){
  
  out <- list()
  
  for (ct in unique(df_case$cell_type)) {
    
    df_ct <- df_case %>%
      filter(cell_type == ct)
    
    pos <- df_ct %>%
      filter(Virus_Pos) %>%
      pull(cell)
    
    neg <- df_ct %>%
      filter(!Virus_Pos) %>%
      pull(cell)
    
    display_ok <-
      length(pos) >= MIN_POS_PER_TYPE_DISPLAY &&
      length(neg) >= 1
    
    test_ok <-
      length(pos) >= MIN_POS_PER_TYPE_TEST &&
      length(neg) >= MIN_POS_PER_TYPE_TEST
    
    if (!display_ok)
      next
    
    for (mod in names(MODULES)) {
      
      genes <- intersect(
        MODULES[[mod]],
        rownames(xe_global)
      )
      
      if (length(genes) < 1)
        next
      
      sc_pos <- module_score_vec(
        xe_global,
        pos,
        genes,
        layer = "data"
      )
      
      sc_neg <- module_score_vec(
        xe_global,
        neg,
        genes,
        layer = "data"
      )
      
      p_wc <- if (test_ok) {
        
        x <- c(
          as.numeric(sc_pos),
          as.numeric(sc_neg)
        )
        
        g <- factor(
          c(
            rep("pos", length(sc_pos)),
            rep("neg", length(sc_neg))
          )
        )
        
        wilcox.test(
          x ~ g
        )$p.value
        
      } else {
        
        NA_real_
        
      }
      
      out[[length(out) + 1]] <- tibble(
        cell_type = ct,
        module = mod,
        n_pos = length(sc_pos),
        n_neg = length(sc_neg),
        delta =
          mean(sc_pos, na.rm = TRUE) -
          mean(sc_neg, na.rm = TRUE),
        p_wilcox = p_wc,
        test_ok = test_ok
      )
      
    }
  }
  
  bind_rows(out)
}

# P17-post tumor cells only
module_df <- meta %>%
  filter(
    Patient_Match == "P17",
    Timepoint == "post",
    cell_type %in% TUMOR_LABELS
  )

modules_long <- compute_module_deltas(
  module_df
) %>%
  mutate(
    q_wilcox = p.adjust(
      replace_na(p_wilcox, 1),
      method = "BH"
    ),
    Stars = case_when(
      q_wilcox < 0.001 ~ "***",
      q_wilcox < 0.01  ~ "**",
      q_wilcox < 0.05  ~ "*",
      TRUE             ~ ""
    ),
    delta_plot = ifelse(
      test_ok,
      delta,
      NA_real_
    )
  )

write_csv(
  modules_long,
  "../../output/Fig5_Module_Delta.csv"
)

# Ordering
ct_order <- c(
  "MES-like/Hypoxic Tumor",
  "Astrocytes/AC-like",
  "Oligodendrocytes/OPC-like"
)

modules_long$cell_type <- factor(
  modules_long$cell_type,
  levels = ct_order
)

# Heatmap
p_modules <- ggplot(
  modules_long,
  aes(
    x = module,
    y = cell_type,
    fill = delta_plot
  )
) +
  
  geom_tile(
    color = "grey80",
    linewidth = 0.2
  ) +
  
  geom_text(
    aes(label = Stars),
    size = pt_size / 2.2
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Mean module score\n(Virus+ − Virus−)",
    na.value = "grey95"
  ) +
  
  scale_x_discrete(
    labels = c(
      "Interferon response"    = "Interferon\nresponse",
      "Inflammatory chemokine" = "Inflammatory\nchemokine",
      "Antigen processing"     = "Antigen\nprocessing",
      "Stress response"        = "Stress\nresponse",
      "CXCL16"                 = "CXCL16"
    )
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box.margin = margin(
      t = 20,
      r = 0,
      b = 0,
      l = 0
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  ) +
  
  labs(
    x = NULL,
    y = NULL
  )

save_pdf(
  "../../plots/Fig5_Module_Heatmap.pdf",
  p_modules,
  2*P_Width,
  P_Height
)

# ============================================================
# Representative virus-near triads and virus-distant triads
# ============================================================

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
virus_xy <- p17 %>%
  filter(Virus_Pos) %>%
  select(x, y)

# Triad coordinates
tri_xy <- p17 %>%
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
cx <- 9500   # manual x coordinate
cy <- 2100   # manual y coordinate

WINDOW <- 1000

xmin <- cx - WINDOW
xmax <- cx + WINDOW

ymin <- cy - WINDOW
ymax <- cy + WINDOW


# Data in plotting window
plot_cells <- p17 %>%
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
    title = "P17"
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
  "../../plots/Fig5_Spatial_Virus_Near_Distant_Triads.pdf",
  p_triad_zoom,
  P_Width,
  P_Height
)

# ============================================================
# Virus-near vs virus-distant triad phenotypes
# ============================================================

MODULES_EXT <- list(
  
  # Myeloid
  "Inflammatory chemokine" = c("CXCL9","CXCL10","CXCL11","CCL5"),
  "Antigen presentation"   = c("CD74","CTSS","CD80","CD86"),
  "IFNγ response" = c("STAT1","IRF1","IFNGR1"),
  
  # T/NK Cells
  "Stemness"     = c("TCF7","IL7R","SELL"),
  "Cytotoxicity" = c("PRF1","NKG7","GNLY","GZMB"),
  "Exhaustion"   = c("PDCD1","LAG3","HAVCR2","TOX"),
  "Regulatory T" = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  "Chemokine receptor" = c("CCR5","CXCR6"),
  
  # B/Plasma Cells
  "Lymphoid organization" = c("CXCL13","LTB","CCR7","CCL19"),
  "Regulatory B"          = c("CD274","IL10","TGFB1")
)

module_order <- tibble(
  Module = c(
    "Inflammatory chemokine",
    "Antigen presentation",
    "IFNγ response",
    
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T",
    "Chemokine receptor",
    
    "Lymphoid organization",
    "Regulatory B"
  ),
  Compartment = c(
    rep("Myeloid", 3),
    rep("T/NK", 5),
    rep("B/Plasma", 2)
  )
)

module_info <- module_order

T_all <- p17 %>%
  filter(grepl(T_R, cell_type, ignore.case = TRUE)) %>%
  select(cell, x, y)

B_all <- p17 %>%
  filter(grepl(B_R, cell_type, ignore.case = TRUE)) %>%
  select(cell, x, y)

tri_xy$T_members <- cells_within_radius(
  tri_xy,
  T_all,
  radius = r_tri
)

tri_xy$B_members <- cells_within_radius(
  tri_xy,
  B_all,
  radius = r_tri
)

tri_nf <- tri_xy %>%
  filter(
    Zone %in% c("Virus-near","Virus-distant")
  )

tri_nf$Zone <- factor(
  tri_nf$Zone,
  levels = c("Virus-distant","Virus-near")
)

mean_over <- function(ids, vec){
  
  ids <- intersect(ids, names(vec))
  
  if(!length(ids))
    return(NA_real_)
  
  mean(vec[ids], na.rm = TRUE)
  
}

forest_list <- list()

for(i in seq_len(nrow(module_info))) {
  
  mod  <- module_info$Module[i]
  comp <- module_info$Compartment[i]
  
  genes <- intersect(
    MODULES_EXT[[mod]],
    rownames(xe_global)
  )
  
  if(length(genes) < 2)
    next
  
  if(comp == "Myeloid") {
    
    score_raw <- module_score_vec(
      xe_global,
      tri_nf$anchor,
      genes
    )
    
    score <- score_raw[tri_nf$anchor]
    
  } else if(comp == "T/NK") {
    
    ids <- unique(
      unlist(tri_nf$T_members)
    )
    
    score_vec <- module_score_vec(
      xe_global,
      ids,
      genes
    )
    
    score <- vapply(
      tri_nf$T_members,
      mean_over,
      numeric(1),
      vec = score_vec
    )
    
  } else {
    
    ids <- unique(
      unlist(tri_nf$B_members)
    )
    
    score_vec <- module_score_vec(
      xe_global,
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
  
  tmp <- tri_nf %>%
    mutate(Score = score) %>%
    filter(is.finite(Score))
  
  if(nrow(tmp) < 20)
    next
  
  tmp$Score_Z <- as.numeric(
    scale(tmp$Score)
  )
  
  near_vals <- tmp$Score_Z[
    tmp$Zone == "Virus-near"
  ]
  
  far_vals <- tmp$Score_Z[
    tmp$Zone == "Virus-distant"
  ]
  
  delta <- mean(
    near_vals,
    na.rm = TRUE
  ) -
    mean(
      far_vals,
      na.rm = TRUE
    )
  
  p_wil <- tryCatch(
    
    wilcox.test(
      near_vals,
      far_vals
    )$p.value,
    
    error = function(e) NA_real_
    
  )
  
  forest_list[[length(forest_list) + 1]] <- tibble(
    
    Module = mod,
    
    Compartment = comp,
    
    Delta = delta,
    
    P = p_wil,
    
    N_Near = length(near_vals),
    
    N_Far = length(far_vals)
    
  )
  
}

forest_tbl <- bind_rows(forest_list) %>%
  
  mutate(
    
    FDR = p.adjust(
      P,
      method = "BH"
    ),
    
    Stars = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
    
  ) %>%
  
  arrange(FDR)

write_csv(
  forest_tbl,
  "../../output/Fig5_Virus_Triad_Modules.csv"
)

module_dot <- forest_tbl %>%
  
  dplyr::select(
    Module,
    Compartment,
    Delta,
    FDR,
    Stars
  )

module_dot$Compartment <- factor(
  module_dot$Compartment,
  levels = c(
    "Myeloid",
    "T/NK",
    "B/Plasma"
  )
)

module_dot$Module <- factor(
  module_dot$Module,
  levels = rev(
    module_order$Module
  )
)

delta_lim <- max(
  abs(module_dot$Delta),
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
      fill = Delta
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
      -delta_lim,
      delta_lim
    ),
    name = "Mean module\nscore\n(Virus-near\nvs\nVirus-distant)"
    ) +
  
  scale_x_continuous(
    limits = c(
      0.85,
      1.18
    ),
    breaks = NULL
  ) +
  
  facet_grid(
    Compartment ~ .,
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
  "../../plots/Fig5_Triad_Module_Virus_Distance.pdf",
  p_triad_modules,
  1.2 * P_Width,
  P_Height
)

# ============================================================
# Triad enrichment near virus
# (% myeloid cells that are triad anchors)
# ============================================================

# Distance of each myeloid cell to nearest virus+ cell
myeloid_zone_df <- p17 %>%
  filter(
    grepl(
      MY_R,
      cell_type,
      ignore.case = TRUE
    )
  ) %>%
  mutate(
    DistToVirus = dist_to_set(
      .,
      virus_xy
    ),
    
    Virus_Proximity = case_when(
      DistToVirus <= NEAR_MAX ~ "Virus-near",
      DistToVirus >= FAR_MIN  ~ "Virus-distant",
      TRUE                    ~ NA_character_
    ),
    
    Triad = cell %in% triad_anchors
  ) %>%
  filter(
    !is.na(Virus_Proximity)
  )

# Fisher test
ctable <- table(
  myeloid_zone_df$Virus_Proximity,
  myeloid_zone_df$Triad
)

ft <- fisher.test(ctable)

or_result <- tibble(
  Comparison = "Virus-near",
  OR      = as.numeric(ft$estimate),
  CI_Low  = ft$conf.int[1],
  CI_High = ft$conf.int[2],
  P       = ft$p.value
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

write_csv(
  or_result,
  "../../output/Fig5_Triad_Enrichment_OR.csv"
)

# Prevalence table
triad_prevalence <- myeloid_zone_df %>%
  group_by(Virus_Proximity) %>%
  summarise(
    Fraction_Triad = mean(Triad),
    Triad_Cells    = sum(Triad),
    Total_Cells    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Virus_Proximity = factor(
      Virus_Proximity,
      levels = c(
        "Virus-distant",
        "Virus-near"
      )
    )
  )

print(triad_prevalence)

write_csv(
  triad_prevalence,
  "../../output/Fig5_Triad_Prevalence.csv"
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
      label = scales::percent(
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
      "Virus-near"    = "blue"
    )
  ) +
  
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = max(triad_prevalence$Fraction_Triad) * 1.3,
    
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
    y = max(triad_prevalence$Fraction_Triad) * 1.1,
    label = or_result$Sig,
    size = pt_size / 1.4,
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(
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
    legend.position = "none"
  )

save_pdf(
  "../../plots/Fig5_Triad_Enrichment_Near_Virus.pdf",
  p_triad_enrichment,
  0.8 * P_Width,
  P_Height
)