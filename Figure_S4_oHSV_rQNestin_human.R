# ==============================================================================
# Figure S4 analyses using rQNestin34.5v.2 data
# ==============================================================================

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
  library(ggrastr)
  
})

set.seed(34)
dir.create("../plots", showWarnings = FALSE)
dir.create("../output", showWarnings = FALSE)

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

oHSV_xenium <- readRDS(
  "../NCT03152318_rGBM_oHSV/analysis_output/oHSV_Global_with_Harmonized_Labels.rds"
)

oHSV_xenium <- UpdateSeuratObject(oHSV_xenium)

DefaultAssay(oHSV_xenium) <- "Xenium"

# ============================================================
# Build metadata with coordinates
# ============================================================

raw_coords <- do.call(
  rbind,
  lapply(Images(oHSV_xenium), function(fov) {
    
    GetTissueCoordinates(oHSV_xenium[[fov]]) |>
      dplyr::mutate(FOV = fov)
    
  })
)

meta <- oHSV_xenium@meta.data |>
  tibble::rownames_to_column("cell") |>
  dplyr::left_join(raw_coords, by = "cell")

# Patient/sample identifiers
meta$sample_id <- stringr::str_extract(
  meta$orig.ident,
  "P[0-9]+_(POST|PRE)"
)

meta$Patient_Match <- stringr::str_extract(
  meta$sample_id,
  "P[0-9]+"
)

meta$Timepoint <- ifelse(
  grepl("POST", meta$sample_id, ignore.case = TRUE),
  "POST",
  "PRE"
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
    grepl(MY_R, final_harmonized_label, ignore.case = TRUE)
  )
  
  B <- subset(
    df,
    grepl(B_R, final_harmonized_label, ignore.case = TRUE)
  )
  
  T <- subset(
    df,
    grepl(T_R, final_harmonized_label, ignore.case = TRUE)
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
      final_harmonized_label,
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
# PRE VS POST TRIAD FREQUENCY
# ============================================================

meta$Timepoint <- dplyr::case_when(
  grepl("_PRE$",  meta$sample_id) ~ "PRE",
  grepl("_POST$", meta$sample_id) ~ "POST",
  TRUE ~ NA_character_
)

# Sample-level triad fraction
triad_summary <- meta %>%
  dplyr::filter(
    grepl(MY_R, final_harmonized_label, ignore.case = TRUE)
  ) %>%
  dplyr::group_by(
    Patient_Match,
    sample_id,
    Timepoint
  ) %>%
  dplyr::summarise(
    n_myeloid = n(),
    n_triad   = sum(Triad),
    triad_frac = n_triad / n_myeloid,
    .groups = "drop"
  )

write.csv(
  triad_summary,
  "../output/Fig3_TriadFraction_BySample.csv",
  row.names = FALSE
)

# Matched PRE/POST pairs
triad_wide <- triad_summary %>%
  dplyr::select(
    Patient_Match,
    Timepoint,
    triad_frac
  ) %>%
  tidyr::pivot_wider(
    names_from  = Timepoint,
    values_from = triad_frac
  ) %>%
  dplyr::filter(
    !is.na(PRE),
    !is.na(POST)
  ) %>%
  dplyr::mutate(
    Delta = POST - PRE
  )

cat(
  "\nMatched PRE/POST pairs:",
  nrow(triad_wide),
  "\n"
)

print(triad_wide)

# Paired Wilcoxon test
p_triad <- wilcox.test(
  triad_wide$PRE,
  triad_wide$POST,
  paired = TRUE,
  exact = FALSE
)$p.value

cat(
  "\nPaired Wilcoxon P = ",
  signif(p_triad, 3),
  "\n"
)

# Plot
plot_df <- triad_wide %>%
  tidyr::pivot_longer(
    cols = c(PRE, POST),
    names_to = "Timepoint",
    values_to = "TriadFraction"
  )

plot_df$Timepoint <- factor(
  plot_df$Timepoint,
  levels = c("PRE", "POST")
)

p_triad_pair <- ggplot(
  plot_df,
  aes(
    x = Timepoint,
    y = TriadFraction
  )
) +
  geom_boxplot(
    aes(fill = Timepoint),
    outlier.shape = NA,
    width = 0.5,
    alpha = 0.85
  ) +
  geom_line(
    aes(group = Patient_Match),
    color = "grey50",
    alpha = 0.6
  ) +
  geom_point(
    aes(fill = Timepoint),
    shape = 21,
    size = 1.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "PRE"  = "#E41A1C",
      "POST" = "#377EB8"
    )
  ) +
  theme_fig() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "none"
  ) +
  labs(
    title = "Meylan et al. oHSV dataset",
    x = NULL,
    y = "Triad fraction among myeloid cells"
  )

y_annot <- max(
  plot_df$TriadFraction,
  na.rm = TRUE
) * 1.05

p_triad_pair <- p_triad_pair +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = y_annot,
    label = format_p_plotmath(p_triad),
    parse = TRUE,
    size = pt_size / 2.2
  )

save_pdf(
  "../plots/FigS4_oHSV_PRE_POST_TriadFraction.pdf",
  p_triad_pair,
  P_Width,
  P_Height
)

# ============================================================
# Distance of myeloid cells to vasculature/stroma
# ============================================================

df_my_dist <- meta |>
  dplyr::filter(
    grepl(MY_R, final_harmonized_label, ignore.case = TRUE)
  ) |>
  dplyr::group_by(sample_id, FOV) |>
  dplyr::group_modify(~{
    
    vas <- meta |>
      dplyr::filter(
        sample_id == .y$sample_id,
        FOV == .y$FOV,
        grepl(VASC_R, final_harmonized_label, ignore.case = TRUE)
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "none"
  ) +
  labs(
    title = "Meylan et al. oHSV dataset",
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
  "../plots/FigS4_oHSV_Myeloid_Distance_to_Vasculature.pdf",
  p_struct,
  P_Width,
  P_Height
)

# ============================================================
# Triad-associated immune programs
# ============================================================

# Use linear mixed-effects models (LMMs) to evaluate triad
# associations, with triad status as a fixed effect and
# FOV nested within sample as random intercepts.
#
# Null hypothesis (H0):
# β = 0
#
# Cells inside triads have the same average module score as
# cells outside triads.
#
# Statistical testing is performed on z-scored module scores.
# Percent change is calculated from raw module scores for
# biological interpretation and supplementary reporting.

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
  rownames(oHSV_xenium)
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
        final_harmonized_label,
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
    oHSV_xenium,
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
  "../output/FigS4_oHSV_Triad_Module_Statistics.csv",
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
    title = "Meylan et al. oHSV dataset",
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),
    
    axis.text.y = element_text(
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
  "../plots/FigS4_oHSV_Triad_Immune_Modules.pdf",
  p_module_summary,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Tumor cells
# ============================================================

# Use linear mixed-effects models (LMMs) to evaluate
# triad-associated tumor programs.
#
# Null hypothesis (H0):
# β = 0
#
# Tumor cells in regions with higher triad density have the
# same average module score as tumor cells in regions with
# lower triad density.
#
# Statistical testing is performed on z-scored module scores.
# Triad density is calculated using a Gaussian kernel density
# estimate and standardized prior to modeling.

tumor_meta <- meta %>%
  dplyr::filter(
    final_harmonized_label %in% TUMOR_LABELS
  )

tumor_meta$Tumor_State <- factor(
  tumor_meta$final_harmonized_label,
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
  rownames(oHSV_xenium)
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
  "../output/FigS4_oHSV_TriadDensity_Summary.csv",
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
    oHSV_xenium,
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
  "../output/FigS4_oHSV_Tumor_Program_Statistics.csv",
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
    title = "Meylan et al. oHSV dataset",
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),
    
    axis.text.y = element_text(
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
  "../plots/FigS4_oHSV_Triad_Tumor_Modules.pdf",
  p_tumor,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# UMAP
# ============================================================

red_name <- intersect(
  c("full.umap", "umap"),
  Reductions(oHSV_xenium)
)[1]

stopifnot(!is.na(red_name))

emb <- Embeddings(oHSV_xenium, red_name)

umap_df <- data.frame(
  UMAP_1        = emb[, 1],
  UMAP_2        = emb[, 2],
  final_harmonized_label     = oHSV_xenium$final_harmonized_label
)

p_umap <- ggplot(
  umap_df,
  aes(UMAP_1, UMAP_2, color = final_harmonized_label)
) +
  geom_point_rast(
    size = 0.32,
    raster.dpi = 600
  ) +
  scale_color_manual(values = pal_best) +
  theme_fig() +
  guides(color = "none") +
  labs(
    title = "Meylan et al. oHSV dataset",
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    )
  )

save_pdf(
  "../plots/FigS4_oHSV_UMAP.pdf",
  p_umap,
  P_Width,
  P_Height
)

# ============================================================
# Cell-type marker dot plot
# ============================================================

markers <- list(
  "Myeloid" = c(
    "CD68", "ITGAM", "CSF1R", "TREM2"
  ),
  "T/NK Cells" = c(
    "CD3E", "CD8A", "NKG7", "PRF1"
  ),
  "B/Plasma Cells" = c(
    "MS4A1", "CD79A", "MZB1", "JCHAIN", "IGKC"
  ),
  "Vasculature/Stroma" = c(
    "FLT1", "PLVAP", "RGS5", "DCN"
  ),
  "Neurons" = c(
    "DCX", "STMN2", "ENO2", "NELL2"
  ),
  "MES-like/Hypoxic Tumor" = c(
    "CHI3L1", "VCAN", "CD44", "VEGFA"
  ),
  "Astrocytes/AC-like" = c(
    "GFAP", "SLC1A3", "S100B", "SPARCL1"
  ),
  "Oligodendrocytes/OPC-like" = c(
    "OLIG1", "OLIG2", "PDGFRA", "BCAN"
  )
)

cell_order <- c(
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
cells_use <- rownames(oHSV_xenium@meta.data)[
  oHSV_xenium$final_harmonized_label %in% cell_order
]

expr_df <- FetchData(
  oHSV_xenium,
  vars = gene_order,
  cells = cells_use
)

expr_df$final_harmonized_label <- oHSV_xenium$final_harmonized_label[
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text  = element_text(size = 4.5),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm")
  ) +
  labs(
    title = "Meylan et al. oHSV dataset",
    x = NULL,
    y = NULL
  )

save_pdf(
  "../plots/FigS4_oHSV_DotPlot_CellTypeMarkers.pdf",
  p_dot,
  2 * P_Width,
  P_Height
)

# ============================================================
# Cell-type composition across all samples
# ============================================================

composition_facet_levels <- c("PRE","POST") 

meta$Timepoint <- factor(
  meta$Timepoint,
  levels = composition_facet_levels
)

comp_all <- meta |>
  count(
    Patient_Match,
    Timepoint,
    final_harmonized_label,
    name = "n"
  ) |>
  group_by(
    Patient_Match,
    Timepoint
  ) |>
  mutate(
    Freq = n / sum(n)
  ) |>
  ungroup()

comp_all$final_harmonized_label <- factor(
  comp_all$final_harmonized_label,
  levels = rev(cell_order)
)

p_comp_all <- ggplot(
  comp_all,
  aes(Patient_Match, Freq, fill = final_harmonized_label)
) +
  geom_col(
    width = 0.8,
    position = "fill"
  ) +
  facet_grid(
    ~ Timepoint,
    scales = "free_x",
    space = "free"
  ) +
  scale_fill_manual(
    values = pal_best,
    name = "Cell type"
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    legend.position = "right",
    legend.key.height = unit(0.18, "cm"),
    legend.key.width  = unit(0.18, "cm")
  ) +
  labs(
    title = "Meylan et al. oHSV dataset",
    x = NULL,
    y = "Proportion"
  )

save_pdf(
  "../plots/FigS4_oHSV_Composition_AllSamples.pdf",
  p_comp_all,
  2 * P_Width,
  P_Height
)

# ============================================================
# PRE vs POST myeloid transcriptional programs
# ============================================================

PREPOST_MODULES <- list(
  
  "Inflammatory chemokine" = c(
    "CXCL9","CXCL10","CXCL11","CCL5"
  ),
  
  "Antigen presentation" = c(
    "CD74","CTSS","CD80","CD86"
  ),
  
  "IFNγ response" = c(
    "STAT1","IRF1","IFNGR1"
  )
  
)

PREPOST_MODULES <- lapply(
  PREPOST_MODULES,
  intersect,
  rownames(oHSV_xenium)
)

PREPOST_MODULES <- PREPOST_MODULES[
  sapply(PREPOST_MODULES, length) >= 2
]

# Myeloid cells only
myeloid_df <- meta %>%
  dplyr::filter(
    grepl(
      MY_R,
      final_harmonized_label,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::filter(
    !is.na(Timepoint)
  ) %>%
  dplyr::select(
    cell,
    sample_id,
    Patient_Match,
    FOV,
    Timepoint
  )

myeloid_df$Timepoint <- factor(
  myeloid_df$Timepoint,
  levels = c(
    "PRE",
    "POST"
  )
)

# Mixed models
prepost_stats <- list()

prepost_plot_data <- list()

for(mod_name in names(PREPOST_MODULES)) {
  
  message("Processing: ", mod_name)
  
  genes <- PREPOST_MODULES[[mod_name]]
  
  expr_mat <- FetchData(
    oHSV_xenium,
    vars = genes,
    cells = myeloid_df$cell,
    layer = "data"
  )
  
  module_score <- rowMeans(
    expr_mat,
    na.rm = TRUE
  )
  
  if(sd(module_score, na.rm = TRUE) == 0)
    next
  
  df_lmm <- myeloid_df
  
  df_lmm$Module_Score <- module_score
  
  df_lmm$Score_Z <- as.numeric(
    scale(module_score)
  )
  
  prepost_plot_data[[mod_name]] <- df_lmm
  
  fit <- tryCatch(
    
    lmer(
      Score_Z ~
        Timepoint +
        (1 | Patient_Match) +
        (1 | sample_id/FOV),
      data = df_lmm
    ),
    
    error = function(e) NULL
    
  )
  
  if(is.null(fit))
    next
  
  coef_tab <- summary(fit)$coefficients
  
  beta <- coef_tab[
    "TimepointPOST",
    "Estimate"
  ]
  
  se <- coef_tab[
    "TimepointPOST",
    "Std. Error"
  ]
  
  pval <- coef_tab[
    "TimepointPOST",
    "Pr(>|t|)"
  ]
  
  prepost_stats[[mod_name]] <- data.frame(
    
    Module = mod_name,
    
    Beta = beta,
    
    CI_low  = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    
    P = pval,
    
    stringsAsFactors = FALSE
  )
}

# Combine
prepost_stats <- bind_rows(
  prepost_stats
)

prepost_stats$FDR <- p.adjust(
  prepost_stats$P,
  method = "BH"
)

prepost_stats <- prepost_stats %>%
  mutate(
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Percent positive
pct_df <- lapply(
  names(prepost_plot_data),
  function(mod_name){
    
    df <- prepost_plot_data[[mod_name]]
    
    data.frame(
      
      Module = mod_name,
      
      Percent_Positive =
        mean(
          df$Module_Score > 0,
          na.rm = TRUE
        ) * 100
    )
  }
) %>%
  bind_rows()

prepost_dot <- prepost_stats %>%
  left_join(
    pct_df,
    by = "Module"
  )

prepost_dot$Module <- factor(
  prepost_dot$Module,
  levels = rev(c(
    "Inflammatory chemokine",
    "Antigen presentation",
    "IFNγ response"
  ))
)

prepost_dot$Group <- "Myeloid"

# Plot
beta_lim <- max(
  abs(prepost_dot$Beta),
  na.rm = TRUE
)

p_prepost_modules <- ggplot(
  prepost_dot,
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
      -beta_lim,
      beta_lim
    ),
    name = "POST vs PRE β"
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
    title = "Meylan et al. oHSV dataset",
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
    panel.grid = element_blank(),
    
    strip.placement = "outside",
    
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    
    strip.text.y.left = element_text(
      angle = 0
    ),
    
    axis.text.y = element_text(
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
  "../plots/FigS4_oHSV_PRE_POST_Myeloid_Modules.pdf",
  p_prepost_modules,
  1.4 * P_Width,
  P_Height
)

# ============================================================
# Myeloid subtype analysis for triads
# ============================================================

# Use mixed-effects logistic regression (GLMM) to evaluate
# myeloid subtype enrichment within triads, with triad status
# as a fixed effect and FOV nested within sample as random
# intercepts.
#
# Null hypothesis (H0):
# β = 0
#
# The odds of a cell belonging to a given myeloid subtype are
# the same inside and outside triads.
#
# For each subtype, a binary outcome is modeled
# (subtype vs. all other assigned myeloid cells):
#
# logit(P(Subtype)) ~ Triad + (1 | sample_id/FOV)
#
# Odds ratios (ORs) are calculated as exp(β):
#
# OR > 1 = enrichment within triads
# OR < 1 = depletion within triads
# OR = 1 = no association
#
# Statistical significance is assessed using Wald tests and
# corrected for multiple comparisons using the Benjamini-
# Hochberg false discovery rate (FDR).
#
# To assess reproducibility across specimens, sample-level
# consistency is additionally reported as the number of
# samples showing enrichment or depletion in the same
# direction (e.g. 8/11+).

# Define myeloid subtype signatures
SIG_LIST <- list(
  
  Microglia = c(
    "CX3CR1",
    "GPR34",
    "CST3"
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
  oHSV_xenium,
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
        final_harmonized_label,
        Triad,
        sample_id,
        FOV
      ),
    by = "cell"
  ) %>%
  
  filter(
    grepl(
      "Myeloid",
      final_harmonized_label
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
  "../output/FigS4_oHSV_Myeloid_OR_LMM.csv",
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
      x = CI_high + 0.3,
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
    title = "Meylan et al. oHSV dataset",
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
    legend.position = "none",
    plot.margin = margin(
      5,
      35,
      5,
      5
    )
  )

save_pdf(
  "../plots/FigS4_oHSV_Myeloid_OR_Forest.pdf",
  p_myeloid_or,
  P_Width,
  P_Height
)

# Myeloid subtype marker validation
marker_order <- c(
  
  # Microglia
  "CX3CR1",
  "GPR34",
  "CST3",
  
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
  "../output/FigS4_oHSV_Myeloid_Marker_Dotplot_Data.csv",
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
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size + 1
    ),
    
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
    title = "Meylan et al. oHSV dataset",
    x = NULL,
    y = NULL
  )

save_pdf(
  "../plots/FigS4_oHSV_Myeloid_Marker_DotPlot.pdf",
  p_dot,
  1.4 * P_Width,
  P_Height
)
