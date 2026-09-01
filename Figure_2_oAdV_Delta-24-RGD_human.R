# ==============================================================================
# Figures 2 and S2 analyses using Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts
#
# Required input files: 
# ../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds
# ../../NCT00805376_rGBM_oAdV/analysis_output/TCR_Diversity_Groups.rds
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
  library(SeuratObject)
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(ggrepel)
  library(ggrastr)
  
  library(readr)
  library(RANN)
  library(survival)
  
})

set.seed(42)

dir.create("../../plots",  showWarnings = FALSE)
dir.create("../../output", showWarnings = FALSE)

# ============================================================
# Figure settings
# ============================================================

mm_to_in <- function(mm) mm / 25.4

P_Width  <- 58
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

# ============================================================
# Parameters
# ============================================================

r_tri <- 30

MY_R   <- "Myeloid"
B_R    <- "B/Plasma"
T_R    <- "T/NK"
VASC_R <- "Vasculature/Stroma"

immune_cell_types <- c(
  "Myeloid",
  "T/NK Cells",
  "B/Plasma Cells"
)

# ============================================================
# Helper functions
# ============================================================

triad_flags_myeloid <- function(df, r = r_tri) {
  
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
  
  if (
    nrow(A) == 0 ||
    nrow(B) == 0 ||
    nrow(T) == 0
  ) {
    return(
      data.frame(
        cell  = character(),
        Triad = logical()
      )
    )
  }
  
  dB <- nn2(
    B[, c("x", "y")],
    A[, c("x", "y")],
    k = 1
  )$nn.dists[, 1]
  
  dT <- nn2(
    T[, c("x", "y")],
    A[, c("x", "y")],
    k = 1
  )$nn.dists[, 1]
  
  data.frame(
    cell  = A$cell,
    Triad = (dB <= r) & (dT <= r)
  )
  
}

triad_fraction_by_radius <- function(meta_df, r) {
  
  meta_df |>
    group_by(Patient_Match, FOV) |>
    group_modify(function(df, key) {
      
      tri <- triad_flags_myeloid(
        df,
        r = r
      )
      
      data.frame(
        n_myeloid = nrow(tri),
        n_triad   = sum(tri$Triad)
      )
      
    }) |>
    ungroup() |>
    group_by(Patient_Match) |>
    summarise(
      n_myeloid = sum(n_myeloid),
      n_triad   = sum(n_triad),
      triad_frac = n_triad / n_myeloid,
      .groups = "drop"
    ) |>
    filter(
      is.finite(triad_frac)
    )
  
}

# ============================================================
# Colors
# ============================================================

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

cols_div <- c(
  "Low"  = "#9ECAE1",
  "High" = "#08306B"
)

km_cols <- c(
  "Triad Low"  = "#9ECAE1",
  "Triad High" = "#08306B"
)

# ============================================================
# Data loading
# ============================================================

xe_global <- readRDS(
  "../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds"
)

xe_global <- UpdateSeuratObject(xe_global)

div_data <- readRDS(
  "../../NCT00805376_rGBM_oAdV/analysis_output/TCR_Diversity_Groups.rds"
)

clinical_df <- readr::read_csv(
  "../../NCT00805376_rGBM_oAdV/clinical_data.csv",
  show_col_types = FALSE
)

DefaultAssay(xe_global) <- "Xenium"

# ============================================================
# Patient annotations
# ============================================================

xe_global$sample_id <- str_extract(
  xe_global$orig.ident,
  "P[0-9]+-(post|pre)"
)

xe_global$Patient_Match <- str_extract(
  xe_global$sample_id,
  "P[0-9]+"
)

xe_global$Timepoint <- ifelse(
  grepl("post", xe_global$sample_id, ignore.case = TRUE),
  "post",
  "pre"
)

atlas_nums <- str_extract(
  xe_global$Patient_Match,
  "[0-9]+"
)

pid_col <- names(div_data)[
  sapply(div_data, function(x)
    any(grepl("[0-9]+", as.character(x))))
][1]

div_nums <- str_extract(
  div_data[[pid_col]],
  "[0-9]+"
)

xe_global$Div_Group <- div_data$Div_Group[
  match(atlas_nums, div_nums)
]

xe_post <- subset(
  xe_global,
  subset = Timepoint == "post" & !is.na(Div_Group)
)

# ============================================================
# Survival metadata
# ============================================================
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

# ============================================================
# Build metadata with coordinates
# ============================================================

raw_coords <- do.call(
  rbind,
  lapply(Images(xe_post), function(fov) {
    
    GetTissueCoordinates(xe_post[[fov]]) |>
      mutate(FOV = fov)
    
  })
)

meta <- xe_post@meta.data |>
  rownames_to_column("cell") |>
  left_join(raw_coords, by = "cell")


# ============================================================
# Metadata with coordinates (all samples)
# ============================================================

raw_coords_all <- do.call(
  rbind,
  lapply(Images(xe_global), function(fov) {
    
    GetTissueCoordinates(xe_global[[fov]]) |>
      mutate(FOV = fov)
    
  })
)

meta_all <- xe_global@meta.data |>
  rownames_to_column("cell") |>
  left_join(raw_coords_all, by = "cell")

# ============================================================
# Triad anchors
# ============================================================

triads <- meta |>
  group_by(sample_id, FOV) |>
  group_modify(
    ~triad_flags_myeloid(.x)
  ) |>
  ungroup()

meta <- meta |>
  left_join(
    select(triads, cell, Triad),
    by = "cell"
  )

meta$Triad[is.na(meta$Triad)] <- FALSE

# ============================================================
# UMAP
# ============================================================

red_name <- intersect(
  c("full.umap", "umap"),
  Reductions(xe_global)
)[1]

stopifnot(!is.na(red_name))

emb <- Embeddings(xe_global, red_name)

umap_df <- data.frame(
  UMAP_1        = emb[, 1],
  UMAP_2        = emb[, 2],
  cell_type     = xe_global$cell_type,
  Timepoint     = xe_global$Timepoint,
  Div_Group     = xe_global$Div_Group,
  Patient_Match = xe_global$Patient_Match
) |>
  filter(
    Timepoint == "post",
    !is.na(Div_Group)
  )

p_umap <- ggplot(
  umap_df,
  aes(UMAP_1, UMAP_2, color = cell_type)
) +
  geom_point_rast(
    size = 0.32,
    raster.dpi = 600
  ) +
  scale_color_manual(values = pal_best) +
  theme_fig() +
  guides(color = "none") +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  )

save_pdf(
  "../../plots/Fig2_UMAP.pdf",
  p_umap,
  P_Width,
  P_Height
)

# ============================================================
# Cell-type composition
# ============================================================

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

comp_df <- meta |>
  count(Patient_Match, cell_type, name = "n") |>
  group_by(Patient_Match) |>
  mutate(Freq = n / sum(n)) |>
  ungroup()

comp_df$cell_type <- factor(
  comp_df$cell_type,
  levels = rev(cell_order)
)

p_comp <- ggplot(
  comp_df,
  aes(Patient_Match, Freq, fill = cell_type)
) +
  geom_col(width = 0.8, position = "fill") +
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
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "Proportion"
  ) +
  theme(
    legend.key.height = unit(0.18, "cm"),
    legend.key.width  = unit(0.18, "cm")
  )

save_pdf(
  "../../plots/Fig2_Composition.pdf",
  p_comp,
  1.3*P_Width,
  P_Height
)

# ============================================================
# Immune abundance vs TCR diversity
# ============================================================

# Immune cell-type frequencies
freq_by_type <- meta |>
  count(Patient_Match, cell_type, name = "n") |>
  group_by(Patient_Match) |>
  mutate(frac = n / sum(n)) |>
  ungroup() |>
  select(Patient_Match, cell_type, frac)

# TCR diversity metadata
div_clean <- div_data |>
  mutate(
    Patient_Match = paste0(
      "P",
      str_extract(.data[[pid_col]], "[0-9]+")
    )
  ) |>
  select(Patient_Match, TCR_Div, Div_Group) |>
  filter(is.finite(TCR_Div)) |>
  distinct(Patient_Match, .keep_all = TRUE)

# Join abundance and diversity data
freq_div <- freq_by_type |>
  filter(cell_type %in% immune_cell_types) |>
  inner_join(div_clean, by = "Patient_Match") |>
  distinct(Patient_Match, cell_type, .keep_all = TRUE)

# Spearman correlation sweep
cor_results_div <- split(freq_div, freq_div$cell_type) |>
  lapply(function(df) {
    
    df2 <- df |>
      filter(
        is.finite(frac),
        is.finite(TCR_Div)
      ) |>
      distinct(Patient_Match, .keep_all = TRUE)
    
    n_pat <- n_distinct(df2$Patient_Match)
    
    if (
      n_pat < 5 ||
      length(unique(df2$frac)) < 2 ||
      length(unique(df2$TCR_Div)) < 2
    ) {
      return(
        data.frame(
          cell_type = df$cell_type[1],
          n         = n_pat,
          rho       = NA_real_,
          p         = NA_real_
        )
      )
    }
    
    ct <- suppressWarnings(
      cor.test(
        df2$frac,
        df2$TCR_Div,
        method = "spearman",
        exact  = FALSE
      )
    )
    
    data.frame(
      cell_type = df$cell_type[1],
      n         = n_pat,
      rho       = unname(ct$estimate),
      p         = unname(ct$p.value)
    )
    
  }) |>
  bind_rows() |>
  filter(!is.na(rho)) |>
  arrange(desc(abs(rho)))

# Plot data
df_bar_div <- cor_results_div |>
  mutate(
    cell_type = factor(
      cell_type,
      levels = rev(cor_results_div$cell_type)
    ),
    rho_lab = sprintf("%.2f", rho)
  )

# Plot
p_immune_div <- ggplot(
  df_bar_div,
  aes(rho, cell_type, fill = cell_type)
) +
  geom_col(
    width = 0.7,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = rho_lab),
    hjust = -0.1,
    size = 1.8
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.1, 0.3))
  ) +
  scale_fill_manual(values = pal_best) +
  theme_fig() +
  theme(
    plot.margin = margin(5, 20, 5, 5)
  ) +
  labs(
    x = expression(
        paste(
          "Correlation with TCR Diversity (",
          rho,
          ")"
        )
    ),
    y = NULL
  )

save_pdf(
  "../../plots/Fig2_ImmuneAbundance_vs_Diversity.pdf",
  p_immune_div,
  0.8*P_Width,
  P_Height
)

# ============================================================
# Triad fraction vs OS
# ============================================================

df_tri_M <- meta |>
  filter(grepl(MY_R, cell_type, ignore.case = TRUE)) |>
  group_by(Patient_Match) |>
  summarise(
    n_myeloid  = n(),
    n_triad    = sum(Triad),
    triad_frac = mean(Triad),
    .groups = "drop"
  ) |>
  inner_join(
    distinct(meta, Patient_Match, Div_Group),
    by = "Patient_Match"
  ) |>
  inner_join(
    os_df,
    by = "Patient_Match"
  )

ct_tri <- suppressWarnings(
  cor.test(
    df_tri_M$triad_frac,
    df_tri_M$OS_months,
    method = "spearman"
  )
)

x_lab <- min(
  df_tri_M$triad_frac,
  na.rm = TRUE
)

y_lo <- min(
  df_tri_M$OS_months,
  na.rm = TRUE
)

y_hi <- max(
  df_tri_M$OS_months,
  na.rm = TRUE
)

y_lab <- y_hi +
  (y_hi - y_lo) * 0.15

lab_df_sc <- data.frame(
  x = x_lab,
  y = y_lab,
  lab = paste0(
    "atop(rho == ",
    sprintf("%.2f", ct_tri$estimate),
    ", ",
    format_p_plotmath(ct_tri$p.value),
    ")"
  )
)

p_tri_scat <- ggplot(
  df_tri_M,
  aes(triad_frac, OS_months)
) +
  geom_smooth(
    method = "lm",
    color = "black",
    linetype = "dashed",
    se = FALSE,
    linewidth = 0.25
  ) +
  geom_point(
    aes(color = Div_Group),
    size = 1.8
  ) +
  geom_text(
    data = lab_df_sc,
    inherit.aes = FALSE,
    aes(
      x = x,
      y = y,
      label = lab
    ),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = pt_size / 2.2
  ) +
  scale_color_manual(
    values = cols_div,
    name = "TCR diversity"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_fig() +
  labs(
    x = "Triad fraction (Myeloid-T/NK-B/Plasma)",
    y = "OS (months)"
  )

save_pdf(
  "../../plots/Fig2_Triad_OS_scatter.pdf",
  p_tri_scat,
  1.3 * P_Width,
  P_Height
)

# ============================================================
# Triad fraction vs TCR diversity
# ============================================================

df_triad_div <- df_tri_M |>
  select(Patient_Match, triad_frac) |>
  inner_join(div_clean, by = "Patient_Match")

if (nrow(df_triad_div) >= 3) {
  
  ct_tr <- suppressWarnings(
    cor.test(
      df_triad_div$triad_frac,
      df_triad_div$TCR_Div,
      method = "spearman"
    )
  )
  
  x_lab_tr <- min(
    df_triad_div$triad_frac,
    na.rm = TRUE
  )
  
  y_lo_tr <- min(
    df_triad_div$TCR_Div,
    na.rm = TRUE
  )
  
  y_hi_tr <- max(
    df_triad_div$TCR_Div,
    na.rm = TRUE
  )
  
  y_lab_tr <- y_hi_tr +
    (y_hi_tr - y_lo_tr) * 0.15
  
  lab_tr_sc <- data.frame(
    x = x_lab_tr,
    y = y_lab_tr,
    lab = paste0(
      "atop(rho == ",
      sprintf("%.2f", ct_tr$estimate),
      ", ",
      format_p_plotmath(ct_tr$p.value),
      ")"
    )
  )
  
  p_triad_div <- ggplot(
    df_triad_div,
    aes(triad_frac, TCR_Div)
  ) +
    geom_smooth(
      method = "lm",
      color = "black",
      linetype = "dashed",
      se = FALSE,
      linewidth = 0.25
    ) +
    geom_point(
      aes(color = Div_Group),
      size = 1.8
    ) +
    geom_text(
      data = lab_tr_sc,
      inherit.aes = FALSE,
      aes(
        x = x,
        y = y,
        label = lab
      ),
      parse = TRUE,
      hjust = 0,
      vjust = 1,
      size = pt_size / 2.2
    ) +
    scale_color_manual(
      values = cols_div,
      name = "TCR diversity"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.25))
    ) +
    theme_fig() +
    labs(
      x = "Triad fraction (Myeloid-T/NK-B/Plasma)",
      y = "TCR Diversity (Inverse Simpson)"
    )
  
  save_pdf(
    "../../plots/Fig2_Triad_Diversity_scatter.pdf",
    p_triad_div,
    1.3 * P_Width,
    P_Height
  )
  
}

# ============================================================
# Kaplan-Meier analysis
# ============================================================

df_tri_km <- df_tri_M |>
  mutate(
    Group = ifelse(
      triad_frac >= median(triad_frac, na.rm = TRUE),
      "Triad High",
      "Triad Low"
    )
  ) |>
  mutate(
    Group = factor(
      Group,
      levels = c("Triad Low", "Triad High")
    )
  )

sf <- survfit(
  Surv(OS_months, rep(1L, nrow(df_tri_km))) ~ Group,
  data = df_tri_km
)

# KM curve data
make_km_tables <- function(sf_obj, group_levels) {
  
  ss <- summary(sf_obj)
  
  if (!is.null(ss$strata)) {
    
    grp_chr <- sub(
      "^Group=",
      "",
      as.character(ss$strata)
    )
    
  } else if (!is.null(sf_obj$strata)) {
    
    grp_chr <- rep(
      sub("^Group=", "", names(sf_obj$strata)[1]),
      length(ss$time)
    )
    
  } else {
    
    grp_chr <- rep(
      group_levels[1],
      length(ss$time)
    )
    
  }
  
  km_df <- data.frame(
    time  = ss$time,
    surv  = ss$surv,
    lower = pmax(pmin(ss$lower, 1), 0),
    upper = pmax(pmin(ss$upper, 1), 0),
    group = factor(
      grp_chr,
      levels = group_levels
    )
  )
  
  do.call(
    rbind,
    lapply(
      split(km_df, km_df$group, drop = TRUE),
      function(d) {
        
        d <- d[order(d$time), , drop = FALSE]
        
        rbind(
          data.frame(
            time  = 0,
            surv  = 1,
            lower = 1,
            upper = 1,
            group = d$group[1]
          ),
          d
        )
        
      }
    )
  )
  
}

km_pts <- make_km_tables(
  sf,
  levels(df_tri_km$Group)
)

# Stepwise CI polygons
build_ci_polygons <- function(km_df) {
  
  polys <- lapply(
    split(km_df, km_df$group, drop = TRUE),
    function(d) {
      
      d <- d[order(d$time), , drop = FALSE]
      
      if (nrow(d) < 2) {
        return(NULL)
      }
      
      t  <- d$time
      up <- d$upper
      lo <- d$lower
      
      x_top <- as.numeric(
        rbind(
          t[-length(t)],
          t[-1]
        )
      )
      
      y_top <- rep(
        up[-length(up)],
        each = 2
      )
      
      x_bot <- rev(x_top)
      
      y_bot <- rev(
        rep(
          lo[-length(lo)],
          each = 2
        )
      )
      
      data.frame(
        x = c(x_top, x_bot),
        y = c(y_top, y_bot),
        group = d$group[1]
      )
      
    }
  )
  
  do.call(rbind, polys)
  
}

ci_poly <- build_ci_polygons(km_pts)

# Log-rank test
sdif <- survdiff(
  Surv(OS_months, rep(1L, nrow(df_tri_km))) ~ Group,
  data = df_tri_km
)

p_lr <- 1 - pchisq(
  sdif$chisq,
  length(sdif$n) - 1
)

lbl_km <- format_p_plotmath(p_lr)

max_t <- max(
  km_pts$time,
  na.rm = TRUE
)

# Plot
p_km <- ggplot() +
  geom_polygon(
    data = ci_poly,
    aes(
      x = x,
      y = y,
      fill = group,
      group = group
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_step(
    data = km_pts,
    aes(
      x = time,
      y = surv,
      color = group,
      group = group
    ),
    linewidth = 0.6
  ) +
  ggplot2::annotate(
    "text",
    x = max_t * 0.02,
    y = 0.05,
    label = lbl_km,
    parse = TRUE,
    hjust = 0,
    vjust = 0,
    size = pt_size / 2
  ) +
  scale_color_manual(
    values = km_cols,
    name = "Triad group",
    drop = FALSE
  ) +
  scale_fill_manual(
    values = km_cols,
    name = "Triad group",
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = c(0, max_t * 1.05),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  theme_fig() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = "Time (months)",
    y = "Survival probability"
  )

save_pdf(
  "../../plots/Fig2_Triad_KM.pdf",
  p_km,
  1.3 * P_Width,
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
    "SLC17A7"
  ),
  "MES-like/Hypoxic Tumor" = c(
    "LOX", "VEGFA", "CD44", "VCAN"
  ),
  "Astrocytes/AC-like" = c(
    "AQP4", "SLC1A3"
  ),
  "Oligodendrocytes/OPC-like" = c(
    "MOG", "OLIG1", "OLIG2"
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
cells_use <- rownames(xe_global@meta.data)[
  xe_global$cell_type %in% cell_order
]

expr_df <- FetchData(
  xe_global,
  vars = gene_order,
  cells = cells_use
)

expr_df$cell_type <- xe_global$cell_type[
  rownames(expr_df)
]

# Per-cell-type summaries
dp <- bind_rows(
  
  lapply(
    cell_order,
    function(ct) {
      
      df_ct <- expr_df[
        expr_df$cell_type == ct,
        ,
        drop = FALSE
      ]
      
      bind_rows(
        
        lapply(
          gene_order,
          function(g) {
            
            data.frame(
              cell_type = ct,
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

dp$cell_type <- factor(
  dp$cell_type,
  levels = cell_order
)

# Plot
p_dot <- ggplot(
  dp,
  aes(gene, cell_type)
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
  "../../plots/FigS2_DotPlot_CellTypeMarkers.pdf",
  p_dot,
  2 * P_Width,
  P_Height
)

# ============================================================
# Cell-type composition across all samples
# ============================================================

composition_facet_levels <- c("pre","post") 

meta_all$Timepoint <- factor(
  meta_all$Timepoint,
  levels = composition_facet_levels
)

comp_all <- meta_all |>
  count(
    Patient_Match,
    Timepoint,
    cell_type,
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

comp_all$cell_type <- factor(
  comp_all$cell_type,
  levels = rev(cell_order)
)

p_comp_all <- ggplot(
  comp_all,
  aes(Patient_Match, Freq, fill = cell_type)
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
    legend.position = "right",
    legend.key.height = unit(0.18, "cm"),
    legend.key.width  = unit(0.18, "cm")
  ) +
  labs(
    x = NULL,
    y = "Proportion"
  )

save_pdf(
  "../../plots/FigS2_Composition_AllSamples.pdf",
  p_comp_all,
  2 * P_Width,
  P_Height
)

# ============================================================
# Triad radius sensitivity analysis
# ============================================================

radii <- seq(10, 100, by = 10)

tri_df <- bind_rows(
  
  lapply(
    radii,
    function(rr) {
      
      df_rr <- triad_fraction_by_radius(
        meta,
        r = rr
      ) |>
        inner_join(
          os_df,
          by = "Patient_Match"
        )
      
      data.frame(
        radius = rr,
        rho = suppressWarnings(
          cor(
            df_rr$triad_frac,
            df_rr$OS_months,
            method = "spearman",
            use = "complete.obs"
          )
        )
      )
      
    }
  )
  
)

p_radius <- ggplot(
  tri_df,
  aes(radius, rho)
) +
  geom_line(
    linewidth = 0.6,
    color = "#08306B"
  ) +
  geom_point(
    size = 1.5,
    color = "#08306B"
  ) +
  geom_vline(
    xintercept = r_tri,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  theme_fig() +
  labs(
    x = "Triad radius (µm)",
    y = expression(
      paste(
        "Spearman ",
        rho,
        " with OS"
      )
    )
  )

save_pdf(
  "../../plots/FigS2_Triad_Radius_Sweep.pdf",
  p_radius,
  P_Width,
  P_Height
)

# ============================================================
# Triad fraction across radii
# ============================================================

tri_frac_df <- bind_rows(
  
  lapply(
    radii,
    function(rr) {
      
      triad_fraction_by_radius(
        meta,
        r = rr
      ) |>
        summarise(
          radius = rr,
          triad_frac = median(
            triad_frac,
            na.rm = TRUE
          )
        )
      
    }
  )
  
)

p_frac <- ggplot(
  tri_frac_df,
  aes(radius, triad_frac)
) +
  geom_line(
    linewidth = 0.6,
    color = "#D55E00"
  ) +
  geom_point(
    size = 1.5,
    color = "#D55E00"
  ) +
  geom_vline(
    xintercept = r_tri,
    linetype = "dashed",
    color = "red",
    linewidth = 0.4
  ) +
  theme_fig() +
  labs(
    x = "Triad radius (µm)",
    y = "Median triad fraction"
  )

save_pdf(
  "../../plots/FigS2_Triad_Fraction_Radius_Sweep.pdf",
  p_frac,
  P_Width,
  P_Height
)



