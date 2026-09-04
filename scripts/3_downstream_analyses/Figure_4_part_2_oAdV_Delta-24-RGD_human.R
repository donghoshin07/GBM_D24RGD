# ==============================================================================
# Figure 4 analysis using Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts/3_downstream_analyses
#
# Required input files: 
# ../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds
# ../../output/probe_table_all_samples_with_SNR.csv
# ../../objects/per_sample_annotated/*.rds
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
  library(patchwork)
  
})

set.seed(42)

dir.create("../../plots", showWarnings = FALSE)
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

save_pdf <- function(
    filename,
    plot_obj,
    w_mm,
    h_mm
) {
  
  ggsave(
    filename,
    plot_obj,
    width = mm_to_in(w_mm),
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

r_tri <- 30

# ============================================================
# Config
# ============================================================

SEURAT_DIR <- "../../objects/per_sample_annotated"

VALID_ALL_CSV <-
  "../../output/probe_table_all_samples_with_SNR.csv"

CACHE_FILE <-
  "../../output/Fig4_Clonotype_APC_Cache.rds"

FORCE_RECOMPUTE <- FALSE
SAVE_CACHE <- TRUE

PROBE_REGEX <- "P\\d+-post-clonotype\\d+"

files <- list.files(
  SEURAT_DIR,
  pattern = "\\.rds$",
  full.names = TRUE
)

# ============================================================
# Data loading
# ============================================================

xe_global <- readRDS(
  "../../NCT00805376_rGBM_oAdV/analysis_output/Global_Atlas_Res0.1.rds"
)

xe_global <- UpdateSeuratObject(
  xe_global
)

DefaultAssay(
  xe_global
) <- "Xenium"

valid_all <- read.csv(
  VALID_ALL_CSV,
  stringsAsFactors = FALSE
)

# ============================================================
# Helper functions
# ============================================================

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

# Clean probes
extract_core_sample <- function(s) {
  
  s2 <- gsub("_", "-", s)
  
  m <- regexpr(
    "P[0-9]+-(pre|post)",
    s2,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  if (m[1] == -1)
    return(NA_character_)
  
  parts <- strsplit(
    regmatches(s2, m)[1],
    "-",
    fixed = TRUE
  )[[1]]
  
  paste0(
    toupper(parts[1]),
    "-",
    tolower(parts[2])
  )
}

# ============================================================
# Build metadata with coordinates
# ============================================================

raw_coords <- do.call(
  
  rbind,
  
  lapply(
    Images(xe_global),
    
    function(fov) {
      
      GetTissueCoordinates(
        xe_global[[fov]]
      ) |>
        
        dplyr::mutate(
          FOV = fov
        )
      
    }
  )
)

meta <- xe_global@meta.data |>
  
  tibble::rownames_to_column(
    "cell"
  ) |>
  
  dplyr::left_join(
    raw_coords,
    by = "cell"
  )

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
    Triad = (
      dB <= r_tri &
        dT <= r_tri
    )
  )
}

# ============================================================
# Identify triad anchors
# ============================================================

triads <- meta |>
  
  dplyr::group_by(
    sample_id,
    FOV
  ) |>
  
  dplyr::group_modify(
    ~triad_flags_myeloid(.x)
  ) |>
  
  dplyr::ungroup()

meta <- meta |>
  
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

triad_anchors <- triads |>
  
  dplyr::filter(
    Triad
  ) |>
  
  dplyr::pull(
    cell
  )

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
# Clonotype probe validation
# ============================================================

p_clonotype_volcano <- ggplot() +
  theme_void()

if (nrow(valid_all) > 0) {
  
  # ----------------------------------------------------------
  # Create one representative row per probe
  # Use the most significant sample/section for each probe
  # ----------------------------------------------------------
  
  volcano_df <- valid_all %>%
    dplyr::group_by(Probe) %>%
    dplyr::slice_min(
      Fisher_q,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    
    dplyr::mutate(
      
      log2SNR =
        log2(
          pmax(SNR, 1e-6)
        ),
      
      neglog10FDR =
        -log10(
          pmax(Fisher_q, 1e-300)
        ),
      
      Validation_Status =
        dplyr::case_when(
          
          Fisher_q > 0.05 ~
            "Not T/NK enriched",
          
          SNR <= 1 ~
            "Low patient specificity",
          
          TRUE ~
            "High confidence TCR probe"
        )
      
    )
  
  volcano_df$Validation_Status <- factor(
    
    volcano_df$Validation_Status,
    
    levels = c(
      "Not T/NK enriched",
      "Low patient specificity",
      "High confidence TCR probe"
    )
  )
  
  message("\n=== Clonotype validation summary ===")
  print(table(volcano_df$Validation_Status))
  
  message(
    "Unique probes plotted: ",
    nrow(volcano_df)
  )
  
  p_clonotype_volcano <-
    
    ggplot() +
    
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      colour = "black",
      linewidth = 0.3
    ) +
    
    # SNR = 1 cutoff
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "grey50",
      linewidth = 0.3
    ) +
    
  # Background probes
  geom_point(
    data =
      volcano_df %>%
      dplyr::filter(
        Validation_Status ==
          "Not T/NK enriched"
      ),
    aes(
      x = log2SNR,
      y = neglog10FDR,
      colour = Validation_Status
    ),
    alpha = 0.8,
    size = 1.4
  ) +
    
  # Significant but low-specificity probes
  geom_point(
    data =
      volcano_df %>%
      dplyr::filter(
        Validation_Status ==
          "Low patient specificity"
      ),
    aes(
      x = log2SNR,
      y = neglog10FDR,
      colour = Validation_Status
    ),
    alpha = 1,
    size = 1.5
  ) +
    
  # High confidence probes
  geom_point(
    data =
      volcano_df %>%
      dplyr::filter(
        Validation_Status ==
          "High confidence TCR probe"
      ),
    aes(
      x = log2SNR,
      y = neglog10FDR,
      colour = Validation_Status
    ),
    alpha = 1,
    size = 1.5
  ) +
    
    scale_colour_manual(
      
      values = c(
        
        "Not T/NK enriched" =
          "grey80",
        
        "Low patient specificity" =
          "grey30",
        
        "High confidence TCR probe" =
          "#009E73"
        
      ),
      
      name = NULL
      
    ) +
    
    labs(
      x = "log2(Specificity Ratio, SNR)",
      y = "-log10(FDR)"
    ) +
    
    theme_fig() +
    
    theme(
      
      legend.position = "bottom",
      
      legend.direction = "vertical",
      
      legend.box = "vertical",
      
      legend.margin = margin(
        0, 0, 0, 0
      ),
      
      legend.box.margin = margin(
        -4, 0, -4, 0
      ),
      
      legend.spacing.y = unit(
        0,
        "mm"
      ),
      
      legend.key.height = unit(
        2,
        "mm"
      ),
      
      plot.margin = margin(
        2, 2, 10, 2
      )
      
    ) +
    
    guides(
      colour =
        guide_legend(
          nrow = 3,
          byrow = TRUE,
          override.aes = list(
            size = 2.5
          )
        )
    )
}

save_pdf(
  "../../plots/Fig4_Clonotype_Validation.pdf",
  p_clonotype_volcano,
  P_Width,
  P_Height
)

write.csv(
  volcano_df,
  "../../output/Fig4_Clonotype_Validation_Volcano.csv",
  row.names = FALSE
)

valid_high <- volcano_df %>%
  filter(
    Fisher_q < 0.05,
    SNR > 1
  )

# ============================================================
# Clonotype tracked cell recovery
# ============================================================

# Build lookup from integrated atlas
# barcode = original Xenium barcode without merge suffix
cell_lookup <- meta %>%
  dplyr::select(
    cell,
    sample_id,
    Patient_Match
  ) %>%
  dplyr::mutate(
    barcode = sub(
      "_[0-9]+$",
      "",
      cell
    )
  )

tracked_df <- list()

for (f in files) {
  
  message("Processing: ", basename(f))
  
  obj <- readRDS(f)
  
  sample_core <- extract_core_sample(
    basename(f)
  )
  
  if (is.na(sample_core))
    next
  
  probes_here <- valid_high %>%
    dplyr::filter(
      startsWith(
        Probe,
        sample_core
      )
    )
  
  if (nrow(probes_here) == 0)
    next
  
  md <- obj@meta.data
  
  t_cells <- rownames(md)[
    grepl(
      T_R,
      md$cell_type,
      ignore.case = TRUE
    )
  ]
  
  if (length(t_cells) == 0)
    next
  
  mat <- GetAssayData(
    obj,
    assay = "Xenium",
    layer = "counts"
  )
  
  for (i in seq_len(nrow(probes_here))) {
    
    pr <- probes_here$Probe[i]
    
    if (!pr %in% rownames(mat))
      next
    
    thr <- 1
    
    vals <- as.numeric(
      mat[pr, t_cells]
    )
    
    names(vals) <- t_cells
    
    pos_cells <- names(vals)[
      vals >= thr
    ]
    
    if (length(pos_cells) == 0)
      next
    
    tmp <- data.frame(
      
      barcode = pos_cells,
      
      sample_id = sample_core,
      
      Probe = pr,
      
      stringsAsFactors = FALSE
      
    )
    
    tracked_df[[length(tracked_df) + 1]] <- tmp
  }
}

tracked_df <- bind_rows(
  tracked_df
)

# Map back into atlas cell IDs
tracked_df <- tracked_df %>%
  
  dplyr::left_join(
    
    cell_lookup,
    
    by = c(
      "sample_id",
      "barcode"
    )
    
  )

tracked_df <- tracked_df %>%
  dplyr::filter(
    !is.na(cell)
  )

tracked_df <- unique(
  tracked_df
)

# Print results
message(
  "\nValidated probes: ",
  n_distinct(valid_high$Probe)
)

message(
  "\nTracked clonotype cells recovered: ",
  nrow(tracked_df)
)

message(
  "\nPatients represented: ",
  n_distinct(
    sub(
      "-post-clonotype.*",
      "",
      tracked_df$Probe
    )
  )
)

# Summary table
probe_summary <- tracked_df %>%
  
  dplyr::count(
    Probe,
    name = "N_Cells"
  ) %>%
  
  dplyr::arrange(
    desc(N_Cells)
  )

write.csv(
  probe_summary,
  "../../output/Fig4_Tracked_Cells_Per_Probe_Summary.csv",
  row.names = FALSE
)

# Define tracked cells
meta$Tracked_Clonotype <-
  meta$cell %in%
  tracked_df$cell

# Rebuild T/NK metadata
t_meta <- meta %>%
  
  dplyr::filter(
    grepl(
      T_R,
      cell_type,
      ignore.case = TRUE
    )
  )

# QC
message(
  "\nTracked cells by patient:"
)

print(
  t_meta %>%
    
    dplyr::filter(
      Tracked_Clonotype
    ) %>%
    
    dplyr::count(
      Patient_Match,
      sort = TRUE
    )
)

n_tracked <- sum(
  t_meta$Tracked_Clonotype
)

message(
  "\nTotal tracked cells:"
)

print(
  n_tracked
)

# Tracked-cell contribution by validated clonotype
tracked_bar_df <- tracked_df %>%
  
  dplyr::mutate(
    Patient = sub(
      "-post-clonotype.*",
      "",
      Probe
    ),
    
    Probe_Short = gsub(
      "-post-clonotype",
      "_c",
      Probe
    )
  ) %>%
  
  dplyr::count(
    Patient,
    Probe_Short,
    name = "N_Cells"
  )

patient_totals <- tracked_bar_df %>%
  
  dplyr::group_by(Patient) %>%
  
  dplyr::summarise(
    Total_Cells = sum(N_Cells),
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(Total_Cells)

tracked_bar_df$Patient <- factor(
  tracked_bar_df$Patient,
  levels = patient_totals$Patient
)

patient_totals$Patient <- factor(
  patient_totals$Patient,
  levels = patient_totals$Patient
)

write.csv(
  tracked_bar_df,
  "../../output/Fig4_Tracked_Cells_Per_Probe.csv",
  row.names = FALSE
)

p_patient_bar <-
  
  ggplot(
    tracked_bar_df,
    aes(
      x = N_Cells,
      y = Patient,
      fill = Probe_Short
    )
  ) +
  
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.25
  ) +
  
  geom_text(
    data = patient_totals,
    aes(
      x = Total_Cells,
      y = Patient,
      label = Total_Cells
    ),
    inherit.aes = FALSE,
    hjust = -0.2,
    size = 2
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.25)
    )
  ) +
  
  labs(
    x = paste0(
      "Tracked hyperexpanded clonotype cells (n = ",
      n_tracked,
      ")"
    ),
    y = NULL,
    fill = "High confidence clonotype"
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "none"
  )

save_pdf(
  "../../plots/Fig4_Hyperexpanded_clonotypes_tracked.pdf",
  p_patient_bar,
  P_Width,
  P_Height
)

# ============================================================
# Tracked clonotype lineage composition
# ============================================================

expr_cd <- FetchData(
  xe_global,
  vars = c("CD4", "CD8A", "CD8B"),
  layer = "data"
)

expr_cd$cell <- rownames(expr_cd)

tracked_cd <- t_meta %>%
  
  dplyr::filter(
    Tracked_Clonotype
  ) %>%
  
  dplyr::select(
    cell,
    sample_id
  ) %>%
  
  dplyr::left_join(
    expr_cd,
    by = "cell"
  ) %>%
  
  dplyr::mutate(
    
    Lineage = dplyr::case_when(
      
      CD4 > 0 &
        CD8A == 0 &
        CD8B == 0 ~ "CD4",
      
      (CD8A > 0 | CD8B > 0) &
        CD4 == 0 ~ "CD8",
      
      TRUE ~ "Undefined"
      
    )
    
  )

# Summary table
lineage_summary <- tracked_cd %>%
  
  dplyr::count(
    sample_id,
    Lineage
  ) %>%
  
  dplyr::group_by(
    sample_id
  ) %>%
  
  dplyr::mutate(
    Fraction = n / sum(n)
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    sample_id = gsub(
      "-post$",
      "",
      sample_id
    )
  )

write.csv(
  lineage_summary,
  "../../output/Fig4_Tracked_Clonotype_Lineage.csv",
  row.names = FALSE
)

print(lineage_summary)

p_tracked_lineage <-
  
  ggplot(
    lineage_summary,
    aes(
      x = sample_id,
      y = Fraction,
      fill = Lineage
    )
  ) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "CD4" = "#56B4E9",
      
      "CD8" = "#E31A1C",
      
      "Undefined" = "grey70"
      
    )
    
  ) +
  
  scale_y_continuous(
    labels =
      scales::percent_format()
  ) +
  
  labs(
    x = NULL,
    y = "Fraction of tracked clonotype cells",
    fill = NULL
  ) +
  
  theme_fig() +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5
      ),
    
    legend.position =
      "bottom"
    
  )

save_pdf(
  "../../plots/Fig4_Tracked_Clonotype_Lineage.pdf",
  p_tracked_lineage,
  P_Width,
  P_Height
)

# ============================================================
# Tracked clonotype T cell module analysis
# ============================================================

T_MODULES <- list(

  "Stemness"     = c("TCF7","IL7R","SELL"),
  "Cytotoxicity" = c("PRF1","NKG7","GNLY","GZMB"),
  "Exhaustion"   = c("PDCD1","LAG3","HAVCR2","TOX"),
  "Regulatory T" = c("FOXP3","IL2RA","CTLA4","TIGIT")
  
)

# T cell metadata
t_prog <- t_meta %>%
  
  dplyr::select(
    cell,
    sample_id,
    FOV,
    Tracked_Clonotype
  )

# Expression matrix
genes_use <- intersect(
  
  unique(
    unlist(T_MODULES)
  ),
  
  rownames(xe_global)
  
)

expr_df <- FetchData(
  
  xe_global,
  
  vars = genes_use,
  
  cells = t_prog$cell,
  
  layer = "data"
  
) %>%
  
  tibble::rownames_to_column(
    "cell"
  )

t_prog <- t_prog %>%
  
  dplyr::left_join(
    expr_df,
    by = "cell"
  )

# Module scores
for (mod in names(T_MODULES)) {
  
  g <- intersect(
    T_MODULES[[mod]],
    colnames(t_prog)
  )
  
  if (length(g) < 2)
    next
  
  t_prog[[mod]] <-
    
    as.numeric(
      
      scale(
        
        rowMeans(
          
          t_prog[
            ,
            g,
            drop = FALSE
          ],
          
          na.rm = TRUE
          
        )
        
      )
      
    )
}

# Mixed-effects testing
tcell_module_results <- list()

for (mod in names(T_MODULES)) {
  
  if (!mod %in% colnames(t_prog))
    next
  
  dat <- t_prog %>%
    
    dplyr::filter(
      !is.na(.data[[mod]])
    ) %>%
    
    dplyr::mutate(
      
      Module_Score =
        .data[[mod]],
      
      Tracked_Clonotype =
        factor(
          Tracked_Clonotype,
          levels = c(
            FALSE,
            TRUE
          )
        )
      
    )
  
  fit <- tryCatch(
    
    lmer(
      
      Module_Score ~
        
        Tracked_Clonotype +
        
        (1 | sample_id/FOV),
      
      data = dat
      
    ),
    
    error = function(e)
      NULL
    
  )
  
  if (is.null(fit))
    next
  
  coef_tab <- summary(
    fit
  )$coefficients
  
  tcell_module_results[[mod]] <-
    
    data.frame(
      
      Module = mod,
      
      Beta =
        coef_tab[
          "Tracked_ClonotypeTRUE",
          "Estimate"
        ],
      
      SE =
        coef_tab[
          "Tracked_ClonotypeTRUE",
          "Std. Error"
        ],
      
      P =
        coef_tab[
          "Tracked_ClonotypeTRUE",
          "Pr(>|t|)"
        ],
      
      stringsAsFactors = FALSE
      
    )
}

tcell_module_results <-
  
  bind_rows(
    tcell_module_results
  )

tcell_module_results$FDR <-
  
  p.adjust(
    tcell_module_results$P,
    method = "BH"
  )

tcell_module_results <-
  
  tcell_module_results %>%
  
  dplyr::mutate(
    
    CI_Low =
      Beta - 1.96 * SE,
    
    CI_High =
      Beta + 1.96 * SE,
    
    Sig = case_when(
  
      FDR < 0.001 ~ "***",
      
      FDR < 0.01 ~ "**",
      
      FDR < 0.05 ~ "*",
      
      TRUE ~ ""
      
    )
    
  )

write.csv(
  tcell_module_results,
  "../../output/Fig4_Tracked_Clonotype_T_NK_Modules.csv",
  row.names = FALSE
)

print(
  tcell_module_results
)

# Plot formatting
tcell_module_results$Module <- factor(
  
  tcell_module_results$Module,
  
  levels = rev(c(
    "Stemness",
    "Cytotoxicity",
    "Exhaustion",
    "Regulatory T"
  ))
  
)

# position significance labels

star_x <- max(
  tcell_module_results$CI_High,
  na.rm = TRUE
) + 0.10

# Forest plot
p_tracked_modules <-
  
  ggplot(
    tcell_module_results,
    aes(
      x = Beta,
      y = Module
    )
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey70"
  ) +
  
  geom_errorbarh(
    
    aes(
      xmin = CI_Low,
      xmax = CI_High
    ),
    
    height = 0.15,
    
    linewidth = 0.7
    
  ) +
  
  geom_point(
    
    shape = 23,
    
    size = 3,
    
    fill = "#009E73",
    
    colour = "black"
    
  ) +
  
  geom_text(
    
    aes(
      x = star_x,
      label = Sig
    ),
    
    size = pt_size / 2
    
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    
    x =
      "T/NK module score difference\n(tracked clonotypes vs other T/NK)",
    
    y = NULL
    
  ) +
  
  theme_fig() +
  
  theme(
    
    plot.margin =
      margin(
        5,
        12,
        5,
        5
      )
    
  )

save_pdf(
  "../../plots/Fig4_Tracked_Clonotype_T_NK_Modules.pdf",
  p_tracked_modules,
  P_Width,
  P_Height
)

# ============================================================
# Distance to nearest triad anchor: tracked clonotypes vs other
# ============================================================

DENSITY_RADIUS <- 50

t_test <- t_meta %>%
  
  dplyr::select(
    cell,
    sample_id,
    FOV,
    x,
    y,
    Tracked_Clonotype
  ) %>%
  
  dplyr::distinct(
    cell,
    .keep_all = TRUE
  )

# Distance to triad anchors
dist_list <- list()

for (sid in unique(t_test$sample_id)) {
  
  for (fov in unique(
    t_test$FOV[
      t_test$sample_id == sid
    ]
  )) {
    
    t_sub <- t_test %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    triad_sub <- triad_xy %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    if (
      nrow(t_sub) == 0 ||
      nrow(triad_sub) == 0
    ) next
    
    d <- nn_dists(
      t_sub[, c("x","y")],
      triad_sub[, c("x","y")]
    )
    
    dist_list[[length(dist_list)+1]] <-
      
      data.frame(
        cell = t_sub$cell,
        Dist_Triad = d
      )
  }
}

dist_df <- bind_rows(dist_list)

# Local densities
my_df <- meta %>%
  
  dplyr::filter(
    cell_type == "Myeloid"
  ) %>%
  
  dplyr::select(
    sample_id,
    FOV,
    x,
    y
  )

density_list <- list()

for (sid in unique(t_test$sample_id)) {
  
  for (fov in unique(
    t_test$FOV[
      t_test$sample_id == sid
    ]
  )) {
    
    t_sub <- t_test %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    my_sub <- my_df %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    if (nrow(t_sub) == 0)
      next
    
    local_t <- numeric(nrow(t_sub))
    local_my <- numeric(nrow(t_sub))
    
    for (i in seq_len(nrow(t_sub))) {
      
      local_t[i] <-
        
        sum(
          sqrt(
            (t_sub$x - t_sub$x[i])^2 +
              (t_sub$y - t_sub$y[i])^2
          ) <= DENSITY_RADIUS
        ) - 1
      
      local_my[i] <-
        
        sum(
          sqrt(
            (my_sub$x - t_sub$x[i])^2 +
              (my_sub$y - t_sub$y[i])^2
          ) <= DENSITY_RADIUS
        )
    }
    
    density_list[[length(density_list)+1]] <-
      
      data.frame(
        cell = t_sub$cell,
        Local_T_Density = local_t,
        Local_Myeloid_Density = local_my
      )
  }
}

density_df <- bind_rows(
  density_list
)

# Model dataframe
triad_model_df <- t_test %>%
  
  dplyr::left_join(
    dist_df,
    by = "cell"
  ) %>%
  
  dplyr::left_join(
    density_df,
    by = "cell"
  ) %>%
  
  dplyr::filter(
    !is.na(Dist_Triad)
  )

fit_triad_distance <- lmer(
  
  log1p(Dist_Triad) ~
    
    Tracked_Clonotype +
    
    scale(Local_T_Density) +
    
    scale(Local_Myeloid_Density) +
    
    (1 | sample_id/FOV),
  
  data = triad_model_df
  
)

print(
  summary(
    fit_triad_distance
  )
)

write.csv(
  triad_model_df,
  "../../output/Fig4_Triad_Proximity_Model_Data.csv",
  row.names = FALSE
)


# ============================================================
# APC-high proximity: Tracked clonotypes vs other T/NK cells
# ============================================================

APC_GENES <- c(
  "CD74",
  "CTSS",
  "CD80",
  "CD86"
)

DENSITY_RADIUS <- 50

# APC score
myeloid_df <- meta %>%
  
  dplyr::filter(
    cell_type == "Myeloid"
  ) %>%
  
  dplyr::distinct(
    cell,
    .keep_all = TRUE
  )

expr_mat <- FetchData(
  
  xe_global,
  
  vars = APC_GENES,
  
  cells = myeloid_df$cell,
  
  layer = "data"
  
)

apc_scores <- data.frame(
  cell = rownames(expr_mat),
  APC_Score = rowMeans(expr_mat),
  stringsAsFactors = FALSE
)

my_apc <- myeloid_df %>%
  
  dplyr::left_join(
    apc_scores,
    by = "cell"
  )

APC_CUTOFF <- quantile(
  my_apc$APC_Score,
  0.75,
  na.rm = TRUE
)

apc_high_df <- my_apc %>%
  
  dplyr::filter(
    APC_Score >= APC_CUTOFF
  ) %>%
  
  dplyr::select(
    sample_id,
    FOV,
    x,
    y
  )

# Distance to APC-high myeloid cells
t_cells <- t_meta %>%
  
  dplyr::distinct(
    cell,
    .keep_all = TRUE
  )

dist_list <- list()

for (sid in unique(t_cells$sample_id)) {
  
  for (fov in unique(
    t_cells$FOV[
      t_cells$sample_id == sid
    ]
  )) {
    
    t_sub <- t_cells %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    apc_sub <- apc_high_df %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    if (
      nrow(t_sub) == 0 ||
      nrow(apc_sub) == 0
    ) next
    
    d <- nn_dists(
      t_sub[, c("x","y")],
      apc_sub[, c("x","y")]
    )
    
    dist_list[[length(dist_list)+1]] <-
      
      data.frame(
        cell = t_sub$cell,
        Dist_APC_High = d
      )
  }
}

apc_dist_df <- bind_rows(
  dist_list
)

# Local densities
my_df <- meta %>%
  
  dplyr::filter(
    cell_type == "Myeloid"
  ) %>%
  
  dplyr::select(
    sample_id,
    FOV,
    x,
    y
  )

density_list <- list()

for (sid in unique(t_cells$sample_id)) {
  
  for (fov in unique(
    t_cells$FOV[
      t_cells$sample_id == sid
    ]
  )) {
    
    t_sub <- t_cells %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    my_sub <- my_df %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    if (nrow(t_sub) == 0)
      next
    
    local_t <- numeric(nrow(t_sub))
    local_my <- numeric(nrow(t_sub))
    
    for (i in seq_len(nrow(t_sub))) {
      
      local_t[i] <-
        
        sum(
          sqrt(
            (t_sub$x - t_sub$x[i])^2 +
              (t_sub$y - t_sub$y[i])^2
          ) <= DENSITY_RADIUS
        ) - 1
      
      local_my[i] <-
        
        sum(
          sqrt(
            (my_sub$x - t_sub$x[i])^2 +
              (my_sub$y - t_sub$y[i])^2
          ) <= DENSITY_RADIUS
        )
    }
    
    density_list[[length(density_list)+1]] <-
      
      data.frame(
        cell = t_sub$cell,
        Local_T_Density = local_t,
        Local_Myeloid_Density = local_my
      )
  }
}

density_df <- bind_rows(
  density_list
)

# Model
t_apc_adj <- t_cells %>%
  
  dplyr::left_join(
    apc_dist_df,
    by = "cell"
  ) %>%
  
  dplyr::left_join(
    density_df,
    by = "cell"
  ) %>%
  
  dplyr::filter(
    !is.na(Dist_APC_High)
  )

fit_apc_density <- lmer(
  
  log1p(Dist_APC_High) ~
    
    Tracked_Clonotype +
    
    scale(Local_T_Density) +
    
    scale(Local_Myeloid_Density) +
    
    (1 | sample_id/FOV),
  
  data = t_apc_adj
  
)

print(
  summary(
    fit_apc_density
  )
)

# Plot
coef_tab <- summary(
  fit_apc_density
)$coefficients

p_c <- coef_tab[
  "Tracked_ClonotypeTRUE",
  "Pr(>|t|)"
]

plot_df <- t_apc_adj %>%
  
  dplyr::mutate(
    
    Group = ifelse(
      Tracked_Clonotype,
      "Tracked clonotypes",
      "Other T/NK cells"
    )
    
  )


# Global medians
med_tracked <- median(
  plot_df$Dist_APC_High[
    plot_df$Group == "Tracked clonotypes"
  ],
  na.rm = TRUE
)

med_other <- median(
  plot_df$Dist_APC_High[
    plot_df$Group == "Other T/NK cells"
  ],
  na.rm = TRUE
)

write.csv(
  plot_df,
  "../../output/Fig4_APC_High_Distance_Distribution.csv",
  row.names = FALSE
)

p_apc_density <-
  
  ggplot(
    plot_df,
    aes(
      x = Dist_APC_High,
      fill = Group
    )
  ) +
  
  geom_density(
    alpha = 0.4,
    linewidth = 0.8,
    colour = "black"
  ) +
  
  scale_fill_manual(
    
    values = c(
      "Tracked clonotypes" = "#f5f5f5",
      "Other T/NK cells" = "#009E73"
    )
    
  ) +
  
  coord_cartesian(
    xlim = c(0,150)
  ) +
  
  ggplot2::annotate(
    
    "text",
    
    x = 100,
    y = Inf,
    
    vjust = 1.5,
    hjust = 1,
    
    parse = TRUE,
    
    size = pt_size / 2,
    
    label = format_p_plotmath(
      p_c
    )
    
  ) +
  
  geom_vline(
    xintercept = med_tracked,
    colour = "#4D4D4D",
    linetype = 2,
    linewidth = 0.35
  ) +
  
  geom_vline(
    xintercept = med_other,
    colour = "#009E73",
    linetype = 2,
    linewidth = 0.35
  ) +
  
  labs(
    
    x =
      "Distance to APC-high myeloid cell (µm)",
    
    y =
      "Density",
    
    fill = NULL
    
  ) +
  
  theme_fig() +
  
  theme(
    legend.position = "bottom"
  )

save_pdf(
  "../../plots/Fig4_APC_High_Distance_Density.pdf",
  p_apc_density,
  P_Width,
  P_Height
)

# ============================================================
# APC module score by distance to tracked clonotype
# ============================================================

# Distance from myeloid cells to tracked clonotypes
tracked_cells <- t_meta %>%
  
  dplyr::filter(
    Tracked_Clonotype
  ) %>%
  
  dplyr::select(
    cell,
    sample_id,
    FOV,
    x,
    y
  ) %>%
  
  dplyr::distinct(
    cell,
    .keep_all = TRUE
  )

dist_list <- list()

for (sid in unique(my_apc$sample_id)) {
  
  for (fov in unique(
    my_apc$FOV[
      my_apc$sample_id == sid
    ]
  )) {
    
    my_sub <- my_apc %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    t_sub <- tracked_cells %>%
      
      dplyr::filter(
        sample_id == sid,
        FOV == fov
      )
    
    if (
      nrow(my_sub) == 0 ||
      nrow(t_sub) == 0
    ) next
    
    d <- nn_dists(
      my_sub[, c("x", "y")],
      t_sub[, c("x", "y")]
    )
    
    dist_list[[length(dist_list) + 1]] <-
      
      data.frame(
        cell = my_sub$cell,
        Dist_Tracked = d
      )
  }
}

dist_df <- bind_rows(
  dist_list
)

print(head(dist_df))


my_apc_plot <- my_apc %>%
  
  dplyr::left_join(
    dist_df,
    by = "cell"
  ) %>%
  
  dplyr::filter(
    !is.na(Dist_Tracked)
  ) %>%
  
  dplyr::mutate(
    
    Distance_Bin = cut(
      
      Dist_Tracked,
      
      breaks = c(
        0,
        15,
        30,
        50,
        100,
        Inf
      ),
      
      labels = c(
        "0-15",
        "15-30",
        "30-50",
        "50-100",
        ">100"
      ),
      
      include.lowest = TRUE
      
    )
    
  )

# Distance bin model
my_apc_plot$Distance_Bin <- factor(
  
  my_apc_plot$Distance_Bin,
  
  levels = c(
    ">100",
    "50-100",
    "30-50",
    "15-30",
    "0-15"
  )
  
)

fit_apc_bin <- lmer(
  
  APC_Score ~
    
    Distance_Bin +
    
    (1 | sample_id/FOV),
  
  data = my_apc_plot
  
)

print(
  summary(
    fit_apc_bin
  )
)

# Mean APC score
mean_df <- my_apc_plot %>%
  
  dplyr::group_by(
    Distance_Bin
  ) %>%
  
  dplyr::summarise(
    
    Mean_APC =
      mean(
        APC_Score,
        na.rm = TRUE
      ),
    
    N = n(),
    
    .groups = "drop"
    
  )

# Extract LMM coefficients
coef_tab <- summary(
  fit_apc_bin
)$coefficients

coef_names <- rownames(
  coef_tab
)

coef_names <- coef_names[
  grepl(
    "^Distance_Bin",
    coef_names
  )
]

beta_df <- list()

for (cn in coef_names) {
  
  bin_name <- sub(
    "^Distance_Bin",
    "",
    cn
  )
  
  beta_df[[length(beta_df) + 1]] <-
    
    data.frame(
      
      Distance_Bin = bin_name,
      
      Beta =
        coef_tab[
          cn,
          "Estimate"
        ],
      
      P =
        coef_tab[
          cn,
          "Pr(>|t|)"
        ],
      
      stringsAsFactors = FALSE
      
    )
  
}

# Reference bin

beta_df[[length(beta_df) + 1]] <-
  
  data.frame(
    
    Distance_Bin = ">100",
    
    Beta = 0,
    
    P = NA,
    
    stringsAsFactors = FALSE
    
  )

beta_df <- bind_rows(
  beta_df
)

beta_df$FDR <- p.adjust(
  beta_df$P,
  method = "BH"
)

beta_df$Sig <- dplyr::case_when(
  
  beta_df$FDR < 0.001 ~ "***",
  
  beta_df$FDR < 0.01 ~ "**",
  
  beta_df$FDR < 0.05 ~ "*",
  
  TRUE ~ ""
  
)

# Plot dataframe
plot_df <- beta_df %>%
  
  dplyr::left_join(
    mean_df,
    by = "Distance_Bin"
  )

plot_df$Distance_Bin <- factor(
  
  plot_df$Distance_Bin,
  
  levels = c(
    "0-15",
    "15-30",
    "30-50",
    "50-100",
    ">100"
  )
  
)

write.csv(
  plot_df,
  "../../output/Fig4_APC_Module_DistanceBins.csv",
  row.names = FALSE
)

beta_lim <- max(
  abs(plot_df$Beta),
  na.rm = TRUE
)

# Dot plot
plot_df$Group <- "Myeloid"

p_apc_dot <-
  
  ggplot(
    plot_df,
    aes(
      x = Distance_Bin,
      y = "APC score"
    )
  ) +
  
  
  geom_point(
    aes(
      fill = Beta
    ),
    size = 8,
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +
  
  geom_text(
    
    aes(
      label = Sig
    ),
    
    nudge_y = 0.18,
    
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
  
  facet_grid(
    Group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  
  labs(
    x = "Distance to tracked clonotype (µm)",
    y = NULL
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
  "../../plots/Fig4_APC_Module_Distance_DotPlot.pdf",
  p_apc_dot,
  1.4*P_Width,
  P_Height
)

# ============================================================
# Triad participation vs distance to tracked clonotype
# ============================================================

# Reuse my_apc_plot from previous section
my_apc_plot$Triad <-
  my_apc_plot$cell %in% triad_anchors

# Mixed-effects logistic model
fit_triad <- glmer(
  
  Triad ~
    
    scale(
      log1p(Dist_Tracked)
    ) +
    
    (1 | sample_id/FOV),
  
  family = binomial,
  
  data = my_apc_plot
  
)

print(
  summary(
    fit_triad
  )
)

coef_tab <- summary(
  fit_triad
)$coefficients

beta_e <- coef_tab[
  "scale(log1p(Dist_Tracked))",
  "Estimate"
]

se_e <- coef_tab[
  "scale(log1p(Dist_Tracked))",
  "Std. Error"
]

p_e <- coef_tab[
  "scale(log1p(Dist_Tracked))",
  "Pr(>|z|)"
]

sig_e <- dplyr::case_when(
  p_e < 0.001 ~ "***",
  p_e < 0.01 ~ "**",
  p_e < 0.05 ~ "*",
  TRUE ~ "ns"
)

or_e <- exp(beta_e)

ci_low_e <- exp(
  beta_e - 1.96 * se_e
)

ci_high_e <- exp(
  beta_e + 1.96 * se_e
)

or_table <- data.frame(
  OR = or_e,
  CI_Low = ci_low_e,
  CI_High = ci_high_e,
  P = p_e
)

write.csv(
  or_table,
  "../../output/Fig4_Triad_Distance_OR.csv",
  row.names = FALSE
)

# Visualization summary
triad_summary <- my_apc_plot %>%
  
  dplyr::group_by(
    Distance_Bin
  ) %>%
  
  dplyr::summarise(
    
    Fraction_Triad =
      mean(Triad),
    
    N =
      n(),
    
    .groups = "drop"
    
  )

triad_summary$Distance_Bin <- factor(
  
  triad_summary$Distance_Bin,
  
  levels = c(
    "0-15",
    "15-30",
    "30-50",
    "50-100",
    ">100"
  )
  
)

write.csv(
  triad_summary,
  "../../output/Fig4_Triad_Prevalence_ByDistance.csv",
  row.names = FALSE
)

# Plot
p_triad_distance <-
  
  ggplot(
    triad_summary,
    aes(
      x = Distance_Bin,
      y = Fraction_Triad
    )
  ) +
  
  geom_col(
    
    fill = unname(
      pal_best["Myeloid"]
    ),
    
    colour = "black",
    
    width = 0.7,
    
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
  
  ggplot2::annotate(
    "text",
    x = 3.5,
    y = max(triad_summary$Fraction_Triad) * 1.30,
    hjust = 1,
    size = pt_size / 2.5,
    lineheight = 0.9,
    label = paste0(
      "OR = ",
      round(or_e, 2),
      "\n95% CI ",
      round(ci_low_e, 2),
      "–",
      round(ci_high_e, 2)
    )
  ) +
  
  ggplot2::annotate(
    "text",
    x = 3,
    y = max(triad_summary$Fraction_Triad) * 1.10,
    label = sig_e,
    size = pt_size / 1.4,
  ) +
  
  scale_y_continuous(
    
    labels =
      scales::percent_format(
        accuracy = 1
      ),
    
    expand =
      expansion(
        mult = c(
          0,
          0.35
        )
      )
    
  ) +
  
  labs(
    
    x =
      "Distance to tracked clonotype (µm)",
    
    y =
      "% of myeloid cells participating in triads"
    
  ) +
  
  theme_fig()

save_pdf(
  "../../plots/Fig4_Triad_Distance_Enrichment.pdf",
  p_triad_distance,
  1.4 * P_Width,
  P_Height
)
