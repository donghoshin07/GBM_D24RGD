# ==============================================================================
# Figures 1 and S1 analyses using Delta-24-RGD data
# ==============================================================================

# Expects this script to be located in: /scripts/3_downstream_analyses
#
# Required input files: 
# ../../input/TCR/patient/*.txt
#
# Output directories:
# ../../plots
# ../../output

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)), silent = TRUE)
}

suppressPackageStartupMessages({
  
  library(immunarch)
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  
  library(readr)
  
  library(ggplot2)
  library(scales)
  library(treemapify)
  
  library(reshape2)
  library(survival)
  library(vegan)
  
  library(patchwork)
  library(cowplot)
  
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
# Colors
# ============================================================

cols_time <- c(
  "Pre" = "#E41A1C",
  "Post" = "#377EB8",
  "Shared" = "purple",
  "Unique Pre" = "#E41A1C",
  "Unique Post" = "#377EB8"
)

clone_levels <- c("Small", "Medium", "Large", "Hyperexpanded")

# ============================================================
# Data loading and cacheing
# ============================================================

TRB <- repLoad("../input/TCR/patient")

TRB$meta <- TRB$meta %>%
  mutate(
    OS = as.numeric(OS),
    Status = 1L,
    Group = factor(
      ifelse(grepl("pre", Sample, ignore.case = TRUE), "Pre", "Post"),
      levels = c("Pre", "Post")
    )
  )

pre_samples <- c(
  "P15_pre", "P21_pre", "P22_pre",
  "P27_pre", "P31_pre", "P42_pre"
)

post_samples <- c(
  "P13_post", "P15_post", "P17_post",
  "P21_post", "P22_post", "P27_post",
  "P39_post", "P41_post", "P42_post"
)

target_samples <- c(pre_samples, post_samples)

pairs <- c("P15", "P21", "P22", "P27", "P42")

# Inverse Simpson diversity estimated across 100 downsampling iterations

cache_file <- "../../output/TRB_Diversity_N100.rds"

if (file.exists(cache_file)) {
  
  div_robust <- readRDS(cache_file)
  
} else {
  
  div_list <- vector("list", 100)
  
  for (i in seq_len(100)) {
    
    ds_i <- repSample(
      TRB$data[target_samples],
      "downsample"
    )
    
    div_i <- repDiversity(ds_i, "inv.simp")
    
    div_list[[i]] <-
      as.data.frame(div_i) %>%
      setNames(c("Sample", "Div"))
    
  }
  
  div_robust <- bind_rows(div_list) %>%
    group_by(Sample) %>%
    summarise(
      TCR_Div = mean(Div, na.rm = TRUE),
      .groups = "drop"
    )
  
  saveRDS(div_robust, cache_file)
  
}

meta_clean <- TRB$meta %>%
  left_join(div_robust, by = "Sample") %>%
  filter(Sample %in% target_samples)

ds_viz <- repSample(
  TRB$data[target_samples],
  "downsample"
)

# ============================================================
# TCR diversity (pre vs post)
# ============================================================

w_p <- wilcox.test(
  TCR_Div ~ Group,
  data = meta_clean,
  exact = FALSE
)$p.value

lab_df <- data.frame(
  x   = 1.5,
  y   = max(meta_clean$TCR_Div, na.rm = TRUE) * 1.10,
  lab = format_p_plotmath(w_p)
)

p_TCR_diversity <- ggplot(
  meta_clean,
  aes(Group, TCR_Div, color = Group)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.25
  ) +
  geom_jitter(
    width = 0.2,
    size = 1.6,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = cols_time[c("Pre", "Post")]
  ) +
  geom_text(
    data = lab_df,
    aes(x, y, label = lab),
    inherit.aes = FALSE,
    parse = TRUE,
    size = pt_size / 2.1
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_fig() +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "TCR Diversity (Inverse Simpson)"
  )

save_pdf(
  "../../plots/Fig1_TCR_Diversity_pre_vs_post.pdf",
  p_TCR_diversity,
  0.8*P_Width,
  P_Height
)

# ============================================================
# Clone-size distributions
# ============================================================

homeo <- repClonality(
  ds_viz,
  .method = "homeo"
)

clon_clean <- as.data.frame(homeo) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    -Sample,
    names_to = "Clone_Size",
    values_to = "Proportion"
  ) %>%
  mutate(
    Clone_Size = recode(
      Clone_Size,
      "Rare (0 < X <= 1e-05)" = "Rare",
      "Small (1e-05 < X <= 1e-04)" = "Small",
      "Medium (1e-04 < X <= 0.001)" = "Medium",
      "Large (0.001 < X <= 0.01)" = "Large",
      "Hyperexpanded (0.01 < X <= 1)" = "Hyperexpanded"
    ),
    Group = factor(
      ifelse(grepl("pre", Sample, ignore.case = TRUE), "Pre", "Post"),
      levels = c("Pre", "Post")
    ),
    Patient_ID = str_extract(Sample, "P[0-9]+")
  ) %>%
  filter(Clone_Size != "Rare")

# PERMANOVA

perm_df <- as.data.frame(homeo)

perm_meta <- data.frame(
  Group = factor(
    ifelse(
      grepl("pre", rownames(perm_df), ignore.case = TRUE),
      "Pre",
      "Post"
    ),
    levels = c("Pre", "Post")
  )
)

p_clone <- adonis2(
  perm_df ~ Group,
  data = perm_meta
)$`Pr(>F)`[1]

reds  <- c("#FEE0D2", "#FC9272", "#DE2D26", "#A50F15")
blues <- c("#DEEBF7", "#9ECAE1", "#3182BD", "#08519C")

color_mapping <- c(
  "Pre Small"         = reds[1],
  "Pre Medium"        = reds[2],
  "Pre Large"         = reds[3],
  "Pre Hyperexpanded" = reds[4],
  "Post Small"         = blues[1],
  "Post Medium"        = blues[2],
  "Post Large"         = blues[3],
  "Post Hyperexpanded" = blues[4]
)

clon_clean <- clon_clean %>%
  mutate(
    Fill_Group = paste(Group, Clone_Size)
  )

p_clonal_expansion <- ggplot(
  clon_clean,
  aes(
    Patient_ID,
    Proportion,
    fill = forcats::fct_rev(Fill_Group)
  )
) +
  geom_col(
    position = "fill",
    width = 0.8
  ) +
  facet_grid(
    ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(
    values = color_mapping,
    name = "Clone Size"
  ) +
  geom_text(
    data = data.frame(
      Group = factor("Post", levels = c("Pre", "Post")),
      x = Inf,
      y = Inf,
      label = format_p_plotmath(p_clone)
    ),
    aes(x, y, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 1.1,
    vjust = 1.5,
    size = pt_size / 2.2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_fig() +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.18, "cm"),
    legend.key.width  = unit(0.18, "cm"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = pt_size * 0.8
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction"
  )

save_pdf(
  "../../plots/Fig1_Clonal_Expansion.pdf",
  p_clonal_expansion,
  1.2 * P_Width,
  P_Height
)

# ============================================================
# Repertoire Occupancy
# ============================================================

process_patient_data <- function(pre_data, post_data, patient_id) {
  
  joined_data <- full_join(
    pre_data %>% select(CDR3.aa, Value_pre = Proportion),
    post_data %>% select(CDR3.aa, Value_post = Proportion),
    by = "CDR3.aa",
    relationship = "many-to-many"
  ) %>%
    mutate(
      Value_pre = coalesce(Value_pre, 0),
      Value_post = coalesce(Value_post, 0),
      Category = case_when(
        Value_pre > 0 & Value_post == 0 ~ "Unique Pre",
        Value_pre == 0 & Value_post > 0 ~ "Unique Post",
        Value_pre > 0 & Value_post > 0 ~ "Shared"
      )
    )
  
  bind_rows(
    joined_data %>%
      group_by(Category) %>%
      summarise(
        Cumulative_Frequency = sum(Value_pre),
        .groups = "drop"
      ) %>%
      mutate(Sample = paste0(patient_id, "_pre")),
    
    joined_data %>%
      group_by(Category) %>%
      summarise(
        Cumulative_Frequency = sum(Value_post),
        .groups = "drop"
      ) %>%
      mutate(Sample = paste0(patient_id, "_post"))
  ) %>%
    complete(
      Category = c("Unique Pre", "Unique Post", "Shared"),
      fill = list(Cumulative_Frequency = 0)
    ) %>%
    group_by(Sample) %>%
    mutate(
      Cumulative_Frequency =
        Cumulative_Frequency / sum(Cumulative_Frequency)
    ) %>%
    ungroup() %>%
    mutate(Patient = patient_id)
  
}

patient_data_list <- lapply(
  pairs,
  function(id) {
    process_patient_data(
      ds_viz[[paste0(id, "_pre")]],
      ds_viz[[paste0(id, "_post")]],
      id
    )
  }
)

TRB_data_pre_post <- bind_rows(patient_data_list) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c(
        paste0(pairs, "_pre"),
        paste0(pairs, "_post")
      )
    ),
    Category = factor(
      Category,
      levels = c("Unique Pre", "Unique Post", "Shared")
    ),
    Timepoint = factor(
      ifelse(grepl("pre", Sample, ignore.case = TRUE), "Pre", "Post"),
      levels = c("Pre", "Post")
    )
  )

p_repertoire_space <- ggplot(
  TRB_data_pre_post,
  aes(
    Timepoint,
    Cumulative_Frequency,
    fill = Category
  )
) +
  geom_col(
    width = 0.85,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Unique Pre" = "#E41A1C",
      "Unique Post" = "#377EB8",
      "Shared" = "purple"
    )
  ) +
  facet_wrap(
    ~ Patient,
    scales = "free_x",
    nrow = 1
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    legend.position = "bottom"
  ) +
  labs(
    x = NULL,
    y = "Cumulative Frequency",
    fill = NULL
  )

save_pdf(
  "../../plots/Fig1_Repertoire_Occupancy.pdf",
  p_repertoire_space,
  1.6 * P_Width,
  P_Height
)

# ============================================================
# Morisita repertoire overlap
# ============================================================

ov_mat <- repOverlap(
  ds_viz,
  .method = "morisita"
)

ov_long <- as.data.frame(ov_mat) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(
    -Sample1,
    names_to = "Sample2",
    values_to = "Morisita"
  ) %>%
  mutate(
    Sample1 = factor(Sample1, levels = target_samples),
    Sample2 = factor(Sample2, levels = target_samples)
  )

p_overlap <- ggplot(
  ov_long,
  aes(Sample1, Sample2, fill = Morisita)
) +
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_gradient(
    low = "#F2F0F7",
    high = "#54278F",
    name = "Morisita"
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = NULL
  )

save_pdf(
  "../../plots/Fig1_Overlap_Morisita_Heatmap.pdf",
  p_overlap,
  P_Width,
  P_Height
)

# ============================================================
# TCR diversity vs overall survival
# ============================================================

surv_df <- meta_clean %>%
  filter(Group == "Post") %>%
  drop_na(TCR_Div, OS)

sp_res <- cor.test(
  surv_df$TCR_Div,
  surv_df$OS,
  method = "spearman"
)

lab_df_sc <- data.frame(
  x = min(surv_df$TCR_Div, na.rm = TRUE),
  y = max(surv_df$OS, na.rm = TRUE) * 1.15,
  lab = sprintf(
    "atop(rho == %.2f, italic(P) == %.3f)",
    unname(sp_res$estimate),
    sp_res$p.value
  )
)

p_diversity_os <- ggplot(
  surv_df,
  aes(TCR_Div, OS)
) +
  geom_smooth(
    method = "lm",
    color = "black",
    linetype = "dashed",
    linewidth = 0.25,
    se = FALSE
  ) +
  geom_point(
    size = 1.6,
    color = cols_time["Post"]
  ) +
  geom_text(
    data = lab_df_sc,
    aes(x, y, label = lab),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = pt_size / 2.2
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  theme_fig() +
  labs(
    x = "Post-virotherapy TCR diversity (Inverse Simpson)",
    y = "OS (months)"
  )

save_pdf(
  "../../plots/Fig1_Diversity_OS_Scatter.pdf",
  p_diversity_os,
  1.3 * P_Width,
  P_Height
)

# ============================================================
# Kaplan-Meier survival by TCR diversity
# ============================================================

surv_df <- surv_df %>%
  mutate(
    Div_Group = factor(
      ifelse(
        TCR_Div >= median(TCR_Div, na.rm = TRUE),
        "High",
        "Low"
      ),
      levels = c("Low", "High")
    )
  )

km_cols <- c(
  "Low" = "#9ECAE1",
  "High" = "#08306B"
)

sf <- survfit(
  Surv(OS, Status) ~ Div_Group,
  data = surv_df
)

make_km_tables <- function(sf_obj, group_levels) {
  
  ss <- summary(sf_obj)
  
  km_df <- data.frame(
    time  = ss$time,
    surv  = ss$surv,
    lower = pmax(pmin(ss$lower, 1), 0),
    upper = pmax(pmin(ss$upper, 1), 0),
    group = factor(
      sub("^Div_Group=", "", ss$strata),
      levels = group_levels
    )
  )
  
  bind_rows(
    lapply(split(km_df, km_df$group), function(df) {
      
      bind_rows(
        data.frame(
          time = 0,
          surv = 1,
          lower = 1,
          upper = 1,
          group = df$group[1]
        ),
        df
      )
      
    })
  )
  
}

km_pts <- make_km_tables(
  sf,
  levels(surv_df$Div_Group)
)

build_ci_polygons <- function(km_df) {
  
  bind_rows(
    lapply(split(km_df, km_df$group), function(d) {
      
      d <- d[order(d$time), ]
      
      t  <- d$time
      up <- d$upper
      lo <- d$lower
      
      x_top <- as.numeric(
        rbind(t[-length(t)], t[-1])
      )
      
      data.frame(
        x = c(x_top, rev(x_top)),
        y = c(
          rep(up[-length(up)], each = 2),
          rev(rep(lo[-length(lo)], each = 2))
        ),
        group = d$group[1]
      )
      
    })
  )
  
}

ci_poly <- build_ci_polygons(km_pts)

p_lr <- 1 - pchisq(
  survdiff(
    Surv(OS, Status) ~ Div_Group,
    data = surv_df
  )$chisq,
  1
)

max_t <- max(km_pts$time)

p_km <- ggplot() +
  geom_polygon(
    data = ci_poly,
    aes(x, y, fill = group, group = group),
    alpha = 0.18,
    color = NA
  ) +
  geom_step(
    data = km_pts,
    aes(time, surv, color = group, group = group),
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = km_cols,
    name = "Diversity group"
  ) +
  scale_fill_manual(
    values = km_cols,
    name = "Diversity group"
  ) +
  coord_cartesian(
    xlim = c(0, max_t * 1.05),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = max_t * 0.02,
    y = 0.05,
    label = format_p_plotmath(p_lr),
    parse = TRUE,
    hjust = 0,
    vjust = 0,
    size = pt_size / 2
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
  "../../plots/Fig1_Diversity_OS_KM.pdf",
  p_km,
  1.3 * P_Width,
  P_Height
)

# ============================================================
# Clonotype treemap
# ============================================================

build_patient_treemap_post <- function(immunarch_data, patient_id) {
  
  df_post <- immunarch_data[[paste0(patient_id, "_post")]] %>%
    arrange(desc(Proportion))
  
  ggplot(
    df_post,
    aes(area = Proportion, fill = CDR3.aa)
  ) +
    geom_treemap(
      color = "white",
      size = 0.2
    ) +
    scale_fill_viridis_d(
      option = "turbo",
      guide = "none"
    ) +
    labs(
      title = patient_id
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 10
      ),
      plot.margin = margin(5, 5, 5, 5),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
  
}

p_treemap_P15 <- build_patient_treemap_post(
  ds_viz,
  "P15"
)

p_treemap_P21 <- build_patient_treemap_post(
  ds_viz,
  "P21"
)

p_treemap <- p_treemap_P15 + p_treemap_P21 +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  )

save_pdf(
  "../../plots/Fig1_Clone_Treemap_Post.pdf",
  p_treemap,
  1.6 * P_Width,
  P_Height
)

# ==============================================================================
# Save diversity mapping for downstream figures
# ==============================================================================

diversity_mapping <- TRB$meta %>%
  left_join(div_robust, by = "Sample") %>%
  filter(Sample %in% post_samples) %>%
  transmute(Sample,
            Div_Group = ifelse(TCR_Div >= median(TCR_Div, na.rm=TRUE), "High","Low"),
            TCR_Div, OS = as.numeric(OS), Status = 1L,
            Patient_ID = str_extract(Sample, "P[0-9]+"))
saveRDS(diversity_mapping, "../../output/TCR_Diversity_Groups.rds")
message(">>> Diversity mapping saved for downstream figures.")

# ============================================================
# Repertoire size before and after downsampling
# ============================================================

raw_volume <- repExplore(
  TRB$data[target_samples],
  .method = "volume"
)

raw_clones <- repExplore(
  TRB$data[target_samples],
  .method = "clones"
)

ds_volume <- repExplore(
  ds_viz,
  .method = "volume"
)

ds_clones <- repExplore(
  ds_viz,
  .method = "clones"
)

raw_volume$Sample <- factor(
  raw_volume$Sample,
  levels = target_samples
)

raw_clones$Sample <- factor(
  raw_clones$Sample,
  levels = target_samples
)

ds_volume$Sample <- factor(
  ds_volume$Sample,
  levels = target_samples
)

ds_clones$Sample <- factor(
  ds_clones$Sample,
  levels = target_samples
)

p_raw_clonotypes <- ggplot(
  raw_volume,
  aes(Sample, Volume)
) +
  geom_col(
    fill = "grey40",
    width = 0.8
  ) +
  theme_fig() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title = "Raw Clonotypes",
    x = NULL,
    y = "Clonotypes"
  )

p_raw_clones <- ggplot(
  raw_clones,
  aes(Sample, Clones)
) +
  geom_col(
    fill = "grey40",
    width = 0.8
  ) +
  theme_fig() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title = "Raw Clones",
    x = NULL,
    y = "Clones"
  )

p_ds_clonotypes <- ggplot(
  ds_volume,
  aes(Sample, Volume)
) +
  geom_col(
    fill = "grey40",
    width = 0.8
  ) +
  theme_fig() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title = "Downsampled Clonotypes",
    x = NULL,
    y = "Clonotypes"
  )

p_ds_clones <- ggplot(
  ds_clones,
  aes(Sample, Clones)
) +
  geom_col(
    fill = "grey40",
    width = 0.8
  ) +
  theme_fig() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = pt_size
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(
    title = "Downsampled Clones",
    x = NULL,
    y = "Clones"
  )

p_downsampling <- (
  p_raw_clones | p_ds_clones
) /
  (
    p_raw_clonotypes | p_ds_clonotypes
  )

save_pdf(
  "../../plots/FigS1_Downsampling_Effects.pdf",
  p_downsampling,
  2 * P_Width,
  1.5 * P_Height
)