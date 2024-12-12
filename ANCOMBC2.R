## ----setup, message = FALSE, warning = FALSE, comment = NA--------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
  "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

# It appears to be a package compatibility issue between the release version of 
# phyloseq and lme4, a fresh installation of phyloseq might be needed
# See this post: https://github.com/lme4/lme4/issues/743
# remotes::install_github("joey711/phyloseq", force = TRUE)

## ----getPackage, eval=FALSE---------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("ANCOMBC")

## ----load, eval=FALSE---------------------------------------------------------
#  library(ANCOMBC)

## -----------------------------------------------------------------------------
data(QMP, package = "ANCOMBC")
set.seed(123)
n = 150
d = ncol(QMP)
diff_prop = 0.1
lfc_cont = 1
lfc_cat2_vs_1 = -2
lfc_cat3_vs_1 = 1

# Generate the true abundances
abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05, 
                    n = n, lib_mean = 1e8, disp = 0.5)
log_abn_data = log(abn_data + 1e-5)
rownames(log_abn_data) = paste0("T", seq_len(d))
colnames(log_abn_data) = paste0("S", seq_len(n))

# Generate the sample and feature meta data
# Sampling fractions are set to differ by batches
smd = data.frame(samp_frac = log(c(runif(n/3, min = 1e-4, max = 1e-3),
                                   runif(n/3, min = 1e-3, max = 1e-2),
                                   runif(n/3, min = 1e-2, max = 1e-1))),
                 cont_cov = rnorm(n),
                 cat_cov = as.factor(rep(seq_len(3), each = n/3)))
rownames(smd) = paste0("S", seq_len(n))
                      
fmd = data.frame(taxon = paste0("T", seq_len(d)),
                 seq_eff = log(runif(d, min = 0.1, max = 1)),
                 lfc_cont = sample(c(0, lfc_cont), 
                                   size = d,
                                   replace = TRUE,
                                   prob = c(1 - diff_prop, diff_prop)),
                 lfc_cat2_vs_1 = sample(c(0, lfc_cat2_vs_1), 
                                        size = d,
                                        replace = TRUE,
                                        prob = c(1 - diff_prop, diff_prop)),
                 lfc_cat3_vs_1 = sample(c(0, lfc_cat3_vs_1), 
                                        size = d,
                                        replace = TRUE,
                                        prob = c(1 - diff_prop, diff_prop))) %>%
    mutate(lfc_cat3_vs_2 = lfc_cat3_vs_1 - lfc_cat2_vs_1)

# Add effect sizes of covariates to the true abundances
smd_dmy = model.matrix(~ 0 + cont_cov + cat_cov, data = smd)
log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"] )
log_abn_data = log_abn_data + outer(fmd$lfc_cat2_vs_1, smd_dmy[, "cat_cov2"])
log_abn_data = log_abn_data + outer(fmd$lfc_cat3_vs_1, smd_dmy[, "cat_cov3"])

# Add sample- and taxon-specific biases
log_otu_data = t(t(log_abn_data) + smd$samp_frac)
log_otu_data = log_otu_data + fmd$seq_eff
otu_data = round(exp(log_otu_data))

# Create the tse object
assays = S4Vectors::SimpleList(counts = otu_data)
smd = S4Vectors::DataFrame(smd)
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)

## -----------------------------------------------------------------------------
set.seed(123)
output = ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                  fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)

res_prim = output$res
res_pair = output$res_pair

## -----------------------------------------------------------------------------
res_merge1 = res_pair %>%
  dplyr::transmute(taxon, 
                   lfc_est1 = lfc_cat_cov2 * diff_cat_cov2,
                   lfc_est2 = lfc_cat_cov3 * diff_cat_cov3,
                   lfc_est3 = lfc_cat_cov3_cat_cov2 * diff_cat_cov3_cat_cov2) %>%
  dplyr::left_join(fmd %>%
                     dplyr::transmute(taxon, 
                                      lfc_true1 = lfc_cat2_vs_1,
                                      lfc_true2 = lfc_cat3_vs_1,
                                      lfc_true3 = lfc_cat3_vs_2),
                   by = "taxon") %>%
  dplyr::transmute(taxon, 
                   lfc_est1 = case_when(lfc_est1 > 0 ~ 1,
                                        lfc_est1 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_est2 = case_when(lfc_est2 > 0 ~ 1,
                                        lfc_est2 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_est3 = case_when(lfc_est3 > 0 ~ 1,
                                        lfc_est3 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_true1 = case_when(lfc_true1 > 0 ~ 1,
                                         lfc_true1 < 0 ~ -1,
                                         TRUE ~ 0),
                   lfc_true2 = case_when(lfc_true2 > 0 ~ 1,
                                         lfc_true2 < 0 ~ -1,
                                         TRUE ~ 0),
                   lfc_true3 = case_when(lfc_true3 > 0 ~ 1,
                                         lfc_true3 < 0 ~ -1,
                                         TRUE ~ 0))
lfc_est1 = res_merge1$lfc_est1
lfc_true1 = res_merge1$lfc_true1
lfc_est2 = res_merge1$lfc_est2
lfc_true2 = res_merge1$lfc_true2
lfc_est3 = res_merge1$lfc_est3
lfc_true3 = res_merge1$lfc_true3

tp1 = sum(lfc_true1 == 1 & lfc_est1 == 1) +
  sum(lfc_true1 == -1 & lfc_est1 == -1)
fp1 = sum(lfc_true1 == 0 & lfc_est1 != 0) +
  sum(lfc_true1 == 1 & lfc_est1 == -1) +
  sum(lfc_true1 == -1 & lfc_est1 == 1)
fn1 = sum(lfc_true1 != 0 & lfc_est1 == 0)

tp2 = sum(lfc_true2 == 1 & lfc_est2 == 1) +
  sum(lfc_true2 == -1 & lfc_est2 == -1)
fp2 = sum(lfc_true2 == 0 & lfc_est2 != 0) +
  sum(lfc_true2 == 1 & lfc_est2 == -1) +
  sum(lfc_true2 == -1 & lfc_est2 == 1)
fn2 = sum(lfc_true2 != 0 & lfc_est2 == 0)

tp3 = sum(lfc_true3 == 1 & lfc_est3 == 1) +
  sum(lfc_true3 == -1 & lfc_est3 == -1)
fp3 = sum(lfc_true3 == 0 & lfc_est3 != 0) +
  sum(lfc_true3 == 1 & lfc_est3 == -1) +
  sum(lfc_true3 == -1 & lfc_est3 == 1)
fn3 = sum(lfc_true3 != 0 & lfc_est3 == 0)

tp = tp1 + tp2 + tp3
fp = fp1 + fp2 + fp3
fn = fn1 + fn2 + fn3

power1 = tp/(tp + fn)
fdr1 = fp/(tp + fp)

res_merge2 = res_pair %>%
  dplyr::transmute(taxon, 
                   lfc_est1 = lfc_cat_cov2 * diff_cat_cov2 * passed_ss_cat_cov2,
                   lfc_est2 = lfc_cat_cov3 * diff_cat_cov3 * passed_ss_cat_cov3,
                   lfc_est3 = lfc_cat_cov3_cat_cov2 * diff_cat_cov3_cat_cov2* passed_ss_cat_cov3_cat_cov2) %>%
  dplyr::left_join(fmd %>%
                     dplyr::transmute(taxon, 
                                      lfc_true1 = lfc_cat2_vs_1,
                                      lfc_true2 = lfc_cat3_vs_1,
                                      lfc_true3 = lfc_cat3_vs_2),
                   by = "taxon") %>%
  dplyr::transmute(taxon, 
                   lfc_est1 = case_when(lfc_est1 > 0 ~ 1,
                                        lfc_est1 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_est2 = case_when(lfc_est2 > 0 ~ 1,
                                        lfc_est2 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_est3 = case_when(lfc_est3 > 0 ~ 1,
                                        lfc_est3 < 0 ~ -1,
                                        TRUE ~ 0),
                   lfc_true1 = case_when(lfc_true1 > 0 ~ 1,
                                         lfc_true1 < 0 ~ -1,
                                         TRUE ~ 0),
                   lfc_true2 = case_when(lfc_true2 > 0 ~ 1,
                                         lfc_true2 < 0 ~ -1,
                                         TRUE ~ 0),
                   lfc_true3 = case_when(lfc_true3 > 0 ~ 1,
                                         lfc_true3 < 0 ~ -1,
                                         TRUE ~ 0))
lfc_est1 = res_merge2$lfc_est1
lfc_true1 = res_merge2$lfc_true1
lfc_est2 = res_merge2$lfc_est2
lfc_true2 = res_merge2$lfc_true2
lfc_est3 = res_merge2$lfc_est3
lfc_true3 = res_merge2$lfc_true3

tp1 = sum(lfc_true1 == 1 & lfc_est1 == 1) +
  sum(lfc_true1 == -1 & lfc_est1 == -1)
fp1 = sum(lfc_true1 == 0 & lfc_est1 != 0) +
  sum(lfc_true1 == 1 & lfc_est1 == -1) +
  sum(lfc_true1 == -1 & lfc_est1 == 1)
fn1 = sum(lfc_true1 != 0 & lfc_est1 == 0)

tp2 = sum(lfc_true2 == 1 & lfc_est2 == 1) +
  sum(lfc_true2 == -1 & lfc_est2 == -1)
fp2 = sum(lfc_true2 == 0 & lfc_est2 != 0) +
  sum(lfc_true2 == 1 & lfc_est2 == -1) +
  sum(lfc_true2 == -1 & lfc_est2 == 1)
fn2 = sum(lfc_true2 != 0 & lfc_est2 == 0)

tp3 = sum(lfc_true3 == 1 & lfc_est3 == 1) +
  sum(lfc_true3 == -1 & lfc_est3 == -1)
fp3 = sum(lfc_true3 == 0 & lfc_est3 != 0) +
  sum(lfc_true3 == 1 & lfc_est3 == -1) +
  sum(lfc_true3 == -1 & lfc_est3 == 1)
fn3 = sum(lfc_true3 != 0 & lfc_est3 == 0)

tp = tp1 + tp2 + tp3
fp = fp1 + fp2 + fp3
fn = fn1 + fn2 + fn3

power2 = tp/(tp + fn)
fdr2 = fp/(tp + fp)

tab_summ = data.frame(Comparison = c("Without sensitivity score filter", 
                                     "With sensitivity score filter"),
                      Power = round(c(power1, power2), 2),
                      FDR = round(c(fdr1, fdr2), 2))
tab_summ %>%
    datatable(caption = "Power/FDR Comparison")

## -----------------------------------------------------------------------------
data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(sample_data(tse)$bmi)

# Create the region variable
tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")

# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]

print(tse)

## -----------------------------------------------------------------------------
set.seed(123)
# It should be noted that we have set the number of bootstrap samples (B) equal 
# to 10 in the 'trend_control' function for computational expediency. 
# However, it is recommended that users utilize the default value of B, 
# which is 100, or larger values for optimal performance.
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = "age + region + bmi", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "bmi", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2, 1),
                                       solver = "ECOS",
                                       B = 10))

## -----------------------------------------------------------------------------
tab_zero = output$zero_ind
tab_zero %>%
    datatable(caption = "The detection of structural zeros")

## -----------------------------------------------------------------------------
res_prim = output$res

## -----------------------------------------------------------------------------
df_age = res_prim %>%
    dplyr::select(taxon, ends_with("age")) 
df_fig_age = df_age %>%
    dplyr::filter(diff_age == 1) %>% 
    dplyr::arrange(desc(lfc_age)) %>%
    dplyr::mutate(direct = ifelse(lfc_age > 0, "Positive LFC", "Negative LFC"),
                  color = ifelse(passed_ss_age == 1, "aquamarine3", "black"))
df_fig_age$taxon = factor(df_fig_age$taxon, levels = df_fig_age$taxon)
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC"))
  
fig_age = df_fig_age %>%
    ggplot(aes(x = taxon, y = lfc_age, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_age - se_age, ymax = lfc_age + se_age), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = "Log fold changes as one unit increase of age") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,
                                     color = df_fig_age$color))
fig_age

## -----------------------------------------------------------------------------
df_bmi = res_prim %>%
    dplyr::select(taxon, contains("bmi")) 

df_fig_bmi1 = df_bmi %>%
    dplyr::filter(diff_bmilean == 1 | 
                    diff_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(diff_bmioverweight == 1, 
                                round(lfc_bmioverweight, 2), 0),
                  lfc2 = ifelse(diff_bmilean == 1, 
                                round(lfc_bmilean, 2), 0)) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_bmi2 = df_bmi %>%
    dplyr::filter(diff_bmilean == 1 | 
                    diff_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(passed_ss_bmioverweight == 1 & diff_bmioverweight == 1, 
                                "aquamarine3", "black"),
                  lfc2 = ifelse(passed_ss_bmilean == 1 & diff_bmilean == 1, 
                                "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "color") %>%
    dplyr::arrange(taxon)

df_fig_bmi = df_fig_bmi1 %>%
    dplyr::left_join(df_fig_bmi2, by = c("taxon", "group"))

df_fig_bmi$group = recode(df_fig_bmi$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese")
df_fig_bmi$group = factor(df_fig_bmi$group, 
                          levels = c("Overweight - Obese",
                                     "Lean - Obese"))
  
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
fig_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_bmi

## -----------------------------------------------------------------------------
res_global = output$res_global
df_bmi = res_prim %>%
    dplyr::select(taxon, contains("bmi")) 
df_fig_global = df_bmi %>%
    dplyr::left_join(res_global %>%
                       dplyr::transmute(taxon, 
                                        diff_bmi = diff_abn, 
                                        passed_ss = passed_ss)) %>%
    dplyr::filter(diff_bmi == 1) %>%
    dplyr::mutate(lfc_overweight = lfc_bmioverweight,
                  lfc_lean = lfc_bmilean,
                  color = ifelse(passed_ss == 1, "aquamarine3", "black")) %>%
    dplyr::transmute(taxon,
                     `Overweight - Obese` = round(lfc_overweight, 2),
                     `Lean - Obese` = round(lfc_lean, 2), 
                     color = color) %>%
    tidyr::pivot_longer(cols = `Overweight - Obese`:`Lean - Obese`, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_global$group = factor(df_fig_global$group, 
                             levels = c("Overweight - Obese",
                                        "Lean - Obese"))
  
lo = floor(min(df_fig_global$value))
up = ceiling(max(df_fig_global$value))
mid = (lo + up)/2
fig_global = df_fig_global %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_global %>%
                                       dplyr::distinct(taxon, color) %>%
                                       .$color))
fig_global

## -----------------------------------------------------------------------------
res_pair = output$res_pair

df_fig_pair1 = res_pair %>%
    dplyr::filter(diff_bmioverweight == 1 |
                      diff_bmilean == 1 | 
                      diff_bmilean_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(diff_bmioverweight == 1, 
                                round(lfc_bmioverweight, 2), 0),
                  lfc2 = ifelse(diff_bmilean == 1, 
                                round(lfc_bmilean, 2), 0),
                  lfc3 = ifelse(diff_bmilean_bmioverweight == 1, 
                                round(lfc_bmilean_bmioverweight, 2), 0)) %>%
    tidyr::pivot_longer(cols = lfc1:lfc3, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_pair2 = res_pair %>%
    dplyr::filter(diff_bmioverweight == 1 |
                      diff_bmilean == 1 | 
                      diff_bmilean_bmioverweight == 1) %>%
    dplyr::mutate(lfc1 = ifelse(passed_ss_bmioverweight == 1 & diff_bmioverweight == 1, 
                                "aquamarine3", "black"),
                  lfc2 = ifelse(passed_ss_bmilean == 1 & diff_bmilean == 1, 
                                "aquamarine3", "black"),
                  lfc3 = ifelse(passed_ss_bmilean_bmioverweight == 1 & diff_bmilean_bmioverweight == 1, 
                                "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc3, 
                        names_to = "group", values_to = "color") %>%
    dplyr::arrange(taxon)

df_fig_pair = df_fig_pair1 %>%
    dplyr::left_join(df_fig_pair2, by = c("taxon", "group"))

df_fig_pair$group = recode(df_fig_pair$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese",
                          `lfc3` = "Lean - Overweight")

df_fig_pair$group = factor(df_fig_pair$group, 
                          levels = c("Overweight - Obese",
                                     "Lean - Obese", 
                                     "Lean - Overweight"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
fig_pair = df_fig_pair %>%
    ggplot(aes(x = group, y = taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value, color = color), size = 4) +
    scale_color_identity(guide = FALSE) +
    labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
fig_pair

## -----------------------------------------------------------------------------
res_dunn = output$res_dunn

df_fig_dunn1 = res_dunn %>%
    dplyr::filter(diff_bmioverweight == 1 | 
                      diff_bmilean == 1) %>%
    dplyr::mutate(lfc1 = ifelse(diff_bmioverweight == 1, 
                                round(lfc_bmioverweight, 2), 0),
                  lfc2 = ifelse(diff_bmilean == 1, 
                                round(lfc_bmilean, 2), 0)) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_dunn2 = res_dunn %>%
    dplyr::filter(diff_bmioverweight == 1 | 
                      diff_bmilean == 1) %>%
    dplyr::mutate(lfc1 = ifelse(passed_ss_bmioverweight == 1 & diff_bmioverweight == 1, 
                                "aquamarine3", "black"),
                  lfc2 = ifelse(passed_ss_bmilean == 1 & diff_bmilean == 1, 
                                "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "color") %>%
    dplyr::arrange(taxon)

df_fig_dunn = df_fig_dunn1 %>%
    dplyr::left_join(df_fig_dunn2, by = c("taxon", "group"))

df_fig_dunn$group = recode(df_fig_dunn$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese")
df_fig_dunn$group = factor(df_fig_dunn$group, 
                          levels = c("Overweight - Obese",
                                     "Lean - Obese"))

lo = floor(min(df_fig_dunn$value))
up = ceiling(max(df_fig_dunn$value))
mid = (lo + up)/2
fig_dunn = df_fig_dunn %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
    scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_dunn

## -----------------------------------------------------------------------------
res_trend = output$res_trend

df_fig_trend = res_trend %>%
    dplyr::filter(diff_abn == 1) %>%
    dplyr::mutate(lfc1 = round(lfc_bmioverweight, 2),
                  lfc2 = round(lfc_bmilean, 2),
                  color = ifelse(passed_ss == 1, "aquamarine3", "black")) %>%
    tidyr::pivot_longer(cols = lfc1:lfc2, 
                        names_to = "group", values_to = "value") %>%
    dplyr::arrange(taxon)

df_fig_trend$group = recode(df_fig_trend$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese")

df_fig_trend$group = factor(df_fig_trend$group, 
                          levels = c("Overweight - Obese", 
                                     "Lean - Obese"))

lo = floor(min(df_fig_trend$value))
up = ceiling(max(df_fig_trend$value))
mid = (lo + up)/2
fig_trend = df_fig_trend %>%
    ggplot(aes(x = group, y = taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
    labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(color = df_fig_trend %>%
                                         dplyr::distinct(taxon, color) %>%
                                         .$color))
fig_trend

## -----------------------------------------------------------------------------
data(dietswap, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)
print(tse)

## -----------------------------------------------------------------------------
set.seed(123)
# It should be noted that we have set the number of bootstrap samples (B) equal 
# to 10 in the 'trend_control' function for computational expediency. 
# However, it is recommended that users utilize the default value of B, 
# which is 100, or larger values for optimal performance.
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = "nationality + timepoint + group",
                  rand_formula = "(timepoint | subject)",
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "group", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2),
                                       solver = "ECOS",
                                       B = 10))

## -----------------------------------------------------------------------------
res_prim = output$res %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_prim[1:6, ] %>%
    datatable(caption = "ANCOM-BC2 Primary Analysis")

## -----------------------------------------------------------------------------
res_global = output$res_global %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_global[1:6, ] %>%
    datatable(caption = "ANCOM-BC2 Global Test")

## -----------------------------------------------------------------------------
res_pair = output$res_pair %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_pair[1:6, ] %>%
    datatable(caption = "ANCOM-BC2 Multiple Pairwise Comparisons")

## -----------------------------------------------------------------------------
res_dunn = output$res_dunn %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_dunn[1:6, ] %>%
    datatable(caption = "ANCOM-BC2 Dunnett's Type of Test")

## -----------------------------------------------------------------------------
res_trend = output$res_trend %>%
    mutate_if(is.numeric, function(x) round(x, 2))
res_trend[1:6, ] %>%
    datatable(caption = "ANCOM-BC2 Pattern Analysis")

## -----------------------------------------------------------------------------
bias_correct_log_table = output$bias_correct_log_table
# By default, ANCOM-BC2 does not add pseudo-counts to zero counts, which can 
# result in NAs in the bias-corrected log abundances. Users have the option to 
# either leave the NAs as they are or replace them with zeros. 
# This replacement is equivalent to adding pseudo-counts of ones to the zero counts. 
bias_correct_log_table[is.na(bias_correct_log_table)] = 0
# Show the first 6 samples
round(bias_correct_log_table[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log abundances")

## ----sessionInfo, message = FALSE, warning = FALSE, comment = NA--------------
sessionInfo()

