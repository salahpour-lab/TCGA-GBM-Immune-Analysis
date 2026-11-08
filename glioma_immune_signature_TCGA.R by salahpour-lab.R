# ==============================================================================
# Project: TCGA-GBM Immune Gene Analysis
# Description: Survival analysis and correlation of immune-fingerprint genes 
#              in Glioblastoma using TCGA PanCancer Atlas data.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. DATA DOWNLOAD & PREPARATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Step 1: Download the dataset directly from cBioPortal:
#         https://www.cbioportal.org/datasets
#         Search for "gbm_tcga_pan_can_atlas_2018" and download the .tar.gz file.
#
# Step 2: Place the downloaded 'gbm_tcga_pan_can_atlas_2018.tar.gz' file in your 
#         working directory. DO NOT extract (unzip) the file manually. 
#         R can read the data directly from the compressed .tar.gz archive.

# ------------------------------------------------------------------------------

# --Packages ----
# install.packages(c("survival", "survminer", "pheatmap", "ggplot2"))
library(survival)
library(survminer)
library(pheatmap)
library(ggplot2)

# 1. Set your working directory to the folder where the .tar.gz file is located.
#         (e.g., setwd("C:/Users/YourName/Desktop/glioma_R"))
setwd("C:/Users/YourName/Desktop/glioma_R")

# Define the name of the compressed file
study_file <- "gbm_tcga_pan_can_atlas_2018.tar.gz" 
out_dir   <- "results_GBM"
dir.create(out_dir, showWarnings = FALSE)

genes <- c("HMGB1", "CALR", "CD80", "CD86", "TNF", "ITGAX", "CD4", "CD8A", "CD8B")

# ---- 2. Expression matrix ----
# Create a temporary directory
temp_dir <- tempdir()

# IMPORTANT: Include the folder name "gbm_tcga_pan_can_atlas_2018/" inside the archive
untar(study_file, files = "gbm_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt", exdir = temp_dir)
untar(study_file, files = "gbm_tcga_pan_can_atlas_2018/data_clinical_patient.txt", exdir = temp_dir)

# Read the expression file (path now includes the subfolder)
expr_path <- file.path(temp_dir, "gbm_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem.txt")
expr <- read.delim(expr_path, check.names = FALSE, stringsAsFactors = FALSE)

expr <- expr[!is.na(expr$Hugo_Symbol) & expr$Hugo_Symbol != "", ]
expr <- expr[!duplicated(expr$Hugo_Symbol), ]
rownames(expr) <- expr$Hugo_Symbol
mat <- as.matrix(expr[, !(colnames(expr) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))])
mode(mat) <- "numeric"

missing_genes <- setdiff(genes, rownames(mat))
if (length(missing_genes) > 0) {
  warning("Genes not found: ", paste(missing_genes, collapse = ", "))
  genes <- setdiff(genes, missing_genes)
}

is_primary <- substr(colnames(mat), 14, 15) == "01"
mat <- mat[, is_primary]
cat("Primary tumour samples:", ncol(mat), "\n")

lmat <- log2(mat[genes, , drop = FALSE] + 1)

# ---- 3. Signature score ----
z     <- t(scale(t(lmat)))
score <- colMeans(z, na.rm = TRUE)

# ---- 4. Clinical data ----
# Read the clinical file from the temp directory (also inside the subfolder)
clin_path <- file.path(temp_dir, "gbm_tcga_pan_can_atlas_2018", "data_clinical_patient.txt")
clin <- read.delim(clin_path, skip = 4, stringsAsFactors = FALSE, check.names = FALSE)

clin$time   <- as.numeric(clin$OS_MONTHS)
clin$status <- as.integer(sub(":.*", "", clin$OS_STATUS))
keep <- intersect(c("PATIENT_ID", "time", "status", "AGE"), colnames(clin))
clin <- clin[, keep]

samples <- colnames(lmat)
sc <- data.frame(PATIENT_ID = substr(samples, 1, 12),
                 score      = as.numeric(score),
                 t(z),
                 check.names = FALSE)
rownames(sc) <- NULL
sc <- sc[!duplicated(sc$PATIENT_ID), ]

df <- merge(sc, clin, by = "PATIENT_ID")
df <- df[!is.na(df$time) & !is.na(df$status) & df$time > 0, ]
cat("Patients with survival data:", nrow(df), " | events:", sum(df$status), "\n")

# ---- 5. Correlation heatmap ----
cm <- cor(t(lmat), method = "spearman")
png(file.path(out_dir, "fig1_gene_correlation.png"), width = 1600, height = 1400, res = 250)
pheatmap(cm, display_numbers = TRUE, number_format = "%.2f",
         main = "Spearman correlation of immune-fingerprint genes (TCGA-GBM)")
dev.off()

# ---- 6. Kaplan-Meier (median split only) ----
df$group <- factor(ifelse(df$score > median(df$score), "High", "Low"),
                   levels = c("Low", "High"))
fit <- survfit(Surv(time, status) ~ group, data = df)
p_km <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                   xlab = "Months", ylab = "Overall survival",
                   legend.title = "Signature", legend.labs = c("Low", "High"))
png(file.path(out_dir, "fig2_KM_signature.png"), width = 1800, height = 1600, res = 250)
print(p_km)
dev.off()

# ---- 7. Cox regression ----
df$score_sd <- as.numeric(scale(df$score))
cox1 <- coxph(Surv(time, status) ~ score_sd, data = df)
print(summary(cox1))
if ("AGE" %in% colnames(df)) {
  df$AGE <- as.numeric(df$AGE)
  cox2 <- coxph(Surv(time, status) ~ score_sd + AGE, data = df)
  print(summary(cox2))
}

# ---- 8. Per-gene Cox + forest plot ----
res <- do.call(rbind, lapply(genes, function(g) {
  d <- df
  d$x <- as.numeric(scale(d[[g]]))
  s <- summary(coxph(Surv(time, status) ~ x, data = d))
  data.frame(gene  = g,
             HR    = s$coefficients[1, "exp(coef)"],
             lower = s$conf.int[1, "lower .95"],
             upper = s$conf.int[1, "upper .95"],
             p     = s$coefficients[1, "Pr(>|z|)"])
}))
res$padj <- p.adjust(res$p, method = "BH")
print(res)
write.csv(res, file.path(out_dir, "per_gene_cox.csv"), row.names = FALSE)

p_forest <- ggplot(res, aes(x = HR, y = reorder(gene, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  scale_x_log10() +
  labs(x = "Hazard ratio per 1 SD (95% CI)", y = NULL,
       title = "Gene-level association with overall survival (TCGA-GBM)") +
  theme_classic(base_size = 14)
ggsave(file.path(out_dir, "fig3_forest.png"), p_forest, width = 7, height = 5, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Done. Outputs in:", out_dir, "\n")
