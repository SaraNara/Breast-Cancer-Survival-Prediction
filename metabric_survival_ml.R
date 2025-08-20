# ============================================================
# Project: Breast Cancer Survival Prediction (Kaggle - METABRIC)
# Model:   Cox Proportional Hazards + LASSO (glmnet)
# Author:  Saranya Narayana
# ============================================================

# ---- Packages ----
# install.packages(c("glmnet","survival","survminer","dplyr","readr","stringr","ggplot2"))
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(reshape2)
})

set.seed(123)


# ---- Paths ----
#data_path    <- "data/metabric.csv"     # <--- place your Kaggle CSV here
# results_dir  <- "results"
# if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# STEP 1 ---- Load data ----
df <- read.csv("/home/sara/Proj_2_BC/METABRIC_RNA_Mutation.csv", sep=",", header=TRUE)
dim(df)#1904  693
str(df)

#chnge column names to lower case and replace spaces with underscores
names(df) <- tolower(gsub("\\s+", "_", names(df)))


# STEP 2 ---- Exploratory Data Analysis (EDA) of Clinical attributes ----
#extract clinical attributes
df_clinical <- select(df, (1:31))
dim(df_clinical) # 1904  31
head(df_clinical)
str(df_clinical)
colnames(df_clinical)
summary(df_clinical[, 2:31])
df_clinical$overall_survival <- factor(df_clinical$overall_survival)

# Check for missing values
missing_values <- colSums(is.na(df_clinical))
missing_values

# select certain biological relavant clinical attributes
df_clinical_1 <- select(df_clinical, c(1,2,7,11,14,16,20,21,22,24,25,27,29))
dim(df_clinical_1) # 1904  13
head(df_clinical_1)
colnames(df_clinical_1)

#overal prevalence of overall survival
table(df_clinical_1$overall_survival)

#    0    1 
# 1103  801 
prop.table(table(df_clinical_1$overall_survival))
#         0         1 
# 0.5793067 0.4206933 

prop.table(table(df_clinical_1$overall_survival,df_clinical_1$chemotherapy))

  #            0          1
  # 0 0.46796218 0.11134454
  # 1 0.32405462 0.09663866

#survival by all therapies
df_clinical_1 %>%
  group_by(overall_survival, chemotherapy, hormone_therapy, radio_therapy) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

#check the distribution of some numerical attributes:age at diagnosis,mutation_count, nottingham_prognostic_index,overall_survival_months, tumor size, lymph_nodes_examined_positive
head(df_clinical)
str(df_clinical)

# Reshape to long format for plotting

df_clinical_long <- melt(df_clinical,
                id.vars = "overall_survival", 
                measure.vars = c("age_at_diagnosis", "lymph_nodes_examined_positive", "mutation_count", "nottingham_prognostic_index",
                "overall_survival_months", "tumor_size"),
                variable.name = "Attribute",
                value.name = "Value")

# Density plot faceted by attribute
ggplot(df_clinical_long, aes(x = Value, fill = factor(overall_survival))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Attribute, scales = "free") +
  labs(fill = "Overall Survival",
       x = "Value",
       y = "Density",
       title = "Density plots of clinical attributes by Overall Survival status") +
  theme_minimal()

#check correlation of clinical attributes with overall survival
# Select numeric clinical attributes
num_vars <- sapply(df_clinical_1, is.numeric)
num_vars[c("overall_survival", "patient_id")] <- FALSE  # exclude OS itself

df_numeric <- df_clinical_1[, num_vars]

# Compute correlation with overall_survival
cor_with_os <- sapply(df_numeric, function(x) cor(x, df_clinical_1$overall_survival, use="complete.obs"))
cor_with_os





# STEP 3 ---- Multivariable Survival analysis with Clinical attributes only ----
# convert categorical variables to factors
str(df_clinical_1)
cat_vars <- c("stage", "chemotherapy", "hormone_therapy", "radiation")  # replace with your columns
df_clinical_1[cat_vars] <- lapply(df_clinical_1[cat_vars], factor)










# ---- Helper: find survival time and event columns across variants ----
find_first <- function(cands, nm) {
  cand <- cands[cands %in% nm]
  if (length(cand) == 0) NA_character_ else cand[1]
}

nm <- names(df)

# Common time columns seen in METABRIC-like sets
time_col <- find_first(c("overall_survival_months","os_months","survival_months",
                         "months","time","survival_time"), nm)

# Common event/status columns (1/0, Dead/Alive, Deceased/Alive, True/False)
event_col <- find_first(c("overall_survival","os_status","event","status","vital_status"), nm)

if (is.na(time_col) || is.na(event_col)) {
  stop("Could not automatically detect survival time/event columns. ",
       "Please rename your time column to 'overall_survival_months' and event column to 'overall_survival'.")
}

message(sprintf("Using time column:  %s", time_col))
message(sprintf("Using event column: %s", event_col))

# ---- Clean survival outcome ----
df <- df %>%
  mutate(
    .time  = suppressWarnings(as.numeric(.data[[time_col]])),
    .event_raw = .data[[event_col]]
  )

# Convert event to numeric {0,1}
if (is.character(df$.event_raw) || is.factor(df$.event_raw)) {
  ev <- tolower(as.character(df$.event_raw))
  # map common labels to 0/1
  df$.event <- case_when(
    ev %in% c("dead","deceased","1","yes","true","event","died") ~ 1,
    ev %in% c("alive","0","no","false","censored","living") ~ 0,
    TRUE ~ NA_real_
  )
} else {
  # assume already numeric/binary
  df$.event <- as.numeric(df$.event_raw)
}

# Basic filtering
df <- df %>%
  filter(!is.na(.time), !is.na(.event), .time > 0)

# ---- Feature matrix (genes + numeric covariates) ----
# Exclude id-like and obvious non-feature columns
exclude_cols <- unique(c(
  time_col, event_col, ".time", ".event", ".event_raw",
  # common non-feature columns you might see; add/remove as needed
  "patient_id","sample_id","id","tumor_id","er_status","pr_status","her2_status"
))

# Keep numeric columns only
num_cols <- nm[ sapply(df, is.numeric) ]
feat_cols <- setdiff(num_cols, exclude_cols)

if (length(feat_cols) < 10) {
  stop("Not enough numeric feature columns detected. ",
       "Ensure your CSV includes gene/numeric features.")
}

# Optional: add a few clinical covariates if present and numeric
clin_cands <- c("age_at_diagnosis","age","tumor_stage","stage","tumor_size","n_lymph_nodes")
clin_keep <- intersect(clin_cands, feat_cols)  # already numeric subset





# Feature selection: keep top-variance features to speed training
# (Keep up to 1000, but include any clinical numeric covariates)
gene_only <- setdiff(feat_cols, clin_keep)
keep_g <- min(length(gene_only), 1000)
var_rank <- order(apply(df[, gene_only, drop=FALSE], 2, stats::var, na.rm=TRUE), decreasing = TRUE)
top_genes <- gene_only[var_rank[seq_len(keep_g)]]
final_feats <- unique(c(clin_keep, top_genes))

# Impute NAs with column medians (quick + simple)
median_impute <- function(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v }
X <- as.matrix( dplyr::mutate(df[, final_feats, drop=FALSE], across(everything(), median_impute)) )

# ---- Train/Test split (80/20) ----
n <- nrow(df)
idx <- sample(seq_len(n), size = floor(0.8*n))
X_train <- X[idx, , drop=FALSE]
X_test  <- X[-idx, , drop=FALSE]
y_train <- Surv(df$.time[idx], df$.event[idx])
y_test  <- Surv(df$.time[-idx], df$.event[-idx])

# ---- LASSO-Cox with cross-validation ----
set.seed(123)
cvfit <- cv.glmnet(X_train, y_train, family = "cox", alpha = 1, nfolds = 5)
png(file.path(results_dir, "cv_glmnet_plot.png"), width = 900, height = 700, res = 130)
plot(cvfit)
dev.off()

best_lambda <- cvfit$lambda.min
fit <- glmnet(X_train, y_train, family = "cox", alpha = 1, lambda = best_lambda)

# Coefficients (non-zero)
coef_df <- data.frame(
  feature = rownames(coef(fit)),
  coef = as.numeric(coef(fit))
) %>% filter(coef != 0) %>% arrange(desc(abs(coef)))

write_csv(coef_df, file.path(results_dir, "lasso_nonzero_coefficients.csv"))

# ---- Predict risk and evaluate ----
risk_train <- as.numeric(predict(fit, newx = X_train, type = "link"))
risk_test  <- as.numeric(predict(fit, newx = X_test,  type = "link"))

# Concordance (C-index)
c_train <- survConcordance(y_train ~ risk_train)$concordance
c_test  <- survConcordance(y_test  ~ risk_test )$concordance
message(sprintf("C-index (train): %.3f", c_train))
message(sprintf("C-index (test):  %.3f", c_test))

# Save C-index
write_lines(sprintf("C-index (train): %.3f\nC-index (test):  %.3f", c_train, c_test),
            file.path(results_dir, "cindex.txt"))

# ---- KM curves on Test set: split by median risk ----
test_df <- data.frame(
  time  = y_test[, "time"],
  event = y_test[, "status"],
  risk  = risk_test
)
test_df$risk_group <- ifelse(test_df$risk > median(test_df$risk, na.rm = TRUE), "High", "Low")

km_fit <- survfit(Surv(time, event) ~ risk_group, data = test_df)

p <- ggsurvplot(
  km_fit,
  data = test_df,
  pval = TRUE,
  risk.table = TRUE,
  title = "METABRIC (Test Set): Survival by ML-based Risk Group",
  legend.title = "Risk Group"
)

ggsave(filename = file.path(results_dir, "km_curve_test.png"),
       plot = p$plot, width = 8, height = 6, dpi = 300)

# Also save risk table as an image
ggsave(filename = file.path(results_dir, "km_curve_test_risktable.png"),
       plot = p$table, width = 8, height = 4, dpi = 300)

# ---- Optional: top coefficients bar plot ----
top_coef <- coef_df %>% slice_max(order_by = abs(coef), n = min(20, n()))
if (nrow(top_coef) > 0) {
  g <- ggplot(top_coef, aes(x = reorder(feature, coef), y = coef)) +
    geom_col() + coord_flip() +
    labs(title = "Top LASSO Cox Coefficients (Train)", x = "Feature", y = "Coefficient")
  ggsave(file.path(results_dir, "top_coefficients.png"), g, width = 8, height = 6, dpi = 300)
}

message("Done. Results saved in 'results/'")
