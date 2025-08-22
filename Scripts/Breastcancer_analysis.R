# ============================================================
# Author:  Saranya Narayana
# Project: Breast Cancer Survival & ML Analysis:Integrative Survival Prediction in Breast Cancer Using Clinical and Gene Expression Data with Random Forest Machine Learning (Kaggle - METABRIC)
# Date:   2025-08-19
# Purpose: Perform survival analysis using clinical and gene expression data from METABRIC dataset
#         using Random Forest classifier for prediction.
#         Includes EDA, survival analysis, and machine learning model training.
#         Visualizations include density plots, Kaplan-Meier curves, and ROC curves.
#         Top 20 important features are identified.
# Dataset: METABRIC (Molecular Taxonomy of Breast Cancer International Consortium)
# Source: https://www.kaggle.com/datasets/uciml/metabric
# License: CC0: Public Domain
# Note: Ensure to have the required packages installed and the dataset available in the specified path.
# 
# ============================================================

#-----------------------------------
# Logs all output and saves plots
#-----------------------------------

# Install required packages if not already installed
# install.packages(c("glmnet","survival","survminer","dplyr","readr","stringr","ggplot2","reshape2","randomForest","caret","pROC"))

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(reshape2)
  library(randomForest)
  library(caret)
  library(pROC)
})

set.seed(123)

#-----------------------------------
# Create results folder
#-----------------------------------
if(!dir.exists("results")) dir.create("results")

#-----------------------------------
# Logging
#-----------------------------------
log_file <- "/home/sara/Proj_2_BC/results/run_log.txt"

# Open connection for logging both stdout and messages
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)         # redirect stdout
sink(log_con, type = "message")     # redirect messages

cat("===== Analysis Started =====\n")

#-----------------------------------
# STEP 1: Load Data
#-----------------------------------
df <- read.csv("/home/sara/Proj_2_BC/data/METABRIC_RNA_Mutation.csv", sep=",", header=TRUE)
cat("Data dimensions:", dim(df), "\n")
cat("column preview:\n")
head(df[, 1:10])  # Display first 10 columns

# Clean column names
names(df) <- tolower(gsub("\\s+", "_", names(df)))

#-----------------------------------
# STEP 2: Exploratory Data Analysis (Clinical)
#-----------------------------------
df_clinical <- df %>% select(1:31) #extract clinical attributes only
cat("Clinical data dimensions:", dim(df_clinical), "\n")

df_clinical$overall_survival <- factor(df_clinical$overall_survival) # Convert to factor

# Select certain biologically relevant clinical features
df_clinical_1 <- df_clinical %>% select(patient_id, overall_survival, age_at_diagnosis,
                                        tumor_size,lymph_nodes_examined_positive,
                                        chemotherapy, hormone_therapy, radio_therapy,
                                        mutation_count, nottingham_prognostic_index,
                                        overall_survival_months)
cat("Clinical data dimensions:", dim(df_clinical_1), "\n")
cat("\nMissing values per column:\n")
print(colSums(is.na(df_clinical_1)))


#remove rows with missing values for simplicity
df_clinical_2 <- na.omit(df_clinical_1)
cat("Clinical data dimensions after removing NAs:", dim(df_clinical_2), "\n")

# Overall survival distribution
cat("\nOverall survival distribution:\n")
print(table(df_clinical_2$overall_survival))

cat("\nOverall survival proportions:\n")
print(prop.table(table(df_clinical_2$overall_survival)))

cat("\nOverall survival by chemotherapy:\n")
print(prop.table(table(df_clinical_2$overall_survival, df_clinical_2$chemotherapy)))    

cat("\nOverall survival by all therapies:\n")
print(df_clinical_2 %>%
        group_by(overall_survival, chemotherapy, hormone_therapy, radio_therapy) %>%
        summarise(count = n()) %>%
        mutate(prop = count / sum(count)))



# Density plots of clinical features
df_long <- melt(df_clinical_2, id.vars = "overall_survival",
                measure.vars = c("age_at_diagnosis", "lymph_nodes_examined_positive", 
                                 "mutation_count", "nottingham_prognostic_index",
                                 "overall_survival_months", "tumor_size"),
                variable.name = "Attribute",
                value.name = "Value")

plot_density <- ggplot(df_long, aes(x = Value, fill = factor(overall_survival))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Attribute, scales = "free") +
  labs(fill = "Overall Survival", x = "Value", y = "Density",
       title = "Density plots of clinical attributes by Overall Survival") +
    theme_classic() 

ggsave("results/density_clinical_survival.png", plot_density, width = 10, height = 8, dpi = 300)
cat("Density plot saved.\n")


#check correlation of clinical attributes with overall survival
num_vars <- sapply(df_clinical_2, is.numeric)
num_vars[c("overall_survival", "patient_id")] <- FALSE  # exclude OS itself

df_numeric <- df_clinical_2[, num_vars]
df_clinical_2$overall_survival <- as.numeric(df_clinical_2$overall_survival)

# Compute correlation with overall_survival
cor_with_os <- sapply(df_numeric, function(x) cor(x, df_clinical_2$overall_survival, use="complete.obs"))
cat("\ncorrelation of clinical attributes with overall_survival:\n")
print(cor_with_os)




#-----------------------------------
# STEP 3: Cox Survival Analysis (Clinical Only)
#-----------------------------------
df_complete <- df_clinical_2 %>%
  select(patient_id, age_at_diagnosis, tumor_size, lymph_nodes_examined_positive,
         chemotherapy, hormone_therapy, radio_therapy, overall_survival_months, overall_survival) %>%
  na.omit()
cat("df_complete data dimensions:", dim(df_complete), "\n")


# Convert categorical vars to factors
cat_vars <- c("chemotherapy", "hormone_therapy", "radio_therapy")
df_complete[cat_vars] <- lapply(df_complete[cat_vars], factor)

# Survival object
df_complete$overall_survival <- ifelse(df_complete$overall_survival == 1, 1, 0)
surv_obj <- Surv(time = df_complete$overall_survival_months, event = df_complete$overall_survival)

# Fit Cox model
cox_model <- coxph(surv_obj ~ age_at_diagnosis + tumor_size + lymph_nodes_examined_positive +
                     chemotherapy + hormone_therapy + radio_therapy,
                   data = df_complete)
cat("\nCox Model Summary:\n")
print(summary(cox_model))

# Kaplan-Meier plot by chemotherapy
fit_km <- survfit(surv_obj ~ chemotherapy, data = df_complete)
km_plot <- ggsurvplot(fit_km, data = df_complete,
                      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
                      title = "Survival by Chemotherapy Status")
ggsave("results/km_clinical_chemotherapy_full.png", km_plot$plot, width = 8, height = 6, dpi = 300)
cat("KM plot saved.\n")

# Risk groups by median linear predictor
df_complete$risk_score <- predict(cox_model, type = "lp")
df_complete$risk_group <- ifelse(df_complete$risk_score > median(df_complete$risk_score), "High", "Low")

fit_risk <- survfit(surv_obj ~ risk_group, data = df_complete)
risk_plot <- ggsurvplot(fit_risk, data = df_complete, risk.table = TRUE, pval = TRUE,
                        title = "High vs Low Risk Survival")
ggsave("results/km_risk_group_full2.png", risk_plot$plot, width = 8, height = 6, dpi = 300)
cat("KM risk group plot saved.\n")



#-----------------------------------
# STEP 4: Random Forest (Clinical + Gene Expression)
#-----------------------------------
# Extract gene expression features
df_gene <- df %>% filter(patient_id %in% df_complete$patient_id) %>% select(32:ncol(df))
df_gene <- as.data.frame(df_gene)
cat("Gene expression data dimensions:", dim(df_complete), "\n")


# Combine clinical numeric features + gene expression
df_clinical_num <- df_complete %>% select(age_at_diagnosis, tumor_size, lymph_nodes_examined_positive, chemotherapy, hormone_therapy, radio_therapy)
X <- cbind(df_clinical_num, df_gene)
y <- factor(df_complete$overall_survival, levels = c(0,1))

# Train/test split
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]
X_test <- X[-train_idx, ]
y_train <- y[train_idx]
y_test <- y[-train_idx]


# Fit Random Forest
rf_model <- randomForest(x = X_train, y = y_train, ntree = 500, importance = TRUE)
cat("\nRandom Forest Model:\n")
print(rf_model)

# Predictions
y_pred_class <- predict(rf_model, X_test)
y_pred_prob  <- predict(rf_model, X_test, type="prob")[,2]

# Confusion matrix
cm <- confusionMatrix(y_pred_class, y_test, positive="1")
cat("\nConfusion Matrix:\n")
print(cm)

# ROC & AUC
roc_obj <- roc(y_test, y_pred_prob, levels=c("0","1"), direction="<")
auc_val <- auc(roc_obj)
cat("\nAUC =", auc_val, "\n")

# ROC plot
roc_plot <- ggroc(roc_obj, color = "blue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  ggtitle("ROC Curve - RF Clinical + Gene Expression") +
  theme_classic() +
  coord_equal() +  # Make x and y scale equal
  xlim(0, 1) + ylim(0, 1)  # Force 0-1 scale for both axes

# Save the plot
ggsave("results/ROC.png", roc_plot, width = 8, height = 6, dpi = 300)
cat("ROC plot saved.\n")


# Top 20 important variables
var_imp <- importance(rf_model)

var_imp_df <- data.frame(Variable = rownames(var_imp), MDA = var_imp[,3], Gini = var_imp[,4])

top20 <- var_imp_df %>% arrange(desc(MDA)) %>% head(20)
cat("\nTop 20 Important Variables:\n")
print(top20)

top20_plot <- ggplot(top20, aes(x = reorder(Variable, MDA), y = MDA, fill = MDA)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Top 20 Important Variables (RF)", x="Variable", y="MeanDecreaseAccuracy") +
  theme_classic()
ggsave("results/rf_top20_variables.png", top20_plot, width = 8, height = 6, dpi = 300)
cat("Top 20 variable plot saved.\n")

#-----------------------------------
# Finish logging
#-----------------------------------
cat("===== Analysis Finished =====\n")
sink(type = "message")
sink()
close(log_con)
