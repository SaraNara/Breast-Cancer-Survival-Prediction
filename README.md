## Title: Breast Cancer Survival Analaysis & ML prediction

## Subtitle: Integrative Survival Prediction in Breast Cancer Using Clinical and mRNA sequence Gene Expression Data with Random Forest Machine Learning (Kaggle - METABRIC)

## Project Overview
This project aims to predict overall survival in breast cancer patients by integrating **clinical attributes** and **gene expression data** using statistical modeling and machine learning approaches. The analysis uses the **METABRIC dataset** (Kaggle) to explore risk factors and develop predictive models.

---

## Objectives
- Perform **exploratory data analysis (EDA)** on clinical variables.
- Conduct **survival analysis** using:
- Cox proportional hazards models
- Kaplan-Meier curves
- Combine **clinical and gene expression data** for prediction using **Random Forest classifier**.
- Evaluate model performance using **ROC curves**, **confusion matrix**, and **feature importance**.

---

## Data
- **Source:** [METABRIC Breast Cancer Dataset on Kaggle](https://www.kaggle.com/datasets)
- **Clinical Attributes:** Age at diagnosis, tumor size, lymph nodes, treatment variables, survival outcome.
- **Gene Expression:** RNA-seq data with multiple gene features.
---
## Project Structure

```
Breast-Cancer-Survival-Prediction/
├── data/           # Raw and processed data files (e.g., METABRIC dataset)
├── results/        # Output plots, model results, and summaries
├── scripts/        # R scripts for analysis and modeling
├── README.md       # Project documentation
```


## Requirements
- R >= 4.0
- Packages:
 glmnet, survival, survminer, dplyr, readr, stringr, ggplot2, reshape2,
  randomForest, caret, pROC

## How to Run

Clone the repository:

git clone https://github.com/SaranyaNarayana/Breast-Cancer-Survival-Analysis-and-ML-Prediction.git

cd Breast-Cancer-Survival-Analysis-and-ML-Prediction


Place the METABRIC dataset CSV file in the data/ folder.

Run the analysis:

source("Breastcancer_analysis.R")


All plots and outputs will be saved in the results/ folder.

## Outputs

EDA Plots: Density plots of clinical attributes by survival.

Survival Analysis: Kaplan-Meier curves, Cox model summaries.

Machine Learning: Random Forest classification, ROC curves, top 20 important features.

## Notes

The code uses reproducible random seeds for train/test splits.

Binary variables are converted to factors where necessary for modeling.

Ensure all required R packages are installed before running the script.

# References

https://pubmed.ncbi.nlm.nih.gov/37370847/

https://pubmed.ncbi.nlm.nih.gov/33681976/

https://www.researchgate.net/publication/377025287_Breast_Cancer_Detection_using_Machine_Learning_Algorithms

https://pmc.ncbi.nlm.nih.gov/articles/PMC11152379/




# Author
Saranya Narayana
(Ph.D. in Quantitative Genetics and Epidemiology)
