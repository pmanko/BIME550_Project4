library(data.table)

p <- function(x, name) {
  x[]
}

my_s <- function(levels, probs = NULL) {
  factor(sample(levels, 1000, replace = TRUE, prob=probs))
}

# ID

# Pre-populate from recurrance dataset
indeces <- sample(1:nrow(recurrance), 1000, replace = TRUE)

sample_data <- recurrance[indeces,]
sample_data[,index:=1:1000]

## Demographics and Patient Health

# Family History
c("yes", "no")

# Overall Health
c("poor", "medium", "good")

# BMI
c("low", "medium", "high")


## DIAGNOSIS
# Histology
# [ductal, lobular, mixed, metaplastic] [tubular, mucinous]

# ER  Status
# [+/-]


# PR Status
#[+/-]

# HER2 Status
#[+/-]

# Genetic Recurrence Score (21-gene panel) (page 17)
c("none", "low", "intermediate", "high")

## STAGING (pages 74 - 77)
# Primary Tumor Classification
c("T1a", "T1b", "T1c", "T2", "T3", "T4")

# Regional Lymph Nodes
c("N0", "N1", "N2", "N3")

# Pathologic Lymph Nodes
c("pNO", "pN1", "pN2", "pN3")

# Distal Metastasis
c("M0", "M1", "M")

# Stages
c("0", "IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")

# Histological Grade
c("G1", "G2", "G3")


## TREATMENT (page 13)
# Masectomy
c("none", "partial", "total")

# Lumpectomy
c("yes", "no")

# Breast Reconstruction
c("yes", "no")

# Radiation therapy
c("none", "whole_breast", "regional", "APBI")

# Endocrine Therapy
c("yes", "no")

# Chemotherapy
c("yes", "no")


## OUTCOMES

# Recurrence/remission
c("none", "local", "regional")


# Fertility

# Bayes Network
modelstr <- "[tumor_grade][metastasis][molecular_subtype][chemotherapy][age][tumor_size][inv_nodes][breast][menopause|age][node_caps|inv_nodes][breast_quad|breast][deg_malig|node_caps][irradiat|node_caps]
[recurrance|deg_malig]"