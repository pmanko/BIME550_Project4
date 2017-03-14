library(bnlearn)
library(graph)
library(Rgraphviz)
library(RBGL)
library(gRain)
library(mice)
library(data.table)

# Load Datasets:

## Wisconsin

# https://shiring.github.io/machine_learning/2017/01/15/rfe_ga_post

### breast-cancer-wisconsin
#  Attribute                     Domain
# -- -----------------------------------------
#   1. Sample code number            id number
# 2. Clump Thickness               1 - 10
# 3. Uniformity of Cell Size       1 - 10
# 4. Uniformity of Cell Shape      1 - 10
# 5. Marginal Adhesion             1 - 10
# 6. Single Epithelial Cell Size   1 - 10
# 7. Bare Nuclei                   1 - 10
# 8. Bland Chromatin               1 - 10
# 9. Normal Nucleoli               1 - 10
# 10. Mitoses                       1 - 10
# 11. Class:                        (2 for benign, 4 for malignant)

wisconsin <- read.table("data/wisconsin_data/breast-cancer-wisconsin.data", header = FALSE, sep = ",")
colnames(wisconsin) <- c("sample_code_number", "clump_thickness", "uniformity_of_cell_size", "uniformity_of_cell_shape", "marginal_adhesion", "single_epithelial_cell_size", "bare_nuclei", "bland_chromatin", "normal_nucleoli", "mitosis", "classes")
wisconsin$classes <- ifelse(wisconsin$classes == "2", "benign", ifelse(wisconsin$classes == "4", "malignant", NA))
wisconsin[wisconsin == "?"] <- NA

# impute missing data
wisconsin[,2:10] <- apply(wisconsin[, 2:10], 2, function(x) as.numeric(as.character(x)))
dataset_impute <- mice(wisconsin[, 2:10],  print = FALSE)
wisconsin <- cbind(wisconsin[, 11, drop = FALSE], mice::complete(dataset_impute, 1))
wisconsin$classes <- as.factor(wisconsin$classes)
wisconsin <- as.data.table(wisconsin)

# 1) ID number
# 2) Outcome (R = recur, N = nonrecur)
# 3) Time (recurrence time if field 2 = R, disease-free time if 
#          field 2	= N)
# 4-33) Ten real-valued features are computed for each cell nucleus:
#   
#   a) radius (mean of distances from center to points on the perimeter)
# b) texture (standard deviation of gray-scale values)
# c) perimeter
# d) area
# e) smoothness (local variation in radius lengths)
# f) compactness (perimeter^2 / area - 1.0)
# g) concavity (severity of concave portions of the contour)
# h) concave points (number of concave portions of the contour)
# i) symmetry 
# j) fractal dimension ("coastline approximation" - 1)
# 
# Several of the papers listed above contain detailed descriptions of
# how these features are computed. 
# 
# The mean, standard error, and "worst" or largest (mean of the three
#                                                   largest values) of these features were computed for each image,
# resulting in 30 features.  For instance, field 4 is Mean Radius, field
# 14 is Radius SE, field 24 is Worst Radius.
# 
# Values for features 4-33 are recoded with four significant digits.
# 
# 34) Tumor size - diameter of the excised tumor in centimeters
# 35) Lymph node status - number of positive axillary lymph nodes
# observed at time of surgery


### wdbc
wdbc <- as.data.table(read.table("data/wisconsin_data/wdbc.data", header = FALSE, sep = ","))
phenotypes <- rep(c("radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity", "concave_points", "symmetry", "fractal_dimension"), 3)
types <- rep(c("mean", "se", "largest_worst"), each = 10)
colnames(wdbc) <- c("ID", "diagnosis", paste(phenotypes, types, sep = "_"))

### wpbc

wpbc <- read.table("data/wisconsin_data/wpbc.data", header = FALSE, sep = ",")
colnames(wpbc) <- c("ID", "outcome", "time", paste(phenotypes, types, sep = "_"), "tumor_size", "lymph_node_status")
wpbc[wpbc == "?"] <- NA

# impute missing data
wpbc[,3:35] <- apply(wpbc[,3:35], 2, function(x) as.numeric(as.character(x)))
dataset_impute <- mice(wpbc[,3:35],  print = FALSE)
wpbc <- cbind(wpbc[, 2, drop = FALSE], mice::complete(dataset_impute, 1))
wpbc <- as.data.table(wpbc)

## UCI

# Class: no-recurrence-events, recurrence-events
# 
# age: 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80-89, 90-99.
# 
# menopause: lt40, ge40, premeno.
# 
# tumor-size: 0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59.
# 
# inv-nodes: 0-2, 3-5, 6-8, 9-11, 12-14, 15-17, 18-20, 21-23, 24-26, 27-29, 30-32, 33-35, 36-39.
# 
# node-caps: yes, no.
# 
# deg-malig: 1, 2, 3.
# 
# breast: left, right.
# 
# breast-quad: left-up, left-low, right-up, right-low, central.
# 
# irradiat: yes, no


recurrance <- as.data.table(read.csv("data/datasets-uci-breast-cancer.csv", header=FALSE, quote="\'"))
colnames(recurrance) <- c("age", "menopause", "tumor_size", "inv_nodes", "node_caps", "deg_malig", "breast", "breast_quad", 'irradiat', 'recurrance')

# Fix nans manually
recurrance[node_caps=='nan']
recurrance[node_caps=='nan' & deg_malig == "3", node_caps:="yes"]
recurrance[node_caps=='nan' & irradiat == "yes", node_caps:="yes"]
recurrance[node_caps=='nan', node_caps:="no"]
recurrance$node_caps <- factor(recurrance$node_caps)

recurrance[breast_quad=='nan', breast_quad:="left_up"]
recurrance$breast_quad <- factor(recurrance$breast_quad)

recurrance$deg_malig <- factor(recurrance$deg_malig)



## BCSC
### Risk Estimate

#### Single (risk_rand.txt)
# Variable Name Columns Coding
# 1 menopaus 1 0 = premenopausal; 1 = postmenopausal or age>=55 ; 9 = unknown
# 2 agegrp 3-4 1 = 35-39; 2 = 40-44; 3 = 45-49; 4 = 50-54; 5 = 55-59; 6 = 60-64; 7 = 65-69; 8 = 70-74; 9 = 75-79; 10 = 80-84
# 3 density 6
# BI-RADS breast density codes 1 = Almost entirely fat; 2 = Scattered fibroglandular densities; 3 = Heterogeneously dense; 4 =
#   Extremely dense; 9 = Unknown or different measurement system
# 4 race 8 1 = white; 2 = Asian/Pacific Islander; 3 = black; 4 = Native American; 5 = other/mixed; 9 = unknown
# 5 Hispanic 10 0 = no; 1 = yes; 9 = unknown
# 6 bmi 12 Body mass index: 1 = 10-24.99; 2 = 25-29.99; 3 = 30-34.99; 4 = 35 or more; 9 = unknown
# 7 agefirst 14 Age at first birth: 0 = Age < 30; 1 = Age 30 or greater; 2 = Nulliparous; 9 = unknown
# 8 nrelbc 16 Number of first degree relatives with breast cancer: 0 = zero; 1= one; 2 = 2 or more; 9 = unknown
# 9 brstproc 18 Previous breast procedure: 0 = no; 1 = yes; 9 = unknown
# 10 lastmamm 20 Result of last mammogram before the index mammogram: 0 = negative; 1 = false positive; 9 = unknown
# 11 surgmeno 22 Surgical menopause: 0 = natural; 1 = surgical; 9 = unknown or not menopausal (menopaus=0 or menopaus=9)
# 12 hrt 24 Current hormone therapy: 0 = no; 1 = yes; 9 = unknown or not menopausal (menopaus=0 or menopaus=9)
# 13 invasive 26 Diagnosis of invasive breast cancer within one year of the index screening mammogram: 0 = no; 1 = yes
# 14 cancer 28 Diagnosis of invasive or ductal carcinoma in situ breast cancer within one year of the index screening mammogram: 0 = no; 1 = yes
# 15 training 30 Training data: 0 = no (validation); 1 = yes (training)
# 16 count 32-37 Frequency count of this combination of covariates and outcomes (all variables 1 to 15)


risk <- as.data.table(read.table("data/risk_dataset_v2/risk_rand.txt", header = FALSE, sep=''))
colnames(risk) <- c("menopaus", "agegrp", "density", "race", "Hispanic", "bmi", "agefirst", "nrelbc", "brstproc", 'lastmamm', 'surgmeno', 'hrt', 'invasive', 'cancer', 'training', 'count')

### Risk Factor
# year	Calendar year of observation	Numerical, 2000-2009
# age_group_5_years	Age (years) in 5 year groups	1 = Age 18-29; 2 = Age 30-34; 3 = Age 35-39; 4 = Age 40-44; 5 = Age 45-49; 6 = Age 50-54; 7 = Age 55-59; 8 = Age 60-64; 9 = Age 65-69; 10 = Age 70-74; 11 = Age 75-79; 12 = Age 80-84; 13 = Age ≥85
# race_eth	Race/ethnicity	1 = Non-Hispanic white; 2 = Non-Hispanic black; 3 = Asian/Pacific Islander; 4 = Native
# American; 5 = Hispanic; 6 = Other/mixed; 9 = Unknown
# first_degree_hx	History of breast cancer in a first degree relative	0 = No; 1 = Yes; 9 = Unknown
# age_menarche	Age (years) at menarche	0 = Age ≥14; 1 = Age 12-13; 2 = Age <12; 9 = Unknown
# age_first_birth	Age (years) at first birth	0 = Age < 20; 1 = Age 20-24; 2 = Age 25-29; 3 = Age ≥30; 4 = Nulliparous; 9 = Unknown
# BIRADS_breast_density	BI-RADS breast density	1 = Almost entirely fat; 2 =
#   Scattered fibroglandular densities; 3 = Heterogeneously
# dense; 4 = Extremely dense; 9 = Unknown or different
# measurement system
# current_hrt	Use of hormone replacement therapy	0 = No; 1 = Yes; 9 = Unknown
# menopaus	Menopausal status	1 = Pre- or peri-menopausal; 2 = Post-menopausal; 3 = Surgical menopause; 9 = Unknown
# bmi_group	Body mass index	1 = 10-24.99; 2 = 25-29.99; 3 = 30-34.99; 4 = 35 or more; 9 = Unknown
# biophx	Previous breast biopsy or aspiration	0 = No; 1 = Yes; 9 = Unknown
# breast_cancer_history	Prior breast cancer diagnosis	0 = No; 1 = Yes; 9 = Unknown
# count	Frequency count of this combination of covariates	Numerical

risk_factor <- as.data.table(read.table("data/bcsc_rf_data/BCSC_risk_factors_summarized.csv", header = TRUE, sep=','))




