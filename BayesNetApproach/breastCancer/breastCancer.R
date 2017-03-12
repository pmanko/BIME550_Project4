library(bnlearn)
library(graph)
library(Rgraphviz)
library(RBGL)
library(gRain)

# Load Datasets:

## Wisconsin


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

### wpbc
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


## UCI

## BCSC

### Risk Estimate

#### Multiple (risk.txt)

#### Single (risk_rand.txt)

### Risk Factor

### Mamogram


