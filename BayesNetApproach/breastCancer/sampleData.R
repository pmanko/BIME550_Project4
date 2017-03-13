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


# age
#sample_data[,age:=round(rnorm(1000, mean=50, sd=10))]
# blood pressure
sample_data[,systolic_bp:=round(rnorm(1000, mean=130, sd=10))]
sample_data[,diastolic_bp:=systolic_bp - round(rnorm(1000, mean=40, sd=5))]

# ethicity
ethnicityL <- c("nh_white", "nh_black", "asian_pi", "native_american", "hispanic", "other")
sample_data[,ethnicity:=sample(ethnicityL, 1000, replace=TRUE, prob=c(.4, .3, .1, .05, .1, .05))]

# menopause


# Targets
# recurrance
# 


# Bayes Network
modelstr <- "[tumor_grade][metastasis][molecular_subtype][chemotherapy][age][tumor_size][inv_nodes][breast][menopause|age][node_caps|inv_nodes][breast_quad|breast][deg_malig|node_caps][irradiat|node_caps]
[recurrance|deg_malig]"