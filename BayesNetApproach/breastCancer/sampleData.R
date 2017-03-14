## Points to Make
# 1. In this analysis, we focused on discrete Bayes Networks. Given more time, we would have liked to dive into Gausian Bayesian Networks to account
#    for continuous variables, and hybrid Bayesian Networks for combos.
# 2. Having access to a dataset containing all of these variables would be incredibly useful for learning or evaluating graph structure, determining
#    determining posterior probabilities, etc. We struggled with figuring out what was the best approach - creating dummy data, or dummy conditional probabilites. 
#    Since we wanted to practice setting 
# 3. Due to the concepts of the markov blanket and conditional probability, the relationships between nodes become much more managable. Studies focues on
#    a ssubset of relationships can be used as basis for prior evidence to estimate the conditional probability tables. I understand much more clearly
#    how BNs are created in the real world
# 4. Using R to play around with these concepts was very rewarding. I read through some tutorials on BN software packages (Bayesia), but
#    the user interface obfuscates the concepts and makes understanding what's going on way more complicated. However, if I were comfortable
#    using BNs, I would might switch to a dedicated environment to facilitate actual analyses.
# 5. The problem of over-fitting came up with the risk and risk_factor datasets. Learning structure without tweaking the thresholds for, for example,
#    significance generated incredibly dense graphs whose structure was not very informative.
# 6. There's a sense of authenticity about the representation of the world in terms of BNs - real life seldom has absolute certainty, and many 
#    exceptions exist to the rule. However, when using BNs for decision support on critical issues, the uncertainty might be disconcerning in some cases.
# 7. The state space of possible graphs representing given concepts is incredibly large, and multiple viable representations are likely exist. It is diffuclt
#    to know which are optimal, and algorithms are not guaranteed to find the global maximum. Also, the directions of the arcs are not always intuitive - the meaning
#    of the direction is not entirely clear, since inference can work either way. I need to develop a better understanding of the underlying probability theory
#    in order to gain a true appreciation of BNs. 



library(data.table)
library(pathRender)

my_s <- function(levels, probs = NULL) {
  factor(sample(levels, 1000, replace = TRUE, prob=probs))
}

add_node <- function(name, levels, prob=NULL) {
  sample_data[[name]] <<- sample(levels, 1000, replace = TRUE)
  nodes <<- c(nodes, name)
}

# Pre-populate from recurrance dataset
indeces <- sample(1:nrow(recurrance), 1000, replace = TRUE)
sample_data <- recurrance[indeces,]

# Get initial nodes
nodes <- colnames(sample_data)[c(1:6, 9, 10)]

# ID
# sample_data[,index:=1:1000]

## Demographics and Patient Health

# Family History of Breast Cancer
add_node('family_history_bc', c('yes', 'no'))

# Overall Health
add_node('overall_health', c("poor", "medium", "good"))

# BMI
add_node('bmi', c("low", "medium", "high"))

## DIAGNOSIS
# Histology
add_node('tumor_histology', c('ductal', 'lobular', 'mixed', 'metaplastic', 'tubular', 'mucinous'))

# ER  Status
add_node('er_status', c('positive', 'negative'))

# PR Status
add_node('pr_status', c('positive', 'negative'))

# HER2 Status
add_node('her2_status', c('positive', 'negative'))

# Genetic Recurrence Score (21-gene panel) (page 17)
add_node('21_gene_assay_rec_risk', c("none", "low", "intermediate", "high"))

## STAGING (pages 74 - 77)
# Primary Tumor Classification
add_node('primary_tumor_class', c("T1a", "T1b", "T1c", "T2", "T3", "T4"))

# Regional Lymph Nodes
add_node('regional_lymph_class', c("N0", "N1", "N2", "N3"))

# Pathologic Lymph Nodes
add_node('pathologic_lymph_class', c("pNO", "pN1", "pN2", "pN3"))

# Distal Metastasis
add_node('metastasis_class', c("M0", "M1"))

# Stages
nodes <- c(nodes, 'cancer_stage')
sample_data[primary_tumor_class == 'T1a' | primary_tumor_class == "T1b", cancer_stage:='IA',]
sample_data[primary_tumor_class == 'T1c', cancer_stage:='IB',]
sample_data[primary_tumor_class %in% c("T1a", "T1b", "T1c") & regional_lymph_class == 'N1', cancer_stage:='IIA',]
sample_data[primary_tumor_class == 'T2', cancer_stage:='IIB',]
sample_data[regional_lymph_class == 'T3', cancer_stage:='IIIA',]
sample_data[regional_lymph_class == 'N2', cancer_stage:='IIIB',]
sample_data[regional_lymph_class == 'N3', cancer_stage:='IIIC',]
sample_data[metastasis_class == 'M1', cancer_stage:='IV',]

# Histological Grade
add_node('histological_grade', c("G1", "G2", "G3"))


## TREATMENT (page 13)
# Masectomy
add_node('masectomy', c("none", "partial", "total"))

# Lumpectomy
add_node('lumpectomy', c("yes", "no"))

# Breast Reconstruction
add_node('breast_reconstruction', c("yes", "no"))

# Radiation therapy
add_node('radiation_therapy', c("none", "whole_breast", "regional", "APBI"))

# Endocrine Therapy
add_node('endocrine_therapy', c("yes", "no"))

# Chemotherapy
add_node('chemotherapy', c("yes", "no"))


## OUTCOMES
# Primary Treatment Survival

# Recurrence/remission
# Covered by sample data
#c("none", "local", "regional")

## Long-term quality of Life  

# Physical Function
# Body Image
# Sexual Function
# Coping
# Cognitive Function
# Social Support
# Anxiety
add_node('physical_function', c('low', 'medium', 'high'))
add_node('body_image', c('low', 'medium', 'high'))
add_node('sexual_function', c('low', 'medium', 'high'))
add_node('coping', c('low', 'medium', 'high'))
add_node('cognitive_function', c('low', 'medium', 'high'))
add_node('social_support', c('low', 'medium', 'high'))
add_node('anxiety', c('low', 'medium', 'high'))

# Fertility
add_node('fertility', c('lost', 'kept', 'unknown'))

## Bayes Network Construction

# Get arcs from computed recurrance DAG
original_arcs <- arcs(dag2)

# Remove breast location arc
original_arcs <- original_arcs[1:5,]

# Set up empty BN
final_dag <- empty.graph(nodes=nodes)
arcs(final_dag) <- original_arcs

# Add arcs that make sense
final_dag <- set.arc(final_dag, from='age', to='breast_reconstruction')
final_dag <- set.arc(final_dag, from='age', to='fertility')
final_dag <- set.arc(final_dag, from='menopause', to='fertility')
final_dag <- set.arc(final_dag, from='masectomy', to='body_image')
final_dag <- set.arc(final_dag, from='social_support', to='anxiety')
final_dag <- set.arc(final_dag, from='tumor_size', to='lumpectomy')
final_dag <- set.arc(final_dag, from='tumor_size', to='masectomy')



# Use modelstr to add arcs:
"[age][tumor_size][inv_nodes][family_history_bc][overall_health|age:bmi][bmi|age][tumor_histology][er_status][pr_status][her2_status][21_gene_assay_rec_risk|er_status:pr_status:her2_status][primary_tumor_class|tumor_size:tumor_histology][regional_lymph_class|inv_nodes:node_caps:deg_malig][pathologic_lymph_class|inv_nodes:node_caps:deg_malig][metastasis_class|histological_grade:tumor_histology:pathologic_lymph_class:primary_tumor_class:deg_malig][cancer_stage][histological_grade|tumor_histology][radiation_therapy][endocrine_therapy][chemotherapy][physical_function][sexual_function][coping][cognitive_function][social_support][menopause|age][node_caps|inv_nodes][masectomy|tumor_size][lumpectomy|tumor_size][breast_reconstruction|age][anxiety|social_support][deg_malig|node_caps][irradiat|node_caps][body_image|masectomy][fertility|age:menopause][recurrance|deg_malig:cancer_stage:21_gene_assay_rec_risk:endocrine_therapy:irradiat:chemotherapy:er_status:pr_status:her2_status:cancer_stage:histological_grade]"

"[age][tumor_size][inv_nodes][family_history_bc][overall_health|age:bmi][bmi|age][tumor_histology][er_status][pr_status][her2_status][21_gene_assay_rec_risk|er_status:pr_status:her2_status:family_history_bc][primary_tumor_class|tumor_size][regional_lymph_class|inv_nodes:node_caps][pathologic_lymph_class|inv_nodes:node_caps][metastasis_class|pathologic_lymph_class:tumor_histology:regional_lymph_class]
   [cancer_stage][histological_grade][radiation_therapy][endocrine_therapy][chemotherapy][physical_function][sexual_function][coping][cognitive_function][social_support][menopause|age][node_caps|inv_nodes][masectomy|tumor_size]
[lumpectomy|tumor_size][breast_reconstruction|age][anxiety|social_support][deg_malig|node_caps][irradiat|node_caps][body_image|masectomy][fertility|age:menopause][recurrance|deg_malig]"


"
[age][cancer_stage][chemotherapy][cognitive_function][coping|anxiety:social_support:recurrance:body_image][endocrine_therapy][er_status][family_history_bc][her2_status]
   [histological_grade][inv_nodes][physical_function][pr_status][radiation_therapy][sexual_function|menopause:masectomy:breast_reconstruction][social_support|masectomy:breast_reconstruction:recurrance:fertility]
[tumor_histology][tumor_size][21_gene_assay_rec_risk|er_status:family_history_bc:her2_status:pr_status][anxiety|social_support]
[bmi|age][breast_reconstruction|age][lumpectomy|tumor_size][masectomy|tumor_size][menopause|age][node_caps|inv_nodes]
[primary_tumor_class|tumor_size][body_image|masectomy][deg_malig|node_caps][fertility|age:menopause][irradiat|node_caps]
[overall_health|age:bmi][pathologic_lymph_class|inv_nodes:node_caps][regional_lymph_class|inv_nodes:node_caps]
[metastasis_class|pathologic_lymph_class:regional_lymph_class:histological_grade][recurrance|deg_malig]
"


"[age][cancer_stage|histological_grade:metastasis_class:primary_tumor_class:regional_lymph_class:pathologic_lymph_class][chemotherapy][cognitive_function][endocrine_therapy][er_status][family_history_bc][her2_status]
[histological_grade][inv_nodes][physical_function][pr_status][tumor_histology][tumor_size]
[21_gene_assay_rec_risk|er_status:family_history_bc:her2_status:pr_status][bmi|age][breast_reconstruction|age]
[lumpectomy|tumor_size][masectomy|tumor_size][menopause|age][node_caps|inv_nodes][primary_tumor_class|tumor_size]
[body_image|masectomy][deg_malig|node_caps][fertility|age:menopause][irradiat|node_caps][overall_health|age:bmi]
[pathologic_lymph_class|inv_nodes:node_caps][regional_lymph_class|inv_nodes:node_caps]
[sexual_function|breast_reconstruction:masectomy:menopause]
[metastasis_class|histological_grade:pathologic_lymph_class:regional_lymph_class][recurrance|irradiat:deg_malig:cancer_stage:chemotherapy:endocrine_therapy:21_gene_assay_rec_risk:histological_grade]
[social_support|breast_reconstruction:fertility:masectomy:recurrance][anxiety|social_support]
[coping|anxiety:body_image:recurrance:social_support]"

"[age][chemotherapy|overall_health:cancer_stage][endocrine_therapy|overall_health:cancer_stage:er_status:her2_status:pr_status][er_status][family_history_bc][her2_status][histological_grade|tumor_histology]
[inv_nodes][physical_function|recurrance:chemotherapy:lumpectomy:masectomy][pr_status][tumor_histology][tumor_size]
[21_gene_assay_rec_risk|er_status:family_history_bc:her2_status:pr_status][bmi|age][breast_reconstruction|age]
[lumpectomy|tumor_size:age:overall_health][masectomy|tumor_size:age:overall_health][menopause|age][node_caps|inv_nodes][primary_tumor_class|tumor_size]
[body_image|masectomy][deg_malig|node_caps][fertility|age:menopause][irradiat|node_caps:overall_health:cancer_stage][overall_health|age:bmi]
[pathologic_lymph_class|inv_nodes:node_caps][regional_lymph_class|inv_nodes:node_caps]
[sexual_function|breast_reconstruction:masectomy:menopause:physical_function]
[metastasis_class|histological_grade:pathologic_lymph_class:regional_lymph_class]
[cancer_stage|histological_grade:metastasis_class:pathologic_lymph_class:primary_tumor_class:regional_lymph_class]
[recurrance|21_gene_assay_rec_risk:cancer_stage:chemotherapy:deg_malig:endocrine_therapy:histological_grade:irradiat:lumpectomy:masectomy]
[social_support|breast_reconstruction:fertility:masectomy:recurrance][anxiety|social_support]
[coping|anxiety:body_image:recurrance:social_support:sexual_function]"



final_dag <- model2network("[age][er_status][family_history_bc][her2_status][inv_nodes][pr_status][tumor_histology][tumor_size]
[21_gene_assay_rec_risk|er_status:family_history_bc:her2_status:pr_status][bmi|age][breast_reconstruction|age]
[histological_grade|tumor_histology][menopause|age][node_caps|inv_nodes][primary_tumor_class|tumor_size][deg_malig|node_caps]
[fertility|age:menopause:chemotherapy][overall_health|age:bmi][pathologic_lymph_class|inv_nodes:node_caps]
[regional_lymph_class|inv_nodes:node_caps][lumpectomy|age:overall_health:tumor_size][masectomy|age:overall_health:tumor_size]
[metastasis_class|histological_grade:pathologic_lymph_class:regional_lymph_class][body_image|masectomy]
[cancer_stage|histological_grade:metastasis_class:pathologic_lymph_class:primary_tumor_class:regional_lymph_class]
[chemotherapy|cancer_stage:overall_health][endocrine_therapy|cancer_stage:er_status:her2_status:overall_health:pr_status]
[irradiat|cancer_stage:node_caps:overall_health]
[recurrance|21_gene_assay_rec_risk:cancer_stage:chemotherapy:deg_malig:endocrine_therapy:histological_grade:irradiat:lumpectomy:masectomy]
[physical_function|chemotherapy:lumpectomy:masectomy:recurrance]
[social_support|breast_reconstruction:fertility:masectomy:recurrance][anxiety|social_support]
[sexual_function|breast_reconstruction:masectomy:menopause:physical_function]
[coping|anxiety:body_image:recurrance:sexual_function:social_support]")

sample_data[is.na(cancer_stage), cancer_stage:='IIA']
for (j in colnames(sample_data)) set(sample_data, j = j, value = factor(sample_data[[j]]))
remove_cols <- setdiff(colnames(sample_data), nodes(final_dag))
sample_data[,(remove_cols):=NULL]

final.bn.bayes1 <- bn.fit(final_dag, data=sample_data, method='bayes', iss=10)




# Visualize
g <- graphviz.plot(final_dag, layout='neato') ; plot(g, attrs = list(graph = list(rankdir="LR"), node = list(fillcolor = "lightblue", fontsize=20)))

modelstr <- "[tumor_grade][metastasis][molecular_subtype][chemotherapy][age][tumor_size][inv_nodes][breast][menopause|age][node_caps|inv_nodes][breast_quad|breast][deg_malig|node_caps][irradiat|node_caps]
[recurrance|deg_malig]"


