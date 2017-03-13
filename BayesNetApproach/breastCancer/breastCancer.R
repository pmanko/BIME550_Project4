
library(bnlearn)
library(graph)
library(Rgraphviz)
library(RBGL)
library(gRain)


# Bayes Network Exploration for Project 4
# Breast Cancer Datsaets
mplot <- function(x) {
  graphviz.plot(x, shape='ellipse')
}

## Learning Bayes Network Structure from Data using Hill-Climbing Algorithm and TABU search

# Wisconsis BN
cols <- colnames(wisconsin)[2:3]
wisconsin_d <- as.data.table(discretize(wisconsin[,2:10], method='interval', breaks=rep(3,9), ordered=FALSE))
wisconsin_d[,class:=wisconsin$classes]

learned_dag <- hc(wisconsin_d, restart = 100, perturb = 100)
mplot(learned_dag)

learned_dag_rsmax2 <- rsmax2(wisconsin_d, restrict = "si.hiton.pc", test = "x2", maximize = "tabu", score = "bde", maximize.args = list(iss = 5), alpha = 0.25)
mplot(learned_dag_rsmax2)


learned_dag_tabu <- tabu(wisconsin_d, debug=TRUE)
graphviz.plot(learned_dag_tabu, layout='circo', shape='ellipse')



# wpbc
wpbc
ncol(wpbc)
wpbc_d <- as.data.table(discretize(wpbc[,2:34], method = 'quantile', breaks=c(rep(3, 32), 2)))

learned_wpbc_dag_tabu <- tabu(wpbc, debug=TRUE)
mplot(learned_wpbc_dag_tabu)

# recurrance dataset
recurrance
learned_recurrance_dag_rsmax2 <- rsmax2(recurrance, restrict = "si.hiton.pc", test = "x2", maximize = "tabu", score = "bde", maximize.args = list(iss = 5))
mplot(learned_recurrance_dag_rsmax2)

learned_recurrance_dag <- hc(recurrance, perturb = 100, restart=100)
mplot(learned_recurrance_dag)

learned_recurrance_dag_tabu <- tabu(recurrance, debug= TRUE)
mplot(learned_recurrance_dag_tabu)

# risk
ncol(risk)

# Remove meta-data (unrelated)
risk_subset <- risk[,1:14]

# We're working w/ discrete dags, so change ints to factors
for (j in colnames(risk_subset)) set(risk_subset, j = j, value = as.factor(risk_subset[[j]]))
learned_dag_risk <- rsmax2(risk_subset, restrict = "si.hiton.pc", test = "x2", maximize = "tabu", score = "bde", maximize.args = list(iss = 5), alpha = 0.0001)
mplot(learned_dag_risk)

# risk_factor
ncol(risk_factor)

risk_factor_subset <- risk_factor[,2:12]
for (j in colnames(risk_factor_subset)) set(risk_factor_subset, j = j, value = as.factor(risk_factor_subset[[j]]))

learned_dag_risk_factor <- rsmax2(risk_factor_subset, restrict = "si.hiton.pc", test = "x2", maximize = "tabu", score = "bde", maximize.args = list(iss = 5), alpha=0.0001)
graphviz.plot(learned_dag_risk_factor, shape='ellipse')


## Parameter Estimation: Using Wisconsin/Recurrance BNs to determine conditional probabilities

# Network Graphs:
learned_dag_tabu <- tabu(wisconsin_d, debug=TRUE)
learned_recurrance_dag <- hc(recurrance, perturb = 100, restart=100)

dag1 <- learned_dag_tabu
dag2 <- learned_recurrance_dag

graphviz.plot(dag1, layout='circo', shape='ellipse')
graphviz.plot(dag2, layout='dot', shape='ellipse')

# Clean up the dataset to match what was used in determening dag structure
# I wanted to go back to non-discretized data for the rest of analysis. 
evidence1 <- wisconsin
for (j in colnames(evidence1)) set(evidence1, j = j, value = as.factor(evidence1[[j]]))
setnames(evidence1, 'classes', 'class')

evidence2 <- recurrance

# Graph Properties:
print(dag1)
print(dag2)

# Arc Strengths (using conditional independence tests)
print(arc.strength(dag1, data = evidence1, criterion='x2')) # Pearson's X^2 test for contigency tables (p-value)
print(arc.strength(dag2, data = evidence2, criterion='x2')) # Pearson's X^2 test for contigency tables (p-value)

# Adding an edge: how it affects evidence:
testdag <- dag2
set.arc(testdag, from='tumor_size', to='recurrance', debug=TRUE)

# (For some reason, set.arc was not working so had to use modelstring workaround -_-)
test_string <- modelstring(dag2)
test_string <- "[age][tumor_size][inv_nodes][breast][menopause|age][node_caps|inv_nodes][breast_quad|breast][deg_malig|node_caps][irradiat|node_caps][recurrance|deg_malig:tumor_size]"
testdag <- model2network(test_string)

arcs(testdag)

print(arc.strength(testdag, data = evidence2, criterion='x2')) # Pearson's X^2 test for contigency tables (p-value)
# The p-val for the new edge is .45

# Using Maximum Likelihood/Bayesian posterior distributions to estimate conditional probabilities:
bn.mle1 <- bn.fit(dag1, data=evidence1, method='mle')
bn.bayes1 <- bn.fit(dag1, data=evidence1, method='bayes', iss=10)

bn.mle2 <- bn.fit(dag2, data=evidence2, method='mle')
bn.bayes2 <- bn.fit(dag2, data=evidence2, method='bayes', iss=10)

## Inference!

# Conditional Probability Queries
#  We can check if two vars are related to each other with the concept of D-separation (test for conditional independence)

dsep(dag1, x = 'uniformity_of_cell_shape', y='clump_thickness')
# Returns false - no Dsep

dsep(dag1, x = 'uniformity_of_cell_shape', y='clump_thickness', z='class')
# Returns true - there is dsep since we condition on class

dsep(dag2, y = 'recurrance', x='irradiat')


# Exact Inference
#  gRain package - requires building a junction tree

# 1. build junction tree
junction1 <- compile(as.grain(bn.mle1))
junction2 <- compile(as.grain(bn.bayes2))

plot(junction1)
plot(junction2)

# 2. Asking questions using querygrain

# Q: If we know someone had recurrance, what is the probability they had irradiation:
querygrain(junction2, nodes='irradiat') # general prob dist around 75/25%

jreccurance <- setEvidence(junction2, nodes='recurrance', states='recurrence-events') # Set evidence: recurrance happened!
querygrain(jreccurance, nodes='irradiat') # a bit less likely to have had irradiations

# Q: If we know someone had node caps, what are the likely numbers of invasive nodes, likely malignancy, and recurrance?
a1 <- querygrain(junction2, nodes=c('inv_nodes', 'irradiat', 'deg_malig', 'recurrance'))

jnodecaps <- setEvidence(junction2, nodes='node_caps', states='yes')
a2 <- querygrain(jnodecaps, nodes=c('inv_nodes', 'irradiat', 'deg_malig', 'recurrance')) # all of them change...

# For example, here's inv_nodes:
plot(a1$inv_nodes)
plot(a2$inv_nodes)

bn.fit.barchart(bn.mle1$class, main = "Class", xlab = "", ylab = "")
bn.fit.dotplot(bn.mle1$class, main = "Class", xlab = "", ylab = "")


# Approximate Infrerence
# Q: if the cell is very non-uniform in shape, how likely is the tumor to be malignant?

cpquery(bn.mle1, event= (class=='malignant'), evidence=(uniformity_of_cell_shape=='3')) # .4
cpquery(bn.mle1, event= (class=='malignant'), evidence=(uniformity_of_cell_shape=='6')) # .9
cpquery(bn.mle1, event= (class=='malignant'), evidence=(uniformity_of_cell_shape=='1'), n=10^6) # .006

# Conditional Independence Queries

# Most Likely Explanation Queries
