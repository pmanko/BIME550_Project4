# source("http://bioconductor.org/biocLite.R")
# biocLite(c("graph", "Rgraphviz", "RBGL"))
# install.packages("gRain", 'bnlearn')

library(bnlearn)
library(graph)
library(Rgraphviz)
library(RBGL)
library(gRain)


# Set up: Age, Sex, Education, Occupation, Residence, Travel
dag <- empty.graph(nodes = c("A", "S", "E", "O", "R", "T"))

# Set up arcs
dag <- set.arc(dag, from = "A", to = "E")
dag <- set.arc(dag, from="S", to = "E")

dag <- set.arc(dag, from ='E', to="O")
dag <- set.arc(dag, from="E", to="R")

dag <- set.arc(dag, from='O', to='T')
dag <- set.arc(dag, from="R", to="T")

modelstring(dag)
graphviz.plot(dag, shape='ellipse')

# Use matrix to set arcs
dag2 <- empty.graph(nodes = c("A", "S", "E", "O", "R", "T"))
arc_matrix <- arcs(dag)
arcs(dag2) <- arc_matrix

# Ensures DAG:
try(set.arc(dag, from='T', to="E"))

# Set Variable Levels
A.lv <- c("young", "adult", "old")
S.lv <- c("M", "F")
E.lv <- c("high", "uni")
O.lv <- c("emp", "self")
R.lv <- c("small", "big")
T.lv <- c("car", "train", "other")

# Variables not linked by arc are conditionally independent:
# Pr(A, S, E, O, R, T) = Pr(A) Pr(S) Pr(E | A, S) Pr(O | E) Pr(R | E) Pr(T | O, R).

# Given this, set probabilites:
A.prob <- array(c(0.30, 0.50, 0.20), dim = 3, dimnames = list(A = A.lv))
S.prob <- array(c(0.60, 0.40), dim = 2, dimnames = list(S = S.lv))

O.prob <- array(c(0.96, 0.04, 0.92, 0.08), dim = c(2, 2), dimnames = list(O = O.lv, E = E.lv))
R.prob <- array(c(0.25, 0.75, 0.20, 0.80), dim = c(2, 2), dimnames = list(R = R.lv, E = E.lv))

E.prob <- array(c(0.75, 0.25, 0.72, 0.28, 0.88, 0.12, 0.64, 0.36, 0.70, 0.30, 0.90, 0.10), dim = c(2, 3, 2), dimnames = list(E = E.lv, A = A.lv, S = S.lv))
T.prob <- array(c(0.48, 0.42, 0.10, 0.56, 0.36, 0.08, 0.58, 0.24, 0.18, 0.70, 0.21, 0.09), dim = c(3, 2, 2), dimnames = list(T = T.lv, O = O.lv, R = R.lv))

# Recreate DAG using modelstring
dag3 <- model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")
all.equal(dag, dag3)

# Local Probability Distributions
cpt <- list(A=A.prob, S = S.prob, E = E.prob, O = O.prob, R = R.prob, T = T.prob)

# Combine dag and prob dists into a bayes network:
bn <- custom.fit(dag, cpt)

# Some useful fns
arcs(bn)
nparams(bn)
bn$R
bn
plot(bn)

# The code above assumes we know the graph structure AND conditional probabilites


