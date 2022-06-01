# This file contains code to produce Figs. 14 and 15
# in the supplementary material of Wu et al (2021).

# Set working directory to files pane location
# setwd("~/Desktop/ADA/code")

# Import utility files and packages
source("wiki_utils.R")
source("bt_utils.R")
suppressPackageStartupMessages(library(ggplot2))

prefix <- 'planyc'

if (prefix == "planyc") {
  question.text <- "Which do you think is a better idea for creating a greener, greater New York City?"
  cex.main.val <- 0.8
  question.text.2line <- "Which do you think is a better idea for creating\na greener, greater New York City?"
  wiki.survey.num <- 608
  votes.filename <- paste("data/wikisurvey_608_votes","_cleaned.csv", sep="")
  ideas.filename <- paste("data/wikisurvey_608_ideas", "_cleaned.csv", sep="")
  nonvotes.filename <- paste("data_precleaned/wikisurvey_608_nonvotes", "_cleaned.csv", sep="")
  start.date <- as.Date("2010-10-07")
  stop.date <- as.Date("2011-01-30")
}

## Load data
votes <- read.csv(votes.filename, header=TRUE, sep = ",", dec=".")
ideas <- read.csv(ideas.filename, header=TRUE, sep=",", dec=".")

## Preprocessing
ideas = preprocess.ideas(ideas,thres = 0)
votes = preprocess.votes(votes, ideas)

#Bradley-Terry model fit
W = make.cont.matrix(votes,ideas)
bt.scores = fit.bt.model(W)

n = nrow(W) # Number of objects
V = W + t(W)
A = matrix(data = 0, nrow = n, ncol = n)
for(i in 1:n) A[,i] = bt.scores
for(i in 1:n) A[i,] = A[i,] - bt.scores
B = sigmoid(A)

H.beta = V * B * (1-B)
tmp = rowSums(H.beta)
diag(H.beta) = -tmp
#H.beta is the Hessian matrix of \beta
H.gamma = H.beta[2:n,2:n]
#H.gamma is the Hessian matrix of \gamma
Var.gamma = solve(-H.gamma)
#Var.gamma is the variance of \gamma
P = matrix(data = 0, nrow = n, ncol = n-1) #\beta = P\gamma
for(j in 1:(n-1)) P[j+1,j] = 1
P = P - 1 / n
Var.beta = P %*% Var.gamma %*% t(P)

df = matrix(data = 0, nrow = n, ncol = 3)
colnames(df) = c('name','num.app','score')
df = as.data.frame(df)
df$name = names(bt.scores)
df$score = bt.scores
df$num.app = rowSums(V)

df$variance = diag(Var.beta)
df$appxvar = df$variance * df$num.app

# Relationship between the estimated variance and the
# number of appearances 
# (Fig.1 on page 3 of the supplementary material)
ggplot(df,
       aes(x = Num.Appearances, y = variance)) + 
  geom_point(aes(color = abs(bt.scores))) +
  xlab(TeX('N_i')) + 
  ylab(TeX('$\\widehat{Var}(\\hat{\\beta}_i)$')) +
  scale_colour_continuous(TeX('$|\\hat{\\beta}_i|$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  theme(legend.position = c(0.9,0.9))

# Relationship bewteen the estimated Bradley-Terry scores and
# the product of number of appearances and the estimated variance
# (Fig.2 on page 3 of the supplementary material)
ggplot(df,
       aes(x = bt.scores, y = appxvar)) +
  geom_point() +
  xlab(TeX('$\\hat{\\beta}_i$')) + 
  ylab(TeX('$N_i\\widehat{Var}(\\hat{\\beta}_i)$'))+ 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_y_continuous(breaks = c(4,8,12,16))

