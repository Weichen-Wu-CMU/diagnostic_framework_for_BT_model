# This file contains code for experiments 1 and 2 in the paper
# "A Diagnostic Framework for the Bradley-Terry model", by Wu et al (2022)

## Load packages
source("bt_utils.R")
suppressPackageStartupMessages(library(BradleyTerry2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(latex2exp))

# Real scores
real.scores = seq(from = -10, to = 10, by = 1)
real.scores = real.scores / 10
tmp = matrix(data = 0, nrow = 21, ncol = 21)
for(i in 1:21) tmp[,i] = real.scores
for(i in 1:21) tmp[i,] = tmp[i,] - real.scores

# To replicate the results in the paper, set the random seed to 0
set.seed(0)

# Example 1

P.star = sigmoid(tmp) 
#P.star[i,j] = P(the current subject chooses i over j)
#In example 1, this needs to be updated after each subject
P = P.star
W = matrix(data = 0, nrow = 21, ncol = 21) # the contingency matrix
rownames(W) = seq(1,21)
colnames(W) = seq(1,21)

# The contingency matrix for the current subject
cur = matrix(data = 0, nrow = 21, ncol = 21) 
for(t in 1:100) # iterate through 100 subjects
{
  for(i in 1:20)
  {
    for(j in (i+1):21)
    {
      cur[i,j] = rbinom(1,1,P.star[i,j])
      cur[j,i] = 1 - cur[i,j]
    }
  }
  W = W + cur
  P.star = (P + W/t) / 2 # Update P.star
}

# Fit the Bradley-Terry model
bt.scores = fit.bt.model(W)

# Relationship between real and estimated ranks
est.rank = order(bt.scores$bt.scores)
rankings = matrix(data = NA, nrow = 21, ncol = 2)
colnames(rankings) = c('real','est')
rankings = as.data.frame(rankings)
rankings$real = seq(1,21)
rankings$est = est.rank
# Figure 2 on page 10 of the paper
ggplot(rankings,aes(x = real,y=est)) + 
  geom_point()+
  xlab(TeX('Rank in  $\\beta^*_j$')) +
  ylab(TeX('Rank in $\\hat{\\beta}_j$')) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'none')
  

# Object diagnostic plots
obj.diag.plots = bt.obj.diagnostics(W, bt.scores)

# QQ plot of residuals, Figure 1 on page 10 in the paper
obj.diag.plots$r.qqplot +
  theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))

# Example 2
P.star = sigmoid(tmp**3)
W = matrix(data = 0, nrow = 21, ncol = 21)
rownames(W) = seq(1,21)
colnames(W) = seq(1,21)

# Simulate data
for(i in 1:20)
{
  for (j in (i+1):21)
  {
    W[i,j] = rbinom(1, 100, P.star[i,j])
    W[j,i] = 100 - W[i,j]
  }
}

# Fit the Bradley-Terry model
bt.results = fit.bt.model(W)

# Uncomment the following line to check the relationship bewteen the 
# fitted score and the real scores
plot(bt.results$bt.scores,real.scores)

# Object diagnostic plots
obj.diag.plots = bt.obj.diagnostics(W, bt.results)


# Relationship bewteen residuals and Bradley-Terry scores
# (Figure 3 on page 10 of the paper)
obj.diag.plots$r.vs.scores +
  theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 15))

# Closer look at objects 1 and 11
probs = matrix(data = 0, nrow = 21, ncol = 5)
colnames(probs) = c('Index','real1','est1','real11','est11')
probs = as.data.frame(probs)
probs$Index = seq(1,21)

# probs$obs1[i] = Observed proportion that 1 beats i
probs$obs1 = W[1,]/100
probs$obs1[1] = 0.5

# probs$obs11[i] = Observed proportion that 11 beats i
probs$obs11 = W[11,]/100
probs$obs11[11] = 0.5

probs$scores = bt.results$bt.scores
bt.scores = bt.results$bt.scores

for(i in 1:21)
{
  # probs$real1[i] = Real probability that 1 beats i
  probs$real1[i] = sigmoid(-(real.scores[i]+1)**3)
  
  # probs$est1[i] = Estimated probability that 1 beats i
  probs$est1[i] = sigmoid(-bt.scores[i] + bt.scores[1])
  
  # probs$real11[i] = Real probability that 11 beats i
  probs$real11[i] = sigmoid(-real.scores[i]**3)
  
  # probs$est11[i] = Estimated probability that 11 beats i
  probs$est11[i] = sigmoid(-bt.scores[i])
}

# Exploration of object 1 (Left plot of Figure 4 in the paper)
ggplot(probs, aes(x = bt.scores)) +
  geom_point(aes(y = real1, 
                 shape = 'Real')) +
  geom_line(aes(y = real1)) +
  geom_point(aes(y = est1, 
                 shape = 'Estimated')) +
  geom_line(aes(y = est1))+
  xlab(TeX('Estimated Bradley-Terry scores $\\hat{\\beta}_j$')) +
  ylab('P(1 beats j)') +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'none')

# Exploration of object 11 (Right plot of Figure 4 in the paper)
ggplot(probs, aes(x = bt.scores)) +
  geom_point(aes(y = real11, shape = 'Real')) +
  geom_line(aes(y = real11)) +
  geom_point(aes(y = est11, shape = 'Estimated')) +
  geom_line(aes(y = est11)) +
  xlab(TeX('Estimated Bradley-Terry scores $\\hat{\\beta}_j$')) +
  ylab('P(11 beats j)') + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = 'none')

