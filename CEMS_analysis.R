# This file contains code to analyzing CEMS data 
# using the framework proposed by Wu et al (2022).

# Set working directory to files pane location

# Import utility files and packages
source("bt_utils.R")
library(QTLRel)

# Load data
data(CEMS)

# Delete NAs
votes = CEMS$preferences[!is.na(CEMS$preferences$win1),]
votes = votes[!is.na(votes$win2),]
# Delete ties
# votes = votes[votes$tied==0, ]
# Change attribute names
votes$object1 = votes$school1
votes$object2 = votes$school2
votes$Subject.ID = votes$student

# Make contingency matrix
votes.to.W <- function(votes, obj.names = NULL)
{
  n = length(obj.names)
  W = matrix(data = 0,nrow = n, ncol = n)
  rownames(W) = obj.names
  colnames(W) = obj.names
  for (i in 1:nrow(votes))
  {
    obj1 = votes$object1[i]
    obj2 = votes$object2[i]
    W[obj1,obj2] = W[obj1,obj2] + votes$win1.adj[i] 
    W[obj2,obj1] = W[obj2,obj1] + votes$win2.adj[i] 
  }
  return(W)
}

school.names = c('London', 'Paris', 'Milano', 'St.Gallen', 
                 'Barcelona','Stockholm')

W = votes.to.W(votes,obj.names = school.names)

# Fit the Bradley-Terry model
bt.results = fit.bt.model(W,verbose = TRUE)

# Visualize the results
bt.visual = bt.visualize(W, bt.results, top = 6, point = F,
                         obj.names = school.names)

# Ranking graph, Figure 12 (left) on page 25 of the paper
bt.visual$rank.graph + theme(axis.title = element_text(size = 15),
      axis.text = element_text(size = 15))

# z-score heatmap, Figure 12 (right) on page 25 of the paper
bt.visual$zscore.heatmap+ theme(plot.margin = unit(c(1,0,1,0),'cm'))

# Object diagnostics plots
# Figure 13 on page 26 of the paper
diag.plots = bt.obj.diagnostics(W, bt.results)

# Subject diagnostics
new.sbj.diagnostics <- function(W, bt.results, votes)
{
  bt.scores = bt.results$bt.scores
  obj.list = names(bt.scores)
  n = length(obj.list)
  B = bt.results$P
  R = bt.results$R
  Var.beta = bt.results$Var.beta
  
  total.votes = sum(W)
  sbj.list = unique(votes$Subject.ID)
  sbj.diag = matrix(data = 0, nrow = length(sbj.list),ncol = 5)
  rownames(sbj.diag) = sbj.list
  colnames(sbj.diag) = c('num.votes','difference','influence','dif.q','inf.q')
  sbj.diag = as.data.frame(sbj.diag)
  
  for(sbj in sbj.list)
  {
    delta.votes = votes[votes$Subject.ID == sbj,]
    delta.W = votes.to.W(delta.votes,obj.names = obj.list)
    delta.V = delta.W + t(delta.W)
    delta.w = rowSums(delta.W)
    E.delta.w = rowSums(delta.V * B)
    diff = delta.w - E.delta.w
    infl = Var.beta %*% diff
    sbj.diag[sbj,'num.votes'] = nrow(delta.votes)
    sbj.diag[sbj,'difference'] = sum(diff**2)
    sbj.diag[sbj,'influence'] = sum(infl**2)
  }
  
  # Reference for difference
  
  mus = eigen(Var.beta,symmetric = TRUE, 
              only.values = TRUE)$values[1:(n-1)]
  lambdas = 1/mus
  
  sbj.diag$avg.diff = sbj.diag$difference / sbj.diag$num.votes
  sbj.diag$avg.infl = sbj.diag$influence / sbj.diag$num.votes
  num.voters = nrow(sbj.diag)
  
  # Q-Q plot for the generalized chi-square distribution
  sbj.diag = sbj.diag[order(sbj.diag$avg.diff),]
  # Theoretical Quantiles
  sbj.diag$dif.q = 0
  for(i in 1:num.voters)
  {
    sbj.diag$dif.q[i] = qgenchisq(p = (i-0.5)/num.voters, 
                                  lambda = lambdas, max_iter = 100)
  }

  # Difference plot
  #sbj.diag = sbj.diag[order(sbj.diag$avg.diff),]
  sbj.diag$dif.q = sbj.diag$dif.q / total.votes
  
  # Influence 
  sbj.diag = sbj.diag[order(sbj.diag$avg.infl),]
  #sbj.diag$inf.q = 0
  for(i in 1:num.voters)
    sbj.diag$inf.q[i] = qgenchisq(p = (i-0.5)/num.voters, 
                                  lambda = mus, max_iter = 100)
  sbj.diag$inf.q = sbj.diag$inf.q / total.votes
  
  return(sbj.diag)
}

sbj.diag = new.sbj.diagnostics(W, bt.results, votes)

# Difference plot, Figure 14 (left) on page 27 of the paper
qqPlot(sbj.diag$avg.diff,sbj.diag$dif.q,
       xlab = 'Theoretical Quantiles',
       ylab = 'Empirical Quantiles',
       main = '',
       qqline = 'expected')

# Influence plot, Figure 14 (right) on page 27 of the paper
qqPlot(sbj.diag$avg.infl,sbj.diag$inf.q,
       xlab = 'Theoretical Quantiles',
       ylab = 'Empirical Quantiles',
       main = '',
       qqline = 'expected')


# Parametric bootstrap experiment, on page 26-27 of the paper
new.votes = votes
set.seed(0)

for(i in 1:nrow(votes))
{
  obj1 = new.votes$object1[i]
  obj2 = new.votes$object2[i]
  p = bt.results$P[obj1,obj2]
  new.votes$win1.adj[i] = rbinom(1,1,p)
  new.votes$win2.adj[i] = 1 - new.votes$win1.adj[i]
}
new.W = votes.to.W(new.votes,obj.names = school.names)

# Fit the Bradley-Terry model
new.bt.results = fit.bt.model(new.W,verbose = TRUE)

# Subject Diagnostics
new.sbj.diag = new.sbj.diagnostics(new.W, new.bt.results, new.votes)

# Difference plot, Figure 15 (left) on page 28 of the paper
qqPlot(new.sbj.diag$avg.diff,new.sbj.diag$dif.q,
       xlab = 'Theoretical Quantiles',
       ylab = 'Empirical Quantiles',
       main = '',
       qqline = 'expected')

# Influence plot, Figure 15 (right) on page 28 of the paper
qqPlot(new.sbj.diag$avg.infl,new.sbj.diag$inf.q,
       xlab = 'Theoretical Quantiles',
       ylab = 'Empirical Quantiles',
       main = '',
       qqline = 'expected')



