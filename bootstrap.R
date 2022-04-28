# This file contains code to perform the bootstrap experiments
# for the PlaNYC wiki-survey in by Wu et al (2022).

# Set working directory to files pane location

# Import utility files
source("wiki_utils.R")
source("bt_utils.R")

votes.filename <- paste("data/wikisurvey_608_votes","_cleaned.csv", sep="")
ideas.filename <- paste("data/wikisurvey_608_ideas", "_cleaned.csv", sep="")

## Load data
votes <- read.csv(votes.filename, header=TRUE, sep = ",", dec=".")
ideas <- read.csv(ideas.filename, header=TRUE, sep=",", dec=".")

## Preprocessing
ideas = preprocess.ideas(ideas,thres = 0)
votes = preprocess.votes(votes, ideas)

#Bradley-Terry model fit
W = make.cont.matrix(votes,ideas)
bt.results = fit.bt.model(W)
print('Model fitted.')

bt.scores = bt.results$bt.scores
obj.list = names(bt.scores)
n = length(obj.list)
B = bt.results$P
R = bt.results$R
Var.beta = bt.results$Var.beta

# Reference for difference
mus = eigen(Var.beta,symmetric = TRUE, 
            only.values = TRUE)$values[1:(n-1)]
lambdas = 1/mus

upper.dif = qgenchisq(p = 0.975, lambda = lambdas, max_iter = 100)
lower.dif = qgenchisq(p = 0.025, lambda = lambdas, max_iter = 100)
middle.dif = qgenchisq(p = 0.5, lambda = lambdas, max_iter = 100)

# Reference for influence
upper.inf = qgenchisq(p = 0.975, lambda = mus, max_iter = 100)
lower.inf = qgenchisq(p = 0.025, lambda = mus, max_iter = 100)
middle.inf = qgenchisq(p = 0.5, lambda = mus, max_iter = 100)

rownames(B) = obj.list
colnames(B) = obj.list
total.votes = nrow(votes)
voters.list= unique(votes$Subject.ID)

# Experiment 1: Every voter is faced with the same choices
sim.voter.diag = matrix(nrow = length(voters.list),
                        ncol = 3)
rownames(sim.voter.diag) = voters.list
colnames(sim.voter.diag) = c('num.votes','difference','influence')
sim.voter.diag = as.data.frame(sim.voter.diag)
for(voter in voters.list)
{
  delta.votes = votes[votes$Session.ID == voter,]
  delta.W = matrix(data = 0, nrow = n, ncol = n)
  rownames(delta.W) = obj.list
  colnames(delta.W) = obj.list
  for(i in 1:nrow(delta.votes))
  {
      left = delta.votes$Left.Choice.ID[i]
      right = delta.votes$Right.Choice.ID[i]
      x = runif(1)
      if(x < B[left,right]) # left wins
        delta.W[left,right] = delta.W[left,right] + 1
      else
        delta.W[right,left] = delta.W[right,left] + 1
  }
  delta.V = delta.W + t(delta.W)
  delta.w = rowSums(delta.W)
  E.delta.w = rowSums(delta.V * B)
  diff = delta.w - E.delta.w
  infl = Var.beta %*% diff
  sim.voter.diag[voter,'num.votes'] = nrow(delta.votes)
  sim.voter.diag[voter,'difference'] = sum(diff**2)
  sim.voter.diag[voter,'influence'] = sum(infl**2)
}

# Left plot of Figure 9 on page 22 of the paper
diff_bootstrap1 = 
  ggplot(data = sim.voter.diag, aes(x = num.votes, y = difference))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.dif / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.dif / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.dif / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Difference $D^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
diff_bootstrap1

# Right plot of Figure 9 on page 22 of the paper
inf_bootstrap1 = 
  ggplot(data = sim.voter.diag, aes(x = num.votes, y = influence))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.inf / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.inf / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.inf / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Influence $I^{(s)}$'))+ 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
inf_bootstrap1

# Experiment 2: Every voter faces new choices
sim.voter.diag2 = matrix(nrow = length(voters.list), ncol = 3)
colnames(sim.voter.diag2) = c('num.votes','difference','influence')
rownames(sim.voter.diag2) = voters.list
sim.voter.diag2 = as.data.frame(sim.voter.diag2)

for(voter in voters.list)
{
  delta.votes = votes[votes$Session.ID == voter,]
  num.votes = nrow(delta.votes)
  sim.voter.diag2[voter,'num.votes'] = num.votes
  indices = sample(x = total.votes, size = num.votes)
  delta.votes = votes[indices,]
  delta.W = matrix(data = 0, nrow = n, ncol = n)
  rownames(delta.W) = obj.list
  colnames(delta.W) = obj.list
  for(i in 1:nrow(delta.votes))
  {
    left = delta.votes$Left.Choice.ID[i]
    right = delta.votes$Right.Choice.ID[i]
    x = runif(1)
    if(x < B[left,right]) # left wins
      delta.W[left,right] = delta.W[left,right] + 1
    else
      delta.W[right,left] = delta.W[right,left] + 1
  }
  delta.V = delta.W + t(delta.W)
  delta.w = rowSums(delta.W)
  E.delta.w = rowSums(delta.V * B)
  diff = delta.w - E.delta.w
  infl = Var.beta %*% diff
  
  sim.voter.diag2[voter,'difference'] = sum(diff**2)
  sim.voter.diag2[voter,'influence'] = sum(infl**2)
}

# Left plot of Figure 10 on page 23 of the paper
diff_bootstrap2 = 
  ggplot(data = sim.voter.diag2, aes(x = num.votes, y = difference))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.dif / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.dif / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.dif / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Difference $D^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
diff_bootstrap2

# Right plot of Figure 10 on page 23 of the paper
inf_bootstrap2 = 
  ggplot(data = sim.voter.diag2, aes(x = num.votes, y = influence))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.inf / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.inf / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.inf / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Influence $I^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
inf_bootstrap2
