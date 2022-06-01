# This file contains code to analyze the PlaNYC wiki-survey
# using the framework proposed by Wu et al (2022).

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

# Distribution of the number of times each idea
# appeared in comparisons
# (Fig.6 on page 18 of the paper)
show.num.app(ideas)

#Bradley-Terry model fit
W = make.cont.matrix(votes,ideas)
bt.results = fit.bt.model(W)
print('Model fitted.')

# Summary of the top-10 ideas
ideas.Summary = c('Ban fracking',
                  'Invest in transportation',
                  'Plug ships into electricity grids',
                  'Enhance bike lane network',
                  'Year-long greenmarkets',
                  'Local food in public schools',
                  'Protected bike paths',
                  'Transit service outside Manhattan',
                  'Support community gardens',
                  'Upgrade building energy efficiency')

# Visualize the results of the Bradley-Terry model
bt.visual = bt.visualize(W,bt.results, 
                         obj.names = ideas.Summary,
                         rank.graph = T,
                         prob.heatmap = T,
                         zscore.heatmap = T)
# Ranking graph of the top-10 ideas
# (Fig.7 on page 19 of the paper)
bt.visual$rank.graph

# Z-score heatmap of the top-10 ideas
# (Fig.8 on page 20 of the paper)
bt.visual$zscore.heatmap

# Object(idea) diagnostics
# (Fig.9 on page 21 of the paper)
obj.diag.plots = bt.obj.diagnostics(W, bt.results)

# Subject (voter) diagnostics
# (Fig.10 on page 22 of the paper)
sbj.diag.plots = bt.sbj.diagnostics(votes, bt.results)

# Distribution of the number of votes cast by each voter
# (Fig.5 on page 18 of the paper)
sbj.diag.plots$rank.plot

