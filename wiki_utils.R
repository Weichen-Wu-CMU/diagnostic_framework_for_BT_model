# This file contains functions to process
# wiki-survey data

suppressPackageStartupMessages(library(ggplot2))

get.active.ideas <- function(ideas)
# This function gets the active ideas in the idea data frame 
# (Mostly for internal use)
{
  ideas = ideas[ideas$Active,]
  return(ideas)
}

rename.ideas <- function(ideas)
# This function adds "o." to the IDs of the ideas
# (Mostly for internal use)
{
  ideas$Idea.ID = as.character(ideas$Idea.ID)
  ideas$Idea.ID = paste('o',ideas$Idea.ID,sep='.')
  return(ideas)
}

change.type.ideas <- function(ideas)
# This function changes the data type of and text of
# ideas into character
# (Mostly for internal use)
{
  ideas$Idea.Text = as.character(ideas$Idea.Text)
  return(ideas)
}

threshold.ideas <- function(ideas, thres = 0)
# This function gets the ideas that
# win and lose for at least thres times
# Inputs:
# ideas, the data frame representing the ideas in the wiki-survey;
# thres, the threshold
# Output:
# ideas, the data frame representing selected ideas
# (Mostly for internal use)
{
  flag = vector(length = nrow(ideas))
  for(i in 1:length(flag))
  {
    flag[i] = min(ideas$Wins[i],ideas$Losses[i])
  }
  ideas = ideas[flag>thres,]
  return(ideas)
}

preprocess.ideas <- function(ideas, thres = 0)
# This function completes the preprocessing of ideas
# by getting the active ideas, making necessary changes to data types
# of certain attributes, and selecting the ideas that win and lose
# for at least thres times.
# Inputs:
# ideas, the data frame representing the ideas in the wiki-survey;
# thres, the threshold
# Output:
# ideas, the data frame representing preprocessed ideas
# (For external use)
{
  ideas = get.active.ideas(ideas)
  ideas = rename.ideas(ideas)
  ideas = change.type.ideas(ideas)
  ideas = threshold.ideas(ideas, thres)
  return(ideas)
}

rename.votes <- function(votes)
# This function adds "s." at the beginning of the Session IDs
# and "o." at the beginning of the idea IDS of each vote
# Input:
# votes, the data frame representing the votes
# Output:
# votes, the data frame representing the votes with session IDs and
#        idea IDs changed.
{
  votes$Session.ID = as.character(votes$Session.ID)
  votes$Session.ID = paste('s',votes$Session.ID,sep = '.')
  votes$Left.Choice.ID = as.character(votes$Left.Choice.ID)
  votes$Left.Choice.ID = paste('o',votes$Left.Choice.ID,sep = '.')
  votes$Right.Choice.ID = as.character(votes$Right.Choice.ID)
  votes$Right.Choice.ID = paste('o',votes$Right.Choice.ID,sep = '.')
  votes$Winner.ID = as.character(votes$Winner.ID)
  votes$Winner.ID = paste('o',votes$Winner.ID,sep = '.')
  votes$Loser.ID = as.character(votes$Loser.ID)
  votes$Loser.ID = paste('o',votes$Loser.ID,sep = '.')
  votes$Subject.ID = votes$Session.ID
  return(votes)
}

get.valid.votes <- function(votes, act.ideas.list)
# This function gets the valid votes 
# Valid votes must be between two active ideas
# Inputs:
# votes, the data frame representing the votes
#        notice that the names of this data frame must have been
#        pre-processed
# act.ideas.list, the list of the names of active ideas
# Output:
# votes, the data frame representing the valid votes
# (Mostly for internal use)
{
  votes = votes[votes$Valid,]
  votes = votes[votes$Left.Choice.ID %in% act.ideas.list,]
  votes = votes[votes$Right.Choice.ID %in% act.ideas.list,]
  return(votes)
}

give.vote.results <- function(votes)
# This function adds an attribute to the data frame of votes
# called "results". The results are set to be 1 of the idea shown 
# on the left won the comparison.
# Input:
# votes, the data frame representing the votes
# Output:
# votes, the data frame representing the votes, after adding "results"
# (Mostly for internal use)
{
  votes$results = (votes$Left.Choice.ID == votes$Winner.ID)
  return(votes)
}

preprocess.votes <- function(votes,ideas)
# This function completes the preprocessing of votes
# by running the three functions above
# Inputs:
# votes, the data frame representing the votes
# ideas, the data frame representing the preprocessed ideas
# Output:
# votes, the data frame representing the preprocessed votes
# (For external use)
{
  votes = rename.votes(votes)
  act.ideas.list = ideas$Idea.ID
  votes = get.valid.votes(votes,act.ideas.list)
  votes = give.vote.results(votes)
  return(votes)
}

  
show.num.app <- function(ideas)
# This function shows the distribution of the number of times 
# that each idea is involved in comparisons
# (For external use)
{
  ideas$Num.App = ideas$Wins + ideas$Losses
  ideas$Provider = 'Designer'
  ideas$Provider[ideas$User.Submitted] = 'Voter'
  ggplot(ideas,aes(Num.App, fill = Provider)) + 
    geom_histogram() +
    scale_fill_manual(values = c('Voter' = 'blue',
                                 'Designer' = 'orange')) +
    xlab('Number of appearances per idea') + 
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.position = c(0.9,0.9))
}


make.cont.matrix <- function(votes,
                             ideas = NULL,
                             act.ideas.list = NULL)
# This function creates the contingency matrix that can be used
# in the Bradley-Terry and Thurstone models
# Inputs:
# votes, the data frame representing the preprocessed votes
# ideas, the data frame representing the preprocessed ideas
# Output:
# cont.matrix, the contingency matrix
# (For external use)
{
  if(is.null(act.ideas.list))
    act.ideas.list = ideas$Idea.ID
  cont.matrix = matrix(nrow = length(act.ideas.list),
                       ncol = length(act.ideas.list))
  colnames(cont.matrix) = act.ideas.list
  rownames(cont.matrix) = act.ideas.list
  cont.matrix[,] = 0
  for(i in 1:nrow(votes))
  {
    cont.matrix[votes$Winner.ID[i],votes$Loser.ID[i]] = 
      cont.matrix[votes$Winner.ID[i],votes$Loser.ID[i]] + 1
  }
  return(cont.matrix)
}

# The following functions are used to generate data for
# STAN estimation of the hierarchical Thurstone model

get.theta.v.indexes <- function(votes)
# This function creates the indexes for theta_v in the 
# hierarchical Thurstone (Salganik-Levy) model
# (Mostly for internal use)
{
  theta.v.indexes = list()
  for(i in 1:nrow(votes))
  {
    left.object = votes$Left.Choice.ID[i]
    right.object = votes$Right.Choice.ID[i]
    session = votes$Session.ID[i]
    left.index = paste(left.object,session,sep=',')
    right.index = paste(right.object,session,sep=',')
    theta.v.indexes = c(theta.v.indexes,left.index,right.index)
  }
  theta.v.indexes = unique(theta.v.indexes)
  theta.v.indexes = unlist(theta.v.indexes)
  object.indexes = theta.v.indexes
  for(i in 1:length(theta.v.indexes))
  {
    object.indexes[i] = unlist(strsplit(theta.v.indexes[i],','))[1]
  }
  
  theta.v.indexes = theta.v.indexes[order(object.indexes)]
  object.indexes = object.indexes[order(object.indexes)]
  return(list(theta.v.indexes = theta.v.indexes,
              object.indexes = object.indexes))
}


make.stan.data <- function(votes,
                           sigma2 = 1,
                           tau0 = 4)
# This function creates the data for STAN implementation of the 
# hierarchical Thurstone (Salganik-Levy) model
# Inputs:
# votes, the data frame representing preprocessed votes
# Output:
# stan.data, a list that can be used for STAN estimation
{
  print('Preparing indexes......')
  # Prepare labels
  prep.indexes = get.theta.v.indexes(votes)
  theta.v.indexes = prep.indexes$theta.v.indexes
  #print(length(theta.v.labels))
  object.indexes = prep.indexes$object.indexes
  print('Indexes prepared.')
  
  # The design matrix
  print('Preparing the design matrix......')
  left.is = vector(length = nrow(votes))
  right.is = vector(length = nrow(votes))
  win.is = vector(length = nrow(votes))
  
  for(i in 1:nrow(votes))
  {
    left.object = votes$Left.Choice.ID[i]
    right.object = votes$Right.Choice.ID[i]
    session = votes$Session.ID[i]
    left.index = paste(left.object,session,sep=',')
    right.index = paste(right.object,session,sep=',')
    left.i = which(theta.v.indexes == left.index)
    right.i = which(theta.v.indexes == right.index)
    win.i = votes$results[i]
    left.is[i] = left.i
    right.is[i] = right.i
    win.is[i] = win.i
    if(i%%1000 == 0) print(i)
  }
  print('Design matrix prepared.')
  
  # N_k: the number that theta_jk is in theta_v, for every k
  print('Preparing Nk......')
  act.ideas.list = unique(object.indexes)
  Ns = vector(length = length(act.ideas.list) + 1)
  Ns[1] = 1
  cnt = 1
  for(i in 2:length(theta.v.indexes))
  {
    if(object.indexes[i] != object.indexes[i-1])
    {
      cnt = cnt + 1
      Ns[cnt] = i
    }
    if(i%%1000 == 0) print(i)
  }
  Ns[length(act.ideas.list) + 1] = length(theta.v.indexes) + 1
  print('Nk prepared.')
  
  stan.data <- list(V = nrow(votes),
                    K = length(act.ideas.list),
                    N = length(theta.v.indexes),
                    Ns = Ns,
                    left_indices = left.is,
                    right_indices = right.is,
                    y = as.integer(win.is),
                    mu0 = rep(0,length(act.ideas.list)),
                    tau0 = c(1e-3,rep(tau0,length(act.ideas.list)-1)),
                    sigma = sigma2,
                    theta.v.indexes = theta.v.indexes)
  print('STAN data prepared.')
  return(stan.data)
}
