library(mixOmics)
library(tidyverse)

#Question 1.1

mirna = read.csv("mirna.csv")
mrna = read.csv("mrna.csv")
protein = read.csv("protein.csv")
group = read.csv("sample_group.csv")
dim(mirna)
dim(mrna)
dim(protein)
dim(group)

#Question 1.2

f_cv= function(vector){
  cv = sd(vector)/mean(vector)
  return(cv)
}

#Question 1.3

f_cv_data = function(data){
  cv_data = rep(NA,dim(data)[2]-1)
  for(i in 1:length(cv_data)){
    cv_data[i] = f_cv(data[,i+1])
  }
  return(cv_data)
}

cv_mirna = f_cv_data(mirna)
hist(cv_mirna)

cv_mrna = f_cv_data(mrna)
hist(cv_mrna)

cv_protein = f_cv_data(protein)
hist(cv_protein)


