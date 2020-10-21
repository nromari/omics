library(mixOmics)
library(tidyverse)

#Parti I - Preparation

mirna = read.csv("mirna.csv")
mrna = read.csv("mrna.csv")
protein = read.csv("protein.csv")
group = read.csv("sample_group.csv")

  #Question 1.1

dim(mirna) #150 ech et 185 variables
dim(mrna)  #150 ech et 201 variables
dim(protein) #150 ech et 143 variables
dim(group) #150 ech et 2 variables

  #Question 1.2

f_cv= function(vector){
  cv = sd(vector)/mean(vector)
  return(cv)
}

  #Question 1.3

f_cv_bloc = function(bloc){
  cv_bloc = apply(bloc[,2:dim(bloc)[2]], MARGIN=2, FUN=f_cv)
  return(cv_bloc)
}

#f_cv_bloc = function(bloc){
  #cv_bloc = rep(NA,dim(bloc)[2]-1)
  #for(i in 1:length(cv_bloc)){
    #cv_bloc[i] = f_cv(bloc[,i+1])
  #}
  #return(cv_bloc)
#}

  #Question 1.4

cv_mirna = f_cv_bloc(mirna)
hist(cv_mirna)

cv_mrna = f_cv_bloc(mrna)
hist(cv_mrna)

cv_protein = f_cv_bloc(protein)
hist(cv_protein)

  #Visuelement il semble que le bloc protéine présente la distribution de cv la plus variable.
  #Calculons les étendues et les variances pour confirmer

 #Question 2.1

f_variance = function (cv_bloc){
  cv_bloc_var = var(cv_bloc)
  return(cv_bloc_var)
}

f_etendue = function (cv_bloc){
  cv_bloc_etendue = max(cv_bloc)-min(cv_bloc)
  return(cv_bloc_etendue)
}

print(round(f_variance(cv_mirna),4)) #var = 0.0053
print(round(f_etendue(cv_mirna),2)) #etendue =0.33

print(round(f_variance(cv_mrna),4)) #var = 0.0039
print(round(f_etendue(cv_mrna),2)) #etendue = 0.37

print(round(f_variance(cv_protein),4)) #var = 1070.215
print(round(f_etendue(cv_protein),2)) #etendue = 335.01
  #On confirme la grande variabilité dans le bloc protéine,
  #en lien avec la variabilité des évènements traductionels et post traductionnels (epissage, repliement...)

  #Question 2.2

#f_filtrage = function(bloc, cv_bloc){
  #for (i in length(cv_bloc)){
    #if(abs(cv_bloc[i])>=0.15){
     # bloc_f=bloc[,-(i+1)]
   # }else{
      #bloc_f=bloc[,(i+1)]
   # }
  #}
  #return(bloc_f)
#}
#Cette 1ere fonction m'a donné des résultats cohérents avec ceux de la classe (reste 142 prot)
#Mais inchérente avec les vrais valeurs et avec les hist...
#C'est la 2ème qui est cohérente... A voir avec le prof

f_filtrage = function(bloc, cv_bloc){
  filtre = which(cv_bloc>=0.15)
  bloc_f = bloc[,-filtre]
  return(bloc_f)
}

 #Question 3

mirna_f = f_filtrage(mirna, cv_mirna)
dim(mirna_f) #Reste 103
mrna_f = f_filtrage(mrna, cv_mrna)
dim(mrna_f) #Reste 27
protein_f = f_filtrage(protein, cv_protein)
dim(protein_f) #Reste 53

 #Question 4

# j'ai supposé qu'on devait utiliser les résultats non filtrés... à voir avec le prof
which(cv_mrna==max(cv_mrna))
 # Le gene le plus variable est PLCD4 (conordant avec le groupe). il code pour la phopholipase C delta 4

which(colnames(protein)=="PLCD4")
find("Phospholipase C Delta 4",colnames(protein))
# j'ai essayé plusieurs combinaisons... introuvable
