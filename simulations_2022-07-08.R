# Simple Synchrony Simulations Model
# Miao, Dale & Galati (2022)

set.seed(42)
setwd("../") # Swap this with your working directory 
source('functions_2022-07-08.R')

# Set parameters 
allResults = c()
#load("allResults_2022-07-08.Rd") # if loading from prior .Rd file
d=c(-1,0,1)
C = expand.grid(d,d,d,d)
colnames(C) = c('s_1','o_1','o_2','s_2')

# Generate simulations
for (C_i in 1:nrow(C)) {
  print(C_i)
  C_line = as.numeric(C[C_i,])
  for (sim_j in 1:100) {
    res = runModel(matrix(C_line,nrow=2,ncol=2),iterations=500)  
    allResults = rbind(allResults,data.frame(C_i=C_i,
                                             sim_j=sim_j,
                                             C[C_i,],
                                             cond=paste(C_line,collapse=','),
                                             r=cor(res[[1]],res[[2]])))
  }
}

# Visualize results
head(allResults)
hist(allResults$r)

# Visualize behavioral time series plots with specific context matrix C
res = runModel(matrix(c(1,0,1,0),nrow=2,ncol=2),500)
plotResults(res[[1]],res[[2]],'C = (1, 0; 1, 0)','blue','red')

res = runModel(matrix(c(1,1,1,1),nrow=2,ncol=2),500)
plotResults(res[[1]],res[[2]],'C = (1, 1; 1, 1)','blue','red')

res = runModel(matrix(c(1,1,1,0),nrow=2,ncol=2),500)
plotResults(res[[1]],res[[2]],'C = (1, 1; 1, 0)','blue','red')

res = runModel(matrix(c(0,-1,-1,0),nrow=2,ncol=2),500)
plotResults(res[[1]],res[[2]],'C = (0, -1; -1, 0)','blue','red')


# Statistical analysis of turn-taking pattern: with and without inhibition in C matrix 

## View all cases with strong negative r value 
allResults[allResults$r<(-.25),][,3:6]

## Obtain tables counting the number of cases with -1 in the C matrix
negativeCorTable = table(apply(allResults[allResults$r<(-.25),][,3:6],1,function(x){return(-1 %in% x)}))
positiveCorTable = table(apply(allResults[allResults$r>(.25),][,3:6],1,function(x){return(-1 %in% x)}))

negativeCorTable/sum(negativeCorTable)
positiveCorTable/sum(positiveCorTable)

## Chi-square test on uniform probability
chisq.test(negativeCorTable)

## Chi-square test on corresponding positive correlation
chisq.test(negativeCorTable,p=positiveCorTable/sum(positiveCorTable))



# Multicategorical Analysis (must be run after loading data or running simulations)
allResults <- transform(allResults,
                          s_1p = (s_1 == 1)*1,
                          s_1n = (s_1 == -1)*1,
                          o_1p = (o_1 == 1)*1,
                          o_1n = (o_1 == -1)*1,
                          o_2p = (o_2 == 1)*1,
                          o_2n = (o_2 == -1)*1,
                          s_2p = (s_2 == 1)*1,
                          s_2n = (s_2 == -1)*1)


# Linear Interaction Analysis
modelSub01 <- lm(r ~ s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n), data=allResults)
summary(modelSub01)
  #Multiple R-squared:  0.4045,	Adjusted R-squared:  0.4035 

modelSub02 <- lm(r ~ s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n), data=allResults)
summary(modelSub02)
  #Multiple R-squared:  0.4273,	Adjusted R-squared:  0.4264

modelSub03 <- lm(r ~ o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n), data=allResults)
summary(modelSub03)
  #Multiple R-squared:  0.4055,	Adjusted R-squared:  0.4046 

modelSub04 <- lm(r ~ o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n), data=allResults)
summary(modelSub04)
  #Multiple R-squared:  0.4057,	Adjusted R-squared:  0.4048 

modelSub05 <- lm(r ~ o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n), data=allResults)
summary(modelSub05)
  #Multiple R-squared:  0.406,	Adjusted R-squared:  0.405 

modelSub06 <- lm(r ~ o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n), data=allResults)
summary(modelSub06)
  #Multiple R-squared:  0.4047,	Adjusted R-squared:  0.4038 

modelSub07 <- lm(r ~ s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n), data=allResults)
summary(modelSub07)
  #Multiple R-squared:  0.4038,	Adjusted R-squared:  0.4029

modelSub08 <- lm(r ~ s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n), data=allResults)
summary(modelSub08)
  #Multiple R-squared:  0.4266,	Adjusted R-squared:  0.4256

# Overall regression model
modelAllEverything <- lm(r ~ s_1p + o_1p + o_2p + s_2p + s_1n + o_1n + o_2n + s_2n 
                         + (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                         + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                         + (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                         + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                         + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                         + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                         + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                         + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelAllEverything)
  #Multiple R-squared:  0.9451,	Adjusted R-squared:  0.9449


# Further exploration: Removing different parts of the model to see what minimal combinations give us R-square > 0.9 
## A model that only includes interactions  
modelAllInteractionsOnly <- lm(r ~ (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                               + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                               + (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                               + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                               + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                               + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                               + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                               + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelAllInteractionsOnly)
  #Multiple R-squared:  0.9451,	Adjusted R-squared:  0.9449 

## A model that only includes positive predictors 
modelInteractionsP <- lm(r ~ (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                               + (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                               + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                               + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelInteractionsP)
  #Multiple R-squared:  0.8856,	Adjusted R-squared:  0.8853 

## A model that only includes negative predictors 
modelInteractionsN <-lm(r ~ (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                        + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                        + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                        + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelInteractionsN)
  #Multiple R-squared:  0.8964,	Adjusted R-squared:  0.8961 
  #Interesting! All negative components has slightly greater R-square than all positive.

## A model that includes Person A's influences on self and on Person B
modelPersonA1 <- lm(r ~ s_1p + o_2p + s_1n + o_2n  
                   + (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                   + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                   + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                   + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n)), data=allResults)
summary(modelPersonA1)
  #Multiple R-squared:  0.9084,	Adjusted R-squared:  0.9081

## A model that includes Person A's influences by self and by Person B 
modelPersonA2 <- lm(r ~ s_1p + s_2p + s_1n + s_2n 
                    + (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                    + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                    + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                    + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelPersonA2)
  #Multiple R-squared:  0.9405,	Adjusted R-squared:  0.9403 

## A model that includes Person B's influences on self and on Person A
modelPersonB1 <- lm(r ~ o_1p + s_2p + o_1n + s_2n 
                    + (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                    + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                    + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                    + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelPersonB1)
  #Multiple R-squared:  0.9082,	Adjusted R-squared:  0.9079 

## A model that includes Person B's influences by self and by Person A 
modelPersonB2 <- lm(r ~ o_1p + o_2p + o_1n + o_2n 
                    + (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                    + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                    + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                    + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n)), data=allResults)
summary(modelPersonB2)
  #Multiple R-squared:  0.9451,	Adjusted R-squared:  0.9449



## Removing the main effect terms do not impact R-square much 
modelInteractionsPersonA1 <- lm(r ~ (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                                + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                                + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                                + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n)), data=allResults)
summary(modelInteractionsPersonA1)
  #Multiple R-squared:  0.9084,	Adjusted R-squared:  0.9081 


modelInteractionsPersonA2 <- lm(r ~ (s_1p * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n))
                                + (s_1n * (o_1p + o_2p + s_2p + o_1n + o_2n + s_2n)) 
                                + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                                + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelInteractionsPersonA2)
  #Multiple R-squared:  0.9405,	Adjusted R-squared:  0.9403

modelInteractionsPersonB1 <- lm(r ~ (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                                + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                                + (s_2p * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n))
                                + (s_2n * (s_1p + o_1p + o_2p + s_1n + o_1n + o_2n)), data=allResults)
summary(modelInteractionsPersonB1)
  #Multiple R-squared:  0.9082,	Adjusted R-squared:  0.9079 

modelInteractionsPersonB2 <- lm(r ~ (o_1p * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                                + (o_1n * (s_1p + o_2p + s_2p + s_1n + o_2n + s_2n))
                                + (o_2p * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n))
                                + (o_2n * (s_1p + o_1p + s_2p + s_1n + o_1n + s_2n)), data=allResults)
summary(modelInteractionsPersonB2)
  #Multiple R-squared:  0.9451,	Adjusted R-squared:  0.9449
  


  
