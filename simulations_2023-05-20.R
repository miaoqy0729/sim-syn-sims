# Simple Synchrony Simulations Model
# Miao, Dale & Galati (2023)

set.seed(42)
source('functions_2023-05-20.R')

# Set parameters 
allResults = c()
allStivers = c()

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
    #print(res[1])
    
    ccfs = t(as.numeric(res[[3]][[1]]))
    colnames(ccfs) = gsub('-','m',paste('ccf_lag',-20:20,sep=''))
    allStivers = rbind(allStivers, data.frame(row.names=paste(C_i,'-',sim_j,'-',1:length(res[[4]]),sep=''),
                                              C_i=C_i,
                                              sim_j=sim_j,
                                              C[C_i,],
                                              cond=paste(C_line,collapse=','),
                                              mids=as.numeric(res[[4]]),
                                              counts=as.numeric(res[[5]]/sum(res[[5]]))))
    allResults = rbind(allResults,data.frame(C_i=C_i,
                                             sim_j=sim_j,
                                             C[C_i,],
                                             cond=paste(C_line,collapse=','),
                                             r=cor(res[[1]],res[[2]]),
                                             ccfs
                                             #h_mids=res[[4]],
                                             #h_mids[,length(res[[4]])]=res[[4]]
                                             #h_counts[,]=res[[5]]
                                             )
                       )
  }
}

save(allStivers, file = "allStivers_2023-04-18.Rd")
save(allResults, file = "allResults_2023-04-18.Rd")

allStivers$midBins = round(allStivers$mids)
plot(aggregate(counts~midBins,data=allStivers,mean),type='o',xlab='Lag to next turn',ylab='% of turn switches in convo')

pos = 0
diff_pos = c()

for (i in 1:length(allStivers$midBins)){
  print(allStivers$midBins[i])
  print(i)
  
  compare = v1[i] == v2[i]
  if (compare == FALSE){
    diff_pos <- append(diff_pos, i)
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

res = runModel(matrix(c(0,1,1,-1),nrow=2,ncol=2),500)
plotResults(res[[1]][1:100],res[[2]][1:100],'C = (0, 1; 1, -1)','blue','red')
plotResults(res[[1]][1:200],res[[2]][1:200],'C = (0, 1; 1, -1)','blue','red')
plotResults(res[[1]],res[[2]],'C = (0, 1; 1, -1)','blue','red')

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
  



############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
# NEW! for main simulation/analysis section:
#
#

#### CCF ####

# s1 o1 o2 s2
# C = (1, 1, 0, 0)
res1 = allResults[allResults$o_1==1&allResults$o_2==0&allResults$s_1==1&allResults$s_2==0,]
# C = (1, -1, 0, 0)
res2 = allResults[allResults$o_1==-1&allResults$o_2==0&allResults$s_1==1&allResults$s_2==0,]
# C = (1, 0, 0, 1)
res3 = allResults[allResults$o_1==0&allResults$o_2==0&allResults$s_1==1&allResults$s_2==1,]
# C = (-1, 0, 1, 1)
res4 = allResults[allResults$o_1==0&allResults$o_2==1&allResults$s_1==-1&allResults$s_2==1,]

p = plotCcf(res1,'pink')
p = plotCcf(res2,'orange',p_add=p)
p = plotCcf(res3,'green',p_add=p)
p = plotCcf(res4,'purple',p_add=p)
p

# C = (-1, 1, 1, -1)
res4 = allResults[(allResults$o_1==1&allResults$o_2==1)&(allResults$s_1==-1&allResults$s_2==-1),]
p = plotCcf(res4,'red')
p

#### Stivers ####

allStivers$midBin = round(allStivers$mids/10)*10

# First, calculate the sum of counts within each sim_j value
res_sum = allStivers %>%
  group_by(sim_j) %>%
  summarise(total_counts = sum(counts))

res_with_total_counts = allStivers %>%
  left_join(res_sum, by = "sim_j")

allStivers$normalizedCounts = (res_with_total_counts %>%
  mutate(proportion = counts / total_counts))$proportion

# s1 o1 o2 s2
# C = (1, 1, 0, 0)
res1 = allStivers[allStivers$o_1==1&allStivers$o_2==0&allStivers$s_1==1&allStivers$s_2==0,]
# C = (1, -1, 0, 0)
res2 = allStivers[allStivers$o_1==-1&allStivers$o_2==0&allStivers$s_1==1&allStivers$s_2==0,]
# C = (1, 0, 0, 1)
res3 = allStivers[allStivers$o_1==0&allStivers$o_2==0&allStivers$s_1==1&allStivers$s_2==1,]
# C = (-1, 0, 1, 1)
res4 = allStivers[allStivers$o_1==0&allStivers$o_2==1&allStivers$s_1==-1&allStivers$s_2==1,]
# C = (-1, 1, 1, -1)
res5 = allStivers[allStivers$o_1==1&allStivers$o_2==1&allStivers$s_1==-1&allStivers$s_2==-1,]

p = plotStivers(res1,'pink')
p = plotStivers(res2,'orange',p_add=p)
p = plotStivers(res3,'green',p_add=p)
p = plotStivers(res4,'purple',p_add=p)
p = plotStivers(res5,'red')
p


  
