# Author: Felicity Lamrock

# MSM code and data from the MSM package


###########
# Set up  #
###########

install.packages("msm")
library(msm)

# Pre-existing data set
data<-cav

cav[1:21,]

str(cav)

statetable.msm(state, PTNUM, data=cav)




# Fitting the msm model

Q  <-  rbind ( c(0, 0.25, 0, 0.25),
               c(0.166, 0, 0.166, 0.166),
               c(0, 0.25, 0, 0.25),
               c(0, 0, 0, 0) )

Q.crude  <- crudeinits.msm(state ~ years, PTNUM, data=cav,qmatrix=Q)

cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,qmatrix = Q, deathexact = 4)

# Show the msm model results

cav.msm

# We want to use the covariates in the model

cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav, qmatrix = Q, deathexact = 4, covariates = ~ sex)

cavsex.msm

# Add covariates in one by one using previous Q matrix

model0.msm <- msm( state ~ years, subject=PTNUM, data = data, qmatrix = Q, 
                   control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), 
                   deathexact = 4)

qq<-qmatrix.msm(model0.msm)$estimates

model1.msm <- msm( state ~ years, subject=PTNUM, data = data, qmatrix = qq, 
                  control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), 
                  deathexact = 4, covariates = ~ sex)

qq<-qmatrix.msm(model1.msm)$estimates

model2.msm <- msm( state ~ years, subject=PTNUM, data = data, qmatrix = qq, 
                  control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), 
                  deathexact = 4, covariates = ~ sex + cumrej)

qq<-qmatrix.msm(model2.msm)$estimates

model3.msm <- msm( state ~ years, subject=PTNUM, data = data, qmatrix = qq, 
                   control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), 
                   deathexact = 4, covariates = ~ sex + cumrej + dage)

# Added in control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), for covergence issues and to watch it!

# We could have added all in at once...
model.msm <- msm( state ~ years, subject=PTNUM, data = data, qmatrix = Q, 
                   control=list(trace=1, REPORT=1, maxit=1000, fnscale=100000), 
                   deathexact = 4, covariates = ~ sex + cumrej + dage)

model.msm
model3.msm

# Look at the very slight differences!

# Can loop over different starting Q matrix values and pick the best msm model using model.msm$minus2loglik

# q.list <- boot.msm(model.msm, stat=function(x){qmatrix.msm(x)$estimates})
#q.array <- array(unlist(q.list), dim=c(4,4,1000))
#apply(q.array, c(1,2), sd)
#apply(q.array, c(1,2), function(x)quantile(x, c(0.025, 0.975)))


##################
#  Probabilities #
##################


# What if we want the probabilities? We can use the inbuilt functions...
pmatrix.msm(model.msm, t=10)

# Assumes Q is constant within the desired time interval

# Cut points for when the measurement for the cumulative number of rejection changes - would be better if it was a treatment change example!
times <- c(0, 5)
# Cut points for the change in cumulative number of rejections
covariates <- list( list (cumrej=0), list (cumrej=0), list(cumrej=1))

pmatrix.piecewise.msm(model3.msm, 0, 10, times, covariates)


pmatrix.msm(model.msm, t=10, ci="normal")


# The option ci="normal" computes a confidence interval for P(t) 
# by repeated sampling from the asymptotic normal distribution of the maximum likelihood estimates of the log(qrs). 
# Default 1000 samples

pmatrix.msm(model.msm, t=10, ci="boot")

# Nonparametric bootstrap resampling (ci="boot"). 
# Bootstrap datasets of transitions are drawn with replacement and the model refitted repeatedly to estimate the 
# sampling uncertainty surrounding the estimates. 
# This method is more accurate but much slower due to the need to refit the model for each resample.




# I want the probability for EACH person, using lots of covariates



# Use Q matrices 
# Excuse the bad coding!

Element1<-  rbind(c(0, 1, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0) )


Element2<-  rbind(c(0, 0, 0, 1),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0) )


Element3<-  rbind(c(0, 0, 0, 0),
               c(1, 0, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0) )

Element4<-  rbind(c(0, 0, 0, 0),
               c(0, 0, 1, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0) )

Element5<-  rbind(c(0, 0, 0, 0),
               c(0, 0, 0, 1),
               c(0, 0, 0, 0),
               c(0, 0, 0, 0) )

Element6<-  rbind(c(0, 0, 0, 0),
               c(0, 0, 0, 0),
               c(0, 1, 0, 0),
               c(0, 0, 0, 0) )

Element7<-  rbind(c(0, 0, 0, 0),
               c(0, 0, 0, 0),
               c(0, 0, 0, 1),
               c(0, 0, 0, 0) )

base1<- model.msm$Qmatrices$baseline[1,2]  
base2<- model.msm$Qmatrices$baseline[1,4] 
base3<- model.msm$Qmatrices$baseline[2,1]  
base4<- model.msm$Qmatrices$baseline[2,3]
base5<- model.msm$Qmatrices$baseline[2,4]
base6<- model.msm$Qmatrices$baseline[3,2]
base7<- model.msm$Qmatrices$baseline[3,4]


QBase1<-Element1*base1
QBase2<-Element2*base2
QBase3<-Element3*base3
QBase4<-Element4*base4
QBase5<-Element5*base5
QBase6<-Element6*base6
QBase7<-Element7*base7

QBase<-QBase1+QBase2+QBase3+QBase4+QBase5+QBase6+QBase7


# Long winded way to make model.msm$Qmatrices$baseline QBase without the State information for the loops

model.msm$Qmatrices$baseline
QBase

# Note the diagonals don't need to be calculated as they are sum of the rest of the row set to 0

# Do this for the other covariates now. Take the exponential of the covariate ones though so that you don't have zeros on the diagnoal

sex1<- exp(model.msm$Qmatrices$sex[1,2])
sex2<- exp(model.msm$Qmatrices$sex[1,4]) 
sex3<- exp(model.msm$Qmatrices$sex[2,1])  
sex4<- exp(model.msm$Qmatrices$sex[2,3])
sex5<- exp(model.msm$Qmatrices$sex[2,4])
sex6<- exp(model.msm$Qmatrices$sex[3,2])
sex7<- exp(model.msm$Qmatrices$sex[3,4])


QSex1<-Element1*sex1
QSex2<-Element2*sex2
QSex3<-Element3*sex3
QSex4<-Element4*sex4
QSex5<-Element5*sex5
QSex6<-Element6*sex6
QSex7<-Element7*sex7

QSex<-QSex1+QSex2+QSex3+QSex4+QSex5+QSex6+QSex7

# So QSex is the exponential of model.msm$Qmatrices$sex

QSex
exp(model.msm$Qmatrices$sex)

# Cumrej

cumrej1<- exp(model.msm$Qmatrices$cumrej[1,2])
cumrej2<- exp(model.msm$Qmatrices$cumrej[1,4]) 
cumrej3<- exp(model.msm$Qmatrices$cumrej[2,1])  
cumrej4<- exp(model.msm$Qmatrices$cumrej[2,3])
cumrej5<- exp(model.msm$Qmatrices$cumrej[2,4])
cumrej6<- exp(model.msm$Qmatrices$cumrej[3,2])
cumrej7<- exp(model.msm$Qmatrices$cumrej[3,4])


QCumrej1<-Element1*cumrej1
QCumrej2<-Element2*cumrej2
QCumrej3<-Element3*cumrej3
QCumrej4<-Element4*cumrej4
QCumrej5<-Element5*cumrej5
QCumrej6<-Element6*cumrej6
QCumrej7<-Element7*cumrej7

QCumrej<-QCumrej1+QCumrej2+QCumrej3+QCumrej4+QCumrej5+QCumrej6+QCumrej7

# So QCumrej is the exponential of model.msm$Qmatrices$cumrej

QCumrej 
exp(model.msm$Qmatrices$cumrej)


# Final covariate

dage1<- exp(model.msm$Qmatrices$dage[1,2])
dage2<- exp(model.msm$Qmatrices$dage[1,4]) 
dage3<- exp(model.msm$Qmatrices$dage[2,1])  
dage4<- exp(model.msm$Qmatrices$dage[2,3])
dage5<- exp(model.msm$Qmatrices$dage[2,4])
dage6<- exp(model.msm$Qmatrices$dage[3,2])
dage7<- exp(model.msm$Qmatrices$dage[3,4])


QDage1<-Element1*dage1
QDage2<-Element2*dage2
QDage3<-Element3*dage3
QDage4<-Element4*dage4
QDage5<-Element5*dage5
QDage6<-Element6*dage6
QDage7<-Element7*dage7

QDage<-QDage1+QDage2+QDage3+QDage4+QDage5+QDage6+QDage7

# So QDage is the exponential of model.msm$Qmatrices$dage

QDage
exp(model.msm$Qmatrices$dage)

# Now we have QBase, QSex, Qcumrej and QDage

# The equation for the model is qrs = qoexp(beta z(t))


# Only looking at probability from CAV-free to mild CAV in this example
# Moving from state 1 to state 2
# We want to look at the individuals who are in state 1 currently



# Dimensions of the data are 2846 observations/records/row
# It is a 4 x 4 matrix
# Set all the values to -99 so you can see when they're filled up

install.packages("expm")
library(expm)

Probability_CAV<-array(-99,dim=c(2846))


for (k in 1:2846)
{record<-data[k,]

Qmatrix<-QBase *
  QSex^record$sex*
  QCumrej^record$cumrej * 
  QDage^record$dage


Qmatrix[1,1]<- -sum(Qmatrix[1,2]+Qmatrix[1,3]+Qmatrix[1,4])
Qmatrix[2,2]<- -sum(Qmatrix[2,1]+Qmatrix[2,3]+Qmatrix[2,4])
Qmatrix[3,3]<- -sum(Qmatrix[3,1]+Qmatrix[3,2]+Qmatrix[3,4])
Qmatrix[4,4]<- -sum(Qmatrix[4,1]+Qmatrix[4,2]+Qmatrix[4,3])

Pmatrix<-expm(Qmatrix)

Probability_CAV[k]<-Pmatrix[1,2]


}

data$P_CAV<-Probability_CAV



########################
# Confidence intervals #
########################



# Same calculation as before, but we can sum these probabilities over the cohort after we simulate more information
# Use lower and upper confidence intervals of Q matrix

base1<- model.msm$QmatricesL$baseline[1,2]  
base2<- model.msm$QmatricesL$baseline[1,4] 
base3<- model.msm$QmatricesL$baseline[2,1]  
base4<- model.msm$QmatricesL$baseline[2,3]
base5<- model.msm$QmatricesL$baseline[2,4]
base6<- model.msm$QmatricesL$baseline[3,2]
base7<- model.msm$QmatricesL$baseline[3,4]


QBase1<-Element1*base1
QBase2<-Element2*base2
QBase3<-Element3*base3
QBase4<-Element4*base4
QBase5<-Element5*base5
QBase6<-Element6*base6
QBase7<-Element7*base7

QBase<-QBase1+QBase2+QBase3+QBase4+QBase5+QBase6+QBase7


sex1<- exp(model.msm$QmatricesL$sex[1,2])
sex2<- exp(model.msm$QmatricesL$sex[1,4]) 
sex3<- exp(model.msm$QmatricesL$sex[2,1])  
sex4<- exp(model.msm$QmatricesL$sex[2,3])
sex5<- exp(model.msm$QmatricesL$sex[2,4])
sex6<- exp(model.msm$QmatricesL$sex[3,2])
sex7<- exp(model.msm$QmatricesL$sex[3,4])


QSex1<-Element1*sex1
QSex2<-Element2*sex2
QSex3<-Element3*sex3
QSex4<-Element4*sex4
QSex5<-Element5*sex5
QSex6<-Element6*sex6
QSex7<-Element7*sex7

QSex<-QSex1+QSex2+QSex3+QSex4+QSex5+QSex6+QSex7


cumrej1<- exp(model.msm$QmatricesL$cumrej[1,2])
cumrej2<- exp(model.msm$QmatricesL$cumrej[1,4]) 
cumrej3<- exp(model.msm$QmatricesL$cumrej[2,1])  
cumrej4<- exp(model.msm$QmatricesL$cumrej[2,3])
cumrej5<- exp(model.msm$QmatricesL$cumrej[2,4])
cumrej6<- exp(model.msm$QmatricesL$cumrej[3,2])
cumrej7<- exp(model.msm$QmatricesL$cumrej[3,4])


QCumrej1<-Element1*cumrej1
QCumrej2<-Element2*cumrej2
QCumrej3<-Element3*cumrej3
QCumrej4<-Element4*cumrej4
QCumrej5<-Element5*cumrej5
QCumrej6<-Element6*cumrej6
QCumrej7<-Element7*cumrej7

QCumrej<-QCumrej1+QCumrej2+QCumrej3+QCumrej4+QCumrej5+QCumrej6+QCumrej7


dage1<- exp(model.msm$QmatricesL$dage[1,2])
dage2<- exp(model.msm$QmatricesL$dage[1,4]) 
dage3<- exp(model.msm$QmatricesL$dage[2,1])  
dage4<- exp(model.msm$QmatricesL$dage[2,3])
dage5<- exp(model.msm$QmatricesL$dage[2,4])
dage6<- exp(model.msm$QmatricesL$dage[3,2])
dage7<- exp(model.msm$QmatricesL$dage[3,4])


QDage1<-Element1*dage1
QDage2<-Element2*dage2
QDage3<-Element3*dage3
QDage4<-Element4*dage4
QDage5<-Element5*dage5
QDage6<-Element6*dage6
QDage7<-Element7*dage7

QDage<-QDage1+QDage2+QDage3+QDage4+QDage5+QDage6+QDage7

# Bootstrap the msm models and run the same code
# Take the average of the cohort to get the probability
# Run this probability for time dependent covariates 

# Take the 95% upper and lower values for these simulations to be the 95% confidence intervals

Probability_CAV_CI_L<-array(-99,dim=c(2846))

for (k in 1:2846)
{record<-data[k,]

Qmatrix<-QBase *
  QSex^record$sex*
  QCumrej^record$cumrej * 
  QDage^record$dage


Qmatrix[1,1]<- -sum(Qmatrix[1,2]+Qmatrix[1,3]+Qmatrix[1,4])
Qmatrix[2,2]<- -sum(Qmatrix[2,1]+Qmatrix[2,3]+Qmatrix[2,4])
Qmatrix[3,3]<- -sum(Qmatrix[3,1]+Qmatrix[3,2]+Qmatrix[3,4])
Qmatrix[4,4]<- -sum(Qmatrix[4,1]+Qmatrix[4,2]+Qmatrix[4,3])

Pmatrix<-expm(Qmatrix)

Probability_CAV_CI_L[k]<-Pmatrix[1,2]


}

data$P_CAV_CI_L<-Probability_CAV_CI_L

# Compare the previous probability and this lower confidence interval in the data set

#################
# Final remarks #
#################

# Use only those with the first record for CAV-free to mild-CAV example

data_CAV_free<-data[!duplicated(data$PTNUM), ]

# We have 622 patients where we are predicting the next year if they have mild CAV
# Note we don't use the time dependent variables in the calculations for this illustration
# Assume that all of the variables stay the same

# Run simulations of starting Q matrix
# Add covariates one by one
# Use upper and lower CI if wanting a probability for EACH individual over time
# Can get a probability for EACH individual EACH year - and incorporate those changing Q matrices
# Average over the cohort and take these probabilities in CEM to see what difference that makes



