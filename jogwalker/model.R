
## Josephine Walker 
## https://github.com/jogwalker/InfDis_RforHTA

#Load libraries
library(ggplot2)
library(deSolve)
library(tidyverse)

# Initial conditions

initial_values=c(SP=6000,IP=4000,TIP=0,CP=0,TCP=0,RP=0,RCP=0,SX=0,IX=0,TIX=0,CX=0,TCX=0,RX=0,RCX=0,D=0) 

# Time points

time=seq(from=1,to=100,by=1)

# SIR model function 

sir_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # force of infection assuming only transmission in PWID
    NP=SP+IP+TIP+CP+TCP+RP+RCP
    lambda=beta*((IP+TIP+CP+TCP)/NP) 
    B = muP*NP + muX*(SX+IX+CX+TIX+TCX+RX+RCX)
    dSP= B -lambda*SP - tau*SP - muP*SP # birth (non HCV deaths), - infection, cessation, death
    dIP= lambda*SP + (1-svr)*tl*TIP -treat*treatP*IP  - zeta*IP - muP*IP - tau*IP# +infection, failed treatment, - treatment, -progression, -death, -cessation
    dTIP= treat*treatP*IP - tl*TIP - muP*TIP # treatment, - leaving treatment, - death
    dCP= -(treat*treatP+muC+muP+tau)*CP + zeta*IP + (1-svr)*tl*TCP# +progression, failed treatment, -death, - treatment, - cessation
    dTCP= treat*treatP*CP - tl*TCP - muC*TCP - muP*TCP
    dRP= svr*tl*TIP - muP*RP - tau*RP# recovered, death, cessation
    dRCP = svr*tl*TCP - muP*RCP -muC*RCP - tau*RCP
      
    dSX= tau*SP - muX*SX
    dIX= tau*IP - muX*IX - treat*IX + (1-svr)*tl*TIX
    dTIX= - muX*TIX + treat*IX - tl*TIX
    dCX= tau*CP - muX*CX - muC*CX - treat*CX + (1-svr)*tl*TCX
    dTCX= - muX*TCX - muC*TCX + treat*CX - tl*TCX
    dRX= tau*RP - muX*RX + svr*tl*TIX
    dRCX = tau*RCP - muX*RCX - muC*RCX + svr*tl*TCX

    dD= muC*(CP+RCP+TCP+CX+TCX+RCX) #+B don't track non-HCV deaths
      
    return(list(c(dSP,dIP,dTIP,dCP,dTCP,dRP,dRCP,dSX,dIX,dTIX,dCX,dTCX,dRX,dRCX,dD)))
  }
  )
}

# Baseline parameters
parameters=c(
  treat=0, # baseline rate of treatment 
  treatP=1,# scaling factor for treating PWID
  beta=0.3, # transmission rate per contact per year 
  tl=1/0.25, # rate of leaving treatment 
  svr=0.95, # cure rate
  muC=1/20, # death rate with cirrhosis
  muX=1/70, # death rate for exPWID
  muP=1/40, # death rate for PWID
  zeta=1/20, # progression to cirrhosis
  tau=1/20 # average injecting duration of 20 years
)

# Treatment parameters
parameters_treat=c(
  treat=0.2, # baseline rate of treatment 
  treatP=1,# scaling factor for treating PWID
  beta=0.3, # transmission rate per contact per year 
  tl=1/0.25, # rate of leaving treatment 
  svr=0.95, # cure rate
  muC=1/20, # death rate with cirrhosis
  muX=1/70, # death rate for exPWID
  muP=1/40, # death rate for PWID
  zeta=1/20, # progression to cirrhosis
  tau=1/20 # average injecting duration of 20 years
)

# Treat PWID at a higher rate
parameters_treatPWID=c(
  treat=0.2, # baseline rate of treatment 
  treatP=3,# scaling factor for treating PWID
  beta=0.3, # transmission rate per contact per year 
  tl=1/0.25, # rate of leaving treatment 
  svr=0.95, # cure rate
  muC=1/20, # death rate with cirrhosis
  muX=1/70, # death rate for exPWID
  muP=1/40, # death rate for PWID
  zeta=1/20, # progression to cirrhosis
  tau=1/20 # average injecting duration of 20 years
)

# Don't treat PWID
parameters_noPWID=c(
  treat=0.2, # baseline rate of treatment 
  treatP=0,# scaling factor for treating PWID
  beta=0.3, # transmission rate per contact per year 
  tl=1/0.25, # rate of leaving treatment 
  svr=0.95, # cure rate
  muC=1/20, # death rate with cirrhosis
  muX=1/70, # death rate for exPWID
  muP=1/40, # death rate for PWID
  zeta=1/20, # progression to cirrhosis
  tau=1/20 # average injecting duration of 20 years
)

#Solving the differential equations
output<-as.data.frame(ode(y=initial_values,func = sir_model,parms=parameters,times = time))
# head(output)

output_treat <- as.data.frame(ode(y=initial_values,func = sir_model,parms=parameters_treat,times = time))

output_treatPWID <- as.data.frame(ode(y=initial_values,func = sir_model,parms=parameters_treatPWID,times = time))

output_noPWID <- as.data.frame(ode(y=initial_values,func = sir_model,parms=parameters_noPWID,times = time))


outputAll <- bind_rows(list(Baseline=output,Treat=output_treat,TreatPWID=output_treatPWID,NoPWID=output_noPWID),.id="Scenario")
out_long <- outputAll %>% pivot_longer(3:ncol(.))

out_long$PWID <- ifelse(grepl("P",out_long$name),"PWID","Ex-PWID")
out_long$group <- out_long$name
out_long$group <- gsub("X","",out_long$group)
out_long$group <- gsub("P","",out_long$group)

# To plot number in each compartment over time
p1 <- ggplot(data = out_long,          
       aes(x = time, y = value/10000, colour = group)) +  
  geom_line() +xlab("Time (years)")+ylab("Proportion of the population")+ facet_grid(Scenario~PWID) + theme_minimal() + scale_color_brewer(type="qual",palette=2,name="State")

pdf("~/git/InfDis_RforHTA/plot1.pdf",width=8,height=6)
p1
dev.off()



## Costs and QALY weights

treatcost <- 375
cirrhosiscost <- 200

QoL <- c(0.94,0.77,0.77,0.55,0.55,0.82,0.61)
# Ex-PWID: Uninfected 0.94; Infected 0.77; Cirrhotic 0.55; Recovered 0.82; Recovered cirrhosis 0.61 # treatment as infected
QP <- (0.85/0.94)

outputEE <- outputAll

# costs
outputEE$treatcost <- (outputEE$TIP + outputEE$TIX + outputEE$TCP + outputEE$TCX)*treatcost
outputEE$cirrhosiscost <- (outputEE$CP + outputEE$CX + outputEE$RCP + outputEE$RCX)*cirrhosiscost
outputEE$tot.cost <- outputEE$treatcost + outputEE$cirrhosiscost

# outcomes
outputEE$QALY.X <- rowSums(t(apply(outputEE[,3:9],1,function(x){x*QoL})))
outputEE$QALY.P <- rowSums(t(apply(outputEE[,10:16],1,function(x){x*QoL*QP})))
outputEE$QALYs <- outputEE$QALY.X + outputEE$QALY.P

# but for DALYs deaths count as 1 (for simplicity assume DALY weights are inverse of QALY weights)
outputEE$DALY.X <- rowSums(t(apply(outputEE[,3:9],1,function(x){x*(1-QoL)})))
outputEE$DALY.P <- rowSums(t(apply(outputEE[,10:16],1,function(x){x*(1-QoL*QP)})))
outputEE$DALYs <- outputEE$DALY.X + outputEE$DALY.P + outputEE$D

## discount
discounting <- function(years,value,rate,baseyear) { 
    time <- years - baseyear 
    out <- value/(1+rate)^time
    return(out)
}

outputEE$costD <- discounting(outputEE$time,outputEE$tot.cost,0.03,1)
outputEE$QALYD <- discounting(outputEE$time,outputEE$QALYs,0.03,1)
outputEE$DALYD <- discounting(outputEE$time,outputEE$DALYs,0.03,1)
outputEE$population <- rowSums(outputEE[,3:17])

## calculate mean ICER
pop = 10000
TotalCE <- outputEE %>% group_by(Scenario) %>% summarise(costD=sum(costD)/pop,DALYD=sum(DALYD)/pop,QALYD=sum(QALYD)/pop) %>% arrange(QALYD) # sum over time horizon
TotalCE$diffcost <- c(NA,diff(TotalCE$costD))
TotalCE$diffQALY <-  c(NA,diff(TotalCE$QALYD))
TotalCE$ICER <- TotalCE$diffcost / TotalCE$diffQALY
TotalCE

TotalCEshort <- outputEE %>% filter(time<5) %>% group_by(Scenario) %>% summarise(costD=sum(costD)/pop,DALYD=sum(DALYD)/pop,QALYD=sum(QALYD)/pop) %>% arrange(QALYD)
TotalCEshort$diffcost <- c(NA,diff(TotalCEshort$costD))
TotalCEshort$diffQALY <-  c(NA,diff(TotalCEshort$QALYD))
TotalCEshort$ICER <- TotalCEshort$diffcost / TotalCEshort$diffQALY
TotalCEshort


                        