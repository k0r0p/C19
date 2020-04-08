rm(list = ls())

setwd("~/Dropbox/Research/Corona/C19/")

library(deSolve)
library(tidyr)
library(reshape)
library(magrittr)
library(plyr)
library(ggplot2)
library(gdata)
library(matlib)

#CONTACT MATRIX
n_mat <- read.xls("Data/nepal_cmat.xls", header = F)

n_mat1 <- n_mat[,-1]                             # contact matrix without the labeling column
n_matA <- unname(as.matrix(n_mat1[1:16,]))       # all contact matrix
n_matH <- unname(as.matrix(n_mat1[17:32,]))      # home contact matrix
n_matS <- unname(as.matrix(n_mat1[33:48,]))      # school contact matrix
n_matW <- unname(as.matrix(n_mat1[49:64,]))      # work contact matrix
n_matO <- unname(as.matrix(n_mat1[65:80,]))      # other contact matrix

# matrix of total number of contacts
n.matT <- unname(n_matH + n_matS + n_matW + n_matO) # sum matrices and remove row and col names

# duration per contact
avg.cont <- mean(apply(n.matT, 1, sum))      # average contacts under normal circumstances 15.6 contacts   
avg.dur <- 1/(apply(n.matT, 1, sum))         # vector of contact duration by age cohort
avg.dur.mat <- matrix(rep(avg.dur,16), ncol = 16, byrow=F) #matrix of contact duration by age cohort


#INTERVENTION EFFECTS
# i16 x 16 identity matrix
id.mat <- diag(16)

# matrix of contacts under lockdown /social distancing scenario
ld.f <- (1-0.7)              # lockdown reuces contact down to this proportion
sd.f <- (1-0.35)             # social distancing reduces contact down to this proportion

ld.mat <- n.matT %*% (ld.f * id.mat)
sd.mat <- n.matT %*% (sd.f * id.mat)


# POPULATION SIZE
N0 = 2.6e6                 # the population of kathmandu

E0 = 66                    # Estimated number exposed on day simulatino started
Q0 = 0                     
I0 = 121                   # Estimated number infected on day simulation started
S0 = 2.6e6 - E0 - I0     # subtract exposed and infected from total to get susceptible
J0 = 0
H0 = 0
U0 = 0
R0 = 0
D0 = 0    
ac = 16                                           # number of age cohorts (16 5-year bins 0-75+)
np <-  c(9.7,9.7,10.4,11.2,11.1,8.7,6.8,          # age-cohort specific population proportion
        5.9,5.4,4.9,4.2,3.5,2.8,2.3,1.6,1.8)/100

# population bins based on age cohorts
N <- c(rep(N0, ac))*np
S <- c(rep(S0, ac))*np
E <- c(rep(E0, ac))*np
Q <- c(rep(Q0, ac))*np
I <- c(rep(I0, ac))*np
J <- c(rep(J0, ac))*np
R <- c(rep(R0, ac))*np
H <- c(rep(H0, ac))*np
U <- c(rep(U0, ac))*np
D <- c(rep(D0, ac))*np

#time variable
t <- 365               # time (days) to simulate
dt <- 1                # time step to simualte
times <- seq(0, t, by = dt)

# time controls for lockdown, starting at ld.on and lasting for ld.dur more days
ld.on <- 1
ld.dur <- 30

# time controls for social distancing, starting at sd.on and  lasting for s.dur more days
sd.on <- 1
sd.dur <- 364

# PARAMETERS
Ro    = 2.4              #reproduction number
e.dur = 3.0              # epidemic duration
i.dur = 12.5             #infectius duration

# probaility of infection transmission per contact
## i.prob.avg <- Ro / (avg.cont*i.dur) 


# probaility of infection transmission per contact (poisson distributed)

i.prob.poi <- 1 - (exp)(-(Ro/i.dur)*(avg.dur.mat))  # beta by age cohort

#mean(apply(beta.nl.mat, 1, sum))/mean(apply(b.poiss, 1, sum))

## flow parameters
beta.nl.mat <- i.prob.poi*n.matT      # beta matrix at normal circmstances, scalar multiplication
beta.ld.mat <- i.prob.poi*ld.mat      # beta matrix during lockdown
beta.sd.mat <- i.prob.poi*sd.mat      # beaa matrix during social distancing
sigma <- 1 / e.dur               # exposed to infectious transition rate
gamma <- 1 / i.dur               # recovery rate, recovery includes death and healing
asy.f <- 0.25                   # factor by which asymptomatic individuals are infectious
iso.p <- 0                     # proportion isolated and quarantined
iso.inf <- 0.15                # calc as avg.contIsolation/avg.cont; double this for quarantined people 

#fatality rates

ifr <- c(rep(0.0016, 2), rep(0.00695, 2), rep(0.0309, 2), rep(0.0844, 2),
         rep(0.161, 2), rep(0.595, 2), rep(1.93, 2), 4.28, 7.80)/100      # age specific infn fatality ratio
                                                                          # From Verity et al; Lancet infectious diseases

 
#health service specific parameters
#age-specific hospitalization rates
hos.pac <- c(rep(0, 2), rep(0.0408, 2), rep(1.04, 2), rep(3.43, 2), 
             rep(4.25, 2), rep(8.16,2),  rep(11.8,2), rep(16.6,2))/100    # age cohort based hospialization rate
                                                                          # From Verity et al; Lancet infectious diseases
hos.pac <- hos.pac
# hos.p <- 0.028 #sum(hos.pac*np) # oveall proportion of  infectious people who require hospitalization,
                         # calculated from above
           
icu.p <- 0.3           # percentage of hospitalized patients that require ICU 
dis.r1 <- 1/10         # discharge rate for patients admitted to general ward (LOS 10 days)
dis.r2 <- 1/6          # discharge rate for patients admitted to ICU (LoS 6 days)
g.cap <- 2500          # total bed capacity, in ICU and General ward (for COVID patients)
i.cap <- 200           # total bed capacity, in ICU and General ward (for COVID patients)
dg.p <- ifr[16]        # mortality rate among those who require hospitalization but are not able to get it
di.p <- 0.75           # dg.p roughly estimated assuming to qual th highst mortality of any ag cohort
                       # and 75% among those req icu bed
                       # 

# Creating a list of beta matices 

# baseline scenario
beta.list <- list(length=366)
 for (i in 1:366) {
  beta.list[[i]] <-  beta.nl.mat 
 }

# under lockdown starting at ld.on and ending at ld.off
ld.blf <- function(ld.on, ld.dur) {
  for (i in ld.on:(ld.on+ld.dur)) {
   beta.list[[i]] <-  beta.ld.mat
  }
  return(beta.list)
} 

# under social distancing starting at sd.on and ending at sd.off
sd.blf <- function(sd.on, sd.dur) {
  for (i in sd.on:(sd.on+sd.dur)) {
    beta.list[[i]] <-  beta.sd.mat 
  }
  return(beta.list)
} 

beta.nlist <- beta.list               # list of beta matrices under baseline scenario
beta.llist <- ld.blf(ld.on, ld.dur)   # list of beta matrices during lockdown
beta.slist <- sd.blf(sd.on, sd.dur)   # list of beta matrices during social distancing

# State and parameter values to pass on to LSODA 
state.val <- c(S=S, E=E, Q=Q, I=I, J=J, R=R, H=H, U=U, D=D) 

# parameter values; change the value of beta to test various scenarios
param <- list(beta=beta.nlist, sigma=sigma, gamma=gamma, ifr = ifr,
              iso.p = iso.p, asy.f=asy.f, iso.inf = iso.inf,
              hos.pac = hos.pac, icu.p = icu.p, dis.r1 = dis.r1, dis.r2 = dis.r2,
              g.cap = g.cap, i.cap = i.cap, dg.p = dg.p, di.p = di.p)                                          #parameters 
  
# Transmission model

K.model <- function(times, Y, param){
  
  dY <- numeric(length(Y))
  
  with(param,{

    for(i in 1:ac){ 
      
      dY[i] <-   -(beta[[times+1]][i,]%*%((Y[3*ac + seq(1:ac)] + 
                              (Y[ac + seq(1:ac)]*asy.f) +
                              (Y[2*ac + seq(1:ac)]*asy.f*iso.inf*2) + 
                              (Y[4*ac+seq(1:ac)]*iso.inf))/N)) *  # N is a vector representing Nj compartments
                              Y[i]                                # #Susceptible comp; infectious people contributing to transmission 
                                                                  # includes infectious, exposed, quarantined and isolated individuals
      dY[1*ac+i] <-  (beta[[times+1]][i,]%*%((Y[3*ac + seq(1:ac)] + 
                              (Y[ac + seq(1:ac)]*asy.f) +
                              (Y[2*ac + seq(1:ac)]*asy.f*iso.inf*2) + 
                              (Y[4*ac+seq(1:ac)]*iso.inf))/N)) * 
                              Y[i] - 
                              (sigma + iso.p)*Y[1*ac + i]         #Exposed
                                     
      dY[2*ac+i] <- iso.p * Y[1*ac + i] - sigma * Y[2*ac+i]                              #Quarantined
      dY[3*ac+i] <- sigma * Y[1*ac+i]  - (gamma + iso.p)* Y[3*ac + i]                    #Infectious
      dY[4*ac+i] <- iso.p * Y[3*ac+i] + sigma * Y[2*ac+i] - gamma * Y[4*ac+i]            #Isolated  # iso.p is the prop of people in isolaton / quarantine
      dY[5*ac+i] <- gamma * (1-ifr[i]) * (Y[3*ac+i] + Y[4*ac+i])                         #Recovered
      dY[6*ac+i] <- hos.pac[i] * sigma * (Y[1*ac+i] + Y[2*ac+i]) - dis.r1 * Y[6*ac+i]         # Hospitalized
      dY[7*ac+i] <- icu.p * hos.pac[i] * sigma * (Y[1*ac+i] + Y[2*ac+i]) - dis.r2 * Y[7*ac+i] # IC
      dY[8*ac+i] <- ifelse(((sum(Y[6*ac]) - 2700) > 0),
                        (Y[6*ac+i]) * ifr[i] + gamma * ifr[i] * (Y[3*ac+i] + Y[4*ac+i]),
                           gamma * ifr[i] * (Y[3*ac+i] + Y[4*ac+i]) ) 
      
                    
    }
    list(dY)
    
    })
}
                      
#Solving diff equations
valS = lsoda(state.val, times, K.model, param) 



# Compartment tables
valC <- cbind(valS[,1],
              apply(valS[,2:17], 1, sum),
              apply(valS[,18:33], 1, sum),
              apply(valS[,34:49], 1, sum),
              apply(valS[,50:65], 1, sum),
              apply(valS[,66:81], 1, sum),
              apply(valS[,82:97], 1, sum),
              apply(valS[,98:113], 1, sum),
              apply(valS[,114:129], 1, sum),
              apply(valS[,130:145], 1, sum))
              
             
              
colnames(valC) <- c("time","S", "E","Q","I","J","R","H", "U", "D")

tail(valC)
     

## data tables with poisson distributed beta

# Output tables saved  
# age specific ifr and hos.pac, baseline capacity
# valS.b252 <- valS  # output data by age cohort, 
# valC.b252 <- valC  # output data consolidated 


# valC.b252T <- valS  # Considering excess deaths due to lack of healthcare access


#Table of final compartment sizes
round(tail(valC))

sum(valC[,11])


# Day with maximum infectious
which.max(valC.b252[,"I"])


## Calculating hospital gen. bed need

# Total hospital bed-days
sum(valC.b252[,"H"])

# Maximum hospital bed requirement in a day
max(valC.b252[,"H"])

# Day with maximum hospital bed-requirement
which.max(valC.b252[,"H"])

## Calculating ICU bed need

# Total ICU bed-days
sum(valC.b252[,"U"])

# Maximum ICU bed requirement
max(valC.b252[,"U"])

# Day with maximum ICU bed-requirement
which.max(valC.b252[,"U"])


## Plotting functions

round(head(valC.b252,15L))

# TempA
 
  plot(valCp.b252.h1[,"H"], type = "l", col = "blue", ylim = c(0,2.5e4), xaxt = "n", lwd = 2.0,
       #xlim = c(0, 400),
       ylab = "Number of Individuals", 
       xlab = "Days", 
       main = "Corona virus transmission dynamics for Kathmandu\n Total Infected",
       sub = "R0 = 2.3, IncuP = 3.0, InfDur = 12.5, N = 2.6MM, Initial Exp. = 66, Initial Inf. = 121, InfFR = 0.28%", 
       cex.sub = 0.8, col.sub = "red") 
  axis(1, at = 50*(0:9))
  lines(valCp.b252.h2[,"H"], col = "red2", lwd = 2)
  lines(valCp.b252.h1[,"D"], col = "purple", lwd = 2)  
  #lines(valCp.b252.h2[,"D"], col = "orange3", lwd = 2)
  lines(valC.sdcf10252[,"I"], col = "green4", lwd = 2.5)
  abline(h = 2700, col = "maroon", lwd = 1, lty = 4)
  #text(40, 5000, "SYSTEM CAPACITY \n(2700 BEDS/day)", col = "maroon", cex = 0.9)
  legend(190, 9.0e5, legend = c("Base Case", "Lock down", "Social Distancing", "SD + Test/\nCase Finding(10%)"),
         col = c("red2", "blue4","orange3","green4"), lwd =2, lty=1, cex = 0.8)#, #"blue","green3"), 
  
  # TempB
  
  plot(valC[,"U"], type = "l", col = 1, ylim = c(0,3.2e3), xaxt = "n", lwd = 2.0,
       #xlim = c(0, 400),
       ylab = "Number of patients", 
       xlab = "Days", 
       main = "ICU Hospitalization Need versus System Capacity",
       sub = "R0 = 2.3, IncuP = 3.0, InfDur = 12.5, N = 2.6MM, Initial Exp. = 66, Initial Inf. = 121, InfFR = 0.28%", 
       cex.sub = 0.8, col.sub = "red") 
  axis(1, at = 50*(0:9))
  lines(valC.ld252[,"U"], col = 2, lwd = 2)
  lines(valC.sd252[,"U"], col = "green4", lwd = 2)  
  lines(valC.sdcf10252[,"U"], col = "blue3", lwd = 1.5)
  #lines(t, sol2[,"TtDead"], col = "black", lwd = 1.5)
  abline(h = 200, col = "maroon", lwd = 2, lty = 4)
  text(40, 450, "SYSTEM CAPACITY \n(200 ICU BEDS/day)", col = "maroon", cex = 0.9)
  legend(200, 2.50e3, legend = c("Base Case", "Lock down", "Social Distancing", "SD + Test(10%)"),
         col = c(1:2, "green4","blue3"), lty=1, cex = 0.8, lwd = 2)
  
  
  # TempB
  
  plot(valC.b252[,"S"], type = "l", col = 1, ylim = c(0,2.6e6), xaxt = "n", lwd = 2.0,
       #xlim = c(0, 400),
       ylab = "Number of individuals", 
       xlab = "Days", 
       main = "COVID-19 Transmission Dynamics in Kathmandu",
       sub = "R0 = 2.3, IncuP = 3.0, InfDur = 12.5, N = 2.6MM, Initial Exp. = 66, Initial Inf. = 121, InfFR = 0.28%", 
       cex.sub = 0.8, col.sub = "maroon") 
  axis(1, at = 50*(0:9))
  lines(valC.b252[,"E"], col = "orange", lwd = 2)
  lines(valC.b252[,"I"], col = "red3", lwd = 2)  
  lines(valC.b252[,"R"], col = "blue3", lwd = 1.5)
  lines(valC.b252[,"D"], col = "black", lwd = 1.5)
  abline(h = 200, col = "maroon", lwd = 2, lty = 4)
  #text(40, 450, "SYSTEM CAPACITY \n(200 ICU BEDS/day)", col = "maroon", cex = 0.9)
  legend(180, 2.00e6, legend = c("Susceptible", "Exposed", "Infectious", "Recovered", "Dead"),
         col = c("green4", "orange","red3", "blue3", "black"), lty=1, cex = 0.8, lwd = 2)
  