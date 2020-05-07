rm(list = ls())

setwd("~/Dropbox/Research/Corona/C19/")


library(deSolve)
library(reshape)
library(magrittr)
library(dplyr)
library(ggplot2)
library(gdata)
library(matlib)
library(RColorBrewer)
library(purrr)
library(plot.matrix)


#CONTACT MATRIX
n_mat <- read.xls("Data/nepal_cmat.xls", header = F)

## Intervention matrix
inter.mat <- read.csv("Data/covid_inter.csv", header = TRUE, row.names = 1)

n_mat1 <- n_mat[,-1]                             
n_matA <- unname(as.matrix(n_mat1[1:16,]))       # all contact matrix
n_matH <- unname(as.matrix(n_mat1[17:32,]))      # home contact matrix
n_matS <- unname(as.matrix(n_mat1[33:48,]))      # school contact matrix
n_matW <- unname(as.matrix(n_mat1[49:64,]))      # work contact matrix
n_matO <- unname(as.matrix(n_mat1[65:80,]))      # other contact matrix

# matrix of total number of contacts
n.matT <- unname(n_matH + n_matS + n_matW + n_matO) # sum matrices and remove row and col names

# duration per contact
avg.cont <- mean(apply(n.matT, 1, sum))                    # average number of contacts
avg.dur <- 1/(apply(n.matT, 1, sum))                       # vector of contact duration by age cohort
avg.dur.mat <- matrix(rep(avg.dur,16), ncol = 16, byrow=F) # matrix of contact duration by age cohort


#INTERVENTION EFFECTS
id.mat <- diag(16)

ld.f <- (1-0.7)              # lockdown reuces contact down to this proportion
sd.f <- (1-0.35)             # social distancing reduces contact down to this proportion

ld.mat <- n.matT %*% (ld.f * id.mat)
sd.mat <- n.matT %*% (sd.f * id.mat)


# POPULATION SIZE
N0 = 2.6e6                 # the population of kathmandu

E0 = 65                    # Estimated number exposed on day simulatino started
Q0 = 0                     
I0 = 120                   # Estimated number infected on day simulation started
S0 = N0 - E0 - I0          # subtract exposed and infected from total to get susceptible
J0 = 0
H0 = 0
U0 = 0
D0 = 0
R0 = 0    
ac = 16                                           # number of age cohorts (16 5-year bins 0-75+)
np <-  c(9.7,9.7,10.4,11.2,11.1,8.7,6.8,          # age-cohort specific population proportion
        5.9,5.4,4.9,4.2,3.5,2.8,2.3,1.6,1.8)/100

# population bins based on age cohorts

S <- c(rep(S0, ac))*np
E <- c(rep(E0, ac))*np
Q <- c(rep(Q0, ac))*np
I <- c(rep(I0, ac))*np
J <- c(rep(J0, ac))*np
H <- c(rep(H0, ac))*np
U <- c(rep(U0, ac))*np
D <- c(rep(D0, ac))*np
R <- c(rep(R0, ac))*np
N <- c(rep(N0, ac))*np

#time variable
t <- 365                       # time (days) to simulate
dt <- 1                        # time step to simualte
times <- seq(0, t, by = dt)


#fatality rates

ifr <- c(rep(0.0016, 2), rep(0.00695, 2), rep(0.0309, 2), rep(0.0844, 2),
         rep(0.161, 2), rep(0.595, 2), rep(1.93, 2), 4.28, 7.80)/100      



#age-specific hospitalization rates
hos.pac <- c(rep(0.01, 2), rep(0.0408, 2), rep(1.04, 2), rep(3.43, 2), 
             rep(4.25, 2), rep(8.16,2),  rep(11.8,2), rep(16.6,2))/100 


#nsims
n.sims = 2000


asy.f <- 0.1                # factor by which asymptomatic individuals are infectious
iso.inf <- 0.15             # calc as avg.contIsolation/avg.cont; double this for quarantined people 
i.dur <- 7    
icu.p <- 0.20               # percentage of hospitalized patients that require ICU 
dis.r1 <- 1/8               # discharge rate for patients admitted to general ward (LOS 10 days)
dis.r2 <- 1/6               # discharge rate for patients admitted to ICU (LoS 6 days)
g.cap <- 5400               # total bed capacity, in ICU and General ward (for COVID patients)
hs.f <- 0                   # 0.7 for 1000 additional beds, 0.52 for 2000 additional beds baseline,    # health service capacity factor: (bed ratio china - b.ratio kathamandu)/b.ratio kathmandu)
f <- 0.8                    # proportion of deaths that happen among patients who require ICU

# State and parameter values to pass on to LSODA 
state.val <- c(S=S, E=E, Q=Q, I=I, J=J, H=H, U=U, D=D, R=R, N=N) 


## Simulations
system.time(
for (i in 1:nrow(inter.mat)) {

oput <- replicate(n = n.sims, 
                  
                  OvFun(ld.on = inter.mat[i,1],
                       ld.off = inter.mat[i,2],
                       sd.on = inter.mat[i,3],
                       sd.off = inter.mat[i,4],
                       is.on = inter.mat[i,5],
                       is.off = inter.mat[i,6],
                       iso.p = inter.mat[i,7]), 
                            simplify = F)

assign(paste("Int", i, sep = ""), oput)
}
)


## List of simultion results

IntL <- list(Int1, Int2, Int3, Int4, Int5, Int6, Int7, Int8, 
             Int9, Int10, Int11, Int12, Int13, Int14, Int15)


## Summary variables

purrr::map_dbl(Int1, ~max(.x[,"I"])) %>% summary() 

purrr::map_dbl(Int1, ~which.max(.x[,"I"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~max(.x[,"E"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~which.max(.x[,"E"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~max(.x[,"H"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~which.max(.x[,"H"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~max(.x[,"U"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~which.max(.x[,"U"])) %>% summary() 
                      
purrr::map_dbl(Int1, ~.x[366,"D"]) %>% summary() 
                      

# Variable sweep and ODE function

OvFun <- function(ld.on,
                  ld.off,
                  sd.on,
                  sd.off,
                  is.on,
                  is.off,
                  iso.p) {
                  

  Ro =  runif(1, min = 2.0, max = 2.8)
  
  e.dur = runif(1,  min = 3, max = 7)
  
  i.prob.poi <- 1 - (exp)(-(Ro/i.dur)*(avg.dur.mat))  
  
  
  beta.nl.mat <- i.prob.poi*n.matT      # beta matrix at normal circmstances, scalar multiplication
  beta.ld.mat <- i.prob.poi*ld.mat      # beta matrix during lockdown
  beta.sd.mat <- i.prob.poi*sd.mat      # beta matrix during social distancing
  
 
  
  beta.fn <- function(ld.on, ld.off, sd.on, sd.off) {
    
  beta.list <- list(length=366)  
  
  for (i in 1:366) {
    beta.list[[i]] <-  beta.nl.mat 
  }
  
  if (sd.on > 0) {
  for (i in sd.on:sd.off) {
      beta.list[[i]] <-  beta.sd.mat
  }
  }
  
  if (ld.on > 0) {
    for (i in ld.on:ld.off) {
      beta.list[[i]] <-  beta.ld.mat
    }
  } 
  return(beta.list)
  }
  
  
   beta.hlist <- beta.fn(ld.on, ld.off, sd.on, sd.off)
  
  
  is.t <- rep_len(0, 366) 
    for (i in is.on:(is.on+is.off)){
     is.t[i] <- 1
  }
  
  
  param <- list(beta = beta.hlist, 
                sigma = 1 / e.dur, 
                gamma = 1 / i.dur, 
                ifr = ifr,
                iso.p = iso.p, 
                is.t = is.t, 
                asy.f = asy.f, 
                iso.inf = iso.inf,
                hos.pac = hos.pac, 
                icu.p = icu.p, 
                dis.r1 = dis.r1, 
                dis.r2 = dis.r2,
                g.cap = g.cap, 
                hs.f = hs.f, 
                f = f)  
  
  

# ODE solver function


K.model <- function(times, Y, param){
  
  dY <- numeric(length(Y))
  
  with(param,{

    for(i in 1:ac){ 
     
      dY[i] <-    -beta[[times+1]][i,]%*%((
                                         (Y[1*ac + seq(1:ac)]*asy.f) +
                                         (Y[2*ac + seq(1:ac)]*asy.f*iso.inf*2) + 
                                         (Y[3*ac + seq(1:ac)]) + 
                                         (Y[4*ac + seq(1:ac)]*iso.inf) + 
                                         (Y[5*ac + seq(1:ac)]*iso.inf) +
                                         (Y[6*ac + seq(1:ac)]*iso.inf))/N) * 
                                         Y[i]      
                                                                      # #Susceptible comp; infectious people contributing to transmission 
                                                                      # includes infectious, exposed, quarantined and isolated individuals
      dY[1*ac+i] <-  beta[[times+1]][i,]%*%((
                                           (Y[1*ac + seq(1:ac)]*asy.f) +
                                           (Y[2*ac + seq(1:ac)]*asy.f*iso.inf*2) + 
                                           (Y[3*ac + seq(1:ac)]) + 
                                           (Y[4*ac + seq(1:ac)]*iso.inf) + 
                                           (Y[5*ac + seq(1:ac)]*iso.inf) +
                                           (Y[6*ac + seq(1:ac)]*iso.inf))/N) * 
                                           Y[i]  - 
                                           (sigma + iso.p*is.t[times + 1])*Y[1*ac + i]         #Exposed
                                  
      dY[2*ac+i] <- iso.p*is.t[times + 1] * Y[1*ac+i] - sigma * Y[2*ac+i]                                              #Quarantined
      dY[3*ac+i] <- sigma * Y[1*ac+i]  - (gamma  + iso.p*is.t[times + 1] + hos.pac[i])* Y[3*ac + i]                  #Infectious
      dY[4*ac+i] <- iso.p*is.t[times + 1] * Y[3*ac+i] + sigma * Y[2*ac+i] - (gamma + hos.pac[i])* Y[4*ac+i]             #Isolated  # iso.p is the prop of people in isolaton / quarantine                                                                                                                 #Recovered
      dY[5*ac+i] <- hos.pac[i] * (Y[4*ac+i] + Y[3*ac + i]) - (dis.r1)* Y[5*ac+i]                               # Hospitalized
      dY[6*ac+i] <- icu.p * (dis.r1)* Y[5*ac+i]  - dis.r2 * Y[6*ac+i]                                     # ICU
      
      dY[7*ac+i] <- ifelse(((sum(Y[5*ac+seq(1:ac)] + Y[6*ac+seq(1:ac)]) - g.cap) > 0),
        
                    
                    (1-icu.p) * dis.r1 * Y[5*ac+i] * (ifr[i]/hos.pac[i])*(1-f)*(1 + hs.f) + 
                    dis.r2 * Y[6*ac+i] * (ifr[i]/hos.pac[i])*f*(1 + hs.f),
                    
                    (1-icu.p) * dis.r1 * Y[5*ac+i] * (ifr[i]/hos.pac[i])*(1-f) + 
                    dis.r2 * Y[6*ac+i] * (ifr[i]/hos.pac[i])*f)
                                             
      dY[8*ac+i] <-  ifelse(((sum(Y[5*ac+seq(1:ac)] + Y[6*ac+seq(1:ac)]) - g.cap) > 0),
                            
                    - ((1-icu.p) * dis.r1 * Y[5*ac+i] * (ifr[i]/hos.pac[i])*(1-f)*hs.f + 
                    dis.r2 * Y[6*ac+i] * (ifr[i]/hos.pac[i])*f*hs.f) +
                      
                    (1-icu.p) * dis.r1 * Y[5*ac+i] * (1 - (ifr[i]/hos.pac[i])*(1-f)) +
                    dis.r2 * Y[6*ac+i] * (1 - (ifr[i]/hos.pac[i])*f) +  
                    gamma * (Y[3*ac+i] + Y[4*ac+i]),
                    
                    (1-icu.p) * dis.r1 * Y[5*ac+i] * (1 - (ifr[i]/hos.pac[i])*(1-f)) +
                    dis.r2 * Y[6*ac+i] * (1 - (ifr[i]/hos.pac[i])*f) +  
                    gamma * (Y[3*ac+i] + Y[4*ac+i]))
      
                     
        
      dY[9*ac+i] <- dY[i] + dY[1*ac+i] + dY[2*ac+i] + dY[3*ac+i] + dY[4*ac+i] + 
                       dY[5*ac+i] + dY[6*ac+i] + dY[7*ac+i] + dY[8*ac+i] 
        
      
    }
    list(dY)
    
    })
}
                      

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
              apply(valS[,130:145], 1, sum),
              apply(valS[,146:161], 1, sum))


colnames(valC) <- c("time","S", "E","Q","I","J","H", "U", "D","R", "N")

return(round(valC))

}


