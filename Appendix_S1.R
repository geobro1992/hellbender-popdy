library(reshape2)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyverse)

#########################
# mu of survival function
#########################

#max age across all groups in years
Amax <- c(30, 30) # Taber et al. 1975, Peterson et al. 1983, Browne et al. 2014

# natural mortality (should be size/age specific) - #Taber (1975), Peterson (1983, 1985) and Olson (2013) provided survivorship curves
#natmort = c(0.2, 0.2) 
# ML = Mu * (L ^ b)   
Mu = c(200, 300, 400, 500)
b = -1.38

#mm size at "recruitment"
initialsize = 100 # Taber 1975, Peterson et al. 1983  (when they are ~ 1.5 years old)
# hatchlings are 25mm (Smith 1907)

# VBGE parameters
Linf = 570  
# Taber 1975 - 540mm
# Wiggs 1977 - 600mm
# Peterson et al. 1983 - 470 OZARK???

k = 0.13 # still need estimates
# Smith 1907
# Bodinof 2012
# Unger 2013 - size dependent recruitment?

# Storage vectors
L <- array(dim=c(Amax[1]+1), data=NA)    #length-at-age
ML <- array(dim=c(Amax[1], length(Mu)), data=0)      #mortality-at-age/length

L[1] = initialsize #initial size when entering the population model 

# calculate growth and survival
for (a in 1:(Amax[1])) { # age-specific growth, maturity, and selectivity functions
    
    # VB GROWTH EQUATION		 
  L[a+1]=Linf*(1-exp(-k)) + (L[a])*exp(-k) #length at age
    
  } #a 


for(i in 1:length(Mu)){
  for (a in 1:(Amax[1])) {
    ML[a,i] = Mu[i] * (L[a] ^ b)      
  }
}

ML = gather(as.data.frame(ML))

ML$key = as.factor(ML$key)
levels(ML$key) = c("Mu = 200", "Mu = 300", "Mu = 400", "Mu = 500")
ML$a = rep(1:30, 4)

ggplot(ML, aes(x = a, y = 1-value, color = key)) +
  geom_line() +
  ylim(0,1) + xlab("Age") + ylab("Survival") +
  theme_classic()



#######################
# removing DD processes
#######################

nestbox = "YES"

#-------------------------------
#Population model parameters

Npop <- 1000               #initial pop size (not important)
Tmax <- 200                 #projection years
nsims <- 1                # number of simulations
Tcatastrophe = 100              # time nest failure occurs 

# cannibalism / recruitment failure
scenarios = seq(from = 0, to = 1, by = 0.05)  #range of probability of nest failure
Nscenarios = length(scenarios)

#---------------------------------
# species data

Ngroup <- 2                 #Number of groups - IE females,males
prop=c(0.5, 0.5)            #initial proportions females and males
# Peterson 1988 sex ratio 1 : 1
# Humphries and Pauley 2005 sex ratio of 1.2 : 1
# Hillis and Bellis 1971 sex ratio of 1.6 : 1
# Foster et al. 2009 sex ratio of 1.8 : 1
# Burgmeier et al. 2011 sex ratio of 2.6 : 1
# Smith 1907 sex ratio of 8 : 1


# Makowsky et al 2010 found no sexual dimorphism
#but Peterson et al (1983, 1988) and Humphries and Pauley (2005) found differences in mass:length relationships


#max age across all groups in years
Amax <- c(30, 30) # Taber et al. 1975, Peterson et al. 1983, Browne et al. 2014

# natural mortality (should be size/age specific) - #Taber (1975), Peterson (1983, 1985) and Olson (2013) provided survivorship curves
#natmort = c(0.2, 0.2) 
# ML = Mu * (L ^ b)   
Mu = 500
b = -1.38

####################
## Growth parameters
####################

#mm size at "recruitment"
initialsize = c(100, 100) # Taber 1975, Peterson et al. 1983  (when they are ~ 1.5 years old)
# hatchlings are 25mm (Smith 1907)

# VBGE parameters
Linf = c(570, 570)  
# Taber 1975 - 540mm
# Wiggs 1977 - 600mm
# Peterson et al. 1983 - 470 OZARK???

k = c(0.13, 0.13) # still need estimates
# Smith 1907
# Bodinof 2012
# Unger 2013 - size dependent recruitment?


########################
# age-dependent maturity
########################

#Nickerson and Mays 1973
#Taber et al. 1975 (males 5 years at ~ 300mm; females 7-8 years at ~ 380mm)
# Topping and Ingersol 1981
# Peterson et al. 1983 (males 4-5 years at 250-300mm; females 7-8 years at 330-380mm)

matbeta=c(7, 5) # mean age at maturity
matalpha=c(2, 2) # steepness of maturity curve


#################
# storage vectors 
#################

pmat <- array(dim=c(Amax[1],Ngroup), data=NA)   #probability of maturation
L <- array(dim=c(Amax[1]+1,Ngroup), data=NA)    #length-at-age

ML <- array(dim=c(Amax[1]), data=0)      #mortality-at-age/length

N.F <- array(dim=c(Amax[1],Tmax,Nscenarios), data=0)  #this matrix holds population numbers in each year class over time
N.M <- array(dim=c(Amax[2],Tmax,Nscenarios), data=0)  #this matrix holds population numbers in each year class over time

P <- array(dim=c(Tmax, Nscenarios), data=0)  #holds number of larvae at each time

###############################################################
# growth, maturation, and selectivity - SIZE SPECIFIC MORTALITY COULD GO HERE
###############################################################

for(g in 1:Ngroup) {  
  
  L[1,g] = initialsize[g] #initial size when entering the population model 
  
  
  for (a in 1:(Amax[1])) { # age-specific growth, maturity, and selectivity functions
    
    # VB GROWTH EQUATION		 
    L[a+1,g]=Linf[g]*(1-exp(-k[g])) + (L[a,g])*exp(-k[g]) #length at age
    
    # AGE DEPENDENT MATURATION
    pmat[a,g] <- 1/(1+exp(-matalpha[g]*(a-matbeta[g])))
    
  } #a 
} #g

# make sure animals don't mature too early
pmat[which(pmat < 0.15)] = 0

#######################
## Fec / Rec parameters
#######################

#estimate of length-specific egg production
eggs <- 3.025*L[, 1] - 1140.89  # Topping and Ingersol 1981
eggs[which(eggs < 0)] = 0
eggs = eggs[-1]

#assuming Beverton Holt recruitment, these are general parameters I adjusted by eye
# increase beta to reduce carrying capacity (stochastic beta)
alpha = 0.5 
beta = 0.001

# Nickerson and Mays (1973) estimated density to be one hellbender per 6-7 m2 to one per 13-16 m2
# equating to 428 (CI: 341-573) hellbenders per km stretch of river 
# only 50-60% of the riffle was considered suitable habitat

# Wiggs and Barke in Peterson et al. 1983, estimated 506 (CI: 375-677) hellbenders per km of river

# Peterson et al. 1983 estimated ~5 hellbenders per 100m2


# Humphries and Pauley 2005 estimated ~1 hellbender per 100m2
# Burgmeier et al. 2011 estimated 0.06 hellbenders per 100m2
# see Browne et al. 2014 and Burgmeier et al. 2011 for summary table of densities


#####################
# Population Dynamics
#####################

#For every level of fishing pressure and cannibalism, simulate population dynamics 
#start from an arbitrary population size and let the population reach a stable age distribution, then start fishing. 
#The population will reach a new, disturbed, steady state (stable age dist). 
#*Note:  recruitment is based only on Female abundance, and assuming 50:50 offspring sex ratio.

# nest availability
if (nestbox == "YES") {
  n.nests = 100
} else  n.nests = 10

for(n in 1:nsims){
  
  for (s in 1:Nscenarios) {
    
    N.F[1,1,s] = Npop*prop[1] #initial larvae females
    N.M[1,1,s] = Npop*prop[2] #initial larvae males
    
    N.F[10,1,s] = 50 #initial adult females of age 10
    N.M[10,1,s] = 50 #initial adult males of age 10
    
    for (t in 1:(Tmax-1)) {
      
      Etemp <- rep(eggs, round_any(N.F[ , t, s]*pmat[ , 1], 1))       # temp vector for clutch sizes produced by each female at time t
      
      clutches <- if(length(Etemp) >= n.nests) {sample(Etemp, n.nests, replace = FALSE)} else{clutches <- Etemp}      # temp vector for eggs deposited in nests
      
      emergers <- clutches    # NO DENSITY DEPENDENCE
      emergers = emergers[sample(length(emergers), length(emergers) - round_any(length(emergers) * (ifelse(t < Tcatastrophe, scenarios[s], 0.99)),1))] # levels of cannibalism for each scenario, and increase in nest failure after Tcatastrophe
      
      
      P[t,s] = sum(emergers) # total larvae emerging from all nests in a given year
      
      # NO DD Beverton Holt recruitment
      N.F[1,t+1,s]= P[t,s]*prop[1] #females that emerge and recruit to the population in the next time step 
      N.M[1,t+1,s]= P[t,s]*prop[2] #males that emerge and recruit to the population in the next time step 
      
      for (a in 1:(Amax[1])) {
        # AGE/SIZE DEPENDENT MORTALITY (Lorenzen 1996)
        ML[a] = Mu * (L[a,1] ^ b)      
      }
      
      for(a in 1:(Amax[1]-1)){
        
        N.F[a+1,t+1,s] <- N.F[a,t,s]*(1-ML[a]) #number surviving through the time interval
        N.M[a+1,t+1,s] <- N.M[a,t,s]*(1-ML[a])  
        
      } #a
      
    } #t
    
  } #s
  
} #n




#######
# Plots
#######

# stable equilibria population projections
op <- par(mfrow = c(1,1), oma = c(3,2,2,2) + 0.1, cex.main = 2, mar = c(3, 4, 1, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(colSums(N.F[-(1:7),-1,1])[-(1:50)], type="l", lwd=2,   ylab="", xlab="", col='black', axes = FALSE, ylim=c(0, max(colSums(N.F[-(1:7),,1], na.rm=TRUE))))
lines(colSums(N.F[-(1:7),-1,3])[-(1:50)], lwd=2, col='blue')
# lines(colSums(N.F[-(1:7),-1,5]), lwd=2, col='grey') 
# lines(colSums(N.F[-(1:7),-1,7]), lwd=2, col='grey')
# lines(colSums(N.F[-(1:7),-1,9]), lwd=2, col='grey') 
lines(colSums(N.F[-(1:7),-1,11])[-(1:50)], lwd=2, col='grey') 
# lines(colSums(N.F[-(1:7),-1,13]), lwd=2, col='grey') 
# lines(colSums(N.F[-(1:7),-1,15]), lwd=2, col='grey') 
lines(colSums(N.F[-(1:7),-1,17])[-(1:50)], lwd=2, col='red') 
# lines(colSums(N.F[-(1:7),-1,19]), lwd=2, col='grey') 
# lines(colSums(N.F[-(1:7),-1,21]), lwd=2, col='grey') 

axis(1)
axis(2)
par(las = 0)
mtext("Time", side = 1, line = 2.5, cex = 1.5)
mtext("Abundance", side = 2, line = 3.7, cex = 1.5)


##plot age structure over time
op <- par(mfrow = c(4,1), oma = c(3,2,2,2) + 0.1, cex.main = 2, mar = c(3, 4, 1, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 1.3, bty = "n", las = 1)

barplot(N.F[-(1:5), 99, 11], names.arg=c(6:30), axes = FALSE, xlab = "", ylab= "", ylim = c(0, max(N.F[-(1:5), 99, 11])),   ann=FALSE)
axis(2, cex.axis = 1.5, las = 1)
par(las = 0)

barplot(N.F[-(1:5), 110, 11], names.arg=c(6:30), axes = FALSE, xlab = "", ylab= "", ylim = c(0, max(N.F[-(1:5), 99, 11])),  ann=FALSE)
axis(2, cex.axis = 1.5, las = 1)
par(las = 0)

barplot(N.F[-(1:5), 115, 11],names.arg=c(6:30), axes = FALSE, xlab = "", ylab= "", ylim = c(0, max(N.F[-(1:5), 99, 11])),   ann=FALSE)
axis(2, cex.axis = 1.5, las = 1)
par(las = 0)

barplot(N.F[-(1:5), 120, 11], names.arg=c(6:30),axes = FALSE, xlab = "", ylab= "", ylim = c(0, max(N.F[-(1:5), 99, 11])),    ann=FALSE)
axis(2, cex.axis = 1.5, las = 1)
par(las = 0)


title(xlab = "Age (years)", cex.lab = 3, font.lab= 2,
      outer = TRUE, line = 0.5)
title(ylab = "Abundance", cex.lab = 3, font.lab= 2,
      outer = TRUE, line = -0.5)




#######################
# constant survival
#######################

library(reshape2)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyr)

harvest = "NO"
#-------------------------------
#Population model parameters

Npop <- 1000               #initial pop size (not important)
Tmax <- 200                 #projection years
nsims <- 10000                # number of simulations
Tcatastrophe = 500              # time nest failure occurs 
nests = 30

# scenarios
beta = seq(0.0005, 0.005, length.out = 20) # nest availability
failures = seq(from = 0, to = 1, by = 0.05)  #range of probability of nest failure
scenarios = expand.grid(beta, failures)
Nscenarios = length(scenarios[,1])

#---------------------------------
# species data

Ngroup <- 2                 #Number of groups - IE females,males
prop=c(0.5, 0.5)            #initial proportions females and males
# Peterson 1988 sex ratio 1 : 1
# Humphries and Pauley 2005 sex ratio of 1.2 : 1
# Hillis and Bellis 1971 sex ratio of 1.6 : 1
# Foster et al. 2009 sex ratio of 1.8 : 1
# Burgmeier et al. 2011 sex ratio of 2.6 : 1
# Smith 1907 sex ratio of 8 : 1


# Makowsky et al 2010 found no sexual dimorphism
#but Peterson et al (1983, 1988) and Humphries and Pauley (2005) found differences in mass:length relationships


#max age across all groups in years
Amax <- c(30, 30) # Taber et al. 1975, Peterson et al. 1983, Browne et al. 2014

# natural mortality (should be size/age specific) - #Taber (1975), Peterson (1983, 1985) and Olson (2013) provided survivorship curves
#natmort = c(0.2, 0.2) 
# ML = Mu * (L ^ b)   
Mu = 500
b = -1.38

####################
## Growth parameters
####################

#mm size at "recruitment"
initialsize = c(100, 100) # Taber 1975, Peterson et al. 1983  (when they are ~ 1.5 years old)
# hatchlings are 25mm (Smith 1907)

# VBGE parameters
Linf = c(570, 570)  
# Taber 1975 - 540mm
# Wiggs 1977 - 600mm
# Peterson et al. 1983 - 470 OZARK???

k = c(0.13, 0.13) # still need estimates
# Smith 1907
# Bodinof 2012
# Unger 2013 - size dependent recruitment?


########################
# age-dependent maturity
########################

#Nickerson and Mays 1973
#Taber et al. 1975 (males 5 years at ~ 300mm; females 7-8 years at ~ 380mm)
# Topping and Ingersol 1981
# Peterson et al. 1983 (males 4-5 years at 250-300mm; females 7-8 years at 330-380mm)

matbeta=c(7, 5) # mean age at maturity
matalpha=c(2, 2) # steepness of maturity curve


#################
# storage vectors 
#################

pmat <- array(dim=c(Amax[1],Ngroup), data=NA)   #probability of maturation
L <- array(dim=c(Amax[1]+1,Ngroup), data=NA)    #length-at-age

ML <- array(dim=c(Amax[1]), data=0)      #mortality-at-age/length
ML.temp <- array(dim=c(Amax[1]), data=0) #temp storage for stoachastic survival

N.F <- array(dim=c(Amax[1],Tmax,Nscenarios), data=0)  #this matrix holds population numbers in each year class over time
N.M <- array(dim=c(Amax[2],Tmax,Nscenarios), data=0)  #this matrix holds population numbers in each year class over time

P <- array(dim=c(Tmax, Nscenarios), data=0)  #holds number of larvae at each time

###############################################################
# growth, maturation, and selectivity - SIZE SPECIFIC MORTALITY COULD GO HERE
###############################################################

for(g in 1:Ngroup) {  
  
  L[1,g] = initialsize[g] #initial size when entering the population model 
  
  
  for (a in 1:(Amax[1])) { # age-specific growth, maturity, and selectivity functions
    
    # VB GROWTH EQUATION		 
    L[a+1,g]=Linf[g]*(1-exp(-k[g])) + (L[a,g])*exp(-k[g]) #length at age
    
    # AGE DEPENDENT MATURATION
    pmat[a,g] <- 1/(1+exp(-matalpha[g]*(a-matbeta[g])))
    
  } #a 
} #g


#######################
## Fec / Rec parameters
#######################

#estimate of length-specific egg production
eggs <- 3.025*L[, 1] - 1140.89  # Topping and Ingersol 1981
eggs[which(eggs < 0)] = 0
eggs = eggs[-1]

#assuming Beverton Holt recruitment, these are general parameters I adjusted by eye
# beta is varying for scenarios
alpha = 0.5 

# Nickerson and Mays (1973) estimated density to be one hellbender per 6-7 m2 to one per 13-16 m2
# equating to 428 (CI: 341-573) hellbenders per km stretch of river 
# only 50-60% of the riffle was considered suitable habitat

# Wiggs and Barke in Peterson et al. 1983, estimated 506 (CI: 375-677) hellbenders per km of river

# Peterson et al. 1983 estimated ~5 hellbenders per 100m2


# Humphries and Pauley 2005 estimated ~1 hellbender per 100m2
# Burgmeier et al. 2011 estimated 0.06 hellbenders per 100m2
# see Browne et al. 2014 and Burgmeier et al. 2011 for summary table of densities


###################
## Stochastic survival
#####################

ML = c(rep(0.5, 6), rep(0.15, 24))      


if(harvest == "YES"){
  
  ML.upr = vector()
  ML.upr[which(ML < 0.15)] = ML[which(ML < 0.15)]
  ML.upr[which(ML > 0.15)] = ML[which(ML > 0.15)]-0.13
  
  ML.lwr = vector()
  ML.lwr[which(ML < 0.15)] = ML[which(ML < 0.15)]+0.13
  ML.lwr[which(ML > 0.15)] = ML[which(ML > 0.15)]+0.13
  
  
} else{
  
  ML.upr = vector()
  ML.upr[which(ML < 0.2)] = ML[which(ML < 0.2)]
  ML.upr[which(ML > 0.2)] = ML[which(ML > 0.2)]-0.13
  
  ML.lwr = vector()
  ML.lwr[which(ML < 0.2)] = ML[which(ML < 0.2)]
  ML.lwr[which(ML > 0.2)] = ML[which(ML > 0.2)]+0.13
  
}


#####################
# Population Dynamics
#####################

#For every level of fishing pressure and cannibalism, simulate population dynamics 
#start from an arbitrary population size and let the population reach a stable age distribution, then start fishing. 
#The population will reach a new, disturbed, steady state (stable age dist). 
#*Note:  recruitment is based only on Female abundance, and assuming 50:50 offspring sex ratio.


n.ex=rep(0, Nscenarios)          # number of sims going extinct in timeframe

for(n in 1:nsims){
  
  for (s in 1:Nscenarios) {
    
    N.F[1,1,s] = Npop*prop[1] #initial larvae females
    N.M[1,1,s] = Npop*prop[2] #initial larvae males
    
    N.F[10,1,s] = 50 #initial adult females of age 10
    N.M[10,1,s] = 50 #initial adult males of age 10
    
    for (t in 1:(Tmax-1)) {
      
      Etemp <- rep(eggs, round_any(N.F[ , t, s]*pmat[ , 1], 1))       # temp vector for clutch sizes produced by each female at time t
      
      clutches <- if(length(Etemp) >= nests) {sample(Etemp, nests, replace = FALSE)} else{clutches <- Etemp}      # temp vector for eggs deposited in nests
      
      emergers <- 200 / (1 + exp(-0.01 * (clutches - 200)))    # density dependent egg survival
      emergers = emergers[sample(length(emergers), length(emergers) - round_any(length(emergers) * scenarios[s,2] ,1))] # levels of cannibalism for each scenario, and increase in nest failure after Tcatastrophe
      
      
      P[t,s] = sum(emergers) # total larvae emerging from all nests in a given year
      
      #Beverton Holt recruitment after cannibalism
      N.F[1,t+1,s]= (alpha*P[t,s]/(1+scenarios[s,1]*P[t,s])*prop[1]) #females that emerge and recruit to the population in the next time step 
      N.M[1,t+1,s]= (alpha*P[t,s]/(1+scenarios[s,1]*P[t,s])*prop[2]) #males that emerge and recruit to the population in the next time step 
      
      for (a in 1:(Amax[1])) {
        # AGE/SIZE DEPENDENT MORTALITY (Lorenzen 1996)
        ML.temp[a] = runif(1, ML.upr[a], ML.lwr[a])      
      }
      
      
      for(a in 1:(Amax[1]-1)){
        
        N.F[a+1,t+1,s] <- N.F[a,t,s]*(1-ML.temp[a]) #number surviving through the time interval
        N.M[a+1,t+1,s] <- N.M[a,t,s]*(1-ML.temp[a])  
        
      } #a
      
    } #t
    
    if(any(colSums(N.F[,,s]*pmat[,1]) < 10)){
      n.ex[s] = n.ex[s]+1
    } else{}
    
  } #s
  
} #n

# calculate probability of extinction for each year
p.ex = n.ex/nsims

# save
save(p.ex, file = "pex_deterministic_supp_surv.RData")   # when adult mortality is stochastic (harvest = YES)
#save(p.ex, file = "pex_deterministic_allBH.RData") # when adult mortality is constant (harvest = NO)



#################################
# plotting deterministic contour

load("pex_deterministic_allBH.RData")
p.ex.high = p.ex

z.high = as.data.frame(cbind(scenarios, p.ex.high))
names(z.high) = c("beta", "failure", "p.ex")

z = z.high %>% 
  spread(beta, p.ex)

z = as.matrix(z[,-1])

#########################################
# calculate mean of low and high matrices

x = failures*100
y = beta

op <- par(mfrow = c(1,1), oma = c(3,3,2,2) + 0.1, cex.main = 2, mar = c(5, 8, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 2, bty = "n", las = 1)

# deterministic adult mortality
filled.contour(x = x, y = y, z=z, 
               main = "Deterministic Adult Mortality",
               ylim = c(0.0005,0.004), xlab = "Nest Failure %", ylab = "Strength of \nDD Recruitment", cex.lab = 2,
               plot.axes = {
                 axis(1)
                 axis(2, at = seq(0,0.004, by = 0.001), labels = FALSE, tck = -0.001)
                 contour(x = x, y = y, z=z, 
                         levels = seq(0.1,0.9, by = 0.2), lwd = 2, labcex = 1.2, method = "flattest",
                         add = TRUE)
               }
)




#################################
# plotting stochastic contours

#################
# range of BH parameters
load("pex_stochasticL_allBH.RData")
p.ex.high = p.ex

z.high = as.data.frame(cbind(scenarios, p.ex.high))
names(z.high) = c("beta", "failure", "p.ex")

z2 = z.high %>% 
  spread(beta, p.ex)

z2 = as.matrix(z2[,-1])

#########################################
# calculate mean of low and high matrices

x = failures*100
y = beta

op <- par(mfrow = c(1,1), oma = c(3,3,2,2) + 0.1, cex.main = 2, mar = c(5, 8, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 2, bty = "n", las = 1)

# deterministic adult mortality
filled.contour(x = x, y = y, z=z2, 
               main = "Stochastic Adult Mortality",
               ylim = c(0.0005,0.004), xlab = "Nest Failure %", ylab = "Strength of \nDD Recruitment", cex.lab = 2,
               plot.axes = {
                 axis(1)
                 axis(2, at = seq(0,0.004, by = 0.001), labels = FALSE, tck = -0.001)
                 contour(x = x, y = y, z=z2, 
                         levels = seq(0.1,0.9, by = 0.2), lwd = 2, labcex = 1.2, method = "flattest",
                         add = TRUE)
               }
)


