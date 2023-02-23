library(reshape2)
library(ggplot2)
library(ggpubr)
library(plyr)
library(tidyr)

harvest = "YES"
#-------------------------------
#Population model parameters

Npop <- 1000               #initial pop size (not important)
Tmax <- 200                 #projection years
nsims <- 10000                # number of simulations
Tcatastrophe = 500              # time nest failure occurs 

# scenarios
nests = seq(0, 60, by = 5) # nest availability
failures = seq(from = 0, to = 1, by = 0.05)  #range of probability of nest failure
scenarios = expand.grid(nests, failures)
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
# increase beta to reduce carrying capacity (stochastic beta) 
# highBH - beta = 0.001 
# lowBH - beta = 0.005
alpha = 0.5 
beta = 0.005

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
for (a in 1:(Amax[1])) {
  # AGE/SIZE DEPENDENT MORTALITY (Lorenzen 1996)
  ML[a] = Mu * (L[a,1] ^ b)      
}


if(harvest == "YES"){
  
  ML.upr = vector()
  ML.upr[which(ML < 0.15)] = ML[which(ML < 0.15)]
  ML.upr[which(ML > 0.15)] = ML[which(ML > 0.15)]-0.13
  
  ML.lwr = vector()
  ML.lwr[which(ML < 0.15)] = ML[which(ML < 0.15)]+0.13
  ML.lwr[which(ML > 0.15)] = ML[which(ML > 0.15)]+0.13
  
  
} else{
  
  ML.upr = vector()
  ML.upr[which(ML < 0.15)] = ML[which(ML < 0.15)]
  ML.upr[which(ML > 0.15)] = ML[which(ML > 0.15)]-0.13
  
  ML.lwr = vector()
  ML.lwr[which(ML < 0.15)] = ML[which(ML < 0.15)]
  ML.lwr[which(ML > 0.15)] = ML[which(ML > 0.15)]+0.13
  
}


#####################
# Population Dynamics
#####################

n.ex=rep(0, Nscenarios)          # number of sims going extinct in timeframe

for(n in 1:nsims){
  
  for (s in 1:Nscenarios) {
    
    N.F[1,1,s] = Npop*prop[1] #initial larvae females
    N.M[1,1,s] = Npop*prop[2] #initial larvae males
    
    N.F[10,1,s] = 50 #initial adult females of age 10
    N.M[10,1,s] = 50 #initial adult males of age 10
    
    for (t in 1:(Tmax-1)) {
      
      Etemp <- rep(eggs, round_any(N.F[ , t, s]*pmat[ , 1], 1))       # temp vector for clutch sizes produced by each female at time t
      
      clutches <- if(length(Etemp) >= scenarios[s,1]) {sample(Etemp, scenarios[s,1], replace = FALSE)} else{clutches <- Etemp}      # temp vector for eggs deposited in nests
      
      emergers <- 200 / (1 + exp(-0.01 * (clutches - 200)))    # density dependent egg survival
      emergers = emergers[sample(length(emergers), length(emergers) - round_any(length(emergers) * scenarios[s,2] ,1))] # levels of cannibalism for each scenario, and increase in nest failure after Tcatastrophe
      
      
      P[t,s] = sum(emergers) # total larvae emerging from all nests in a given year
      
      #Beverton Holt recruitment after cannibalism
      N.F[1,t+1,s]= (alpha*P[t,s]/(1+beta*P[t,s])*prop[1]) #females that emerge and recruit to the population in the next time step 
      N.M[1,t+1,s]= (alpha*P[t,s]/(1+beta*P[t,s])*prop[2]) #males that emerge and recruit to the population in the next time step 
      
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
save(p.ex, file = "pex_all_stochasticL_lowBH.RData") # when adult mortality is stochastic (harvest = YES), and density dependence is strong (beta = 0.005 in beverton holt recruitment function)
#save(p.ex, file = "pex_all_stochasticL_highBH.RData") # when adult mortality is stochastic (harvest = YES), and density dependence is weak (beta = 0.001 in beverton holt recruitment function)
#save(p.ex, file = "pex_all_deterministic_lowBH.RData") # when adult mortality is constant (harvest = NO), and density dependence is strong (beta = 0.005 in beverton holt recruitment function)
#save(p.ex, file = "pex_all_deterministic_highBH.RData") # when adult mortality is constant (harvest = NO), and density dependence is weak (beta = 0.001 in beverton holt recruitment function)


#################################
# plotting deterministic contours

#################
# high BH function
load("pex_all_deterministic_highBH.RData")
p.ex.high = p.ex

z.high = as.data.frame(cbind(scenarios, p.ex.high))
names(z.high) = c("nests", "failure", "p.ex")

z1 = z.high %>% 
  spread(nests, p.ex)

z1 = as.matrix(z1[,-1])

#################
# low BH function
load("pex_all_deterministic_lowBH.RData")
p.ex.low = p.ex

z.low = as.data.frame(cbind(scenarios, p.ex.low))
names(z.low) = c("nests", "failure", "p.ex")

z2 = z.low %>% 
  spread(nests, p.ex)

z2 = as.matrix(z2[,-1])

#########################################
# calculate mean of low and high matrices
z.list <- list(z1, z2)
z.array <- do.call(cbind, z.list)
z.array <- array(z.array, dim=c(dim(z.list[[1]]), length(z.list)))

z = apply(z.array, c(1, 2), mean, na.rm = TRUE)


x = failures*100
y = nests

op <- par(mfrow = c(1,1), oma = c(3,3,2,2) + 0.1, cex.main = 2, mar = c(5, 6, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 2, bty = "n", las = 1)

# deterministic adult mortality
filled.contour(x = x, y = y, z=z, 
               main = "Deterministic Adult Mortality",
               ylim = c(0,60), xlab = "Nest Failure %", ylab = "Nest Availability", cex.lab = 2,
               plot.axes = {
                 axis(1)
                 axis(2, at = seq(0,60, by = 10))
                 contour(x = x, y = y, z=z, 
                        levels = seq(0.1,0.9, by = 0.2), lwd = 2, labcex = 1.2, method = "flattest",
                         add = TRUE)
               }
)




#################################
# plotting stochastic contours

#################
# high BH function
load("pex_all_stochasticL_highBH.RData")
p.ex.highS = p.ex

z.highS = as.data.frame(cbind(scenarios, p.ex.highS))
names(z.highS) = c("nests", "failure", "p.ex")

z3 = z.highS %>% 
  spread(nests, p.ex)

z3 = as.matrix(z3[,-1])

#################
# low BH function
load("pex_all_stochasticL_lowBH.RData")
p.ex.lowS = p.ex

z.lowS = as.data.frame(cbind(scenarios, p.ex.lowS))
names(z.lowS) = c("nests", "failure", "p.ex")

z4 = z.lowS %>% 
  spread(nests, p.ex)

z4 = as.matrix(z4[,-1])

#########################################
# calculate mean of low and high matrices
z.listS <- list(z3, z4)
z.arrayS <- do.call(cbind, z.listS)
z.arrayS <- array(z.arrayS, dim=c(dim(z.listS[[1]]), length(z.listS)))

zS = apply(z.arrayS, c(1, 2), mean, na.rm = TRUE)


x = failures*100
y = nests

op <- par(mfrow = c(1,1), oma = c(3,3,2,2) + 0.1, cex.main = 2, mar = c(5, 6, 3, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, font.main = 2, cex.axis = 2, bty = "n", las = 1)

# deterministic adult mortality
filled.contour(x = x, y = y, z=zS, 
               main = "Stochastic Adult Mortality",
               ylim = c(0,60), xlab = "Nest Failure %", ylab = "Nest Availability", cex.lab = 2,
               plot.axes = {
                 axis(1)
                 axis(2, at = seq(0,60, by = 10))
                 contour(x = x, y = y, z=zS, 
                          levels = seq(0.1,0.9, by = 0.2), lwd = 2, labcex = 1.2, method = "flattest",
                          add = TRUE)
               }
)
