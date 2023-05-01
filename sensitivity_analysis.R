rm(list = ls())
library(pse) # for sensitivity analysis
library(tidyverse) # contains ggplot2 for plotting, dplyr and tidyr for data manipulation


############### For the US ###########################

#### Sensitivity analysis varying all parameters ####

# Parameter Space Exploration with Latin Hypercubes. 

# define parameters to vary
factors <- c("etaMI", "ZPF", "ZC", "ZAMin", "ZPM", 
             "etaMinF", "etaFF",  
             "ZAMac", "omega", "ZQ",  "etaMacF")

# define how the values within the range are distributed (e.g. could be qnorm instead)
q <- c("qunif", "qunif", "qunif", "qunif", "qunif", 
       "qunif", "qunif", 
       "qunif", "qunif", "qunif", "qunif")

# define the minimum and maximum boundaries of each parameter
q.arg <- list( list(min=0.02, max=0.6), list(min=41, max=405), list(min=3, max=14), list(min=2, max=18), list(min=1, max=9),
               list(min=-0.3, max=-0.05),  list(min=-0.75, max=-0.21), 
               list(min=1.5, max=15), list(min=0.3, max=0.7), list(min=403, max=999), 
               list(min=0.17, max=0.25)) #


# Create a function and let each parameter within it vary.

oneRun <- function (etaMI, ZPF, ZC, ZAMin, ZPM, 
                    etaMinF, etaFF,  
                    ZAMac, omega, ZQ, etaMacF, gastax1, MinF1, MacF1) { 
  #       Elasticities;
  
  epsiLL = 0.2; epsiCLL = 0.35; 
  
  #       Prices;
  tL = 0.318;
  pf = 186;
  
  #       Baselines;
  fueleff0 = 24;
  gastax0 = 55;
  Mac0 = 56.8;
  Min0 = 2191.8;
  
  betaMin = etaMinF/etaFF; betaMac = etaMacF/etaFF;
  # betaMin = 1; etaMinF/etaFF
  
  tLratio = tL/(1 - tL);
  fuel0 = Min0/fueleff0;
  
  MinF0 = Min0/fuel0;
  MacF0 = Mac0/fuel0; 
  
    #       Output at baseline;
    #  MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 # without physical activity
    MECbase =  ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0;
    MEBbase = tLratio*epsiLL/(1 - tLratio*epsiLL);
    taxOPT=  MECbase/(1 + MEBbase) + tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0;
    taxOPT0 = taxOPT
    #Finding for stable outcome
    
    tmp1 = tLratio*epsiCLL*(1 - etaMI)/(-etaFF);
    
    # The code below then does the following. Let gastax be the implemented gastax and taxOPT the optimal gasoline tax given gastax. \
    # Then if taxOPT and gastax deviate, we have not found the solution. If so, set the gasoline tax midway to the previous taxOPT and gastax. \
    # For this new tax level, compute Mac, Min, fuel, and Q1, and from here taxOPT. Check taxOPT against gastax and repeat if needed. Ideally, \
    # this converges to the interior solution. 

    gastax1 = gastax0; prec = 0.0001; approach = 0.5; dev = 1; i = 0; ALT = 0;
    
    #(* bit that remains constant *)
    
    MEB = tLratio*epsiLL/(1 - tLratio*epsiLL);
    
    while(dev >= prec){
      
      gastax1 = (1 - approach)*gastax1 + approach*taxOPT;
      
      #Endogenous variables;
      
      MinF1 = MinF0*((pf + gastax1)/(pf + gastax0))^(etaMinF - etaFF);
      MacF1 = MacF0*((pf + gastax1)/(pf + gastax0))^(etaMacF - etaFF);
      
      MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1;
      # MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 # rate without physical activity
      
      taxOPT = MEC/(1 + MEB) +tLratio*(pf + gastax1)* epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1;
      
      dev = (taxOPT - gastax1)^2;
      print(dev)
    }
    
    # define what output we are interested in.
  return (taxOPT) 
}


# Run the function above by feeding it an empty dataframe my.data with 11 different parameters that can vary

modelRun <- function (my.data) {
  return(mapply(oneRun, my.data[,1], my.data[,2], my.data[,3], my.data[,4], my.data[,5], my.data[,6], my.data[,7], my.data[,8], my.data[,9], my.data[,10], my.data[,11]))
}

# run the LHS sensitivity command, giving it the function to evaluate in the form of a dataframe, 
# how many times it should run ( generates 200 latin hypercubes), what the parameters are, and the number of boostrap replicates (50).
# Advantage over Monte Carlo is significant time savings in running this over the MonteCarlo package commands (i.e. seconds vs 10+minutes).

myLHS <- LHS(modelRun, factors, 200, q, q.arg, nboot=50)


# Plot results to see where solution is most likely to lie (the steepest part of the curve)
plotecdf(myLHS, main=NULL, ylab="Probability <= X", xlab="Tax, cents/gallon", col="dark blue") # clearly optimum is likely to lie between 1000 and 2000 cents/gallon.  

# Plot which parameter influences the optimal fuel tax the most. 
plotprcc(myLHS, main=NULL, title = NULL, ylab="Partial rank correlation coefficient", col = "royalblue2")

# Plotting prettier, optional
data <- as.data.frame(myLHS[["res"]])
ggplot(data, aes(V1)) + 
  stat_ecdf(geom = "step", col = "dark blue") + 
  geom_vline(xintercept=55, linetype="dashed", color = "red", size=0.5) +
  geom_hline(yintercept=0.001, linetype="dashed", color = "red", size=0.5) +
  theme_minimal() +
  labs(x = "Tax, cents/gallon", y = "Probability <= X")


########## Varying one parameter at a time ##################

# This code requires more editing between each run. 
# 1. Choose which parameter, of the par_set tibble possibilities, you would like to vary, e.g. etaMI. 
# 2. Run the eta_MI par_set, edit line 156 to say eta_MI in par_set, and hashtag eta_MI out of the pre-defined variables in the for loop, and run the for loop code.
# 3. Within the Creating Data Frame section, lines 256-364, find your eta_MI dataframe saving code, and run that. Run the ggplot code to check whether you get results that make sense.
# 4. Run the same thing for the UK, lines 498 onwards. 
# 5. Create a dataframe for the UK results as well, line 620 onwards. 
# 6. Section starting line 775 binds all the dataframes of individual parameter sensitivity together into one and saves it.
# 7. Lines 820 onwards create graphs of the analysis. 

# choose which one to run! 

#par_set <- tibble("omega"= runif(20,0.3,0.7))
#par_set <- tibble("ZPF"= runif(20,41,405))
#par_set <- tibble("ZC"= runif(20,3,14))
#par_set <- tibble("ZAMin"= runif(20,2,18))
#par_set <- tibble("ZPM"= runif(20,2,9))
#par_set <- tibble("ZQ"= runif(20,403,999))
#par_set <- tibble("ZAMac"= runif(20,1.5,15))
#par_set <- tibble("etaFF"= runif(20,-0.75, -0.21))
#par_set <- tibble("etaMacF"= runif(20,0.17,0.25))
#par_set <- tibble("etaMinF"= runif(20,-0.3, -0.05))
par_set <- tibble("etaMI"= runif(20, 0.02, 0.5))


# replace omega with var of interest; and replace the default value! 

for (etaMI in par_set) {
  
  #       Elasticities;
  etaFF = -0.36; 
  etaMinF = -0.25; 
  #etaMI = 0.4; 
  etaMacF = 0.18;
  epsiLL = 0.2; 
  epsiCLL = 0.35;
  
  #       Externalities;
  ZPF = 91; 
  ZPM = 4.5; 
  ZC = 10; 
  ZAMin = 6.4;
  ZAMac = 5.3; 
  ZQ = 691;
  
  #       Other parameters;
  omega = 0.5;
  
  #       Prices;
  tL = 0.318;
  pf = 186;
  
  #       Baselines;
  fueleff0 = 24;
  gastax0 = 55;
  Mac0 = 56.8;
  Min0 = 2191.8;
  
  #       Betas;
  betaMin = etaMinF/etaFF; betaMac = etaMacF/etaFF;
  # betaMin = 1; etaMinF/etaFF
  
  tLratio = tL/(1 - tL);
  fuel0 = Min0/fueleff0;
  
  MinF0 = Min0/fuel0;
  MacF0 = Mac0/fuel0;
  
  
  
  #  MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 # without physical activity
  MECbase =  ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0;
  
  MEBbase = tLratio*epsiLL/(1 - tLratio*epsiLL);
  taxOPT=  MECbase/(1 + MEBbase) + tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0;
  
  
  taxOPT0 = taxOPT;
  
  #Finding for stable outcome
  
  tmp1 = tLratio*epsiCLL*(1 - etaMI)/(-etaFF)
  
  # Finding the interior solution;
  # (*    by implementing the recommended tax, fuel use, inactive travel and active travel will change, which will in turn alter the optimal tax \
  # rate. Hence,instead of using the baseline input values these need  to be updated. Not note that etaZF is the fuel price elasticity of some variable Z \
  # with respect to F: i.e., the percentage increase in Z due to a one percent increase in the fuel price relative to baseline. The baseline \
  # fuel price is pf+gastax0, the percentage increase is (tax-gastax0)/(pf+gastax0). Hence, Z=Z0*(1+etaZF*(tax-gastax0)/(pf+gastax0)) 
  
  # The code below then does the following. Let gastax be the implemented gastax and taxOPT the optimal gasoline tax given gastax. \
  # Then if taxOPT and gastax deviate, we have not found the solution. If so, set the gasoline tax midway to the previous taxOPT and gastax. \
  # For this new tax level, compute Mac, Min, fuel, and Q1, and from here taxOPT. Check taxOPT against gastax and repeat if needed. Ideally, \
  # this converges to the interior solution. 
  
  # An alternative strategy would be the FindRoot command in Mathematica. The advantage of this is that Mathematica figures out \
  # the (probably more efficient) algoritm to find the solution. The disadvantage is that I have to rewrite our solution in a single \
  # lenghty equation (which risks errors). But if the strategy below does not work FindRoot offers a straightforward alternative *)
  
  gastax1 = gastax0; prec = 0.0001; approach = 0.5; dev = 1; i = 0; ALT = 0;
  
  #(* bit that remains constant *)
  
  MEB = tLratio*epsiLL/(1 - tLratio*epsiLL);
  
  while(dev >= prec){
    
    gastax1 = (1 - approach)*gastax1 + approach*taxOPT;
    
    #Endogenous variables;
    
    MinF1 = MinF0*((pf + gastax1)/(pf + gastax0))^(etaMinF - etaFF);
    MacF1 = MacF0*((pf + gastax1)/(pf + gastax0))^(etaMacF - etaFF);
    
    MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1;
    # MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 # rate without physical activity
    
    taxOPT = MEC/(1 + MEB) +tLratio*(pf + gastax1)* epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1;
    
    dev = (taxOPT - gastax1)^2;
    print(dev)
  }
  return (taxOPT)   
}

taxOPT

##### Create save-able data frames out of function output for plotting ######

# for omega

taxomega <- as.data.frame(cbind(taxOPT, omega))
taxomega$country <- "US"

ggplot(taxomega, aes(omega, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for carbon emissions

taxzpf <- as.data.frame(cbind(taxOPT, ZPF))
taxzpf$country <- "US"

ggplot(taxzpf, aes(ZPF, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for congestion

taxzc <- as.data.frame(cbind(taxOPT, ZC))
taxzc$country <- "US"

ggplot(taxzc, aes(ZC, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


# for accidents, inactive

taxzamin <- as.data.frame(cbind(taxOPT, ZAMin))
taxzamin$country <- "US"

ggplot(taxzamin, aes(ZAMin, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


# for air pollution

taxzpm <- as.data.frame(cbind(taxOPT, ZPM))
taxzpm$country <- "US"

ggplot(taxzpm, aes(ZPM, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for physical inactivity

taxzq <- as.data.frame(cbind(taxOPT, ZQ))
taxzq$country <- "US"

ggplot(taxzq, aes(ZQ, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for accidents, active

taxzamac <- as.data.frame(cbind(taxOPT, ZAMac))
taxzamac$country <- "US"

ggplot(taxzamac, aes(ZAMac, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

## PLOTTING ELASTICITIES ##

# fuel price elasticity

taxetaff <- as.data.frame(cbind(taxOPT, etaFF))
taxetaff$country <- "US"

ggplot(taxetaff, aes(etaFF, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# cross price elasticity, eta macf

taxetamacf <- as.data.frame(cbind(taxOPT, etaMacF))
taxetamacf$country <- "US"

ggplot(taxetamacf, aes(etaMacF, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# VMT price elasticity, eta minf

taxetaminf <- as.data.frame(cbind(taxOPT, etaMinF))
taxetaminf$country <- "US"

ggplot(taxetaminf, aes(etaMinF, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# income elasticity, eta mi

taxetami <- as.data.frame(cbind(taxOPT, etaMI))
taxetami$country <- "US"

ggplot(taxetami, aes(etaMI, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()




#################### for the UK ######################################


#### Sensitivity analysis varying all parameters ####


# define parameters to vary
factors <- c("etaMI", "ZPF", "ZC", "ZAMin", "ZPM", 
             "etaMinF", "etaFF",  
             "ZAMac", "omega", "ZQ",  "etaMacF")

# define how the values within the range are distributed (e.g. could be qnorm instead)
q <- c("qunif", "qunif", "qunif", "qunif", "qunif", 
       "qunif", "qunif", 
       "qunif", "qunif", "qunif", "qunif")

# define the minimum and maximum boundaries of each parameter
q.arg <- list( list(min=0.3, max=0.8), list(min=38, max=380), list(min=0.1, max=7.3), list(min=1, max=2.3), list(min=1, max=8),
               list(min=-0.5, max=-0.2),  list(min=-0.9, max=-0.3), 
               list(min=1.2, max=2.6), list(min=0.3, max=0.7), list(min=146, max=683), 
               list(min=0.17, max=0.25)) #


# Create a function and let each parameter within it vary.

oneRun <- function (etaMI, ZPF, ZC, ZAMin, ZPM, 
                    etaMinF, etaFF,  
                    ZAMac, omega, ZQ, etaMacF, gastax1, MinF1, MacF1) { 
  #       Elasticities;
  
  epsiLL = 0.2; epsiCLL = 0.35; 
  
  #       Prices;
  tL = 0.31;
  pf = 186;
  
  #       Baselines;
  fueleff0 = 28;
  gastax0 = 406;
  Mac0 = 33.3;
  Min0 = 410.4;
  
  betaMin = etaMinF/etaFF; betaMac = etaMacF/etaFF;
  # betaMin = 1; etaMinF/etaFF
  
  tLratio = tL/(1 - tL);
  fuel0 = Min0/fueleff0;
  
  MinF0 = Min0/fuel0;
  MacF0 = Mac0/fuel0; 
  
  #       Output at baseline;
  #  MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 # without physical activity
  MECbase =  ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0;
  MEBbase = tLratio*epsiLL/(1 - tLratio*epsiLL);
  taxOPT=  MECbase/(1 + MEBbase) + tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0;
  taxOPT0 = taxOPT
  #Finding for stable outcome
  
  tmp1 = tLratio*epsiCLL*(1 - etaMI)/(-etaFF);
  
  # The code below then does the following. Let gastax be the implemented gastax and taxOPT the optimal gasoline tax given gastax. \
  # Then if taxOPT and gastax deviate, we have not found the solution. If so, set the gasoline tax midway to the previous taxOPT and gastax. \
  # For this new tax level, compute Mac, Min, fuel, and Q1, and from here taxOPT. Check taxOPT against gastax and repeat if needed. Ideally, \
  # this converges to the interior solution. 
  
  gastax1 = gastax0; prec = 0.0001; approach = 0.5; dev = 1; i = 0; ALT = 0;
  
  #(* bit that remains constant *)
  
  MEB = tLratio*epsiLL/(1 - tLratio*epsiLL);
  
  while(dev >= prec){
    
    gastax1 = (1 - approach)*gastax1 + approach*taxOPT;
    
    #Endogenous variables;
    
    MinF1 = MinF0*((pf + gastax1)/(pf + gastax0))^(etaMinF - etaFF);
    MacF1 = MacF0*((pf + gastax1)/(pf + gastax0))^(etaMacF - etaFF);
    
    MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1;
    # MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 # rate without physical activity
    
    taxOPT = MEC/(1 + MEB) +tLratio*(pf + gastax1)* epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1;
    
    dev = (taxOPT - gastax1)^2;
    print(dev)
  }
  
  # define what output we are interested in.
  
  return (taxOPT) 
}


# Run the function above by feeding it an empty dataframe my.data with 11 different parameters that can vary

modelRun <- function (my.data) {
  return(mapply(oneRun, my.data[,1], my.data[,2], my.data[,3], my.data[,4], my.data[,5], my.data[,6], my.data[,7], my.data[,8], my.data[,9], my.data[,10], my.data[,11]))
}

# run the LHS sensitivity command, giving it the function to evaluate in the form of a dataframe, 
# how many times it should run, and what the parameters are. 
# Advantage over Monte Carlo is significant time savings in running this over the MonteCarlo package commands (i.e. seconds vs 10+minutes).

myLHS2 <- LHS(modelRun, factors, 200, q, q.arg, nboot=50)


# Plot results to see where solution is most likely to lie (the steepest part of the curve)
plotecdf(myLHS2, main=NULL, ylab="Probability <= X", xlab="Tax, cents/gallon", col="light blue") # optimum is less defined - could be much higher than we are currently saying!  

# Plot which parameter influences the optimal fuel tax the most. 
plotprcc(myLHS2, main=NULL, ylab="Partial rank correlation coefficient", col = "light blue")

# Plotting prettier, optional
data2 <- as.data.frame(myLHS2[["res"]])
ggplot(data2, aes(V1)) + 
  stat_ecdf(geom = "step", col = "light blue") + 
  geom_vline(xintercept=406, linetype="dashed", color = "red", size=0.5) +
  geom_hline(yintercept=0.15, linetype="dashed", color = "red", size=0.5) +
  theme_minimal() +
  labs(x = "Tax, cents/gallon", y = "Probability <= X")

# Putting them together in one
ggplot(data2, aes(V1)) + 
  stat_ecdf(geom = "step", col = "light blue", size = 1) + 
  geom_vline(xintercept=406, linetype="dashed", color = " light blue", size=0.5) +
 # geom_hline(yintercept=0.15, linetype="dashed", color = "orange", size=0.5) +
  theme_minimal() +
  xlim(0, 2000) + 
  labs(x = "Tax, cents/gallon", y = "Probability <= X") +
  stat_ecdf(data = data, aes(V1), col = "dark blue", size = 1) +
  geom_vline(xintercept=55, linetype="dashed", color = "dark blue", size=0.5) 
# geom_hline(yintercept=0.001, linetype="dashed", color = "dark red", size=0.5)


ggplot(data2) + 
  stat_ecdf(geom = "step", aes(V1, color = "light blue")) + 
  geom_vline(xintercept=406, linetype="dashed", aes(color = "light blue"), size=0.5) +
  # geom_hline(yintercept=0.15, linetype="dashed", color = "orange", size=0.5) +
  theme_minimal() +
  xlim(0, 2000) + 
  labs(x = "Tax, cents/gallon", y = "Probability <= X") +
  stat_ecdf(data = data, aes(V1, color = "blue")) +
  geom_vline(xintercept=55, linetype="dashed", aes(color = "blue"), size=0.5) 
# geom_hline(yintercept=0.001, linetype="dashed", color = "dark red", size=0.5)


############# Varying one parameter at a time ############################

# choose which one to run! 

#par_set <- tibble("omega"= runif(20,0.3,0.7))
par_set <- tibble("ZPF"= runif(20,38,380))
#par_set <- tibble("ZC"= runif(20,0.1,7.3))
#par_set <- tibble("ZAMin"= runif(20,1, 2.3))
#par_set <- tibble("ZPM"= runif(20,1, 8))
#par_set <- tibble("ZQ"= runif(20,146, 683))
#par_set <- tibble("ZAMac"= runif(20,1.2, 2.6))
#par_set <- tibble("etaFF"= runif(20,-0.9, -0.3))
#par_set <- tibble("etaMacF"= runif(20,0.17,0.25))
#par_set <- tibble("etaMinF"= runif(20,-0.5, -0.2))
#par_set <- tibble("etaMI"= runif(20, 0.3, 0.8))


# replace omega with var of interest; and replace the default value! 

for (ZPF in par_set) {
  
  #       Elasticities;
  etaFF = -0.48; 
  etaMinF = -0.35; 
  etaMI = 0.605; 
  etaMacF = 0.18;
  epsiLL = 0.2; epsiCLL = 0.35;
  
  #       Externalities;
  #  ZPF = 86; 
  ZPM = 3.6; 
  ZC = 5; 
  ZAMin = 1.6; 
  ZAMac = 1.6; 
  ZQ = 244;
  
  #       Other parameters;
  omega = 0.5;
  
  #       Prices;
  tL = 0.31;
  pf = 186;
  
  #       Baselines;
  fueleff0 = 28; 
  gastax0 = 406;
  Mac0 = 33.3;
  Min0 = 410.4; 
  
  #       Betas;
  betaMin = etaMinF/etaFF; betaMac = etaMacF/etaFF;
  # betaMin = 1; etaMinF/etaFF
  
  tLratio = tL/(1 - tL);
  fuel0 = Min0/fueleff0;
  
  MinF0 = Min0/fuel0;
  MacF0 = Mac0/fuel0;
  
  
  
  #  MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 # without physical activity
  MECbase =  ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0;
  
  MEBbase = tLratio*epsiLL/(1 - tLratio*epsiLL);
  taxOPT=  MECbase/(1 + MEBbase) + tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0;
  
  
  taxOPT0 = taxOPT;
  
  #Finding for stable outcome
  
  tmp1 = tLratio*epsiCLL*(1 - etaMI)/(-etaFF)
  
  # Finding the interior solution;
  # (*    by implementing the recommended tax, fuel use, inactive travel and active travel will change, which will in turn alter the optimal tax \
  # rate. Hence,instead of using the baseline input values these need  to be updated. Not note that etaZF is the fuel price elasticity of some variable Z \
  # with respect to F: i.e., the percentage increase in Z due to a one percent increase in the fuel price relative to baseline. The baseline \
  # fuel price is pf+gastax0, the percentage increase is (tax-gastax0)/(pf+gastax0). Hence, Z=Z0*(1+etaZF*(tax-gastax0)/(pf+gastax0)) 
  
  # The code below then does the following. Let gastax be the implemented gastax and taxOPT the optimal gasoline tax given gastax. \
  # Then if taxOPT and gastax deviate, we have not found the solution. If so, set the gasoline tax midway to the previous taxOPT and gastax. \
  # For this new tax level, compute Mac, Min, fuel, and Q1, and from here taxOPT. Check taxOPT against gastax and repeat if needed. Ideally, \
  # this converges to the interior solution. 
  
  # An alternative strategy would be the FindRoot command in Mathematica. The advantage of this is that Mathematica figures out \
  # the (probably more efficient) algoritm to find the solution. The disadvantage is that I have to rewrite our solution in a single \
  # lenghty equation (which risks errors). But if the strategy below does not work FindRoot offers a straightforward alternative *)
  
  gastax1 = gastax0; prec = 0.0001; approach = 0.5; dev = 1; i = 0; ALT = 0;
  
  #(* bit that remains constant *)
  
  MEB = tLratio*epsiLL/(1 - tLratio*epsiLL);
  
  while(dev >= prec){
    
    gastax1 = (1 - approach)*gastax1 + approach*taxOPT;
    
    #Endogenous variables;
    
    MinF1 = MinF0*((pf + gastax1)/(pf + gastax0))^(etaMinF - etaFF);
    MacF1 = MacF0*((pf + gastax1)/(pf + gastax0))^(etaMacF - etaFF);
    
    MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1;
    # MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 # rate without physical activity
    
    taxOPT = MEC/(1 + MEB) +tLratio*(pf + gastax1)* epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1;
    
    dev = (taxOPT - gastax1)^2;
    print(dev)
  }
  return (taxOPT)   
}
# ignore warnings

taxOPT




#### Data frames for UK values ######



# plotting omega


taxomega2 <- as.data.frame(cbind(taxOPT, omega))
taxomega2$country <- "UK"
taxomega3 <- rbind(taxomega, taxomega2)
taxomega3$var <- "omega"
names(taxomega3)[names(taxomega3)=="omega"] <- "value"


ggplot(taxomega3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


# for carbon emissions

taxzpf2 <- as.data.frame(cbind(taxOPT, ZPF))
taxzpf2$country <- "UK"
taxzpf3 <- rbind(taxzpf, taxzpf2)
taxzpf3$var <- "ZPF"
names(taxzpf3)[names(taxzpf3)=="ZPF"] <- "value"


ggplot(taxzpf3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for congestion

taxzc2 <- as.data.frame(cbind(taxOPT, ZC))
taxzc2$country <- "UK"
taxzc3 <- rbind(taxzc, taxzc2)
taxzc3$var <- "ZC"
names(taxzc3)[names(taxzc3)=="ZC"] <- "value"


ggplot(taxzc3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


# for accidents, inactive

taxzamin2 <- as.data.frame(cbind(taxOPT, ZAMin))
taxzamin2$country <- "UK"
taxzamin3 <- rbind(taxzamin, taxzamin2)
taxzamin3$var <- "ZAMin"
names(taxzamin3)[names(taxzamin3)=="ZAMin"] <- "value"


ggplot(taxzamin3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


# for air pollution

taxzpm2 <- as.data.frame(cbind(taxOPT, ZPM))
taxzpm2$country <- "UK"
taxzpm3 <- rbind(taxzpm, taxzpm2)
taxzpm3$var <- "ZPM"
names(taxzpm3)[names(taxzpm3)=="ZPM"] <- "value"


ggplot(taxzpm3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for physical inactivity

taxzq2 <- as.data.frame(cbind(taxOPT, ZQ))
taxzq2$country <- "UK"
taxzq3 <- rbind(taxzq, taxzq2)
taxzq3$var <- "ZQ"
names(taxzq3)[names(taxzq3)=="ZQ"] <- "value"


ggplot(taxzq3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# for accidents, active

taxzamac2 <- as.data.frame(cbind(taxOPT, ZAMac))
taxzamac2$country <- "UK"
taxzamac3 <- rbind(taxzamac, taxzamac2)
taxzamac3$var <- "ZAMac"
names(taxzamac3)[names(taxzamac3)=="ZAMac"] <- "value"


ggplot(taxzamac3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

## PLOTTING ELASTICITIES ##

# fuel price elasticity

taxetaff2 <- as.data.frame(cbind(taxOPT, etaFF))
taxetaff2$country <- "UK"
taxetaff3 <- rbind(taxetaff, taxetaff2)
taxetaff3$var <- "eta_ff"
names(taxetaff3)[names(taxetaff3)=="etaFF"] <- "value"



ggplot(taxetaff3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# cross price elasticity, eta macf

taxetamacf2 <- as.data.frame(cbind(taxOPT, etaMacF))
taxetamacf2$country <- "UK"
taxetamacf3 <- rbind(taxetamacf, taxetamacf2)
taxetamacf3$var <- "eta_macf"
names(taxetamacf3)[names(taxetamacf3)=="etaMacF"] <- "value"


ggplot(taxetamacf3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# VMT price elasticity, eta minf

taxetaminf2 <- as.data.frame(cbind(taxOPT, etaMinF))
taxetaminf2$country <- "UK"
taxetaminf3 <- rbind(taxetaminf, taxetaminf2)
taxetaminf3$var <- "eta_minf"
names(taxetaminf3)[names(taxetaminf3)=="etaMinF"] <- "value"


ggplot(taxetaminf3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()

# income elasticity, eta mi

taxetami2 <- as.data.frame(cbind(taxOPT, etaMI))
taxetami2$country <- "UK"
taxetami3 <- rbind(taxetami, taxetami2)
taxetami3$var <- "eta_mi"
names(taxetami3)[names(taxetami3)=="etaMI"] <- "value"


ggplot(taxetami3, aes(value, taxOPT, colour = country)) + 
  geom_line() +# and add an x where relevant
  theme_classic()


#### JOIN THEM ALL TOGETHER ######

allresults <- rbind(taxomega3, taxzpf3, taxzpm3, taxzamac3, taxzamin3, taxzq3, taxzc3, 
                    taxetaff3, taxetamacf3, taxetami3, taxetaminf3)

write.csv(x = allresults, "all_results_sensitivity_proper.csv")


################### Adding annotations #########################

# This code adds Xs to mark the central values in our sensitivity graphs. 
# for the UK
anno_uk <- data.frame(xstar = c(0.5, 86, 3.6, 1.6, 1.6, 244, 5, 
                                -0.48, 0.18, 0.605, -0.35), 
                      ystar = c(454, 454,454, 454,454, 454,454, 454,454, 454, 454),
                      lab = c("x", "x", "x", "x", "x", "x", "x", "x", "x", "x", "x"),
                      var = c("5.Health internalisation rate", "1.Pollution (CO2)", "5.Pollution (air)", "3.Accidents, active", 
                              "2.Accidents, inactive", "4.Inactivity", "6.Congestion", "1.Fuel price elasticity", "3.Cross elasticity of active travel",
                              "4.Income el. of inactive travel", "2.VMT-fuel price elasticity"),
                      country = c("UK", "UK", "UK", "UK", "UK", "UK", "UK", "UK", "UK", "UK", "UK"))

# for the US
anno_us <- data.frame(xstar = c(0.5, 91, 4.5, 5.3, 6.4, 691, 10, 
                                -0.36, 0.18, 0.4, -0.25), 
                      ystar = c(1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013),
                      lab = c("x", "x", "x", "x", "x", "x", "x", "x", "x", "x", "x"),
                      var = c("5.Health internalisation rate", "1.Pollution (CO2)", "5.Pollution (air)", "3.Accidents, active", 
                              "2.Accidents, inactive", "4.Inactivity", "6.Congestion", "1.Fuel price elasticity", "3.Cross elasticity of active travel",
                              "4.Income el. of inactive travel", "2.VMT-fuel price elasticity"), 
                      country = c("US", "US", "US", "US", "US", "US", "US", "US", "US", "US", "US"))

# Join the two
anno <- rbind(anno_uk, anno_us)

# Create subsets to feed into the ggplot faceted graphs. 

anno_elasticities <- anno %>% filter(var %in% c("3.Cross elasticity of active travel", "1.Fuel price elasticity", "4.Income el. of inactive travel", "2.VMT-fuel price elasticity", "5.Health internalisation rate"))
anno_costs <- anno %>% filter(!(var %in% c("5.Health internalisation rate", "3.Cross elasticity of active travel", "1.Fuel price elasticity", "4.Income elasticity of inactive travel", "2.VMT-fuel price elasticity")))



################## PLOTTING FACETED SENSITIVITY GRAPHS ##############################

# load file with all sensitivity results from PSE analysis

#allresults <- read.csv("~/all_results_sensitivity_proper.csv") 


allresults$var <- as.character(allresults$var)
allresults$var[allresults$var=="omega"] <- "5.Health internalisation rate"
allresults$var[allresults$var=="ZPF"] <- "1.Pollution (CO2)"
allresults$var[allresults$var=="ZPM"] <- "5.Pollution (air)"
allresults$var[allresults$var=="ZAMac"] <- "3.Accidents, active"
allresults$var[allresults$var=="ZAMin"] <- "2.Accidents, inactive"
allresults$var[allresults$var=="ZQ"] <- "4.Inactivity"
allresults$var[allresults$var=="ZC"] <- "6.Congestion"
allresults$var[allresults$var=="eta_ff"] <- "1.Fuel price elasticity"
allresults$var[allresults$var=="eta_macf"] <- "3.Cross elasticity of active travel"
allresults$var[allresults$var=="eta_mi"] <- "4.Income el. of inactive travel"
allresults$var[allresults$var=="eta_minf"] <- "2.VMT-fuel price elasticity"

allresults$country <- as.character(allresults$country)
allresults$country[allresults$country=="USA"] <- "US"


# All 11 parameters

ggplot(allresults) +
  aes(x = value, y = taxOPT, colour = country) +
  geom_point(size = 0.94) +
  geom_smooth(span = 0.75) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Value", y = "Optimal tax, cents/gallon", color = "Country") +
  theme_minimal() + 
  theme(panel.spacing = unit(1.5, "lines")) +
  facet_wrap(vars(var), scales = "free_x") +
  theme(strip.text.x = element_text(
    size = 7
  )) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab, group = country), size=7)



# only elasticities

allresults %>%
  filter(var %in% c("5.Health internalisation rate", "1.Fuel price elasticity", "3.Cross elasticity of active travel", "4.Income el. of inactive travel", "2.VMT-fuel price elasticity")) %>%
  ggplot() +
  aes(x = value, y = taxOPT, colour = country) +
  geom_point(size = 0.94) +
  geom_smooth(span = 0.75) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "Value", y = "Optimal tax, cents/gallon", color = "Country") +
  theme_minimal() +
  facet_wrap(vars(var), scales = "free_x")+
  theme(strip.text = element_text(size = 11)) +
  geom_text(data = anno_elasticities, aes(x = xstar,  y = ystar, label = lab, group = country))


# only social costs

allresults %>%
  filter(!(var %in% c("5.Health internalisation rate", "1.Fuel price elasticity", "3.Cross elasticity of active travel", "4.Income el. of inactive travel", "2.VMT-fuel price elasticity"
  ))) %>%
  ggplot() +
  aes(x = value, y = taxOPT, colour = country) +
  geom_point(size = 0.94) +
  geom_smooth(span = 0.75) +
  scale_color_brewer(palette = "Paired") +
  labs(x = "US cents", y = "Optimal tax, cents/gallon", color = "Country") +
  theme_minimal() +
  facet_wrap(vars(var), scales = "free_x")+
  theme(strip.text = element_text(size = 11)) +
  geom_text(data = anno_costs, aes(x = xstar,  y = ystar, label = lab, group = country))

# IF having trouble understanding plotting code, I recommend installing the esquisser package. 
# Then, run esquisse(data.frame) (in our case allresults), and a visual interface will pop up.
# This lets you edit the design, colour scheme, legend, titles, faceting etc of the graphs. 

