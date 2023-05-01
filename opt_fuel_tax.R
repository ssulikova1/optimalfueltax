# clear workspace
rm(list = ls())
switch_UK = T; #T for UK, F for US

# Parameter inputs
if(switch_UK){ #       Parameter inputs - UK;
  
  #       Elasticities;
  etaFF = -0.48; etaMinF = -0.35; etaMI = 0.605; etaMacF = 0.185;
  epsiLL = 0.2; epsiCLL = 0.35;
  
  #       Externalities;
  ZPF = 150.1; ZPM = 4.3; ZC = 6; ZAMin = 1.9; ZAMac = 1.9; ZQ = 291.3;
  
  #       Other parameters;
  omega = 0.5;
  
  #       Prices;
  tL = 0.31;
  pf = 207.3;
  
  #       Baselines;
  fueleff0 = 28; 
  gastax0 = 381.7;
  Mac0 = 33.3;
  Min0 = 410.4;
}
if(!switch_UK){ #       Parameter inputs - US;
  
  #       Elasticities;
  etaFF = -0.36; etaMinF = -0.25; etaMI = 0.4; etaMacF = 0.185;
  epsiLL = 0.2; epsiCLL = 0.35;
  
  #       Externalities;
  ZPF = 160.1; ZPM = 5.4; ZC = 11.9; ZAMin = 7.6; ZAMac = 6.3; ZQ = 825;
  
  #       Other parameters;
  omega = 0.5;
  
  #       Prices;
  tL = 0.318;
  pf = 288.9;
  
  #       Baselines;
  fueleff0 = 24;
  gastax0 = 50.2;
  Mac0 = 56.8;
  Min0 = 2191.8;
}

#       Betas;
betaMin = etaMinF/etaFF; betaMac = etaMacF/etaFF;
# betaMin = 1; etaMinF/etaFF

tLratio = tL/(1 - tL);
fuel0 = Min0/fueleff0;
     
MinF0 = Min0/fuel0;
MacF0 = Mac0/fuel0;
       
       
paste("F=", fuel0, "  Min=", Min0, "  Mac=", Mac0)

#       F=14.6571  Min=410.4  Mac=33.3
       
# # for naive rate, MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0

#       Output at baseline;
    #  MECbase=   ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 # without physical activity
       MECbase =  ZPF + (ZC + ZAMin + ZPM)*betaMin*MinF0 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0;
       
       MEBbase = tLratio*epsiLL/(1 - tLratio*epsiLL);
       taxOPT = MECbase/(1 + MEBbase) + tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0;
       
       
       paste("optimal tax=", taxOPT);
       paste("fuel efficiency=", MinF0);
       
       paste("fuel externality component=", ZPF );
       
       paste("distance externality component=", ZPM *betaMin*MinF0);
       paste("accident in externality component=", ZAMin *betaMin*MinF0);
       paste("congestion  component=", ZC *betaMin*MinF0);
       
       paste("active travel component=", (ZAMac - (1 - omega)*ZQ)*betaMac*MacF0);
       
       paste("active travel accident component=", (ZAMac)*betaMac*MacF0);
       paste("active travel health component=", (-(1 - omega)*ZQ)*betaMac*MacF0);
       paste("MEC=", MECbase);
       paste("MEC/(1+MEBbase)=", MECbase/(1 + MEBbase));
       paste("Excess burden=", MECbase*(1/(1 + MEBbase) - 1));
       
       paste("Ramsey tax component=", tLratio*(pf + gastax0)*epsiCLL*(1 - etaMI)/(-etaFF));

       paste("Congestion feedback=", tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF0);

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
      
      MinF1 = MinF0*((pf + gastax1)/(pf + gastax0))^(etaMinF - etaFF); #inactive miles per gallon of fuel
      MacF1 = MacF0*((pf + gastax1)/(pf + gastax0))^(etaMacF - etaFF); # active miles per gallon of fuel
      
      MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 + (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1; #marginal external cost
   #   MEC = ZPF + (ZC + ZAMin + ZPM)*betaMin * MinF1 # rate without physical activity
      
      taxOPT = MEC/(1 + MEB) +tLratio*(pf + gastax1)* epsiCLL*(1 - etaMI)/(-etaFF) + tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1;

      dev = (taxOPT - gastax1)^2;
      print(dev)
}


print(paste("optimal tax=", taxOPT));
print(paste("fuel efficiency=", MinF1));
print(paste("fuel externality component=", ZPF ));
print(paste("distance externality component=", ZPM *betaMin*MinF1));
print(paste("accident in externality component=", ZAMin *betaMin*MinF1));
print(paste("congestion  component=", ZC *betaMin*MinF1));
print(paste("active travel component=", (ZAMac - (1 - omega)*ZQ)*betaMac*MacF1));
print(paste("active travel accident subcomponent=", (ZAMac)*betaMac*MacF1));
print(paste("active travel health subcomponent=", (-(1 - omega)*ZQ)*betaMac*MacF1));
print(paste("Ramsey tax component=",tLratio*(pf + gastax1)*epsiCLL*(1 - etaMI)/(-etaFF)));
print(paste("Congestion feedback=", tLratio*(epsiLL - epsiCLL*(1 - etaMI))*ZC*betaMin*MinF1));
print(paste("Excess burden=", MEC*(1/(1 + MEB) - 1)))

print(paste("percentage tax increase  due to interior solution= ", 100*(taxOPT/taxOPT0 - 1), " percent"));

print(paste("change in F= ", (((pf + gastax1)/(pf + gastax0))^(etaFF) - 1)*100, " percent"));
print(paste("change in Min= ", (((pf + gastax1)/(pf + gastax0))^(etaMinF) - 1)*100, " percent"));
print(paste("change in Mac= ", (((pf + gastax1)/(pf + gastax0))^(etaMacF) - 1)*100, " percent"))
print(paste("New gas tax revenues as a % proportion old revenues= ",100*(gastax1/gastax0)-1+100*(((pf+gastax1)/(pf+gastax0))^(etaFF)-1), " percent"))

##########################################
# calculating the welfare change
  
  # Note that the specification below endogenizes F, but takes as given the t_L, which is implicit in MEB *)
  
  taxOLD=gastax0;  taxNEW=taxOPT;
  
  # gain=NIntegrate[(-(etaFF/(pf+x))*fuel0*((pf+x)/(pf+gastax0))^(etaFF))*(1+MEB)*(taxOPT-x),{x,taxOLD,taxNEW}];
  integrand <- function(x) {(-(etaFF/(pf+x))*fuel0*((pf+x)/(pf+gastax0))^(etaFF))*(1+MEB)*(taxOPT-x)}
  gain = integrate(integrand, taxOLD, taxNEW)

  print(paste("Consumption-equivalent welfare gain=", gain$value, " cents per year"));# multiply by a billion and divide by population of country to get per person per year value
  
  #(* fuelexp0 is based on the taxOLD level of the tax. Note that the corresponding fuel use at taxOLD is fuel0*((pf+taxOLD)/(pf+gastax0))^(etaFF) *) 
  
  fuelexp0=(pf+gastax0)*fuel0
  
  print(paste("Current tax fuel expenditure=", fuelexp0, " cents per year")); # multiply by a billion and divide by population of country to get per person per year value
  print(paste("Welfare gain as % of pretax fuel expenditure=", 100*gain$value/fuelexp0, " percent"))
  
  #During evaluation of In[1773]:= Consumption-equivalent welfare gain=41.7629 cents per year
  #During evaluation of In[1773]:= Pretax fuel expenditure=8677.03 cents per year
  #During evaluation of In[1773]:= Welfare gain as % of pretax fuel expenditure=0.481304 percent

    
  
  
  