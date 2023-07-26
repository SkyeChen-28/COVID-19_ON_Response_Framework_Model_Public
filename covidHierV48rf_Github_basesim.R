#     BASESIM
#     covidHierV48_Github_basesim.R aggregates parameter sets from fitting algorithm (output from covidHier48_Github_parmfit.R) and runs base simulations
#     Copyright (C) 2021 Skye SG Chen, Kathyrn R Fair, Vadim A Karatayev
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Dependencies
library(abind); 
library(plyr); 
library(dplyr);
library(data.table); 
library(EpiEstim);  
library(tidyverse); 
library(nloptr); 
library(mgcv); 
library(stringr); 
library(ggplot2); 
library(RColorBrewer);      
multiplot <- dget("multiplot.R")

# Suppress summarise info
options(dplyr.summarise.inform = FALSE);

# Directory variables
input_dir = "InputFiles" # Directory where all the input files are stored
output_dir = "OutputFiles/Basesim" # Directory to store output files
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) # Creates output directory if it doesn't exist

#Read in all necessary data for fitting
cmat=readRDS(sprintf("%s/cmat_data_weighted.rds", input_dir)); #Contact matrix data from Prem et al. (2020)
DATA=readRDS(sprintf("%s/covidHierData.rds", input_dir)); #See readme file for description of contents
omgdat<-readRDS(sprintf("%s/mobilitydat_REAL_to2021-02-27_V4.rds", input_dir)) #read in mobility data
regiondat=readRDS(sprintf("%s/ONdat_casesbyPHU_to2021-08-10.rds", input_dir)); ### Contains new region specific known case data for fitting
agedat=readRDS(sprintf("%s/ONdat_newKage_to2021-08-10.rds", input_dir)); ### Contains new age specific known case data for fitting
regionid=readRDS(sprintf("%s/PHUtoREGIONlinker_numeric.rds", input_dir)); #COntains data for linking census divisions to their associated PHU
rfdat = readRDS(sprintf("%s/response_framework_refined_2021-05-01.rds", input_dir)) # read in response framework data
testing_volumes=readRDS(sprintf("%s/tests_by_PHU_aug112021.rds", input_dir)) # Read in testing volumes by PHU data

# Refine the data
agedat$tot<-cumsum(rowSums(agedat[,3:7])) #Add column of total new cases across all ages
regiondat$tot<-cumsum(rowSums(regiondat[,2:35])) #Add column of total new cases across all regions
testing_volumes$ts50 = as.numeric(testing_volumes$Date - as.Date("2020-03-10")) #Add column for ts50, days since the start of the epidemic

#Drop all dates where totals are impacted by reporting lags
end_date<-"2021-02-07"
agedat<-agedat[agedat$Date<=end_date,]
regiondat<-regiondat[regiondat$Date<=end_date,]
rfdat<-rfdat[rfdat$Date<=end_date,]

# Constants to toggle debugging tools
parms_initialized = FALSE

###vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Control panel vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv###
### These variables controls which aspects of the simulation we want active
# When running the script on a supercomputer, turn this variable on! Disables messages and outputs used for debugging that aren't necessary for the supercomputer to have
supercomputer_mode = TRUE
debug_mode = FALSE

### Artifact of previous code, do not modify!
# cftype:        Set closure counterfactual type i.e. schools remain open ("schoolopen"), workplaces remain open ("workopen"), both remain open ("bothopen") or both are shut ("neitheropen", corresponds to what actually occured in the province)
# reopeningtype: Set reopening type; with ("restricted") or without ("unrestricted") NPIs in schools/workplaces
# vdtype:        Set individual NPI adherence counterfactual type i.e. no individual adherence to NPIs ("vdOFF") or individual adherence to NPIs in response to case numbers ("vdON")
cftype<-"neitheropen"; 
reopeningtype<-"restricted"; 
vdtype<-"vdON";


# Testing and debugging variables
printing_timestep = 20            # Print the current day in the simulation if it's a multiple of this int
tepi_message = FALSE              # Prints tepi messages
disable_response_framework<-FALSE # Controls whether or not the response framework is active
init_s_vec <- c(0,30,60,90,300)   # Vector containing initial shutdown lengths (DEPRECATED, but do not remove as parts of the code require this.)
indicators_active = c(            # Controls which indicators are active. Comment out the indicators you want inactive
  ""       # Placeholder, keep "" uncommented
  , "wir"  # 2 (Weekly Incidence Rate)
  , "ppos" # 3 (Percent Positivity)
  , "Rt"   # 4 (Effective Reproduction Number)
)

# Note: The simulation is designed to only calculate indicator values on days when a decision to update the response framework is required!
# This is hard coded in to increase time efficiency of the sim, but you can change when the response framework is activated by modifying rf_start and rf_end.

#### Counterfactual scenarios
zone_change_backup_func = median # The function to use when there is no unique mode from the three indicators. (default: median minimizes the error between data and sim)
###^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Control panel ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^###

tstart<-Sys.time()

# #Function to adjust for census region-level differences in transmission probability based on population (xi_k values in paper):
BetaMod <- function(xi_coefs, xi_PHU = rep(1,34), region_id = regionid, xi_census = rep(1,49)) {
  vals=regionid$phunum[order(regionid$census.region)]
  PHU_populations = aggregate(x = region_id[,"pop2016"], 
                              by = list(region_id$phunum), 
                              FUN = sum)$x # Array to store each PHU's population
  a = xi_coefs[1]
  b = xi_coefs[2]
  c = xi_coefs[3]
  
  xi_PHU = a * exp(-b * PHU_populations) + c
  xi_census = xi_PHU[vals]
  return(xi_census)
}

# vvvvv Constants for fitting vvvvv #
xi_a = 9.315907    
xi_b = 3.240835e-06
xi_c = 4.816625e-01   
xi_coefs = c(xi_a, xi_b, xi_c)
# xi_census = BetaMod(xi_coefs) 


eps_w = 5.095690e-01
alpha_w1 = 8.245924e-01   
alpha_w2 = 7.524312e-01   
alpha_w3 = 4.204813e-01   
alpha_w4 = 3.510599e-01   

eps_s = 6.609320e-01
alpha_s1 = 3.735138e-01 
alpha_s2 = 6.800928e-01   
alpha_s3 = 5.226951e-01 
alpha_s4 = 4.378218e-01   

eps_h = 6.115866e-01
alpha_h1 = 7.044552e-01
alpha_h2 = 4.635093e-01
alpha_h3 = 3.706887e-01   
alpha_h4 = 4.859584e-02   

eps_o = 8.991411e-01
alpha_o1 = 7.326203e-01
alpha_o2 = 6.993630e-01   
alpha_o3 = 5.175376e-01   
alpha_o4 = 6.901249e-01 
# ^^^^^ Constants for fitting ^^^^^ #

# Constant defaults, only modify these variables in the counterfactuals script!
mu = 1
rf_start = 241
rf_end = 290

#Set code version
codeselect<-"v48rf"

path = "ParmfitOutputFiles/OneWeekRun/OutputFiles" #Set to wherever output files from fitting are saved

#Find all output files from parameter fitting
file.names <- dir(path, pattern =sprintf("parmfit_y_%s_", codeselect), full.names=TRUE)
file.names.x <- dir(path, pattern =sprintf("parmfit_x_%s_", codeselect), full.names=TRUE)

#Counts number of fits that meet our criteria for being a reasonably good fit
goodfits<-1;

for (fit in 1:length(file.names))
{
  print(fit)
  y=readRDS(file.names[fit]);
  x.dat=readRDS(file.names.x[fit]);
  
  objective_threshold = 2000
  if (y$objective<=objective_threshold) #Throw out any parm combinations that don't give a reasonably good fit (reduces runtime by throwing out any parameter sets where the cost function value from the fitting is quite high)
  {
    cat("fit = ", fit, "\n")
    print(y$objective)
    
    ###vvvvvvvvvvvvvv Copy and paste simulation code here!!! vvvvvvvvvvvvvv###
    
    #Some parameter name differences from paper:
    #initial testing tauI_t0=cvTl, final testing tauI_tf=tauI, tauA=tauEAf*tauI
    # a1,a2,a3,a4,a5 are age-specific susceptibility modifiers (gamma_i, i=1,2,3,4,5 in the paper)
    # boost is L0 (impact of stay at home orders on NPI adherence)
    #Mfact scales how much travel happens compared to reality; Mfact=1 is just movement from commuting data
    #epsP allows individual NPI adherence efficacy to differ from that of closures (epsP=eps in paper)
    #Tg and Tl are the gobal and local closure thresholds, n0 is fraction initially infected
    #Msave is the travel matrix. Here Msave entries are numbers of commuters; msim3 converts them to proportions
    pops=colSums(DATA$Msave)
    
    ###Adjust to 2020 Q4 pop estimate (from https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000901), as region pop estimates are from 2016 census
    pop2020<-14733119
    popratio<-pop2020/sum(pops)
    popadj2020=round(popratio*pops)
    
    # Parms array for fitting. Comment out for simulation! Uncomment for fitting!
    # Remember to update the simulation parms array if it's modified here!
    parms_initialized = TRUE
    if(debug_mode){
      warning("
    Comment out the hardcoded parms array definition.
    This is for initializing the fitting code.
    Use the fitted parameters instead!")}
    # Note: some parameters named here are artefacts of a previous version of the model. They are retained here to avoid breaking analysis code which requires input files to have a specific dimension but do not impact simulations
    parms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),
                xi_a = xi_a, xi_b = xi_b, xi_c = xi_c,
                # beta=BetaMod(rep(1,34))*exp(0)/mean(popadj2020), # original
                beta=BetaMod(xi_coefs)*exp(0)/mean(popadj2020),
                eps_w=eps_w, eps_s=eps_s,eps_h=eps_h,eps_o=eps_o,omg=exp(9.5),r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
                atravel=0.5, a1=1,a2=1.75,a3=1,a4=5,a5=35, reopen_w=0.5, reopen_s=0.2, B=0.3, phi=80, zeta=0.025,
                tau0=0.035, psi1=0.025,  psi2=0.025,  psi3=0.025,  psi4=0.025,  psi5=0.1, kappa=0.3, taumax=0.5, boost=0.5,
                alpha_w1 = alpha_w1, alpha_w2 = alpha_w2, alpha_w3 = alpha_w3, alpha_w4 = alpha_w4,
                alpha_s1 = alpha_s1, alpha_s2 = alpha_s2, alpha_s3 = alpha_s3, alpha_s4 = alpha_s4,
                alpha_h1 = alpha_h1, alpha_h2 = alpha_h2, alpha_h3 = alpha_h3, alpha_h4 = alpha_h4,
                alpha_o1 = alpha_o1, alpha_o2 = alpha_o2, alpha_o3 = alpha_o3, alpha_o4 = alpha_o4,
                DATA$Msave[-50,]);
    
    # Specification of closure/reopening events
    inCstart=50; #all events occur at a distance from the day total cases >=50 (March 10th)
    inClen_w1=79; #Duration of initial provincial workplace closure (March 25-June 11, Phase 2 commenced June 12)
    inClen_s1=178; #Duration of initial provincial school closure (Mar 14 - Sept 7, Schools reopen Sept 8)
    inClen_w2=46; #Duration of second provincial workplace closure, to first part of staggered reopening (Dec 26th 2020 - Feb 9th 2021)
    inClen_s2=51; #Duration of second provincial school closure, to first part of staggered reopening (Dec 21st 2020 - Feb 9th 2021)
    inBlen_s=70; #Duration of Summer break (June 30th - Sept 8th)
    
    workstaggergap1<-6; #Gap for first 3 reopening stages (Feb 10/16/22)
    workstaggergap2<-20; #Gap between Feb 16th and march 8th
    schoolstaggergap<-8; #Gap between 2 stages of return to in-person learning (Feb 8/16)
    inClen_w2=rep(46+workstaggergap1,49); #days between Dec 26th lockdown and Feb 16th main reopening day (46 days is to the 10th)
    inClen_s2=rep(49,49); #days between Dec 21st first day of Christmas holiday, Feb 8th first reopening day
    
    #Add in mods for workplaces in Hastings Prince Edward (Regions 8,9), Kingston, Frontenac and Lennox & Addington (Regions 6,7), and Renfrew County (Region 38) which reopen 6 days earlier on Feb 10
    inClen_w2[c(6,7,8,9,38)]<-(inClen_w2[c(6,7,8,9,38)]-workstaggergap1)
    #Add in mod for workplaces in York (14) reopening 6 days later on Feb 22
    inClen_w2[c(14)]<-(inClen_w2[c(14)]+workstaggergap1)
    #Add in mods for workplaces in Toronto (15), Peel (16), North Bay - Parry Sound (39,40) reopening 20 days later on March 8th
    inClen_w2[c(15,16,39,40)]<-(inClen_w2[c(15,16,39,40)]+workstaggergap2)
    #Add in mods for schools in Toronto (15), Peel (16), and York (14) reopening 8 days later on feb 16th
    inClen_s2[c(14,15,16)]<-(inClen_s2[c(14,15,16)]+schoolstaggergap)
    
    workclosuregap1=14 #Gap between day total cases >= 50 (Mar 10th) and day workplaces closed (March 25th)
    schoolclosuregap1=3 #Gap between day total cases >= 50 (Mar 10th) and day schools closed for March Break (March 14th)
    workclosuregap2=290 #Gap between day total cases >= 50 (Mar 10th) and day workplaces re-closed (Dec 26th)
    schoolclosuregap2=285 #Gap between day total cases >= 50 (Mar 10th) and day schools re-closed (Dec 21st)
    schoolbreakgap=113; #Gap between day total cases >= 50 (Mar 10th) and day schools would have closed for Summer Holidays (July 1st)
    
    # Specification of response framework events (Look at `agedat[,c("Date", "ts50")]` to find corresponding tepi to each date, each day is 1 less than what's expected in tepi!)
    buffer = 14 # Time between announcement of a tier change and its implementation
    if(disable_response_framework){rf_start <- 1e8}
    
    omggap=144 #Gap between day total cases >= 50 (Mar 10th) and day we switch to wave 2 omega value (Aug 1st)
    Resurge_w=50 #Duration of additional closures 
    Resurge_s=50
    NP=ncol(parms)-nrow(parms)
    
    #Dummy start values for closures s.t. we can do counterfactuals for first wave
    inCstart_w=50;
    inCstart_s=50;
    if (cftype %in% c("workopen", "bothopen")) {inCstart_w=1e8;} #Set so workplaces never close
    
    years<- 1 
    
    #Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds respectively.
    #Nt tracks cumulative # positive cases (including those recovered)
    #C tracks # days since last closure, but in output msim converts all positive C entries into eps
    #VD indicates level of NPI adherence from individuals
    #Sick indicates all individuals who are in any of E, P, A, I (i.e. all those either exposed or currently infectious)
    S=c("S1", "S2", "S3", "S4", "S5"); E=c("E1", "E2","E3", "E4", "E5"); R=c("R1", "R2", "R3", "R4","R5");
    Da=c("A1","Ak1","SA1","SAk1", "A2","Ak2","SA2","SAk2", "A3","Ak3","SA3","SAk3", "A4","Ak4","SA4","SAk4", "A5","Ak5","SA5","SAk5")
    Di=c("I1","Ik1","SI1","SIk1", "I2","Ik2","SI2","SIk2", "I3","Ik3","SI3","SIk3", "I4","Ik4","SI4","SIk4", "I5","Ik5","SI5","SIk5")
    Tn=c("A1","SA1","I1","SI1","A2","SA2","I2","SI2","A3","SA3","I3","SI3","A4","SA4","I4","SI4","A5","SA5","I5","SI5")
    Tk=c("Ak1","SAk1","Ik1","SIk1","Ak2","SAk2","Ik2","SIk2","Ak3","SAk3","Ik3","SIk3","Ak4","SAk4","Ik4","SIk4","Ak5","SAk5","Ik5","SIk5")
    All1<-c("S1", "E1", "R1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); All2<-c("S2", "E2", "R2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); All3<-c("S3", "E3", "R3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); All4<-c("S4", "E4", "R4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); All5<-c("S5", "E5", "R5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
    sick1<-c("E1", "A1","Ak1","SA1","SAk1", "I1","Ik1","SI1","SIk1"); sick2<-c("E2", "A2","Ak2","SA2","SAk2", "I2","Ik2","SI2","SIk2"); sick3<-c("E3", "A3","Ak3","SA3","SAk3", "I3","Ik3","SI3","SIk3"); sick4<-c("E4", "A4","Ak4","SA4","SAk4", "I4","Ik4","SI4","SIk4"); sick5<-c("E5", "A5","Ak5","SA5","SAk5", "I5","Ik5","SI5","SIk5");
    symp1<-c("I1","Ik1","SI1","SIk1"); symp2<-c("I2","Ik2","SI2","SIk2"); symp3<-c("I3","Ik3","SI3","SIk3"); symp4<-c("I4","Ik4","SI4","SIk4"); symp5<-c("I5","Ik5","SI5","SIk5");
    Snm=c(S,E,Da,Di,R,"Nt","Cw", "Cs", "C", "Nt1","Nt2","Nt3","Nt4","Nt5", "Response_Framework", "VD");
    
    #these functions handle all the state transitions. Input x is a matrix where 1st half of columns are source states and 2nd half are destination states
    gpTransB=function(x,Prs,seed,nc=ncol(x)){ 
      xvp=cbind(as.vector(x[,1:(nc/2)]),as.vector(Prs))
      if(max(xvp[,1])==0) return(x); nz=(xvp[,1]*xvp[,2])>0; xvp[!nz,1]=0; 
      set.seed(seed); xvp[nz,1]=apply(matrix(xvp[nz,],ncol=2),1,function(y) rbinom(1,y[1],y[2])); 
      return(x+matrix(c(-xvp[,1],xvp[,1]),ncol=nc))
    }
    
    #gpTrans is a simplified version where one transition probability applies to all states. 
    #If recovery=TRUE have individuals from 1st columns of x all transitioning into the last column
    gpTrans=function(x,Pr,seed,Recovery=FALSE, nc=ncol(x)){
      xv=as.vector(x[,1:c(nc/2,nc-1)[1+Recovery]])
      if(max(xv)==0) return(x); set.seed(seed); xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr));
      if(Recovery){ Tr=matrix(xv,ncol=nc-1); return(x+cbind(-Tr,rowSums(Tr))); }; return(x+matrix(c(-xv,xv),ncol=nc));
    }
    
    #Transition probabilities for travel. Multinomial faster than many binomials. rs=TRUE returns just total # people going to each province
    reshfl2=function(x,M,seed,rs=TRUE,L=ncol(M),xnm=diag(x)){
      set.seed(seed); if(max(x)>0) xnm[,x>0]=apply(matrix(rbind(x,M)[,x>0],nrow=L+1), 2, function(y) rmultinom(1,y[1],y[-1])); 
      if(rs) return(rowSums(xnm)); return(xnm); 
    }
    
    #Modifier of travel matrix as people sick and/or tested positive less likely to travel by a proportion pstay
    Mstay=function(M,pstay,Mod=M*(-(diag(ncol(M))-1)*(1-pstay))) Mod+diag(1-colSums(Mod))
    
    meansd=function(x,wtR=1+x*0,wts=wtR/sum(wtR)){ xm=weighted.mean(x,wts); return(c(xm,sqrt(sum(wts*(x-xm)^2)))); }
    
    #Main function handling all within-day transitions over state space S
    m3iter=function(S,parms,timeinfo,seed,Ns=parms[,"N"]){
      tm <- timeinfo[1]
      tepi <- timeinfo[2]
      
      # Timestep messages
      if((!supercomputer_mode) & (tm %% printing_timestep == 0)){message("------------------ ", "Day: ", tm, " ------------------")}
      if((!supercomputer_mode) & tepi_message){ 
        message("------------------ ", "Date: ", as.Date("2020-03-10") + tepi + 1, " (tepi = ", tepi, ")", " ------------------")}
      
      ## Ramp up testing
      # For epi states P, A
      testPA1<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi1"]*timeinfo[2]));} 
      testPA2<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi2"]*timeinfo[2]));} 
      testPA3<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi3"]*timeinfo[2]));} 
      testPA4<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi4"]*timeinfo[2]));} 
      testPA5<- if(timeinfo[2]<0) {parms[,"taumax"]*0;} else {parms[1,"kappa"]*parms[,"taumax"]*(1-exp(-parms[1,"psi5"]*timeinfo[2]));} 
      # For epi state I
      testI1<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi1"]*timeinfo[2]);}
      testI2<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi2"]*timeinfo[2]);}
      testI3<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi3"]*timeinfo[2]);}
      testI4<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi4"]*timeinfo[2]);}
      testI5<- if(timeinfo[2]<0) {parms[,"tau0"];} else {parms[,"taumax"] - (parms[,"taumax"] - parms[,"tau0"])*exp(-parms[1,"psi5"]*timeinfo[2]);}
      
      if(testPA1 > testI1 || testPA2 > testI2 || testPA3 > testI3 || testPA4 > testI4 || testPA5 > testI5) {print("problem with testing rates"); invisible(readline(prompt="Press [enter] to continue"));}
      
      #Implement testing (test results coming back from previous day)
      S0k=S[,Tk]; S0k_1=S[,Tk[1:4]]; S0k_2=S[,Tk[5:8]]; S0k_3=S[,Tk[9:12]]; S0k_4=S[,Tk[13:16]]; S0k_5=S[,Tk[17:20]];
      S0n_1=S[,Tn[1:4]]; S0n_2=S[,Tn[5:8]]; S0n_3=S[,Tn[9:12]]; S0n_4=S[,Tn[13:16]]; S0n_5=S[,Tn[17:20]];
      S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],	as.vector(cbind(testPA1,testPA1,testI1,testI1,testPA2,testPA2,testI2,testI2,testPA3,testPA3,testI3,testI3,testPA4,testPA4,testI4,testI4,testPA5,testPA5,testI5,testI5)),seed+40); 
      
      #calculate pos, the vector of local active case prevalence
      S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k); pos=rowSums(S[,Tk])/Ns; posglobal=(sum(S[,Tk])/sum(Ns)); 
      S[,"Nt1"]=S[,"Nt1"]+rowSums(S[,Tk[1:4]]-S0k_1);
      S[,"Nt2"]=S[,"Nt2"]+rowSums(S[,Tk[5:8]]-S0k_2);
      S[,"Nt3"]=S[,"Nt3"]+rowSums(S[,Tk[9:12]]-S0k_3);
      S[,"Nt4"]=S[,"Nt4"]+rowSums(S[,Tk[13:16]]-S0k_4);
      S[,"Nt5"]=S[,"Nt5"]+rowSums(S[,Tk[17:20]]-S0k_5);
      
      #Disease progression: zeta is the fraction of people who never show symptoms (modeled implicitly)
      zeta=0.2;
      ### Symptomatic removals by age class
      S[,c(Di[1:4],R[1])]=gpTrans(S[,c(Di[1:4],R[1])],parms[1,"rho"],seed,TRUE); 
      S[,c(Di[5:8],R[2])]=gpTrans(S[,c(Di[5:8],R[2])],parms[1,"rho"],seed+1,TRUE); 
      S[,c(Di[9:12],R[3])]=gpTrans(S[,c(Di[9:12],R[3])],parms[1,"rho"],seed+2,TRUE); 
      S[,c(Di[13:16],R[4])]=gpTrans(S[,c(Di[13:16],R[4])],parms[1,"rho"],seed+3,TRUE); 
      S[,c(Di[17:20],R[5])]=gpTrans(S[,c(Di[17:20],R[5])],parms[1,"rho"],seed+4,TRUE); 
      ### Asymptomatic removals by age class
      S[,c(Da[1:4],R[1])]=gpTrans(S[,c(Da[1:4],R[1])],zeta*prod(parms[1,c("sig","rho")]),seed,TRUE); 
      S[,c(Da[5:8],R[2])]=gpTrans(S[,c(Da[5:8],R[2])],zeta*prod(parms[1,c("sig","rho")]),seed+1,TRUE); 
      S[,c(Da[9:12],R[3])]=gpTrans(S[,c(Da[9:12],R[3])],zeta*prod(parms[1,c("sig","rho")]),seed+2,TRUE); 
      S[,c(Da[13:16],R[4])]=gpTrans(S[,c(Da[13:16],R[4])],zeta*prod(parms[1,c("sig","rho")]),seed+3,TRUE); 
      S[,c(Da[17:20],R[5])]=gpTrans(S[,c(Da[17:20],R[5])],zeta*prod(parms[1,c("sig","rho")]),seed+4,TRUE); 
      ### Transition from pre-symptomatic to symptomatic by age class
      S[,c(Da[1:4],Di[1:4])]=gpTrans(S[,c(Da[1:4],Di[1:4])],parms[1,"sig"]*(1-zeta),seed+5); 
      S[,c(Da[5:8],Di[5:8])]=gpTrans(S[,c(Da[5:8],Di[5:8])],parms[1,"sig"]*(1-zeta),seed+6); 
      S[,c(Da[9:12],Di[9:12])]=gpTrans(S[,c(Da[9:12],Di[9:12])],parms[1,"sig"]*(1-zeta),seed+7); 
      S[,c(Da[13:16],Di[13:16])]=gpTrans(S[,c(Da[13:16],Di[13:16])],parms[1,"sig"]*(1-zeta),seed+8); 
      S[,c(Da[17:20],Di[17:20])]=gpTrans(S[,c(Da[17:20],Di[17:20])],parms[1,"sig"]*(1-zeta),seed+9); 
      ### Shift from exposed to asymptomatic by age class
      S[,c("E1","A1")]=gpTrans(S[,c("E1","A1")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+10); 
      S[,c("E2","A2")]=gpTrans(S[,c("E2","A2")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+11); 
      S[,c("E3","A3")]=gpTrans(S[,c("E3","A3")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+12); 
      S[,c("E4","A4")]=gpTrans(S[,c("E4","A4")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+13); 
      S[,c("E5","A5")]=gpTrans(S[,c("E5","A5")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+14); 
      ### Shift from exposed to superspreader asymptomatic by age class
      S[,c("E1","SA1")]=gpTrans(S[,c("E1","SA1")],parms[1,"alpha"]*parms[1,"s"],seed+15);
      S[,c("E2","SA2")]=gpTrans(S[,c("E2","SA2")],parms[1,"alpha"]*parms[1,"s"],seed+16);
      S[,c("E3","SA3")]=gpTrans(S[,c("E3","SA3")],parms[1,"alpha"]*parms[1,"s"],seed+17);
      S[,c("E4","SA4")]=gpTrans(S[,c("E4","SA4")],parms[1,"alpha"]*parms[1,"s"],seed+18);
      S[,c("E5","SA5")]=gpTrans(S[,c("E5","SA5")],parms[1,"alpha"]*parms[1,"s"],seed+19);
      
      closed_s <- parms[,"eps_s"];
      closed_w <- parms[,"eps_w"];
      closed_h <- parms[,"eps_h"];
      closed_o <- parms[,"eps_o"];
      
      # The hard coded shutdowns activate before the colour coded scheme gets implemented and also after the boxing day provincial shutdown
      if(!between(tepi, rf_start, rf_end)){
        # Conditional statements for opening/closing workplaces and schools
        if(cftype %in% c("workopen", "neitheropen")){
          if(timeinfo[3]==-1 || timeinfo[2]<schoolclosuregap1) closed_s = 0*parms[,"eps_s"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, schools not yet closed (Mar 10 - Mar 13)
          if(timeinfo[3]>=schoolclosuregap1 && timeinfo[3]<schoolclosuregap1+inClen_s1) closed_s = parms[,"eps_s"]; #Initial school closure started (March 14-Sept 7)
          if(timeinfo[3]>=schoolclosuregap1+inClen_s1 && timeinfo[3]< schoolclosuregap2  && reopeningtype=="restricted") closed_s =parms[1,"eps_s"]*parms[,"reopen_s"]; #Schools reopen with covid regs in place
          if(timeinfo[3]>=schoolclosuregap1+inClen_s1 && timeinfo[3]< schoolclosuregap2  && reopeningtype=="unrestricted") closed_s =0*parms[,"eps_s"]; #Schools reopen WITHOUT covid regs in place
          
          if (timeinfo[3]>=schoolclosuregap2) {
            if (reopeningtype=="restricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], parms[1,"eps_s"]*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
            if (reopeningtype=="unrestricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], 0*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
          }
        }
        
        if(cftype %in% c("schoolopen", "bothopen")){
          if(timeinfo[3]==-1 || timeinfo[2]<schoolbreakgap) closed_s = 0*parms[,"eps_s"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, schools not yet closed (Mar 10 - June 30)
          if(timeinfo[3]>=schoolbreakgap && timeinfo[3]<schoolbreakgap+inBlen_s) closed_s = parms[,"eps_s"]; #Initial summer school break started (July 1st-Sept 8th)
          if(timeinfo[3]>=schoolbreakgap+inBlen_s && timeinfo[3]< schoolclosuregap2 && reopeningtype=="restricted") closed_s =parms[1,"eps_s"]*parms[,"reopen_s"]; #Schools reopen with covid regs in place
          if(timeinfo[3]>=schoolbreakgap+inBlen_s && timeinfo[3]< schoolclosuregap2 && reopeningtype=="unrestricted") closed_s =0*parms[,"eps_s"]; #Schools reopen WITHOUT covid regs in place
          
          if (timeinfo[3]>=schoolclosuregap2) {
            if (reopeningtype=="restricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], parms[1,"eps_s"]*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
            if (reopeningtype=="unrestricted") {closed_s <- ifelse(timeinfo[3]<schoolclosuregap2+inClen_s2, parms[,"eps_s"], 0*parms[,"reopen_s"]);} #Christmas holiday to Boxing day lockdown to stay at home, with staggered reopening
          }
          
        }
        
        if(timeinfo[4]==-1 || timeinfo[2]<workclosuregap1) closed_w = 0*parms[,"eps_w"]; #We have not yet reached >=50 cases (prior to Mar 10) or >=50 cases reached, workplaces not yet closed (Mar 10 - Mar 24)
        if(timeinfo[4]>=workclosuregap1 && timeinfo[4]<workclosuregap1+inClen_w1) {closed_w = parms[,"eps_w"];}#Initial workplace closure started (March 25-June 11th)
        if(timeinfo[4]>=workclosuregap1+inClen_w1 && timeinfo[4]< workclosuregap2 && reopeningtype=="restricted") closed_w =parms[1,"eps_w"]*parms[,"reopen_w"]; #Workplaces reopen with covid regs in place
        if(timeinfo[4]>=workclosuregap1+inClen_w1 && timeinfo[4]< workclosuregap2 && reopeningtype=="unrestricted") closed_w = 0*parms[,"eps_w"]; #Workplaces reopen WITHOUT covid regs in place
        
        if (timeinfo[4]>=workclosuregap2) {
          if (reopeningtype=="restricted") {closed_w <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, parms[,"eps_w"], parms[1,"eps_w"]*parms[,"reopen_w"]);} #Boxing day lockdown to stay at home, with staggered reopening
          if (reopeningtype=="unrestricted") {closed_w <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, parms[,"eps_w"], 0*parms[,"reopen_w"]);} #Boxing day lockdown to stay at home, with staggered reopening
        }
      }
      
      S[,"Cs"] = closed_s;
      S[,"Cw"] = closed_w;
      
      #Make modified travel matrices for people feeling sick and/or tested positive, then implement travel
      M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed_w)*(1-closed_s)*(1-closed_h)*(1-closed_o); diag(Mc)=diag(Mc)+1-colSums(Mc);
      ### Revamping McA to include age specific travel rates, first of each pair is for old/young, second for middle
      McA=abind(Mstay(Mc,parms[1,"atravel"]), Mc, Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"eta"])), Mstay(Mc,parms[1,"eta"]), Mstay(Mc,1-(1-parms[1,"atravel"])*(1-parms[1,"r"])), Mstay(Mc,parms[1,"r"]), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])*(1-parms[1,"atravel"])), Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])), along=3);
      #Implement travel
      ### Do Ss for however many age classes you have, added atravel for age-class specific travel
      Ss=abind(reshfl2(S[,"S1"],Mstay(Mc,parms[1,"atravel"]),seed+20,FALSE),reshfl2(S[,"S2"],Mc,seed+21,FALSE),reshfl2(S[,"S3"],Mc,seed+22,FALSE),reshfl2(S[,"S4"],Mstay(Mc,parms[1,"atravel"]),seed+23,FALSE),reshfl2(S[,"S5"],Mstay(Mc,parms[1,"atravel"]),seed+24,FALSE),along=3);
      Rearr=apply(rbind(seed+(25:34), c(c(1,2,2,1,1),c(1,2,2,1,1),c(1,3,1,3,2,4,2,4,2,4,2,4,1,3,1,3,1,3,1,3),c(5,7,5,7,6,8,6,8,6,8,6,8,5,7,5,7,5,7,5,7)), S[,c(R,E,Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
      
      #Age specific contacts, rows are age classes, columns are their contact age classes
      ageSpecifics_w<-cmat[,,1] ## age specific contacts for work
      ageSpecifics_s<-cmat[,,2] ## age specific contacts for school
      ageSpecifics_h<-cmat[,,3] ## age specific contacts for home
      ageSpecifics_o<-cmat[,,4] ## age specific contacts for other
      
      #Seasonal forcing
      scomp<-1+parms[,"B"]*cos((2*pi/365)*(timeinfo[1] + parms[,"phi"]))
      
      S[, "VD"]<-(1-exp(-parms[1,"omg"]*pos));
      if (timeinfo[4]>=workclosuregap2 + 20) {
        S[, "VD"] <- ifelse(timeinfo[4]<workclosuregap2+inClen_w2, (1-exp(-(parms[1,"omg"]*pos + parms[1,"boost"]))), (1-exp(-parms[1,"omg"]*pos))) #Stay at home, with staggered reopening
      }
      
      if (vdtype=="vdOFF") {S[, "VD"]<-0*pos;}
      
      #Base infection probabilities
      Infect1 = scomp*parms[,"a1"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[1,]) + (1-closed_s)%*%t(ageSpecifics_s[1,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[1,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[1,]))
      Infect2 = scomp*parms[,"a2"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[2,]) + (1-closed_s)%*%t(ageSpecifics_s[2,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[2,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[2,]))
      Infect3 = scomp*parms[,"a3"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[3,]) + (1-closed_s)%*%t(ageSpecifics_s[3,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[3,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[3,]))
      Infect4 = scomp*parms[,"a4"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[4,]) + (1-closed_s)%*%t(ageSpecifics_s[4,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[4,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[4,]))
      Infect5 = scomp*parms[,"a5"]*parms[,"beta"]*((1-closed_w)%*%t(ageSpecifics_w[5,]) + (1-closed_s)%*%t(ageSpecifics_s[5,]) + (1-parms[1,"eps_h"]*S[, "VD"])%*%t(ageSpecifics_h[5,]) + (1-parms[1,"eps_o"]*S[, "VD"])%*%t(ageSpecifics_o[5,]))
      
      #Class-specific infection modifiers
      modK=1-parms[1,"eta"]; modS=(1/parms[1,"s"])-1; 
      # modA=cbind(Infect,modK*Infect,modS*Infect,modS*modK*Infect); # Replaced with modA1,2,...
      
      #Introduce age specific modA 3-D arrays for not superspreader or known, not supersrpeader but known, superspreader not known, superspreader and known
      modA1=abind(Infect1,modK*Infect1,modS*Infect1,modS*modK*Infect1,along=3);
      modA2=abind(Infect2,modK*Infect2,modS*Infect2,modS*modK*Infect2,along=3);
      modA3=abind(Infect3,modK*Infect3,modS*Infect3,modS*modK*Infect3,along=3);
      modA4=abind(Infect4,modK*Infect4,modS*Infect4,modS*modK*Infect4,along=3);
      modA5=abind(Infect5,modK*Infect5,modS*Infect5,modS*modK*Infect5,along=3);
      
      #Flatten into 2-D arrays where order is based on age class to match Rearr which is ordered as (all pre/asympt by age class, all sympt by age class)
      modAbind1<-cbind(modA1[,1,],modA1[,2,],modA1[,3,],modA1[,4,],modA1[,5,]);
      modAbind2<-cbind(modA2[,1,],modA2[,2,],modA2[,3,],modA2[,4,],modA2[,5,]);
      modAbind3<-cbind(modA3[,1,],modA3[,2,],modA3[,3,],modA3[,4,],modA3[,5,]);
      modAbind4<-cbind(modA4[,1,],modA4[,2,],modA4[,3,],modA4[,4,],modA4[,5,]);
      modAbind5<-cbind(modA5[,1,],modA5[,2,],modA5[,3,],modA5[,4,],modA5[,5,]);
      
      # Overall infection Pr is then 1-Pr(avoid infection by anyone)
      Infects1=1 - apply((1-cbind(modAbind1, 2*modAbind1))^Rearr[,-(1:10)], 1, prod)
      Infects2=1 - apply((1-cbind(modAbind2, 2*modAbind2))^Rearr[,-(1:10)], 1, prod)
      Infects3=1 - apply((1-cbind(modAbind3, 2*modAbind3))^Rearr[,-(1:10)], 1, prod)
      Infects4=1 - apply((1-cbind(modAbind4, 2*modAbind4))^Rearr[,-(1:10)], 1, prod)
      Infects5=1 - apply((1-cbind(modAbind5, 2*modAbind5))^Rearr[,-(1:10)], 1, prod)
      
      #Implement infection and move susceptibles and newly exposeds back to home county 
      #For each Ss[,,i], entires in each column are individuals from the same region, with the rows showing if/where they travelled (so row sums give current total people in a region including visitors)
      S[,c("S1","E1")]=cbind(0,S[,"E1"]) + colSums(gpTransB(cbind(Ss[,,1],0*Ss[,,1]),round(Infects1,10),seed+35))
      S[,c("S2","E2")]=cbind(0,S[,"E2"]) + colSums(gpTransB(cbind(Ss[,,2],0*Ss[,,2]),round(Infects2,10),seed+36))
      S[,c("S3","E3")]=cbind(0,S[,"E3"]) + colSums(gpTransB(cbind(Ss[,,3],0*Ss[,,3]),round(Infects3,10),seed+37))
      S[,c("S4","E4")]=cbind(0,S[,"E4"]) + colSums(gpTransB(cbind(Ss[,,4],0*Ss[,,4]),round(Infects4,10),seed+38))
      S[,c("S5","E5")]=cbind(0,S[,"E5"]) + colSums(gpTransB(cbind(Ss[,,5],0*Ss[,,5]),round(Infects5,10),seed+39))
      
      return(S)
    }; FUN=m3iter
    
    #Implement initial closure and changes in testing probability
    closeinappl=function(parms,TS,tm=dim(TS)[3],delayInit=10){
      
      Nta=colSums(t(t(TS[,"Nt",])));
      
      omgs=cbind(parms[,"omg"],0);
      
      rampdays1=tm-(which(Nta>=inCstart)[1]+158-1) #Number of days past the beginning of the ramp-down
      parms[,"omg"]=pmax(omgs[,1]*exp(-parms[1,"zeta"]*rampdays1),omgs[,2]) #Ramp-down the omega value from initial to 2nd wave value
      if( (max(Nta)>=inCstart) & (tm-which(Nta>=inCstart)[1])<(158-1) ) { parms[,"omg"]=omgs[,1];} #NPI adherence with risk perception coeff omg_0 before ramp down begins
      if((max(Nta)<inCstart)) parms[,"omg"]=0; #No NPI adherence before start of pandemic (as inCstart is March 10th when total cases>=50)
      
      return(parms);
    }
    
    #Pulls info on what timestep it is
    gettime=function(TS,tm=dim(TS)[3]){
      
      tsNt=colSums(t(t(TS[,"Nt",])));
      if(max(tsNt)>=inCstart) {tepi<-length(tsNt[tsNt>=inCstart]) -1} else {tepi<--1} #Calc days into epidemic (since >=50 cases)
      if(max(tsNt)>=inCstart_s) {tstart_s<-length(tsNt[tsNt>=inCstart_s]) -1} else {tstart_s<--1} #Calc days past trigger date for school closures (==tepi for no counterfactuals)
      if(max(tsNt)>=inCstart_w) {tstart_w<-length(tsNt[tsNt>=inCstart_w]) -1} else {tstart_w<--1} #Calc days past trigger date for work closures (==tepi for no counterfactuals)
      
      timeinfo=c(tm, tepi, tstart_s, tstart_w)
      return(timeinfo);
    }
    
    ### vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Response Framework Code vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ###
    # Code for the color-coded regional restrictions (Ontario's COVID-19 Response Framework)
    # This includes the 5 zones described at https://www.ontario.ca/page/covid-19-response-framework-keeping-ontario-safe-and-open 
    # (above the website has since been taken down, use WayBack Machine to find an archived version)
    # and includes the "Shutdown" zone which is a province wide restriction that is more strict than "Lockdown"
    # zones = c("Prevent", "Protect", "Restrict", "Control", "Lockdown", "Shutdown")
    
    ### Constants
    # Weekly Incident Rate Thresholds
    if("wir" %in% indicators_active){
      if(!supercomputer_mode){message("Weekly incidence rate active")}
      wir_shutdown_threshold <- mu * 100 # Guidelines not clear, determined from eyeballing data
      wir_lockdown_threshold <- mu * 100 # Guidelines not clear, determined from eyeballing data
      wir_control_threshold <- mu * 40
      wir_restrict_threshold <- mu * 25
      wir_protect_threshold <- mu * 10
    } else {
      wir_shutdown_threshold <- 1e6
      wir_lockdown_threshold <- 1e6
      wir_control_threshold <- 1e6
      wir_restrict_threshold <- 1e6
      wir_protect_threshold <- 1e6
    }
    
    # Percent Positivity Thresholds
    if("ppos" %in% indicators_active){
      if(!supercomputer_mode){message("Percent positivity active")}
      ppos_shutdown_threshold <- mu * 6 # Guidelines not clear, determined from eyeballing data
      ppos_lockdown_threshold <- mu * 6 # Guidelines not clear, determined from eyeballing data
      ppos_control_threshold <- mu * 2.5
      ppos_restrict_threshold <- mu * 1.3
      ppos_protect_threshold <- mu * 0.5
    } else {
      ppos_shutdown_threshold <- 1e6
      ppos_lockdown_threshold <- 1e6
      ppos_control_threshold <- 1e6
      ppos_restrict_threshold <- 1e6
      ppos_protect_threshold <- 1e6
    }
    
    # Reproduction Number Thresholds
    if("Rt" %in% indicators_active){
      if(!supercomputer_mode){message("Reproductive number active")}
      Rt_shutdown_threshold <- mu * 1.4  # Guidelines not clear, determined from eyeballing data
      Rt_lockdown_threshold <- mu * 1.4 # Guidelines not clear, determined from eyeballing data
      Rt_control_threshold <- mu * 1.2
      Rt_restrict_threshold <- mu * 1
      Rt_protect_threshold <- mu * 0.95
    } else {
      Rt_shutdown_threshold <- 1e6
      Rt_lockdown_threshold <- 1e6
      Rt_control_threshold <- 1e6
      Rt_restrict_threshold <- 1e6
      Rt_protect_threshold <- 1e6
    }
    
    # Initialize the epsilon and alpha values 
    alpha_w1 <- parms[1, "alpha_w1"]
    alpha_w2 <- parms[1, "alpha_w2"]
    alpha_w3 <- parms[1, "alpha_w3"]
    alpha_w4 <- parms[1, "alpha_w4"]
    eps_w_shutdown <- parms[1, "eps_w"]
    eps_w_lockdown <- eps_w_shutdown
    eps_w_control <- alpha_w1 * eps_w_lockdown
    eps_w_restrict <- alpha_w2 * eps_w_control
    eps_w_protect <- alpha_w3 * eps_w_restrict
    eps_w_prevent <- alpha_w4 * eps_w_protect
    
    alpha_s1 <- parms[1, "alpha_s1"]
    alpha_s2 <- parms[1, "alpha_s2"]
    alpha_s3 <- parms[1, "alpha_s3"]
    alpha_s4 <- parms[1, "alpha_s4"]
    eps_s_shutdown <- parms[1, "eps_s"]
    eps_s_lockdown <- eps_s_shutdown
    eps_s_control <- alpha_s1 * eps_s_lockdown
    eps_s_restrict <- alpha_s2 * eps_s_control
    eps_s_protect <- alpha_s3 * eps_s_restrict
    eps_s_prevent <- alpha_s4 * eps_s_protect
    
    alpha_h1 <- parms[1, "alpha_h1"]
    alpha_h2 <- parms[1, "alpha_h2"]
    alpha_h3 <- parms[1, "alpha_h3"]
    alpha_h4 <- parms[1, "alpha_h4"]
    eps_h_shutdown <- parms[1, "eps_h"]
    eps_h_lockdown <- eps_h_shutdown
    eps_h_control <- alpha_h1 * eps_h_lockdown
    eps_h_restrict <- alpha_h2 * eps_h_control
    eps_h_protect <- alpha_h3 * eps_h_restrict
    eps_h_prevent <- alpha_h4 * eps_h_protect
    
    alpha_o1 <- parms[1, "alpha_o1"]
    alpha_o2 <- parms[1, "alpha_o2"]
    alpha_o3 <- parms[1, "alpha_o3"]
    alpha_o4 <- parms[1, "alpha_o4"]
    eps_o_shutdown <- parms[1, "eps_o"]
    eps_o_lockdown <- eps_o_shutdown
    eps_o_control <- alpha_o1 * eps_o_lockdown
    eps_o_restrict <- alpha_o2 * eps_o_control
    eps_o_protect <- alpha_o3 * eps_o_restrict
    eps_o_prevent <- alpha_o4 * eps_o_protect
    
    eps_df <- data.frame(work = c(eps_w_prevent, eps_w_protect, eps_w_restrict, eps_w_control, eps_w_lockdown, eps_w_shutdown),
                         school = c(eps_s_prevent, eps_s_protect, eps_s_restrict, eps_s_control, eps_s_lockdown, eps_s_shutdown),
                         home = c(eps_h_prevent, eps_h_protect, eps_h_restrict, eps_h_control, eps_h_lockdown, eps_h_shutdown),
                         other = c(eps_o_prevent, eps_o_protect, eps_o_restrict, eps_o_control, eps_o_lockdown, eps_o_shutdown),
                         row.names = c("Prevent", "Protect", "Restrict", "Control", "Lockdown", "Shutdown"))
    
    ### vvvv Functions to implement Ontario's colour-coded lockdown scheme vvvv ###
    # tier_to_number: converts strings of a response_framework to an int.
    tier_to_number <- function(response_framework, timeinfo){
      # The response_framework input for this function is only for one day!!! (A 1-dimensional vector)
      
      today = timeinfo[1]
      tepi = timeinfo[2]
      new_framework <- ifelse(response_framework == "Prevent", 1,
                              ifelse(response_framework == "Protect", 2,
                                     ifelse(response_framework == "Restrict", 3,
                                            ifelse(response_framework == "Control", 4,
                                                   ifelse(response_framework == "Lockdown", 5,
                                                          ifelse(response_framework == "Shutdown", 6,
                                                                 -1
                                                          )
                                                   )
                                            )
                                     )
                              )
      )
      
      # Ensure that the value is 0 for all days before the response framework got implemented
      if (tepi < rf_start){
        new_framework = rep(0, length(new_framework))
      } 
      
      return(new_framework)
    }
    
    
    # number_of_regions: Returns the number of (census) regions in the simulation
    number_of_regions <- function(region_id = regionid){
      return(dim(region_id)[1])
    } 
    
    
    # create_response_framework: Creates an array with rows representing days since start of simulation (tm) and columns representing regions.
    # 2D-array <- (string, int, 2D-array)
    create_response_framework <- function(init_zone, buffer = 14, region_id = regionid){
      # init_zone is the zone that all regions are in at the start of the simulation. Has to be an element of zones.
      # buffer is the time between announcing a switch to a new zone and implementing the zone's restrictions
      
      census_regions = region_id$census.region
      num_of_census <- number_of_regions(region_id)
      response_framework = array(init_zone, dim = c(1+buffer, num_of_census))
      # print(response_framework)
      
      return(response_framework)
    }
    
    # weekly_new_cases: Outputs the new cases/week for a specific day
    # int <- (3D-dataframe, int)
    weekly_new_cases <- function(time, TS){
      one_week_ago = ifelse(time > 7, time - 7 + 1, 1) # slicing is inclusive in R so I have to add one to output a size 7 array.
      result = TS[,"Nt", time] - TS[,"Nt", one_week_ago]
      return(result)
    }
    
    # weekly_incidence_rate_regional: Updates the wir_time array with the weekly incidence rate per 100 000 for each region (columns) over time (rows). 
    # vector <- (3D-dataframe, 2D-array, 2D-array)
    weekly_incidence_rate_regional <- function(wir_time, TS, region_id = regionid){
      # Does nothing if wir was not activated. 
      if(! "wir" %in% indicators_active){
        # print("wir not active")
        return(wir_time)
      }
      
      # Define key variables to check conditions on
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        wir_time = rbind(wir_time, -1)
        return(wir_time) # Modify this so that it's the right size.
      }
      
      # Define the time two weeks ago
      if (today > 14){ # If less than two weeks has passed since the start, then set two_weeks_ago <- 1
        two_weeks_ago <- today - 14 + 1 # slicing is inclusive in R so I have to add one to output a size 14 array.
      } else {
        two_weeks_ago <- 1
      }
      
      # Calculates the weekly new cases for each census region for the past 2 weeks
      weekly_new_cases_2wks = weekly_new_cases(two_weeks_ago:today, TS) 
      
      # Convert census region cases to PHU region cases
      PHU_nums = region_id$phunum # Array to store each census region's PHU number
      PHU_populations = aggregate(x = region_id[,"pop2016"], 
                                  by = list(region_id$phunum), 
                                  FUN = sum)$x # Array to store each PHU's population
      weekly_new_cases_2wks_PHU = aggregate(x = weekly_new_cases_2wks,
                                            by = list(region_id$phunum),
                                            FUN = sum) # Calculate weekly new cases for each PHU for the past 2 weeks
      weekly_new_cases_2wks_PHU = subset(weekly_new_cases_2wks_PHU, select = -c(Group.1)) # Get rid of the unnecessary column that was added from the aggregate function
      weekly_new_cases_PHU = as.matrix(apply(weekly_new_cases_2wks_PHU, 1, mean)) # Calculate the daily average of the weekly new cases for each PHU based on data from the past 2 weeks
      
      # Calculate the weekly incidence rate for each PHU, then expand the values to be the same for each census region with the same PHU.
      weekly_incidence_rate_PHU_today = 100000*weekly_new_cases_PHU/PHU_populations # Calculate the weekly incidence rate for each PHU
      weekly_incidence_rate_regional_today = weekly_incidence_rate_PHU_today[PHU_nums] # Spread the weekly incidence rate to census regions with the same PHU
      
      # Append new weekly cases array to the wir_time array
      wir_time <- rbind(wir_time, weekly_incidence_rate_regional_today)
      return(wir_time)
    }
    
    
    # weekly_incidence_rate_provincial: Outputs the weekly incidence rate of the province. 
    # float <- (3D-array, 2D-array)
    weekly_incidence_rate_provincial <- function(TS, region_id = regionid){
      # Define key variables to check conditions on
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        return(-1) # Modify this so that it's the right size.
      }
      
      if (today > 14){ # If less than two weeks has passed since the start, then set two_weeks_ago <- 1
        two_weeks_ago <- today - 14 + 1 # slicing is inclusive in R so I have to add one to output a size 14 array.
      } else {
        two_weeks_ago <- 1
      }
      
      weekly_new_cases_2wks = colSums(weekly_new_cases(two_weeks_ago:today, TS)) # Calculates the weekly new cases for each census region for the past 2 weeks
      provincial_population = sum(region_id$pop2016) 
      weekly_new_cases = mean(weekly_new_cases_2wks)
      
      weekly_incidence_rate_provincial_today <- 100000*weekly_new_cases/provincial_population
      return(weekly_incidence_rate_provincial_today)
    }
    
    
    # percent_positivitiy_regional: Calculates the % positivity for every census region for each day.
    percent_positivitiy_regional <- function(ppos_time, TS, testing_volumes, region_id = regionid){
      # Does nothing if ppos was not activated.
      if(! "ppos" %in% indicators_active){
        # print("ppos not active")
        return(ppos_time)
      }
      
      # Define key variables to check conditions on
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      epi_sim_diff = today - tepi
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        ppos_time = rbind(ppos_time, -1)
        return(ppos_time) # Modify this so that it's the right size.
      }
      
      # If tepi is before the data records, then set ppos to -1 for that day.
      min_tepi = testing_volumes$ts50[1] 
      if(tepi <= min_tepi){
        ppos_time = rbind(ppos_time, -1)
        return(ppos_time)
      }
      
      if (today > 14){ # If less than two weeks has passed since the start, then set two_weeks_ago <- 1
        two_weeks_ago <- today - 14 + 1 # slicing is inclusive in R so I have to add one to output a size 14 array.
      } else {
        two_weeks_ago <- 1
      }
      
      if (tepi <= min_tepi + 14) { # Condition for days slightly after data starts
        two_weeks_ago = min_tepi + epi_sim_diff + 1
      }
      
      # Define days for indexing
      todays = two_weeks_ago:today
      yesterdays = ifelse(todays > min_tepi + epi_sim_diff, todays - 1, todays)
      new_pos_2wks_census = TS[,"Nt", todays] - TS[,"Nt", yesterdays] # Calculate the new positive cases for the past two weeks for each census region. 
      
      # Convert census region cases to PHU region cases
      PHU_populations = aggregate(x = region_id[,"pop2016"], 
                                  by = list(region_id$phunum), 
                                  FUN = sum)$x # Array to store each PHU's population
      new_pos_2wks_PHU = aggregate(x = new_pos_2wks_census,
                                   by = list(region_id$phunum),
                                   FUN = sum) # Calculate new positive cases for each PHU for the past 2 weeks
      
      # Get rid of the unnecessary column that was added from the aggregate function and transpose
      new_pos_2wks_PHU = subset(new_pos_2wks_PHU, select = -c(Group.1)) 
      new_pos_2wks_PHU = t(new_pos_2wks_PHU)
      
      # Extract testing data from 1 day before relative to the two weeks that we're examining
      yesterday_tests = yesterdays - epi_sim_diff
      num_of_PHUs = length(PHU_populations)
      PHU_names = paste("phu.", c(1:num_of_PHUs), sep = "")
      testing_2wks = testing_volumes[which(testing_volumes$ts50 %in% yesterday_tests), PHU_names]
      
      # Calculate the percent positivity 
      ppos_PHU = 100*new_pos_2wks_PHU/testing_2wks
      ppos_PHU = as.matrix(apply(ppos_PHU, 2, mean)) # Take the average over the two weeks
      
      PHU_nums = region_id$phunum # Array to store each census region's PHU number
      ppos_census = ppos_PHU[PHU_nums] # Convert PHU ppos into census regions
      
      # Append new ppos values to ppos_time
      ppos_time = rbind(ppos_time, ppos_census) 
      return(ppos_time)
    }
    
    
    # percent_positivitiy_provincial: Calculates the % positivity for the province.
    percent_positivitiy_provincial <- function(TS, testing_volumes, region_id = regionid){
      # Define key variables to check conditions on
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      epi_sim_diff = today - tepi
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        return(-1) # Modify this so that it's the right size.
      }
      
      # If tepi is before the data records, then set ppos to -1 for that day.
      min_tepi = testing_volumes$ts50[1] 
      if(tepi <= min_tepi){
        return(-1)
      }
      
      if (today > 14){ # If less than two weeks has passed since the start, then set two_weeks_ago <- 1
        two_weeks_ago <- today - 14 + 1 # slicing is inclusive in R so I have to add one to output a size 14 array.
      } else {
        two_weeks_ago <- 1
      }
      
      if (tepi <= min_tepi + 14) { # Condition for days slightly after data starts
        two_weeks_ago = min_tepi + epi_sim_diff + 1
      }
      
      # Define days for indexing
      todays = two_weeks_ago:today
      yesterdays = ifelse(todays > min_tepi + epi_sim_diff, todays - 1, todays)
      
      # Define new cases for the past two weeks for each census region
      new_pos_2wks_census = (TS[,"Nt", todays] - TS[,"Nt", yesterdays])
      
      if(is.null(dim(new_pos_2wks_census))){ # Base case for when new_pos_2wks_census has only one dimension (one row)
        new_pos_2wks = sum(new_pos_2wks_census)
      } else {
        new_pos_2wks = colSums(new_pos_2wks_census)
      }
      
      # Extract testing data from 1 day before relative to the two weeks that we're examining
      yesterday_tests = yesterdays - epi_sim_diff
      provincial_population = sum(region_id$pop2016)
      num_of_PHUs = length(unique(region_id$phunum))
      PHU_names = paste("phu.", c(1:num_of_PHUs), sep = "")
      testing_2wks = testing_volumes[which(testing_volumes$ts50 %in% yesterday_tests), PHU_names]
      testing_2wks = rowSums(testing_2wks)
      
      # Calculate the percent positivity and it's average over two weeks
      perc_pos = 100*new_pos_2wks/testing_2wks
      perc_pos = mean(perc_pos)
      return(perc_pos)
    }
    
    
    # effective_reproduction_number: Calculates the Rt value given a new_cases array
    effective_reproduction_number <- function(new_cases, 
                                              initial_active_cases,
                                              rolling_window = 7){
      
      # Define parameters for the serial interval
      # https://www.publichealthontario.ca/-/media/data-files/covid-19-data-tool-technical-notes.pdf?la=en
      mu = 4.5
      sig = 2.5
      k = seq(0,20)
      si_distr = discr_si(k,mu,sig)
      
      # Define the rolling window
      t_end = length(new_cases)
      t_start = ifelse(t_end <= rolling_window, 2, t_end - rolling_window + 1)
      
      # Define the incid array
      local = new_cases
      imported = array(0, dim = c(t_end))
      imported[1] = initial_active_cases + local[1]
      local[1] = 0
      incid = data.frame(local = local, imported = imported)
      
      # -1 indicates that there wasn't enough data for the calculation to be accurate.
      not_enough_data = -1
      if(t_end < 5){
        return(not_enough_data)
      }
      
      # Calculate the Rt value, return -1 if it's too early in the epidemic
      output_R = tryCatch(
        {estimate_R(
          incid = incid,
          method = "non_parametric_si",
          config = make_config(list(
            # incid = new_cases,
            # method = method,
            si_distr = si_distr,
            # mean_si = mu,
            # std_si = sig,
            t_start = t_start,
            t_end = t_end)
          )
        )},
        warning = function(w){
          # print(w)
          if(w == warning("You're estimating R too early in the epidemic to get the desired
            posterior CV.")){
            # print(w)
            # print("too early")
            return(not_enough_data)
          }
          else {
            print(w)
            print("another warning has been triggered...")
            return(-2)
          }
        }
      )
      
      # Checks on the output
      if(is.numeric(output_R)){
        Rt = not_enough_data
      } else {
        Rt = output_R$R$`Median(R)`
      }
      
      if(is.na(Rt)){
        print("Rt contains NA")
        print("placeholder")
        print(output_R)
        print("hi")
        # Rt = not_enough_data
      }
      
      return(Rt)
    }
    
    
    # Rt_regional: Calculates Rt for every census region for the latest day.
    Rt_regional <- function(Rt_time, TS, region_id = regionid){
      # Does nothing if Rt was not activated.
      if(! "Rt" %in% indicators_active){
        # print("Rt not active")
        return(Rt_time)
      }
      
      # Define the rolling window and other key variables to check conditions on
      rolling_window = 7
      rw_buffer = 2*rolling_window
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        Rt_time = rbind(Rt_time, -1)
        return(Rt_time) 
      }
      
      if(tepi < 3){ # If not enough time has passed since the start of the epidemic, then output -1
        col_len = dim(Rt_time)[2]
        Rt_time = rbind(Rt_time, rep(-1,col_len))
        return(Rt_time)
      }
      
      # Calculate dates
      epi_sim_diff = today - tepi # Difference between simulation start and start of epidemic (Mar 10, 2020)
      epi_start = as.Date("2020-03-10", format = "%Y-%m-%d") # start of epidemic (Mar 10, 2020)
      end_date = epi_start + tepi # Most recent date of data reporting
      two_weeks_ago = ifelse(tepi >= rw_buffer, today - rw_buffer, epi_sim_diff) 
      start_date = end_date - (today - two_weeks_ago) 
      dates = seq(from = start_date, to = end_date, by = "days")
      
      # Indices used to extract daily cases
      todays = c(two_weeks_ago:today)
      yesterdays = todays - 1
      
      # Calculate daily cases for each census region
      new_cases_census = TS[,"Nt",todays] - TS[,"Nt",yesterdays]
      new_cases_census = as.data.frame(t(new_cases_census))
      num_of_census = number_of_regions(regionid)
      census_nums = paste("census.", c(1:num_of_census), sep = "")
      colnames(new_cases_census) = census_nums
      reported_cases_census = as.data.frame(cbind(date = dates, new_cases_census))
      
      # Convert census region cases to PHU region cases
      PHU_nums = region_id$phunum # Array to store each census region's PHU number
      PHU_populations = aggregate(x = region_id[,"pop2016"],
                                  by = list(region_id$phunum),
                                  FUN = sum)$x # Array to store each PHU's population
      reported_cases_PHU = aggregate(x = t(reported_cases_census[,census_nums]),
                                     by = list(region_id$phunum),
                                     FUN = sum) # Calculate weekly new cases for each PHU for the past 2 weeks
      reported_cases_PHU = subset(reported_cases_PHU, select = -c(Group.1)) # Get rid of the unnecessary column that was added from the aggregate function
      reported_cases_PHU = data.frame((reported_cases_PHU))
      
      # Set initial active cases to be the active cases (Tk) from two weeks ago
      initial_active_cases_census = rowSums(TS[,Tk,two_weeks_ago])
      initial_active_cases_PHU = aggregate(x = initial_active_cases_census,
                                           by = list(region_id$phunum),
                                           FUN = sum)$x
      
      # A wrapper function that makes applying effective_reproduction_number() to the big dataframe easier
      ern_apply <- function(cases){
        cols = length(cases)
        initial_active_cases_PHU = cases[cols]
        reported_cases_PHU = cases[-cols]
        Rts = effective_reproduction_number(reported_cases_PHU, initial_active_cases_PHU)
        return(Rts)
      }
      
      # The big dataframe I was talking about earlier
      cases = cbind(new_cases = reported_cases_PHU, initial_active_cases = initial_active_cases_PHU)
      Rt_PHU = apply(X = cases,
                     FUN = ern_apply,
                     MARGIN = 1
      ) # Calculate Rt for each PHU
      
      
      Rt_census = Rt_PHU[PHU_nums] # Convert Rt back to census regions
      Rt_time = rbind(Rt_time, Rt_census) # Append the new Rt values to the Rt_time array
      colnames(Rt_time) = census_nums
      return(Rt_time)
    }
    
    
    # Rt_provincial: Calculates Rt for the entire province for the latest day.
    Rt_provincial <- function(TS, region_id = regionid){
      # Define key variables to check conditions on
      today = gettime(TS)[1]
      tepi = gettime(TS)[2]
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present 
      
      # If response framework hasn't activated yet, don't go through with the calculation.
      if(!between(implementation_day, rf_start, rf_end)){
        return(-1) # Modify this so that it's the right size.
      }
      
      # If not enough time has passed to allow for accurate estimation of Rt, 
      # then return -1
      not_enough_data = -1
      if(tepi < 3){ 
        return(not_enough_data)
      }
      
      # Set parameters to truncate the data
      rolling_window = 7
      rw_buffer = 2*rolling_window
      
      # Calculate dates starting from Mar 10, 2020
      epi_sim_diff = today - tepi
      epi_start = as.Date("2020-03-10", format = "%Y-%m-%d")
      end_date = epi_start + tepi
      two_weeks_ago = ifelse(tepi >= rw_buffer, today - rw_buffer, epi_sim_diff)
      start_date = end_date - (today - two_weeks_ago)
      dates = seq(from = start_date, to = end_date, by = "days")
      
      # Indices used to extract daily cases starting from Mar 10, 2020
      todays = c(two_weeks_ago:today)
      yesterdays = todays - 1
      
      # Calculate daily cases for each census region
      new_cases_census = TS[,"Nt",todays] - TS[,"Nt",yesterdays]
      new_cases_census = as.data.frame(t(new_cases_census))
      new_cases = rowSums(new_cases_census)
      reported_cases = data.frame(date = dates, new_cases = new_cases)
      
      # Set iniital active cases to be the active cases (Tk) from two weeks ago
      initial_active_cases = sum(TS[,Tk,two_weeks_ago])
      
      # Calculate Rt for the province
      Rt = effective_reproduction_number(new_cases, initial_active_cases)
      return(Rt)
    }
    
    
    # indicator_to_tiernum: Converts a dataframe of indicator values to its corresponding tier number. 
    # indicator must be one of c("wir", "ppos" , "Rt")
    # 2D-array           <-         (2D-array,     string   )
    indicator_to_tiernum <- function(indicator_df, indicator){
      # indicator_df_old = indicator_df
      
      # Check if Date column exists
      dates_exist = FALSE
      if("Date" %in% colnames(indicator_df)){
        indicator_dates = indicator_df$Date
        indicator_df = subset(indicator_df, select = -Date)
        dates_exist = TRUE
      }
      
      # Define boolean variables based on indicator thresholds
      if(indicator == "wir"){
        lockdown_threshold = wir_lockdown_threshold <= indicator_df
        control_threshold  = wir_control_threshold <= indicator_df
        restrict_threshold = between(indicator_df, wir_restrict_threshold, wir_control_threshold)
        protect_threshold  = between(indicator_df, wir_protect_threshold, wir_restrict_threshold) 
        prevent_threshold  = indicator_df < wir_protect_threshold
        
      } else if(indicator == "ppos"){
        lockdown_threshold = ppos_lockdown_threshold <= indicator_df
        control_threshold  = ppos_control_threshold <= indicator_df
        restrict_threshold = between(indicator_df, ppos_restrict_threshold, ppos_control_threshold)
        protect_threshold  = between(indicator_df, ppos_protect_threshold, ppos_restrict_threshold) 
        prevent_threshold  = indicator_df < ppos_protect_threshold
        
      } else if(indicator == "Rt"){
        lockdown_threshold = Rt_lockdown_threshold <= indicator_df
        control_threshold  = Rt_control_threshold <= indicator_df
        restrict_threshold = between(indicator_df, Rt_restrict_threshold, Rt_control_threshold)
        protect_threshold  = between(indicator_df, Rt_protect_threshold, Rt_restrict_threshold) 
        prevent_threshold  = indicator_df < Rt_protect_threshold
        
      } else {
        warning("indicator must be one of c(\"wir\", \"ppos\" , \"Rt\")")
      }
      
      # Convert indicator threshold booleans to their corresponding tier numbers
      indicator_df = ifelse(lockdown_threshold, 5,
                            ifelse(control_threshold, 4,
                                   ifelse(restrict_threshold, 3,
                                          ifelse(protect_threshold, 2,
                                                 ifelse(prevent_threshold, 1, 
                                                        0 # 0 means none of the thresholds were triggered?!
                                                 )
                                          )
                                   )
                            )
      )
      
      if(dates_exist){
        indicator_df = cbind(indicator_dates, indicator_df)
      }
      
      return(indicator_df)
    }
    
    
    # zone_change: Updates the response_framework array with the zones that each region should be in based on the three indicators.
    # 2D-array  <-         (2D-array,     3D-array, 2D-array,  2D-array, 2D-array, int                  , function            , 2D-array)
    zone_change <- function(response_framework, TS, wir_time, ppos_time,  Rt_time, init_shutdown_len = 0, backup_func = median, region_id = regionid){
      # Ensure that the wir, % pos, and R_t values for each census are the same for census regions with identical PHUs. 
      # Some census regions share the same PHU!
      
      # Define key variables to check conditions on
      today <- gettime(TS)[1]
      tepi <- gettime(TS)[2]
      implementation_day = tepi + buffer + 1 # The value of tepi that the zone change will get implemented on. Off-by-one error present
      
      # # First provincial shutdown after implementation of the colour-coded lockdown scheme. 
      # # DEPRECATED
      # if(tepi >= workclosuregap2 - buffer){
      #   region_zone = rep("Shutdown", number_of_regions())
      #   response_framework <- rbind(response_framework, region_zone)
      #   return(response_framework)
      # }
      
      rows_rf = dim(response_framework)[1] # Current number of rows in the response_framework array
      
      # If response_framework is currently not activated, set everything to shutdown
      if(!between(implementation_day, rf_start, rf_end)){
        region_zone = rep("Shutdown", number_of_regions())
      }
      
      # Decisions to change tiers are made on a week-by-week basis
      # If implementation_day is not a Monday, then the new tier is equal to the previous day's tier
      # Exception is if today is the day that the decision is made to activate the response framework in 'buffer' days for the first time
      else if((implementation_day %% 7 != 6) & (implementation_day != rf_start)){ 
        # tepi == 0 is Mar 10th, 2020 which is a Tuesday
        region_zone = response_framework[rows_rf,]
      } 
      
      # Otherwise, do the tier calculations based on the mode of the three indicators
      else {
        
        # print("Monday")
        # This part of the function activates only on Mondays, 
        # or if today is the day that the decision is made to activate the response framework in 'buffer' days for the first time
        # Note: all provincial indicators have been turned off. 
        # Detect which indicators are active in the simulation.
        # If an indicator is active, do the necessary calculation on them,
        # If not, set their value equal to -1
        if("wir" %in% indicators_active){
          wir_reg = wir_time[today,]
          # wir_prov = weekly_incidence_rate_provincial(TS)
        } else {
          wir_reg = rep(-1, num_of_census)
          wir_prov = -1
        }
        
        if("ppos" %in% indicators_active){
          ppos_reg = ppos_time[today,]
        } else {
          ppos_reg = rep(-1, num_of_census)
          ppos_prov = -1
        }
        
        if("Rt" %in% indicators_active){
          Rt_reg = Rt_time[today,]
          # Rt_prov = Rt_provincial(TS)
        } else {
          Rt_reg = rep(-1, num_of_census)
          Rt_prov = -1
        }
        
        wir_reg_tiernum = indicator_to_tiernum(wir_reg, "wir")
        ppos_reg_tiernum = indicator_to_tiernum(ppos_reg, "ppos")
        Rt_reg_tiernum = indicator_to_tiernum(Rt_reg, "Rt")
        
        # Mode function. If there's no unique mode, return the median
        getmode <- function(v, backup_func = median) {
          # print(v)
          uniqv <- unique(v)
          
          if(length(uniqv) == length(v)){ # If there is no mode, return the median
            result = backup_func(v) # Change this line for median/max
            
          } else {
            result = uniqv[which.max(tabulate(match(v, uniqv)))]
          }
          
          # cat("result", result, "\n")
          return(result)
        }
        
        # Determine each region's zone by apply the mode to the three indicators. If there's no unique mode, then apply the median.
        indicator_tiernums = cbind(wir_reg_tiernum, ppos_reg_tiernum, Rt_reg_tiernum)
        region_zone = apply(indicator_tiernums, 1, getmode, backup_func = backup_func)
        region_zone <- ifelse(region_zone == 5, "Lockdown",
                              ifelse(region_zone == 4, "Control",
                                     ifelse(region_zone == 3, "Restrict",
                                            ifelse(region_zone == 2, "Protect",
                                                   ifelse(region_zone == 1, "Prevent",
                                                          response_framework[rows_rf,] # Otherwise, set the zone to what is was determined to be in the previous day
                                                   )
                                            )
                                     )
                              )
        )
      }
      
      # Append new tier levels to the response framework
      response_framework <- rbind(response_framework, region_zone)
      return(response_framework)
    }
    
    
    # Changes the epsilon values eps_w, eps_s, eps_h, eps_o based on what zone each region is in.
    # Outputs a parms array that changes based on the response_framework array
    eps_change <- function(parms, response_framework, eps_df, timeinfo){
      
      # Define key variables to check conditions on
      today = timeinfo[1]
      tepi = timeinfo[2]
      
      if (disable_response_framework){# Set to TRUE to disable the lockdown scheme
        parms[,"eps_w"] <- eps_df["Shutdown", "work"]
        parms[,"eps_s"] <- eps_df["Shutdown", "school"]
        parms[,"eps_h"] <- eps_df["Shutdown", "home"]
        parms[,"eps_o"] <- eps_df["Shutdown", "other"]
        return(parms)
        
      } else { # Response framework is implemented after rf_start
        eps_w <- parms[,"eps_w"]
        eps_s <- parms[,"eps_s"]
        eps_h <- parms[,"eps_h"]
        eps_o <- parms[,"eps_o"]
        
        # Extract today's tier levels
        zones_today <- response_framework[today,]
        
        # Set each eps value based on the tier of each region today
        eps_w <-  ifelse(zones_today == "Prevent", eps_df["Prevent", "work"],
                         ifelse(zones_today == "Protect", eps_df["Protect", "work"],
                                ifelse(zones_today == "Restrict", eps_df["Restrict", "work"],
                                       ifelse(zones_today == "Control", eps_df["Control", "work"],
                                              ifelse(zones_today == "Lockdown", eps_df["Lockdown", "work"], 
                                                     ifelse(zones_today == "Shutdown", eps_df["Shutdown", "work"],
                                                            -1 # if an epsilon value is negative, something went wrong
                                                     )
                                              )
                                       )
                                )
                         )
        )
        eps_s <-  ifelse(zones_today == "Prevent", eps_df["Prevent", "school"],
                         ifelse(zones_today == "Protect", eps_df["Protect", "school"],
                                ifelse(zones_today == "Restrict", eps_df["Restrict", "school"],
                                       ifelse(zones_today == "Control", eps_df["Control", "school"],
                                              ifelse(zones_today == "Lockdown", eps_df["Lockdown", "school"], 
                                                     ifelse(zones_today == "Shutdown", eps_df["Shutdown", "school"],
                                                            -1 # if an epsilon value is negative, something went wrong
                                                     )
                                              )
                                       )
                                )
                         )
        ) 
        eps_h <-  ifelse(zones_today == "Prevent", eps_df["Prevent", "home"],
                         ifelse(zones_today == "Protect", eps_df["Protect", "home"],
                                ifelse(zones_today == "Restrict", eps_df["Restrict", "home"],
                                       ifelse(zones_today == "Control", eps_df["Control", "home"],
                                              ifelse(zones_today == "Lockdown", eps_df["Lockdown", "home"], 
                                                     ifelse(zones_today == "Shutdown", eps_df["Shutdown", "home"],
                                                            -1 # if an epsilon value is negative, something went wrong
                                                     )
                                              )
                                       )
                                )
                         )
        ) 
        eps_o <-  ifelse(zones_today == "Prevent", eps_df["Prevent", "other"],
                         ifelse(zones_today == "Protect", eps_df["Protect", "other"],
                                ifelse(zones_today == "Restrict", eps_df["Restrict", "other"],
                                       ifelse(zones_today == "Control", eps_df["Control", "other"],
                                              ifelse(zones_today == "Lockdown", eps_df["Lockdown", "other"], 
                                                     ifelse(zones_today == "Shutdown", eps_df["Shutdown", "other"],
                                                            -1 # if an epsilon value is negative, something went wrong
                                                     )
                                              )
                                       )
                                )
                         )
        ) 
        
        # Update the parms array to today's eps values. 
        parms[,"eps_w"] <- eps_w
        parms[,"eps_s"] <- eps_s
        parms[,"eps_h"] <- eps_h
        parms[,"eps_o"] <- eps_o
        
        return(parms)
      }
    }
    ### ^^^^ Functions to implement Ontario's colour-coded lockdown scheme ^^^^ ###
    
    
    #Function to implement simulations. InitInf sets which stages the initially infected people are in, InitInf=c(6:10) indicates all initial infections correspond to individuals who are exposed, but split between different age classes
    msim3=function(FUN,parms,Trun=365,seed0=11,plotgive=TRUE,InitInf=c(6:10)){
      #assign initial infections
      L=nrow(parms); nI=InitInf;
      
      set.seed(seed0);
      NsA=sapply(popadj2020, function(x) rmultinom(1,x,(c(3019645,3475990,3855065,2505545,592260)/13448505)));  # pop age demographics, StatCan Census 2016
      
      #Intial infections across age classes
      set.seed(seed0+1);
      Infs1=NI0_1=apply(cbind(NsA[1,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[1]),x[1],prob=x[2]/length(nI[1])));
      set.seed(seed0+2);
      Infs2=NI0_2=apply(cbind(NsA[2,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[2]),x[1],prob=x[2]/length(nI[2])));
      set.seed(seed0+3);
      Infs3=NI0_3=apply(cbind(NsA[3,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[3]),x[1],prob=x[2]/length(nI[3])));
      set.seed(seed0+4);
      Infs4=NI0_4=apply(cbind(NsA[4,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[4]),x[1],prob=x[2]/length(nI[4])));
      set.seed(seed0+5);
      Infs5=NI0_5=apply(cbind(NsA[5,],parms[,"n0"]), 1, function(x) rbinom(n=length(nI[5]),x[1],prob=x[2]/length(nI[5])));
      
      ### vvv Initialize important arrays and dataframes before the simulation vvv ###
      num_of_census = number_of_regions(regionid) # Define number of census regions
      census_nums = paste("census.", c(1:num_of_census), sep = "") # Name every census region by number
      
      # Define objects representing each indicator using tier numbers. rows = time, columns = regions.
      wir_time = array(-1, dim=c(1, num_of_census)) # Object representing the weekly incidence rate, -1 means no data
      ppos_time = array(-1, dim = c(1, num_of_census)) # Object representing the percent positivity, -1 means no data
      Rt_time = array(-1, dim = c(1, num_of_census)) # Object representing the percent positivity, -1 means no data
      
      # Define arrays for checking the provincial values of the indicators
      if(!supercomputer_mode){
        wir_prov = array(-1, dim=c(1, 1))
        ppos_prov = array(-1, dim=c(1, 1))
        Rt_prov = array(-1, dim=c(1, 1))
      }
      
      # Response framework array. rows represent time, columns represent regions.
      response_framework <- create_response_framework("Shutdown", buffer)
      colnames(response_framework) = census_nums
      
      #Create object to store simulation results
      TS=array(0,dim=c(L,length(Snm),1)); 
      TS[,c(1,InitInf[1]),1]=cbind(NsA[1,]-NI0_1,Infs1); TS[,c(2,InitInf[2]),1]=cbind(NsA[2,]-NI0_2,Infs2); TS[,c(3,InitInf[3]),1]=cbind(NsA[3,]-NI0_3,Infs3); TS[,c(4,InitInf[4]),1]=cbind(NsA[4,]-NI0_4,Infs4);	TS[,c(5,InitInf[5]),1]=cbind(NsA[5,]-NI0_5,Infs5);
      TS[,57:59,1]=0; #Set Cw=Cs=C=0
      colnames(TS)=Snm; 
      
      # Input the response framework into TS, represent each tier as an int
      TS[,"Response_Framework",1] <- tier_to_number(response_framework[1,], gettime(TS))
      
      #Modify travel probability as needed
      Mn=round(parms[1,"Mfact"]*parms[,-(1:NP)]*(1-diag(L))); diag(Mn)=parms[,"N"]-colSums(Mn); parms[,-(1:NP)]=Mn%*%diag(1/colSums(Mn));
      ### ^^^ Initialize important arrays and dataframes before the simulation ^^^ ###
      
      #Implement simulation
      set.seed(seed0+6); Seeds=rnorm(Trun,1e6,1e6);
      for(i in 2:Trun) {
        TS=abind(TS,FUN(TS[,,i-1],closeinappl(parms,TS),gettime(TS),seed=Seeds[i]),along=3);
        
        # Arrays for checking the provincial values of the indicators
        if(!supercomputer_mode){
          wir_prov = rbind(wir_prov, weekly_incidence_rate_provincial(TS))
          ppos_prov = rbind(ppos_prov, percent_positivitiy_provincial(TS, testing_volumes))
          Rt_prov = rbind(Rt_prov, Rt_provincial(TS))
        }
        
        # Append new indicator values to their respective xxx_time matrix
        wir_time = weekly_incidence_rate_regional(wir_time, TS) 
        ppos_time = percent_positivitiy_regional(ppos_time, TS, testing_volumes)
        Rt_time = Rt_regional(Rt_time, TS)
        
        # Append new zone data to response_framework
        response_framework <- zone_change(response_framework, TS, wir_time, ppos_time, Rt_time, init_shutdown, backup_func = zone_change_backup_func) 
        
        # Change the parms array to update it with the new epsilon values
        parms <- eps_change(parms, response_framework, eps_df, gettime(TS))
        
        # This is where the response_framework array is appended to TS. 
        # To ensure that every value remains an int or a float, I made the following mapping:
        # Every zone is set to 0 before the response_framework is activated (For least squares fitting), but the epsilon values are set to shutdown
        #   "Prevent" = 1,
        #   "Protect" = 2,
        #   "Restrict" = 3,
        #   "Control" = 4,
        #   "Lockdown" = 5,
        #   "Shutdown" = 6
        TS[,"Response_Framework",i] <- tier_to_number(response_framework[i,], gettime(TS))  
        
        # Save these files at the end of the run
        # These files were used for debugging purposes
        # Comment out this if statement if saving these files is not desired
        if((debug_mode) & (i == Trun) & (!supercomputer_mode)){
          if("wir" %in% indicators_active){
            # Save the weekly incidence rate regional data
            row.names(wir_time) = c(1:Trun)
            colnames(wir_time) = census_nums
            saveRDS(wir_time, file = sprintf("%s/wir_time_sim.rds", output_dir))
            write.csv(wir_time, file = sprintf("%s/wir_time_sim.csv", output_dir))
            
            # Save the weekly incidence rate provincial data
            saveRDS(wir_prov, file = sprintf("%s/wir_prov_sim.rds", output_dir))
            write.csv(wir_prov, file = sprintf("%s/wir_prov_sim.csv", output_dir))
            
            # # # Save the weekly incidence rate tier numbers
            # wir_tiernum_sim = indicator_to_tiernum(wir_time, "wir")
            # write.csv(wir_tiernum_sim, file = sprintf("%s/wir_tiernum_sim.csv", output_dir))
          }
          
          if("ppos" %in% indicators_active){
            # Save the percent positivity regional data
            row.names(ppos_time) = c(1:Trun)
            colnames(ppos_time) = census_nums
            saveRDS(ppos_time, file = sprintf("%s/ppos_time_sim.rds", output_dir))
            write.csv(ppos_time, file = sprintf("%s/ppos_time_sim.csv", output_dir))
            
            # Save the percent positivity provincial data
            saveRDS(ppos_prov, file = sprintf("%s/ppos_prov_sim.rds", output_dir))
            write.csv(ppos_prov, file = sprintf("%s/ppos_prov_sim.csv", output_dir))
            
            # # Save the percent positivity tier numbers
            # ppos_tiernum_sim = indicator_to_tiernum(ppos_time, "ppos")
            # write.csv(ppos_tiernum_sim, file = sprintf("%s/ppos_tiernum_sim.csv", output_dir))
          }
          
          if("Rt" %in% indicators_active){
            # Save the Rt regional data
            row.names(Rt_time) = c(1:Trun)
            colnames(Rt_time) = census_nums
            saveRDS(Rt_time, file = sprintf("%s/Rt_time_sim.rds", output_dir))
            write.csv(Rt_time, file = sprintf("%s/Rt_time_sim.csv", output_dir))
            
            # Save the Rt provincial data
            saveRDS(Rt_prov, file = sprintf("%s/Rt_prov_sim.rds", output_dir))
            write.csv(Rt_prov, file = sprintf("%s/Rt_prov_sim.csv", output_dir))
            
            # # Save the Rt tier numbers
            # Rt_tiernum_sim = indicator_to_tiernum(Rt_time, "Rt")
            # write.csv(Rt_tiernum_sim, file = sprintf("%s/Rt_tiernum_sim.csv", output_dir))
          }
          
          row.names(response_framework) = c(1:(Trun+buffer))
          saveRDS(response_framework, file = sprintf("%s/response_framework_sim.rds", output_dir))
          write.csv(response_framework, file = sprintf("%s/response_framework_sim.csv", output_dir))
          
          saveRDS(TS, file = sprintf("%s/TS_sim.rds", output_dir))
        }
      }
      ### ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Response Framework Code ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ###
      
      
      TS[,"C",]<-apply(TS[,c("Cw","Cs"),],c(1,3),sum);
      
      #Different levels of aggregation in model output. plotgive=3 or 3.5 give shortest output form (tracking only infections and costs) in integer format to reduce output size
      if(plotgive=="TS") return(TS); if(plotgive==TRUE) return(Resagg(TS,parms,plotgive==1));
      if(plotgive%in%c(3,3.5)){ TS[,"C",]=TS[,"C",]*matrix(parms[,"N"],nrow(parms),Trun); TSn=apply(TS,c(2,3),sum); out=matrix(as.integer(TSn),nrow(TSn),ncol(TSn)); if(plotgive==3.5) return(out[c(1:5,59),]); return(rbind(out,colMeans(TS[,"C",]>0))); }
      #In fitting also tracked mean and variance in proportion distancing (omgs) and proportion of infections in Toronto (fracTor)
      if(plotgive=="fit"){ # This is the output for fitting
        omgs=t(apply(TS[,"VD",], 2, function(x) meansd(x,parms[,"N"]))); #print(c(omgs));
        omgs_names = c("meanVD", "sdVD")
        colnames(omgs) = omgs_names
        
        propCits=t(TS[,"Nt",]); states=t(apply(TS,c(2,3),sum)); 
        propCits_names = paste("propCits.", 1:(dim(propCits)[2]), sep = "")
        colnames(propCits) = propCits_names
        
        rf_regions = t(TS[,"Response_Framework",]);
        colnames(rf_regions) = paste("rf.phu.", 1:number_of_regions(regionid), sep = "")
        
        return(cbind(states,omgs,propCits, rf_regions)); 
      } 
    }
    ###^^^^^^^^^^^^^^ Copy and paste simulation code here!!! ^^^^^^^^^^^^^^###
    
    ### Checks to ensure the simulation is in "fitting" mode
    if(!parms_initialized){
      w = warning("parms array not initialized properly!")
      print(w)
    }
    
    ###vvvvvvvvvvvvvv Copy and paste fitting procedure code here!!! vvvvvvvvvvvvvv###
    ### Fitting procedure
    # library(nloptr); library(mgcv);
    agelabels=c("0-19", "20-39", "40-59", "60-79", "80+")
    
    #To aggregate census divisions by PHU
    regionlabels=unique(regionid$phu[order(regionid$phunum)])
    bincase.new=function(x) {aggregate(x, list(PHU=regionid$phunum[order(regionid$census.region)]), FUN=sum)[,2];}
    binrf.new=function(x) {aggregate(x, list(PHU=regionid$phunum[order(regionid$census.region)]), FUN=mean)[,2];}
    
    LLfun=function(parmsTry,parmsFit=parmnames,extras,reps=5){
      
      parmsTry0=parmsTry;
      parmsTry[parmsFit=="zeta"]=parmsTry[parmsFit=="zeta"]/1e3
      parmsTry[parmsFit=="psi1"]=parmsTry[parmsFit=="psi1"]/1e3;
      parmsTry[parmsFit=="psi2"]=parmsTry[parmsFit=="psi2"]/1e3;
      parmsTry[parmsFit=="psi3"]=parmsTry[parmsFit=="psi3"]/1e3;
      parmsTry[parmsFit=="psi4"]=parmsTry[parmsFit=="psi4"]/1e3;
      parmsTry[parmsFit=="psi5"]=parmsTry[parmsFit=="psi5"]/1e3;
      parmsTry[parmsFit=="tau0"]=parmsTry[parmsFit=="tau0"]/1e3;
      parmsTry[parmsFit=="boost"]=parmsTry[parmsFit=="boost"]/1e3;
      
      #Fitting these parameters on a log scale:
      parmsTry[parmsFit=="omg"]=exp(parmsTry[parmsFit=="omg"]);
      parmsTry[parmsFit=="beta"]=exp(parmsTry[parmsFit=="beta"])/mean(popadj2020);
      
      # Add PHU-specific transmission modifiers
      news=parmsFit%in%names(parmsB0); Bmod=rep(1,49); if(length(news)>0) Bmod=BetaMod(parmsTry[news]);
      
      parms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),
                  xi_a = xi_a, xi_b = xi_b, xi_c = xi_c,
                  beta=Bmod*1.654565e-06,
                  eps_w=eps_w, eps_s=eps_s,eps_h=eps_h,eps_o=eps_o,omg=2.28e+05,r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
                  atravel=0.794, a1=1.08,a2=1.05,a3=1,a4=1.89,a5=41.0, reopen_w=0.5, reopen_s=0.5, B=0.3, phi=50, zeta=0.025,
                  tau0=0.5, psi1=0.025,  psi2=0.025,  psi3=0.025,  psi4=0.025,  psi5=0.025, kappa=0.418, taumax=0.8, boost=0.5,
                  alpha_w1 = alpha_w1, alpha_w2 = alpha_w2, alpha_w3 = alpha_w3, alpha_w4 = alpha_w4,
                  alpha_s1 = alpha_s1, alpha_s2 = alpha_s2, alpha_s3 = alpha_s3, alpha_s4 = alpha_s4,
                  alpha_h1 = alpha_h1, alpha_h2 = alpha_h2, alpha_h3 = alpha_h3, alpha_h4 = alpha_h4,
                  alpha_o1 = alpha_o1, alpha_o2 = alpha_o2, alpha_o3 = alpha_o3, alpha_o4 = alpha_o4, 
                  DATA$Msave[-50,]);
      
      #Throw out bad guesses for transmission probabilities (highest probability exceeds 1)
      tcheck<-2*(1+parmsTry[parmsFit=="B"])*((1/parms[1,"s"])-1)*max(Bmod)*parmsTry[parmsFit=="beta"]*max(parmsTry[parmsFit %in% c("a1", "a2", "a3", "a4", "a5")])
      # print("tcheck")
      # print(tcheck)
      if(tcheck>1) {print(sprintf("Dud transmission probability: %f, %f, %f, %f, return %f", parmsTry[parmsFit=="B"], parmsTry0[parmsFit=="beta"], max(parmsTry[parmsFit %in% c("a1", "a2", "a3", "a4", "a5")]), max(Bmod), tcheck)); return(1e8);}
      
      parmsTry=pmax(parmsTry,0); parms[,parmsFit[!news]]=matrix(parmsTry[!news],nrow(parms),length(parmsFit[!news]),byrow=TRUE);
      parms[,"beta"]=Bmod*parms[,"beta"]; #Modify our current beta value by Bmod fcn with region specific transmision modifiers
      if(!"eps_s"%in%parmsFit) parms[,"eps_s"]=parms[,"eps_w"]; #if school distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"eps_h"%in%parmsFit) parms[,"eps_h"]=parms[,"eps_w"]; #if home distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"eps_o"%in%parmsFit) parms[,"eps_o"]=parms[,"eps_w"]; #if other distancing efficacy not specified, assume it equals work distancing efficacy
      if(!"reopen_s"%in%parmsFit) parms[,"reopen_s"]=parms[,"reopen_w"]; #if efficacy of distanced school reopening not specified, assume it equals efficacy of distanced work reopening
      if(is.na(reps)) return(parms); #if only want parameter matrix
      
      #Need to set 2nd dimension of sim to match cols of your msim3 output (115 for geog & age output, 66 for age output)
      sim=array(dim=c(extras["Trun"],166,0)); for(i in 1:reps) sim=abind(sim, msim3(m3iter,parms,extras["Trun"],plotgive="fit",seed0=i*1e3), along=3);
      if(!supercomputer_mode) {saveRDS(sim, file = "sim.rds")}
      
      lenR=extras["dayRfin"];
      FVN=5;
      simp=array(dim=c(lenR+1,5+FVN,0));
      simp_g=array(dim=c(lenR+1,34+FVN,0));
      NTS=simp_g[,1:34,0];
      rf_sim = array(data = 0, dim = c(lenR+1,34,0))
      colnames(rf_sim) = c(paste("rf.phu.", 1:34, sep = ""))
      
      #Counter for skips where simulation does not result in 50 total cases fast enough to have a timeseries of sufficent length for comparison to actual case counts for Ontario
      skipcount=0;
      
      for(i in 1:reps){
        #If with a years sim we don't get to 50 cases fast enough, skip this run
        lcheck=sim[,"Nt",i][which(sim[,"Nt",i]>49)];
        
        if(length(lcheck)<(lenR+1) | max(sim[,"Nt",i])<50 ) {print(sprintf("bad run %i", i)); skipcount=skipcount+1; next;}
        
        trg=sim[,,i][which(sim[,"Nt",i]>49)[1]+(0:lenR),]; # Eliminates rows from sim where known total is less than 50
        if(nrow(trg)<lenR | max(is.na(trg))==1) next;
        
        ###Geog aggregation based on PHU
        nts=t(apply(trg[,c(paste("propCits.", 1:number_of_regions(regionid), sep = ""))],1,bincase.new));
        NTS=abind(NTS,nts,along=3);
        
        rf_sim_temp0=trg[,c(paste("rf.phu.", 1:number_of_regions(regionid), sep = ""))]
        rf_sim_temp = t(apply(rf_sim_temp0, 1, binrf.new))
        rf_sim = abind(rf_sim, rf_sim_temp, along = 3) 
        
        simp=abind(simp, cbind(trg[,"Nt"], rowSums(trg[,c(Da,Di,R)])/trg[,"Nt"], trg[,length(Snm)+1:2], rowSums(trg[,Tk]), trg[, c("Nt1", "Nt2", "Nt3", "Nt4", "Nt5")]),  along=3)
      };
      
      if(reps-skipcount<=1) {print("Too few viable runs"); return(1e8);} #Return the "bad" value if at most 1 sim reached 50 cases quickly enough
      
      ###For geog aggregation based on PHU
      NTS=apply(NTS,2:3,diff);
      
      simp[,2,][simp[,2,]==Inf]=max(simp[,2,][simp[,2,]!=Inf]);
      if(length(dim(simp))<3) return(1e6);
      
      
      #Calculations for cost function based on time series fits
      tots.age=agedat[cumsum(rowSums(agedat[,3:7]))>50,3:7]; to.age=rowSums(tots.age); #Fitting to observed cases by ageclass, no longer drop days as we remove "grey" days in pre-processing
      tots.region=regiondat[cumsum(rowSums(regiondat[,2:35]))>50,2:35]; to.region=rowSums(tots.region); #Fitting to observed cases by ageclass, no longer drop days as we remove "grey" days in pre-processing
      mobility=omgdat$distancing #select residential mobility data, convert from percentages
      tp.age=apply(simp[1:(length(to.age)+1),1,],2,diff); tpts.age=apply(simp[1:(length(to.age)+1),6:10,],2:3,diff); #Fitting modeled new cases by ageclass
      tpts.region=NTS[1:length(to.region),,]; #Fitting modelled new cases by region
      
      tpts_DV=100*simp[1:length(mobility),3,] #Pull out mobility data for days past 50th case (keep only provincial means, not sds)
      tpts_RV=simp[1:length(to.age),2,] #Pull out testing ratio data for days past 50th case
      
      Tlik0.age=apply(abind(tots.age,tpts.age,along=3)[1:length(to.age),,], 1:2, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) { for (j in 1:ncol(tots.age))  {Tlik0.age[i,,j]<-Tlik0.age[i,,j]/(colMeans(tots.age)[j])^2}}
      Tlik.age=sum(apply(Tlik0.age, c(2,3), mean)) #take sum of the mean (for each ts across runs) squared differences
      
      # Region fitting using least squares
      Tlik0.region=apply(abind(tots.region,tpts.region,along=3)[1:length(to.region),,], 1:2, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) { for (j in 1:ncol(tots.region))  {Tlik0.region[i,,j]<-Tlik0.region[i,,j]/(colMeans(tots.region)[j])^2}}
      Tlik.region=sum(apply(Tlik0.region, c(2,3), mean))/50 #take sum of the mean (for each ts across runs) squared differences, /50 to get roughly the same order of mag as Tlik.age
      
      # Colour-Tier system fitting using least squares
      rf_obs = rfdat[-1,-1] # Remove the date column and the first row to align the data. (There's an off-by-one error, but that doesn't matter much)
      rf_to.region=rowSums(rf_obs);
      rf_sim = rf_sim[1:length(rf_to.region),,]
      rf_diff0.region = apply(abind(rf_obs, rf_sim, along = 3)[1:length(rf_to.region),,], 1:2, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) { for (j in 1:ncol(rf_obs)) {rf_diff0.region[i,,j] <- rf_diff0.region[i,,j]/(colMeans(rf_obs)[j])^2}}
      rf_diff.region = sum(apply(rf_diff0.region, c(2,3), mean))/50
      
      if(!supercomputer_mode){
        saveRDS(tots.region, file = sprintf("%s/tots_region.rds", output_dir))
        saveRDS(tpts.region, file = sprintf("%s/tpts_region.rds", output_dir))
        saveRDS(rf_obs, file = sprintf("%s/rf_obs.rds", output_dir))
        saveRDS(rf_sim, file = sprintf("%s/rf_sim.rds", output_dir))
        
        write.csv(rf_obs, file = sprintf("%s/rf_obs.csv", output_dir))
        write.csv(rf_sim, file = sprintf("%s/rf_sim.csv", output_dir))
        
        print("supercomputer mode off")
      }
      
      #New version of omega/individual NPI adherence fitting using least squares, based on google mobility data (replaces Dlik/DV)
      Mlik0=apply(abind(mobility,tpts_DV,along=2)[1:length(mobility),], 1, function(x) (x[1]-x[-1])^2)
      for (i in 1:(reps-skipcount)) {Mlik0[i,]<-Mlik0[i,]/(mean(mobility))^2}
      Mlik=sum(apply(Mlik0, 2, mean)) #take sum of the mean (for each ts across all runs) squared differences
      
      #Calculate likelihood of testing ratio (of positive cases to total infected)
      RV=meansd(tpts_RV[20:40,])*c(1,2) #Check testing ratio over the 20th-40th day after cases exceed 50
      sta=round(RV,2);
      Rlik=dnorm(DATA$testRat,RV[1],RV[2],log=TRUE);
      
      # Check the magnitude of rf_diff.region to see if it needs a boost
      LL=sum(Tlik.age,Tlik.region, 10*rf_diff.region,4*abs(Rlik),5*(RV[2]>16), Mlik);
      stats=round(c(parmsTry0,RV[1],Tlik.age,Tlik.region,rf_diff.region,Rlik,Mlik,LL,min(c(1e8,REPORT[,ncol(REPORT)-1]),na.rm=TRUE)),2); assign("REPORT",rbind(REPORT,stats),.GlobalEnv); print(tail(REPORT,1));
      
      tryCatch({
        if((LL<tail(stats,1) & plotgive==TRUE) | plotgive==99){
          par(mfrow=c(2,4));
          mxplt.age=length(to.age); mxplt_mob=length(mobility); 
          
          # Plot NPI adherence
          matplot(tpts_DV[1:mxplt_mob,]/100,
                  type="l",
                  xlab="Day since 50th positive case",
                  ylab="Individual NPI adherence",
                  lwd=1.5,
                  lty=1,
                  col=1,
                  ylim=range(c(tpts_DV/100, mobility/100)));
          points(mobility/100,pch=16,col=2);
          
          # Plot Daily reported cases by age
          matplot(tp.age[1:mxplt.age,],type="l",xlab="Day since 50th positive case",ylab="Daily reported cases",lwd=1.5,lty=1,col=1,main=paste0(sta,collapse="_"),ylim=range(to.age)); points(to.age,pch=16,col=2);
          for(i in c(1:ncol(tots.age))){ matplot(tpts.age[1:mxplt.age,i,],type="l",lwd=2,lty=1,col=8,ylab="Daily reported cases",main=agelabels[i],ylim=c(0,max(tots.age[,i]))); points(tots.age[,i],col=2,pch=16); };
          
          # Plot Daily reported cases by region
          par(mfrow=c(5,7)); mxplt.region=length(to.region);
          for(i in c(1:ncol(tots.region))){ matplot(tpts.region[1:mxplt.region,i,],type="l",lwd=2,lty=1,col=8,ylab="Daily reported cases",main=regionlabels[i],ylim=c(0,max(tots.region[,i]))); points(tots.region[,i],col=2,pch=16); };
          
          ###vvv Plotting response framework graphs vvv###
          # multiplot <- dget("multiplot.R")
          
          datetag = "may312021"
          input_dir = "InputFiles"
          
          rf_start = 241
          start_date = "2020-11-07"
          end_date = "2021-02-07" # Some regions reopened schools on Feb 8th
          date_seq = seq(as.Date(start_date), as.Date(end_date), by = "day")
          
          rf_obs = rf_obs[-c(1:rf_start),]
          rf_sim = data.frame(rf_sim[-c(1:rf_start),,1])
          colnames(rf_sim) = paste("phu.", 1:34, sep = "")
          
          rf_obs$Date = date_seq
          rf_sim$Date = date_seq
          
          # Define regiondict
          response_framework_raw = read.csv(sprintf("response_framework_%s.csv", datetag))
          regiondict = data.frame(phu = unique(response_framework_raw$Reporting_PHU),
                                  phunum = 1:34)
          
          
          # Pivoted response_framework dataframe to a "long" orientation
          df.long_obs <- rf_obs %>% tidyr::pivot_longer(!Date, names_to = "phunum", values_to = "tiernum", names_prefix="phu.")
          df.long_sim <- rf_sim %>% tidyr::pivot_longer(!Date, names_to = "phunum", values_to = "tiernum", names_prefix="phu.")
          
          #This just adds nice labels to our plots 
          df.tiers<-data.frame(tiers=c("Other", "Prevent", "Protect", "Restrict", "Control", "Lockdown", "Shutdown"), tiernum=0:6)
          
          df.annotated_obs0<-merge(df.long_obs, regiondict, by="phunum")
          df.annotated_obs<-merge(df.annotated_obs0,df.tiers, by="tiernum")
          df.annotated_obs = df.annotated_obs[order(df.annotated_obs$phunum, df.annotated_obs$Date),]
          
          df.annotated_sim0<-merge(df.long_sim, regiondict, by="phunum")
          df.annotated_sim<-merge(df.annotated_sim0,df.tiers, by="tiernum")
          df.annotated_sim = df.annotated_sim[order(df.annotated_sim$phunum, df.annotated_sim$Date),]
          
          
          #Plot ON's response framework
          rf_obs_plot<-ggplot(data=df.annotated_obs, aes(x=Date, y=phu)) +
            geom_tile(aes(group=phu, fill=as.factor(tiernum))) +
            xlab("Date") +
            ylab("Public Health Unit") +
            scale_x_date(expand = c(0,0)) +
            scale_y_discrete(expand=c(0,0)) +
            scale_fill_manual("Tier",
                              limits = factor(0:6),
                              values=c("blue", "green", "yellow", "orange", "red", "grey50", "black"),
                              labels=c("Other", "Prevent", "Protect", "Restrict", "Control", "Lockdown", "Shutdown"), drop=FALSE) +
            theme_bw() +
            ggtitle("Response Framework: Ontario's Data") +
            theme(plot.title = element_text(hjust = 0.5))
          
          #Plot the sim's response framework
          rf_sim_plot<-ggplot(data=df.annotated_sim, aes(x=Date, y=phu)) +
            geom_tile(aes(group=phu, fill=as.factor(tiernum))) +
            xlab("Date") +
            ylab("Public Health Unit") +
            scale_x_date(expand = c(0,0)) +
            scale_y_discrete(expand=c(0,0)) +
            scale_fill_manual("Tier",
                              limits = factor(0:6),
                              values=c("blue", "green", "yellow", "orange", "red", "grey50", "black"),
                              labels=c("Other", "Prevent", "Protect", "Restrict", "Control", "Lockdown", "Shutdown"), drop=FALSE) +
            theme_bw() +
            ggtitle("Response Framework: Simulation Data") +
            theme(plot.title = element_text(hjust = 0.5))
          
          # Custom Colour Palette
          # library(RColorBrewer);
          cust.pal<-brewer.pal(11,"RdBu")
          
          
          #Plot the tier shifts - note that I've set the colours so any case where there was no shift (shift=0), or the shift was only -1,+1 (i.e. green to yellow, or red to orange) are white
          #This choice was made since those are the ones we'd expect, and I wanted to highlight the big shifts where tier(s) were skipped over (i.e. green to red) using the red/blue
          #You'll likely want to tweak that palette for your viz, since unlike with the switch values, a -1,+1 value would still indicate a discrepancy, but it can still be useful to group values
          rf_diff = df.annotated_obs[, c("phunum", "Date", "phu")]
          rf_diff$ON_tiernum  = df.annotated_obs$tiernum
          rf_diff$sim_tiernum = df.annotated_sim$tiernum
          rf_diff <- rf_diff %>% mutate(tiernum_diff = ON_tiernum - sim_tiernum) 
          
          rf_diff_plot<-ggplot(data=rf_diff, aes(x=Date, y=phu)) +
            geom_tile(aes(group=phu, fill=as.factor(tiernum_diff))) +
            xlab("Date") +
            ylab("Public Health Unit") +
            scale_x_date(expand = c(0,0)) +
            scale_y_discrete(expand=c(0,0)) +
            scale_fill_manual("Shift",
                              limits = factor(-5:5),
                              values=c(cust.pal[1],cust.pal[2],cust.pal[3],cust.pal[4],cust.pal[5],cust.pal[6],cust.pal[7],cust.pal[8],cust.pal[9],cust.pal[10],cust.pal[11]),
                              na.value="grey50") +
            theme_bw() +
            ggtitle("Response Framework Differences") +
            theme(plot.title = element_text(hjust = 0.5))
          
          
          multiplot(rf_obs_plot, rf_diff_plot, rf_sim_plot, cols = 2)
          print("")
          ###^^^ Plotting response framework graphs ^^^###
          
        }},
        error = function(e){print(e)}
      )
      
      return(pmin(LL,1e8));
    }
    
    #Fitting implementation:
    finday<-max(nrow(regiondat[cumsum(rowSums(regiondat[,2:35]))>50,2:35]), nrow(agedat[cumsum(rowSums(agedat[,3:7]))>50,]), nrow(omgdat), nrow(rfdat[cumsum(rowSums(rfdat[,-1]))>50,]) #, nrow(rfdat)
    )
    nrep=5; plotgive=FALSE; #Set plotgive to "TRUE" to see fit of model to cases by age classes and cases by PHU
    # extras=c(Trun=400,dayTfin=64,dayRfin=finday,model=3.5); # Original
    extras=c(Trun=375,dayTfin=64,dayRfin=finday,model=3.5); # Modify to be shorter
    
    parmsB0=c(xi_a = xi_a, xi_b = xi_b, xi_c = xi_c);
    
    fitnames=c("omg","eps_w","eps_s","eps_h","eps_o",
               "atravel", "a1","a2","a3", "a4","a5", 
               "reopen_w", "reopen_s", "B", "phi", "beta", "zeta",
               "tau0", "psi1","psi2","psi3", "psi4", "psi5",
               "kappa","taumax", "boost",
               "xi_a", "xi_b", "xi_c",
               "alpha_w1", "alpha_w2", "alpha_w3", "alpha_w4",
               "alpha_s1", "alpha_s2", "alpha_s3", "alpha_s4",
               "alpha_h1", "alpha_h2", "alpha_h3", "alpha_h4",
               "alpha_o1", "alpha_o2", "alpha_o3", "alpha_o4");
    fitnames_len <- length(fitnames)
    pfit=c(1:fitnames_len)
    
    parmsFit=fitnames[pfit];
    
    #Matrix to store fitting results
    REPORT=matrix(nrow=0,ncol=length(pfit)+8); colnames(REPORT)=c(fitnames[pfit],"RV","Tlik.age","Tlik.region", "rf_diff.region", "Rlik","Mlik","LL","pbLL");
    
    #Setting initial values and constraints on parameters
    lowerBex=c(8,0.3,0.3,0.3,0.3,0.25,rep(0.1,5), c(0,0), c(0,1), -5, 0, 0, rep(0,5), 0, 0.3, 0, c(0,0,0), rep(0,16)); 
    upperBex=c(12,0.9,0.9,0.9,0.9,0.75, rep(10,5), c(0.75,0.75), c(0.5,150), 10, 75, 100, rep(125,5), 0.5, 0.8, 750, c(100,1e-5,10), rep(1,16));
    bval.guess<-c(xi_a, xi_b, xi_c)
    parmStr0=c(8.210587, eps_w, eps_s, eps_h, eps_o,
               7.381301e-01   , 7.110524e-01, 1.048301, 1.206083, 2.732836   , 9.892207   ,
               5.804198e-01   , 4.969267e-01, 9.117197e-02, 8.909216e+01, -1.333445, 2.239364e+01,
               3.032116e+01   , 5.806795e+01, 1.336930e+01, 1.099442e+01  , 3.411321e+01 , 1.141062e+02   ,
               3.158605e-01   , 5.839257e-01, 5.441355e+02,
               bval.guess,
               alpha_w1, alpha_w2, alpha_w3, alpha_w4,
               alpha_s1, alpha_s2, alpha_s3, alpha_s4,
               alpha_h1, alpha_h2, alpha_h3, alpha_h4,
               alpha_o1, alpha_o2, alpha_o3, alpha_o4)
    
    names(parmStr0) = fitnames
    ###^^^^^^^^^^^^^^ Copy and paste fitting procedure code here!!! ^^^^^^^^^^^^^^###
    
parmsTry=y$solution
 
if (exists("fullparm.storage")==FALSE) {
  fullparm.storage<-c(y$solution, y$objective)
} else {
  fullparm.storage<-rbind(fullparm.storage, c(y$solution, y$objective))
}
    
#Re-run LLfun to see plot of fit, and to verify all parameter input is working correctly
plotgive="TRUE"
llinfo=LLfun(parmsTry,parmsFit=parmsFit,extras=extras,reps=nrep) 
print("hi")
if (abs(y$objective-llinfo)>1) {
  print("PROBLEM WITH FIT");
  cat("y$objective =", y$objective, "\n")
  cat("llinfo =", llinfo, "\n")
  
  cat("diff =", y$objective-llinfo, "\n")
  } 


# Final parm set, converted from y$solution for input into visualization code
news=parmsFit%in%names(parmsB0);

# Extract the fitted parameters and their variable names
fitted_parms = y$solution
names(fitted_parms) = names(x.dat$x0)
print(fitted_parms)

fitparms=cbind(N=popadj2020,Mfact=popratio*2,s=0.2,Tg_w=1,Tl_w=exp(-9.14),Tg_s=1,Tl_s=exp(-9.14),beta=BetaMod(parmsTry[news])*exp(fitted_parms["beta"])/mean(popadj2020),
               eps_w=fitted_parms["eps_w"],eps_s=fitted_parms["eps_s"],eps_h=fitted_parms["eps_h"],eps_o=fitted_parms["eps_o"],omg=exp(fitted_parms["omg"]),r=0.19,eta=0.8,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,
               atravel=fitted_parms["atravel"], a1=fitted_parms["a1"],a2=fitted_parms["a2"],a3=fitted_parms["a3"],a4=fitted_parms["a4"],a5=fitted_parms["a5"], reopen_w=fitted_parms["reopen_w"], reopen_s=fitted_parms["reopen_s"],
               B=fitted_parms["B"], phi=fitted_parms["phi"], zeta=fitted_parms["zeta"]/1e3,
               tau0=fitted_parms["tau0"]/1e3,
               psi1=fitted_parms["psi1"]/1e3,  psi2=fitted_parms["psi2"]/1e3,  psi3=fitted_parms["psi3"]/1e3,  psi4=fitted_parms["psi4"]/1e3,  psi5=fitted_parms["psi5"]/1e3,
               kappa=fitted_parms["kappa"], taumax=fitted_parms["taumax"], boost=fitted_parms["boost"]/1e3, 
               xi_a = fitted_parms["xi_a"], xi_b = fitted_parms["xi_b"], xi_c = fitted_parms["xi_c"],
               alpha_w1 = fitted_parms["alpha_w1"], alpha_w2 = fitted_parms["alpha_w2"], alpha_w3 = fitted_parms["alpha_w3"], alpha_w4 = fitted_parms["alpha_w4"],
               alpha_s1 = fitted_parms["alpha_s1"], alpha_s2 = fitted_parms["alpha_s2"], alpha_s3 = fitted_parms["alpha_s3"], alpha_s4 = fitted_parms["alpha_s4"],
               alpha_h1 = fitted_parms["alpha_h1"], alpha_h2 = fitted_parms["alpha_h2"], alpha_h3 = fitted_parms["alpha_h3"], alpha_h4 = fitted_parms["alpha_h4"],
               alpha_o1 = fitted_parms["alpha_o1"], alpha_o2 = fitted_parms["alpha_o2"], alpha_o3 = fitted_parms["alpha_o3"], alpha_o4 = fitted_parms["alpha_o4"], 
               DATA$Msave[-50,])

parms=fitparms; # use best fit from fitting
nrun<-5

for (run in 1:nrun) {

  #Run a sim based on a fitted parameter set
  runlength<-375
  x_region=(msim3(m3iter,parms,runlength,plotgive="TS",seed0=fit*run*1e4,InitInf=c(6:10)));

  #Quick plot showing timeseries of total cases to visually verify code is running correctly
  plotgive = FALSE
  tryCatch({
    par(mfrow=c(1,1));
    agedat.plot<-agedat[agedat$Date<"2021-02-17",];
    dat=data.frame(ts50=1:runlength-min(which(colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),]))>50)), Nt=colSums(colSums(x_region[,c("Nt1","Nt2","Nt3","Nt4","Nt5"),])));
    plot(dat$ts50, dat$Nt-dplyr::lag(dat$Nt), col="blue", type="l");
    lines(agedat.plot$ts50, rowSums(agedat.plot[,c( "newK.1", "newK.2", "newK.3", "newK.4", "newK.5")]));},
    error = function(e){print(e)}
  )
  
  ##### Save all output from simulation runs
  #Create all the summary data that would normally be generated with plotgive="TS" to use for the provincial dataset
  omgs=t(apply(x_region[,"VD",], 2, function(x) meansd(x,parms[,"N"])));
  omgs_names = c("meanVD", "sdVD")
  colnames(omgs) = omgs_names

  propCits=t(x_region[,"Nt",]); states=t(apply(x_region,c(2,3),sum));
  propCits_names = paste("propCits.", 1:(dim(propCits)[2]), sep = "")
  colnames(propCits) = propCits_names
  
  x_prov=cbind(states,omgs,propCits);
  
  

  for (i in 1:length(x_prov[,1])) {

    tmp<-as.data.frame(t(x_prov[i,]))
    
    if (i==1) {
      df<-tmp
    } else {
      df<-rbind(df,tmp)
    }
  }

  for (i in 1:length(x_region[,1,1])) {

    tmp_region<-as.data.frame(t(x_region[i,,]))
    tmp_region$ts <- as.numeric(row.names(tmp_region))-1
    tmp_region$region<-i

    if (i==1) {
      df_region<-tmp_region
    } else {
      df_region<-rbind(df_region,tmp_region)
    }
  }

  #Add-ins for regional sims
  df_region$sick<-rowSums(df_region[,colnames(df_region) %in% c(E,Di,Da)]); df_region$know.current<-rowSums(df_region[,colnames(df_region) %in% Tk]); df_region$recovered<-rowSums(df_region[,colnames(df_region) %in% R])

  df_region$known1<-rowSums(df_region[, colnames(df_region) %in% Tk[1:4]]);df_region$known2<-rowSums(df_region[, colnames(df_region) %in% Tk[5:8]]);df_region$known3<-rowSums(df_region[, colnames(df_region) %in% Tk[9:12]]);df_region$known4<-rowSums(df_region[, colnames(df_region) %in% Tk[13:16]]);df_region$known5<-rowSums(df_region[, colnames(df_region) %in% Tk[17:20]]);
  df_region$sick1<-rowSums(df_region[, colnames(df_region) %in% sick1]); df_region$sick2<-rowSums(df_region[, colnames(df_region) %in% sick2]); df_region$sick3<-rowSums(df_region[, colnames(df_region) %in% sick3]); df_region$sick4<-rowSums(df_region[, colnames(df_region) %in% sick4]); df_region$sick5<-rowSums(df_region[, colnames(df_region) %in% sick5]);
  # df_region$fit<-unlist(strsplit(strsplit(file.names[fit], "_")[[1]][5], ".rds"));
  file_names_list = strsplit(file.names[fit], "_")[[1]]
  str_arr_len = length(file_names_list)
  df_region$fit<-unlist(strsplit(file_names_list[str_arr_len], ".rds"));
  df_region$run<-run
  df_region$LL<-y$objective

  df_region <- df_region %>% group_by(fit, run, region) %>% mutate(new.K=Nt - lag(Nt),
                                                                   new.K.1=Nt1 - lag(Nt1), new.K.2=Nt2 - lag(Nt2), new.K.3=Nt3 - lag(Nt3), new.K.4=Nt4 - lag(Nt4), new.K.5=Nt5 - lag(Nt5),
                                                                   new.sick=(sick + recovered) -  lag(sick + recovered),
                                                                   new.sick.1=(sick1 + R1) - lag(sick1  + R1), new.sick.2=(sick2  + R2) - lag(sick2  + R2), new.sick.3=(sick3  + R3) - lag(sick3  + R3), new.sick.4=(sick4  + R4) - lag(sick4  + R4), new.sick.5=(sick5  + R5) - lag(sick5  + R5))


  #Set all first step "new" to equal total at that time since we dont have a previous ts for comparison
  df_region$new.K[df_region$ts==min(df_region$ts)]<-df_region$Nt[df_region$ts==min(df_region$ts)]
  df_region$new.K.1[df_region$ts==min(df_region$ts)]<-df_region$Nt1[df_region$ts==min(df_region$ts)]
  df_region$new.K.2[df_region$ts==min(df_region$ts)]<-df_region$Nt2[df_region$ts==min(df_region$ts)]
  df_region$new.K.3[df_region$ts==min(df_region$ts)]<-df_region$Nt3[df_region$ts==min(df_region$ts)]
  df_region$new.K.4[df_region$ts==min(df_region$ts)]<-df_region$Nt4[df_region$ts==min(df_region$ts)]
  df_region$new.K.5[df_region$ts==min(df_region$ts)]<-df_region$Nt5[df_region$ts==min(df_region$ts)]
  df_region$new.sick[df_region$ts==min(df_region$ts)]<-df_region$sick[df_region$ts==min(df_region$ts)]
  df_region$new.sick.1[df_region$ts==min(df_region$ts)]<-df_region$sick1[df_region$ts==min(df_region$ts)]
  df_region$new.sick.2[df_region$ts==min(df_region$ts)]<-df_region$sick2[df_region$ts==min(df_region$ts)]
  df_region$new.sick.3[df_region$ts==min(df_region$ts)]<-df_region$sick3[df_region$ts==min(df_region$ts)]
  df_region$new.sick.4[df_region$ts==min(df_region$ts)]<-df_region$sick4[df_region$ts==min(df_region$ts)]
  df_region$new.sick.5[df_region$ts==min(df_region$ts)]<-df_region$sick5[df_region$ts==min(df_region$ts)]

  #Add ins for provincial sims
  df$ts<-0:(length(x_prov[,1])-1)

  df$sick<-rowSums(df[,colnames(df) %in% c(E,Di,Da)]); df$know.current<-rowSums(df[,colnames(df) %in% Tk]); df$recovered<-rowSums(df[,colnames(df) %in% R])

  df$known1<-rowSums(df[, colnames(df) %in% Tk[1:4]]);df$known2<-rowSums(df[, colnames(df) %in% Tk[5:8]]);df$known3<-rowSums(df[, colnames(df) %in% Tk[9:12]]);df$known4<-rowSums(df[, colnames(df) %in% Tk[13:16]]);df$known5<-rowSums(df[, colnames(df) %in% Tk[17:20]]);
  df$sick1<-rowSums(df[, colnames(df) %in% sick1]); df$sick2<-rowSums(df[, colnames(df) %in% sick2]); df$sick3<-rowSums(df[, colnames(df) %in% sick3]); df$sick4<-rowSums(df[, colnames(df) %in% sick4]); df$sick5<-rowSums(df[, colnames(df) %in% sick5]);
  file_names_list = strsplit(file.names[fit], "_")[[1]]
  str_arr_len = length(file_names_list)
  df$fit<-unlist(strsplit(file_names_list[str_arr_len], ".rds"));
  df$run<-run
  df$LL<-y$objective
  df$new.K<-df$Nt - lag(df$Nt);
  df$new.K.1<-df$Nt1 - lag(df$Nt1); df$new.K.2<-df$Nt2 - lag(df$Nt2); df$new.K.3<-df$Nt3 - lag(df$Nt3); df$new.K.4<-df$Nt4 - lag(df$Nt4); df$new.K.5<-df$Nt5 - lag(df$Nt5);
  df$new.K[1]<-df$Nt[1]; df$new.K.1[1]<-df$Nt1[1]; df$new.K.2[1]<-df$Nt2[1]; df$new.K.3[1]<-df$Nt3[1]; df$new.K.4[1]<-df$Nt4[1]; df$new.K.5[1]<-df$Nt5[1];
  df$new.sick<-(df$sick + df$recovered) -  lag(df$sick + df$recovered);
  df$new.sick.1<-(df$sick1 + df$R1) - lag(df$sick1  + df$R1); df$new.sick.2<-(df$sick2  + df$R2) - lag(df$sick2  + df$R2); df$new.sick.3<-(df$sick3  + df$R3) - lag(df$sick3  + df$R3); df$new.sick.4<-(df$sick4  + df$R4) - lag(df$sick4  + df$R4); df$new.sick.5<-(df$sick5  + df$R5) - lag(df$sick5  + df$R5);
  df$new.sick[1]<-df$sick[1]; df$new.sick.1[1]<-df$sick1[1]; df$new.sick.2[1]<-df$sick2[1]; df$new.sick.3[1]<-df$sick3[1]; df$new.sick.4[1]<-df$sick4[1]; df$new.sick.5[1]<-df$sick5[1];

  #Tag actual data so we can match based on days since 50th positive case
  model50<-df$ts[min(which(df$Nt >= 50))]
  df$ts50<-df$ts-model50

  #Tag actual data so we can match based on days since 50th positive case
  df_region.sum <- df_region %>%
    group_by(ts, fit, run) %>%
    dplyr::summarize(cases.K=sum(Nt))

  model50_region<-df_region$ts[min(which(df_region.sum$cases.K >= 50))]
  df_region$ts50<-df_region$ts-model50_region

  if (model50_region!=model50) { print("ts50 problem");}
  
  # Reduce the size of dataframes by taking out the state space (SEPAIR) parameters, and the travel matrix for df. 
  df.reduced = df %>% subset(select = !(names(.) %in% c(All1, All2, All3, All4, All5, propCits_names)))
  df_region.reduced = df_region %>% subset(select = !(names(.) %in% c(All1, All2, All3, All4, All5)))

  if (run==1) {
    df.prov.fit<-df.reduced 
    df.region.fit<-df_region.reduced
  } else {
    df.prov.fit<-rbind(df.prov.fit, df.reduced)
    df.region.fit<-rbind(df.region.fit, df_region.reduced)
  }

}

df.fitparms<-as.data.frame(fitparms)
file_names_list = strsplit(file.names[fit], "_")[[1]]
str_arr_len = length(file_names_list)
df.fitparms$fit<-unlist(strsplit(file_names_list[str_arr_len], ".rds"));
df.fitparms$LL<-y$objective

if (goodfits==1){
  df.prov.comb<-df.prov.fit
  df.region.comb<-df.region.fit
  df.fitparms.comb<-df.fitparms
} else {
  df.prov.comb<-rbind(df.prov.comb, df.prov.fit)
  df.region.comb<-rbind(df.region.comb, df.region.fit)
  df.fitparms.comb<-rbind(df.fitparms.comb, df.fitparms)
}

fit_num = length(unique(df.prov.comb$fit))
message("fit = ", sprintf("%s", fit_num), " finished")

goodfits<-goodfits+1;

  } 
}


saveRDS(df.region.comb, file = sprintf("%s/simsREGIONS_%s.rds", output_dir, codeselect));
saveRDS(df.prov.comb, file = sprintf("%s/simsPROVINCE_%s.rds", output_dir, codeselect));
saveRDS(df.fitparms.comb, file = sprintf("%s/simsFITPARMS_%s.rds", output_dir, codeselect));

##Slightly different version of the simsFITPARMS that stores the PHU-specific b1, b2, ..., b33, b34 values
colnames(fullparm.storage)<-c(fitnames, "LL")
fullparm.storage<-as.data.frame(fullparm.storage)
fullparm.storage$beta<-exp(fullparm.storage$beta)/mean(popadj2020)
fullparm.storage$omg<-exp(fullparm.storage$omg)
fullparm.storage[, c("zeta", "tau0", "psi1", "psi2", "psi3", "psi4", "psi5", "boost")]<-fullparm.storage[, c("zeta", "tau0", "psi1", "psi2", "psi3", "psi4", "psi5", "boost")]/1e3
saveRDS(fullparm.storage, file = sprintf("%s/simsFULLFITPARMS_%s.rds", output_dir, codeselect));

tfin<-Sys.time()

print(tfin-tstart)
