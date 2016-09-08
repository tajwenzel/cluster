#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. August 2016

#FITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFIT



#######################################################################################################
library(devtools)
.libPaths('/home/nwenzel/UKflu')
#install.packages('foreach',lib='~/UKflu')
#install.packages('doParallel',lib='~/UKflu')
#install.packages('parallel',lib='~/UKflu')
#install.packages('iterators',lib='~/UKflu')
#install.packages('doSNOW',lib='~/UKflu')
.libPaths('/home/nwenzel/UKflu')
library(fluEvidenceSynthesis)
library(doParallel)
library(foreach)
#library(stringi)
library(plyr)
library(iterators)
library(pander)
library(data.table)
#install.packages('snow',lib='~/UKflu')
library('snow')
library(doSNOW)


if(Sys.info()[7]=='Natasha') {clustertf<-0} else {(clustertf<-1)}

##FILE PATH & CLEAN UP
if(clustertf==0)
    {setwd("/Users/Natasha/Dropbox/UKfluworkGIT")
      } else {
    numCores <- 3
        cl <- makeCluster(numCores)
        registerDoParallel(cl)
        .libPaths('/home/nwenzel/UKflu')
        rm(list = ls())
        setwd('/home/nwenzel/UKflu')
        set.seed(round(runif(n=1, min=1,max=500)))
          }

#######################################################################
#IMPORT DATA
#####################################################################
strain.choice<-function(sti)
 {         

load('/home/nwenzel/UKflu/ili.counts.rda',.GlobalEnv)
load('/home/nwenzel/UKflu/virological.rda',.GlobalEnv)}

  #####################################################################################
  # INPUT
  #####################################################################################
  
  initial.parameters <- dget(file="/home/nwenzel/UKflu/INPUT_UKInitial.R") #unknowns and known
  vstrategy<-dget(file='/home/nwenzel/UKflu/FUNC_cov_strategy.R')
  
  scenario<<-1
  data("age_sizes")

#respecify groups
age.group.limits<-c(1,5,12,15,16,25,45,65,75) #upper limits
  
  
risk.ratios.ce<-matrix(c(0.021,0.055,0.098,0.098,0.098,0.087,0.092,0.183,0.45,0.45,0,0,0,0,0,0,0,0,0,0),ncol=length(age.group.limits)+1, byrow=TRUE)  
  
age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)

polymod<<-fread(input="GBtable10.csv",sep = 'auto');
load('/home/nwenzel/UKflu/coverageH1',.GlobalEnv); coverageH1<-cov.eff;
load('/home/nwenzel/UKflu/coverageH3',.GlobalEnv); coverageH3<-cov.eff;
load('/home/nwenzel/UKflu/coverageB',.GlobalEnv); coverageB<-cov.eff;
cov.eff<-NULL
      
 
  polymod<<-cbind(polymod,polymod$V11)
  colnames(polymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
  current.contact.ids <<- seq(1,nrow(polymod))
  proposed.contact.ids <<- current.contact.ids
  
  
  ########################################################################################
  # Seasonal vaccination plan, using function strain.choice
  ########################################################################################

  cov.eff.data<-list(coverageH1,coverageH3,coverageB)
  cov.eff.in<-cov.eff.data[[sti]]
  
for(sea in 1:19)
 {sea
        vcalendar<-vstrategy(risk.ratios.ce,scenario,cov.eff.in,sea) 
    
    
    ####################################################################################################
    #Likelihood constant
    ####################################################################################################
       buildLL<-dget(file='FUNC_LL.R')
        llikelihood<-buildLL()
    
        llprior <- function(pars) {
        if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
          || pars[9] < log(0.00001) || pars[9] > log(10) )
        return(-Inf)
      
        lprob <- dnorm(pars[5], 0.1653183, 0.02773053, 1) 
        + dlnorm(pars[1], -4.493789, 0.2860455, 1) 
        + dlnorm(pars[2], -4.117028, 0.4751615, 1) 
        + dlnorm(pars[3], -2.977965, 1.331832, 1);
      
        return(lprob)
      }
    
    ####################################################################################################
    #MCMC run iteration 
    burnin<-30; #potatoes
    out<-150; #meat
    saveiteration<-1; #thin size
    
    contact.ids<-list()
    # Run adaptive.mcmc
    ptm <- proc.time() 
    mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, 
                                 nburn=burnin,
                                 initial = initial.parameters,
                                 nbatch = out, blen = saveiteration,
                                 outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids}, acceptfun= function() {current.contact.ids <<- proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce, season=sea, strain=sti, polymod=polymod)
    proc.time() - ptm

    
    strain.name<-c('H1','H3','B')
    date.labels<-format(as.Date(cov.eff.in[[sea]]$V1, origin="1970-01-01"), "%Y")
    
    save(contact.ids,file=paste0('ct.ids',date.labels[1],last(date.labels),strain.name[sti]))
    save(mcmc.result,file=paste0('flu',date.labels[1],last(date.labels),strain.name[sti]))
      }
      #flu.season<-c(1:19)
  #sapply(flu.season,season.choice,strain)
}
  
#select<-c(1:3)
sapply(sti=2, strain.choice)
#foreach(strain=c(1:3),.packages=c('fluEvidenceSynthesis','plyr','pander','data.table', 'Rcpp','tidyr','codetools','curl','gtable','scales','RcppEigen','R6','withr','BH','iterators','chron','codetools','colorspace','digest','labeling','magrittr','munsell','tibble')) %dopar% strain.choice(strain)
