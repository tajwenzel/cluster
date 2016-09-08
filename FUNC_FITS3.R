#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. August 2016

#FITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFIT



#######################################################################################################
.libPaths('/home/nwenzel/UKflu')
#library('devtools')
#install.packages('foreach',lib='~/UKflu')
#install.packages('doParallel',lib='~/UKflu')
#install.packages('parallel',lib='~/UKflu')
#install.packages('iterators',lib='~/UKflu')
#install.packages('doSNOW',lib='~/UKflu')
library('fluEvidenceSynthesis')
library(doParallel)
library(foreach)
#library(stringi)
library(plyr)
library(iterators)
library(pander)
library(data.table)
#install.packages('snow',lib='~/UKflu')
#library('snow')
install.packages('base',lib='~/UKflu')
library('base')
library(doSNOW)
rm(list = ls())


if(Sys.info()[7]=='Natasha') {clustertf<-0} else {(clustertf<-1)}

##FILE PATH & CLEAN UP
if(clustertf==0)
    {setwd("/Users/Natasha/Dropbox/UKfluworkGIT")
      } else {
        numCores <- 8
        cl <- makeCluster(numCores)
        registerDoParallel(cl)
        .libPaths('/home/nwenzel/UKflu')
        setwd('/home/nwenzel/UKflu')
        set.seed(round(runif(n=1, min=1,max=500)))
          }



          load('ili.counts.rda',.GlobalEnv)
          load('virological.rda',.GlobalEnv)
  
  
#load('/home/nwenzel/UKflu/ili.counts.rda',.GlobalEnv)
#load('/home/nwenzel/UKflu/virological.rda',.GlobalEnv)

  #####################################################################################
  # INPUT
  #####################################################################################

data("age_sizes")
age.group.limits<<-c(1,5,12,15,16,25,45,65,75) #upper limits
  
risk.ratios.ce<<-matrix(c(0.021,0.055,0.098,0.098,0.098,0.087,0.092,0.183,0.45,0.45,0,0,0,0,0,0,0,0,0,0),ncol=length(age.group.limits)+1, byrow=TRUE)  
  
age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)

polymod<<-fread(input='/home/nwenzel/UKflu/GBtable10.csv',sep = 'auto');
polymod<<-cbind(polymod,polymod$V11)
colnames(polymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
  current.contact.ids <<- seq(1,nrow(polymod))
  proposed.contact.ids <<- current.contact.ids
  
  vstrategy<<-dget(file='FUNC_cov_strategy.R')  
  load('coverageH1',.GlobalEnv);
  coverageH1<-cov.eff;
  load('coverageH3',.GlobalEnv); 
  coverageH3<-cov.eff;
  load('coverageB',.GlobalEnv); 
  coverageB<-cov.eff;
  cov.eff<-NULL
  cov.eff.data<<-list(coverageH1,coverageH3,coverageB)
  initial.parameters <- dget(file="INPUT_UKInitial.R")
  ########################################################################################
  # Seasonal vaccination plan, using function strain.choice
  ########################################################################################
season.choice<-dget('FUNC_seas_fit.R')

strains<-c(1:3)
seasons<-c(1:19)
jam<-expand.grid(seasons,strains)

mapply(FUN=season.choice,jam[[1]],jam[[2]],program=1)