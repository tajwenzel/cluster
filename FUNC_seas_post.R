season.choice<-function(sea,ii,program)
{  
cov.eff.in<-cov.eff.data[[ii]]
vcalendar<<-vstrategy(risk.ratios.ce,program,cov.eff.in,sea)
master.name<-('flu.t')
strain.name<-c('H1','H3','B')  
date.labels<-format(as.Date(cov.eff.in[[sea]]$V1, origin="1970-01-01"), "%Y") 
ctname<-paste0('ct.ids',date.labels[1],last(date.labels),strain.name[ii])
rname<-paste0('flu',date.labels[1],last(date.labels),strain.name[ii])  
#####################################################################################
# Posterior from previous year
#---must make allowance for year 1, 1995
#strain.pull<-c(glob2rx('flu*H1'),glob2rx('flu*H3'),glob2rx('flu*B'))
#list.files(pattern=strain.pull[ii])
#sea<-11; ii=1; program=1
#load(flu.tablesH1[sea])
#setwd('/Users/Natasha/Dropbox/UKfluworkGIT/cluster')
#####################################################################################

if(sea > 1)
  {load(paste0(master.name,(as.numeric(date.labels[1])-1),(as.numeric(last(date.labels))-1),strain.name[ii]));  sub.year<-(as.numeric(date.labels[1])-1); print('going back 1 year'); print(sub.year)

  if(last(mcmc.result$llikelihoods)==0)
  {load(paste0(master.name,(as.numeric(date.labels[1])-2),(as.numeric(last(date.labels))-2),
               strain.name[ii])); print('zero incidence year; going back 2 years'); 
              sub.year<-(as.numeric(date.labels[1])-2); print(sub.year)}
  
        if(last(mcmc.result$llikelihoods)==0)
        {load(paste0(master.name,(as.numeric(date.labels[1])-3),(as.numeric(last(date.labels))-3),
                     strain.name[ii])); print('zero incidence year; going back 3 years');
          sub.year<-(as.numeric(date.labels[1])-3); print(sub.year)}
  
            if(last(mcmc.result$llikelihoods)==0)
            {load(paste0(master.name,(as.numeric(date.labels[1])-4),(as.numeric(last(date.labels))-4),
              strain.name[ii])); print('zero incidence year; going back 4 years');  
              sub.year<-(as.numeric(date.labels[1])-4); print(sub.year)}
  
              if(last(mcmc.result$llikelihoods)==0)
              {load(paste0(master.name,(as.numeric(date.labels[1])-5),
              (as.numeric(last(date.labels))-5),strain.name[ii])); 
              print('zero incidence year; going back 5 years'); sub.year<-(as.numeric(date.labels[1])-5); print(sub.year)}

                  if(last(mcmc.result$llikelihoods)==0)
                  {load(paste0(master.name,(as.numeric(date.labels[1])-6),
                   (as.numeric(last(date.labels))-6),strain.name[ii])); 
                    print('zero incidence year; going back 6 years'); 
                    ub.year<-(as.numeric(date.labels[1])-6); print(sub.year)}

      med.f<-function(cat) {median(mcmc.result$batch[,cat])}
      initial.parameters<<-sapply(c(1:dim(mcmc.result$batch)[2]),med.f)

      sd.f<-function(dog) {sd(mcmc.result$batch[,dog])} 
      sd.pull<-sapply(c(1:dim(mcmc.result$batch)[2]),sd.f)
       
      mean.f<-function(mouse) {mean(mcmc.result$batch[,mouse])} 
      mean.pull<-sapply(c(1:dim(mcmc.result$batch)[2]),mean.f)
      
      length.f<-dim(mcmc.result$batch)[1]
  
      mcmc.result<-NULL
  }

####################################################################################################
  #Likelihood constant
  ####################################################################################################
  buildLL<-dget(file='FUNC_LL.R')
  llikelihood<-buildLL()
 
  
  if(sea==1){
    llprior <- function(pars) {
    if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
        || pars[9] < log(0.00001) || pars[9] > log(10) )
      return(-Inf)
    
    lprob <- dnorm(pars[5], 0.1653183, 0.02773053, 1) 
    + dlnorm(pars[1], -4.493789, 0.2860455, 1) 
    + dlnorm(pars[2], -4.117028, 0.4751615, 1) 
    + dlnorm(pars[3], -2.977965, 1.331832, 1);
    
    return(lprob)}
  }else{
    options(warn=-1);
  llprior <- function(pars) 
     {
    if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
        || pars[9] < log(0.00001) || pars[9] > log(10) )
      return(-Inf)
    
    parm5<-dnorm(pars[5], mean.pull[5], sd.pull[5], 1) 
    parm1<-dlnorm(pars[1], log(mean.pull[1]), log(sd.pull[1]), 1) 
    parm2<-dlnorm(pars[2], log(mean.pull[2]), log(sd.pull[2]), 1) 
    parm3<-dlnorm(pars[3], log(mean.pull[3]), log(sd.pull[3]), 1)
    
    
    if(is.nan(parm1) || is.nan(parm5) || is.nan(parm2) ||is.nan(parm3))
      {
    lprob <- dnorm(pars[5], mean.pull[5], sd.pull[5], 1) 
    + dt(pars[1], length.f-1, ncp=log(mean.pull[1]), 1) 
    + dt(pars[2], length.f-1, ncp=log(mean.pull[2]), 1) 
    + dt(pars[3], length.f-1,ncp=log(mean,pull[3]), 1)
    print('used students t instead')
      } else {lprob<-parm5+parm1+parm2+parm3}
    
    return(lprob)
    options(warn=0)
    }
  }
  
  ####################################################################################################
 
  burnin<-10
  out<-50
  saveiteration<-1
  
  contact.ids<-list()
  
  mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, nburn=burnin,initial = initial.parameters,nbatch = out, blen = saveiteration,outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids}, acceptfun= function() {current.contact.ids <<-proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce, season=(sea), strain=(ii), polymod=polymod)

    
  save(contact.ids,file=ctname)
  save(mcmc.result,file=rname)

 }

