season.choice<-function(sea,ii,program)
{  
cov.eff.in<-cov.eff.data[[ii]]
vcalendar<<-vstrategy(risk.ratios.ce,program,cov.eff.in,sea)

strain.name<-c('H1','H3','B')  
date.labels<-format(as.Date(cov.eff.in[[sea]]$V1, origin="1970-01-01"), "%Y") 
ctname<-paste0('ct.ids',date.labels[1],last(date.labels),strain.name[ii])
rname<-paste0('flu',date.labels[1],last(date.labels),strain.name[ii])  
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
 
  burnin<-10000
  out<-50000
  saveiteration<-1
  
  contact.ids<-list()
  
  mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, nburn=burnin,initial = initial.parameters,nbatch = out, blen = saveiteration,outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids}, acceptfun= function() {current.contact.ids <<-proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce, season=(sea), strain=(ii), polymod=polymod)

    
  save(contact.ids,file=ctname)
  save(mcmc.result,file=rname)

 }

