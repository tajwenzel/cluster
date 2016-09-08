#Source file for fluEvidenceSynthesis for changing likelihood. Options for varying age group sizes, risk groups, (and susceptibility, and ascertainment through the starting parameters). Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.

data("age_sizes")

build.llikelihood<-function()
  {
llikelihood.f <- function(pars, agegrouplimits, agegroupsizes, riskratios, season, strain,...)
  {
  proposed.contact.ids <<- current.contact.ids
    if (runif(1,0,1) < 1) {
      rs <<- round(runif(2,1,length(proposed.contact.ids)))
      proposed.contact.ids[rs[1]] <<- rs[2]}
    # Resample contact ids.
    contacts<- contact.matrix(as.matrix(polymod[proposed.contact.ids,]),
                              age_sizes[,1], agegrouplimits)
    
    
    age.groups <- stratify_by_age(age_sizes[,1], agegrouplimits)
    
    # Population sizes in each age and risk group
    popv <- stratify_by_risk(age.groups,riskratios) #agegroups*risk groups
    
    
    epsilons <- c(pars[1], pars[1],pars[2],pars[2], pars[3]); #epsilon is ascertainment
    
    initial.risk<-(rep(10^pars[9], length(age.groups)));
    initial.infected <- stratify_by_risk(initial.risk, riskratios) 
    odes <<- infectionODEs(popv, initial.infected,vaccine_calendar=vcalendar,contacts,c(pars[6],pars[6], pars[6], pars[6],pars[6],pars[7], pars[7], pars[7],pars[8], pars[8]),transmissibility = pars[5],infection_delays=c(0.8,1.8), interval=7) 
#interval is in days
    
    
    # Ignore times row
    #relevant.range<-(nrow(riskratios)*length(age.groups))+1
    #converted.odes <- odes[,2:(relevant.range)];
    converted.odes<-matrix(c(rep(0,52*5)),nrow=52,byrow=TRUE)
    dateaxis<-odes[,1]
    
    #Convert age groups and risk groups
    converted.odes[,1] <- rowSums(odes[,c(2,3)])+rowSums(odes[,c(12,13)]) 
    converted.odes[,2] <- rowSums(odes[,c(4,5)])+rowSums(odes[,c(14,15)])
    converted.odes[,3] <- rowSums(odes[,c(6,7,8)])+rowSums(odes[,c(16,17,18)])
    converted.odes[,4] <- odes[,9]+odes[,c(19)]
    converted.odes[,5] <- rowSums(odes[,c(10,11)])+rowSums(odes[,c(20,21)])
    converted.odes <- as.data.frame(converted.odes[,1:5])
    colnames(converted.odes)=c('V1','V2','V3','V4','V5')

    #matplot(dateaxis,converted.odes[,1:length(converted.odes)], type='l')
    #legend('topleft',legend=1:length(converted.odes),col=1:length(converted.odes), pch=2)
    # For each week and each group sum log likelihood
    
    #load in data from global
    positive<-as.matrix(virological$pos.by.strain[[strain]])
    total.sampled<-as.matrix(virological$no.samples)
    colnames(ili.counts$ili)<-c('V1','V2','V3','V4','V5')
    colnames(ili.counts$total.monitored)<-c('V1','V2','V3','V4','V5')
    ili.cases<-as.matrix(ili.counts$ili)
    ili.total<-as.matrix(ili.counts$total.monitored)
    age.groupLL <- stratify_by_age(age_sizes[,1], c(5,15,45,65))
    
    s.index<-(((season-1)*52)+1):(season*52)
     load('zero.compare')
    
    if(identical(positive[s.index,1:5],zero.compare)==TRUE) {ll<-0}
   else {
    ll<-log_likelihood_cases(
      epsilons,pars[4], as.matrix(converted.odes),
      age.groupLL, ili.cases[s.index,],ili.total[s.index,],
     positive[s.index,1:5], total.sampled[s.index,1:5])}
    return(ll)
  }
  return(llikelihood.f)
}
