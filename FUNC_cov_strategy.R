#Source file for fluEvidenceSynthesis with access to vaccine coverage options. Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.
#coveffin<-coverageB; season<-6; scenario<-program; ii<-1

vstrategy<-function(riskratio,scenario,coveffin,season)
{
   cov <- coveffin[[season]][,2:22]
    #USE COV INSTEAD OF COVERAGE SO DON'T HAVE TO DO THE /100 at every step
   non<-matrix(rep(0,length(cov[[1]])),nrow=length(cov[[1]]),byrow=TRUE)
##################################################################################################
#Efficacy
##################################################################################################
   
    eff.pull<-coveffin[[season]][1,23:43]
#pull first row of efficacies from matrix, matrix starts at 23 as 21 risk*age groups in original data +1 column for dates. This encompasses 2 risk*age groups (14 columns)
    
    dates <-as.Date(coveffin[[season]]$V1,origin="1970-01-01")

          calendar <- matrix(rep(0),nrow=length(dates),ncol = length(eff.pull))


############SETTING UP VACCINE PROGRAM, DEFINE WHICH AGE GROUPS ARE VACCINATED 
# Set rate of vaccine uptake for different dates/age groups
          #additional analysis for amount of QALY's potentially lost due to mis-reporting from practices not responding
          ####0) STATUS QUO is vaccination of elderly (65+) and high risk proportion (ages 0-85)

          if(scenario==1)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-cov[[8]] #0-1
            calendar[t,c(2)] <-cov[[9]] #1-4
            calendar[t,c(3)] <-cov[[10]] #5-11
            calendar[t,c(4)] <-cov[[10]] #12-14
            calendar[t,c(5)] <-cov[[10]] #15-16
            calendar[t,c(6)] <-cov[[11]] #17-24
            calendar[t,c(7)] <-cov[[12]] #25-44
            calendar[t,c(8)] <-cov[[13]] #45-64
            calendar[t,c(9)] <-cov[[14]] #65-74
            calendar[t,c(10)] <-cov[[14]] #75+
            calendar[t,c(11)] <-non
            calendar[t,c(12)] <-non
            calendar[t,c(13)] <-non
            calendar[t,c(14)] <-non
            calendar[t,c(15)] <-non
            calendar[t,c(16)] <-non
            calendar[t,c(17)] <-non
            calendar[t,c(18)] <-non
            calendar[t,c(19)] <-cov[[7]]
            calendar[t,c(20)] <-cov[[7]] 
          }
            
          ####1) Preschool Vaccination only; ages 1-4     
          if(scenario==2)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-cov[[2]] #1-4
            calendar[t,c(3)] <-non #5-11
            calendar[t,c(4)] <-non #12-14
            calendar[t,c(5)] <-non #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-cov[[9]] #1-4
            calendar[t,c(13)] <-non #5-11
            calendar[t,c(14)] <-non #12-14
            calendar[t,c(15)] <-non #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          
          ####2) Primary School Vaccination only; ages 5-11            
          if(scenario==3)    
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-non #1-4
            calendar[t,c(3)] <-cov[[3]] #5-11
            calendar[t,c(4)] <-non #12-14
            calendar[t,c(5)] <-non #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-non #1-4
            calendar[t,c(13)] <-cov[[10]] #5-11
            calendar[t,c(14)] <-non #12-14
            calendar[t,c(15)] <-non #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          ####3) Secondary school vaccination only; ages 12-16        
          if(scenario==4)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-non #1-4
            calendar[t,c(3)] <-non #5-11
            calendar[t,c(4)] <-cov[[3]] #12-14
            calendar[t,c(5)] <-cov[[3]] #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-non #1-4
            calendar[t,c(13)] <-non #5-11
            calendar[t,c(14)] <-cov[[10]] #12-14
            calendar[t,c(15)] <-cov[[10]] #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          
          ####4) Preschool+Primary School Vaccination; ages 1-11          
          if(scenario==5)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-cov[[2]] #1-4
            calendar[t,c(3)] <-cov[[3]] #5-11
            calendar[t,c(4)] <-non #12-14
            calendar[t,c(5)] <-non #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-cov[[9]] #1-4
            calendar[t,c(13)] <-cov[[10]] #5-11
            calendar[t,c(14)] <-non #12-14
            calendar[t,c(15)] <-non #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          ####5) Preschool+Primary+Secondary school vaccination; ages 1-16 
          if(scenario==6)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-cov[[2]] #1-4
            calendar[t,c(3)] <-cov[[3]] #5-11
            calendar[t,c(4)] <-cov[[3]] #12-14
            calendar[t,c(5)] <-cov[[3]] #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-cov[[9]] #1-4
            calendar[t,c(13)] <-cov[[10]] #5-11
            calendar[t,c(14)] <-cov[[10]] #12-14
            calendar[t,c(15)] <-cov[[10]] #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          
          ####6) Preschool+Secondary school vaccination; ages 1-4, 12-16                   
          if(scenario==7)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-cov[[2]] #1-4
            calendar[t,c(3)] <-non #5-11
            calendar[t,c(4)] <-cov[[3]] #12-14
            calendar[t,c(5)] <-cov[[3]] #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-cov[[9]] #1-4
            calendar[t,c(13)] <-non #5-11
            calendar[t,c(14)] <-cov[[3]] #12-14
            calendar[t,c(15)] <-cov[[3]] #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          
          ####7) Primary+Secondary school vaccination; ages 5-16
          if(scenario==8)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-non #0-1
            calendar[t,c(2)] <-cov[[2]]#1-4
            calendar[t,c(3)] <-cov[[3]] #5-11
            calendar[t,c(4)] <-cov[[3]] #12-14
            calendar[t,c(5)] <-cov[[3]] #15-16
            calendar[t,c(6)] <-non #17-24
            calendar[t,c(7)] <-non #25-44
            calendar[t,c(8)] <-non #45-64
            calendar[t,c(9)] <-non #65-74
            calendar[t,c(10)] <-non #75+
            calendar[t,c(11)] <-non #0-1
            calendar[t,c(12)] <-cov[[9]] #1-4
            calendar[t,c(13)] <-cov[[10]] #5-11
            calendar[t,c(14)] <-cov[[10]] #12-14
            calendar[t,c(15)] <-cov[[10]] #15-16
            calendar[t,c(16)] <-non #17-24
            calendar[t,c(17)] <-non #25-44
            calendar[t,c(18)] <-non #45-64
            calendar[t,c(19)] <-non #65-74
            calendar[t,c(20)] <-non #75+
          }
          #compile vaccine schedule
        #compile vaccine schedule
                   
vaccine3 <- as.vaccination.calendar(efficacy = eff.pull,dates = dates,coverage = calendar, no_risk_groups=2,no_age_groups=10)

eff.out<-as.vector(unlist(vaccine3$efficacy),mode='numeric')

v.output<-list('efficacy'=eff.out,'dates'=vaccine3$dates,'calendar'=vaccine3$calendar)

return(v.output)
}