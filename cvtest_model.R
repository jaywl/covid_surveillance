########################################################################################
################### model to simulate effect of COVID-19 testing strategies ############
###################         Jay Love, Lindsay Keegan                        ############
########################################################################################

##### README: This model is a hybrid. It includes both population-level (notably, 
##### infections occur via population-level parameters and contacts are not 
##### discretely modeled) and individual-based (recovery, symptom state, etc) dynamics.
##### It isn't very efficient, as it's written with a nested loop structure that iterates
##### across testing proportions, testing strategies, sites, simulations, days, and people.
##### Yes, that's a lot of loops.

##wrap in a function
cvtest_epi<-function(nstrats,nsims,nsites,ntesting_proportions,
                     testprops,
                     quarantine_days,path,
                     pcr.sensitivity,
                     pcr.specificity,
                     antigen.sensitivity,
                     antigen.specificity,
                     theta.pcr,
                     theta.antigen,
                     prob_symptomatic,
                     ILI_incidence,
                     waiting.reduction,
                     r.time,
                     include.rs){ 
  ## load dependencies
  require(data.table)
  nstrats<-nstrats #number of testing strategies
  nsims<-nsims # number of simulations per strategy
  nsites<-nsites # number of sites
  ntestingproportions<-ntesting_proportions # number of testing proportions
  quarantine_days<-quarantine_days
  pcr.sensitivity<-pcr.sensitivity
  pcr.specificity<-pcr.specificity
  antigen.sensitivity<-antigen.sensitivity
  antigen.specificity<-antigen.specificity
  waiting.reduction<-waiting.reduction
  r.time<-r.time ## for ow many days following recovery are infections detectable by PCR? (if include.rs==TRUE)
  incluide.rs<-include.rs ## logical: do you want to include Recovered class detected infections as true positives? 
  t1<-testprops[1]
  t2<-testprops[2]
  t3<-testprops[3]
  t4<-testprops[4]
  t5<-testprops[5]
  ### must copy/paste above 3 lines, replacing "5" with "n" for n testing proportions
  
  step.size<-1 ## this should be 1 unless you want to do something weird
  
  ## to run simulations at different testing proportions
  for(o in 1:ntestingproportions){
    if(o==1){
      daily.testing.proportion<-t1
    }
    if(o==2){
      daily.testing.proportion<-t2
    }
    if(o==3){
      daily.testing.proportion<-t3
    }
    if(o==4){
      daily.testing.proportion<-t4
    }
    if(o==5){
      daily.testing.proportion<-t5
    }
    ##.. etc. ### must copy/paste above 3 lines, replacing "5" with "n" for n testing proportions
    
    
    ## to run simulations at different sites:
    for(u in 1:nsites){
      ## create accounting data frame to be rebuilt for each site:
      megs<-data.frame(rep=1:(nstrats*nsims),test.strategy=NA,location=NA,infectious_period=NA,gamma=NA,R0=NA,prop.testing=NA,infections=NA,missed.infections=NA,hospitalizations=NA,ICU=NA,
                       deaths=NA,days.quarantine=NA,days.false.quarantine=NA, time.waiting=NA,
                       n.reflexes=NA,n.tests=NA,n.alt.tests=NA)
      ##create combo data frame with all simulations to make epi curves. This file will be large if running many simulations.
      sum4plots<-data.frame(t=0,S=0,E=0,I=0,R=0,new.infs=0,tests=0,retests=0,missed=0,caught=0,quarantines=0,
                            false_quarantines=0,p_expose=0,
                            qS=0,qE=0,qI=0,qR=0,positive.tests=0,negative.tests=0,
                            tS=0,tE=0,tI=0,tR=0,site=NA,test.strategy=NA,prop.testing=NA)
      
      
      if(u==1){
        site<-"NH"
      }
      if(u==2){
        site<-"dorm"
      }
      
      for(i in 1:(nstrats*nsims)){
        
        no.testing<-FALSE
        ## establish testing strategy sequence for simulations
        if(i>0 & i<=nsims){
          test.type<-"PCR"
        }
        
        if(i>nsims & i<=2*nsims){
          test.type<-"antigen"
        }
        if(i>2*nsims & i<=3*nsims){
          test.type<-"retest all negatives"
        }
        if(i>3*nsims & i<=4*nsims){
          test.type<-"retest negative symptomatic PCR"
        }
        if(i>4*nsims & i<=5*nsims){
          no.testing<-TRUE
          test.type<-"none"
        }
        
        ## catch in case testing proportion is set to zero
        if(daily.testing.proportion==0){
          no.testing<-TRUE
          test.type<-"none"
          if(i>nsims){ ## end the testing for proportion==0 after 1*nsims
            break
          }
        }
        
        # ############### OVERRIDEEEEEE TESTING STRATEGY ############################################################
        # #########################################################################################################
        # (uncomment to override) 
        #
        # test.type<-"retest all negatives"
        # no.testing<-FALSE
        # 
        # ############### OVERRIDEEEEEE TESTING STRATEGY ############################################################
        # #########################################################################################################
        
        
        pop.size<-ifelse(site=="NH",round(86*1.18),3150) # total population size 
        simdur<-183 ## days (or time steps) to simulate
        community.prevalence<-0.018
        n.seed.events <- round(rnorm(1,mean=ifelse(site=="NH",2*community.prevalence*pop.size,community.prevalence*pop.size),
                                     sd=.3)) # number of infected people at time t=0
        
        # # define SEIR parameters
        sigma <- 1/3  ## we're assuming I to mean infectious rather than infected, so this is "parent period". 
        incubation<-runif(pop.size,min=5,max=6) ## days from exposure to symptom onset (generating a distribution for individuals)
        infectious.period.min<-2.6 ## infectious period low
        infectious.period.max<-6 ## infectious period high
        gamma <- 1/runif(pop.size,infectious.period.min,infectious.period.max)  ## recovery rate
        R0 <- runif(pop.size,1.2,1.5) ##basic reproduction number (distribution)
        beta <- gamma*R0 ## beta is the average number of infection-producing contacts each person makes per unit time. 
        omega <- 1/quarantine_days ## quarantine factor (1/time in quarantine)
        theta.pcr <- theta.pcr ## rate of return of test results for PCR (extra param for special case of test-switch reflexing)
        theta.antigen<-theta.antigen ## for antigen
        theta <- ifelse(test.type=="PCR",theta.pcr,theta.antigen) ## test turn around (1/processing time for test) ## this will change by test type
        tau <- -log(1-daily.testing.proportion) ## -log(proportion of state compartment tested) \/
        phi.p <- ifelse(test.type=="PCR",pcr.specificity,antigen.specificity) ## specificity ## will change by test type: ifelse(test.type=="PCR",something,somethingelse)
        phi.p.pcr <- pcr.specificity # or whatever pcr test specificity is (for use in special case)
        phi.e <- ifelse(test.type=="PCR",pcr.sensitivity,antigen.sensitivity) ## sensitivity ## will change by test type
        phi.e.pcr <- pcr.sensitivity # or whatever pcr test sensitivity is (for use in special case)
        phi.S<-1/10000 ## probability of testing positive when state is S
        q<-waiting.reduction ## reduction in beta due to testing 
      ####### HARDCODE FLAG ##########
      ####### retest.delay sets time steps between retesting in reflex testing strategies.
        retest.delay<-2 ##
      ####### end hardcode flag
        
        if(no.testing==TRUE){
          tau<-0
          test.type<-"none"
        }
        
        ppl<-data.frame(person=1:pop.size,state="S",recovery.day=0,symptom=0,symptom.onset.day=0,ILI=0,tested=0,day.tested=0,new.test=0,waiting=0,time.waiting=0,test.returned=0,test.result=2, needs.retest=0,
                        quarantine=0,new.retest=0,retested=0,retest.day=0,ever.retested=0,ever.infected=0,ever.true.positive=0)
        ### calculate probabilities that will be used consistently throughout
        ### SEIR disease course 
        p.infect <- 1-exp(-step.size*sigma)
        p.recover <- 1-exp(-step.size*gamma)
        p.test <- 1-exp(-step.size*tau) 
        p.return <- 1-exp(-step.size*theta) # probability of a test result bac ###### NOTE: currently just using theta below. I think that makes sense: the rapid testing is rapid... as in same-day test results. with this modified p.return, only two-thirdsish will be returned same day.
        ## probability of leaving/ending quarantine
        p.end.isolate <- 1-exp(-step.size*omega) ## probability of quarantined person coming out of quarantine
        ## probability of being symptomatic 
        p.symptomatic<- prob_symptomatic 
        ## incidence of Influenza-Like Illness in population
        ILI.incidence<-ILI_incidence
        
        
        ## everyone starts as susceptible and untested
        ppl$state<-"S"
        ppl$tested<-0
        
        ##seed with initial infections
        ppl$state[sample(ppl$person,size=n.seed.events)]<-"I"
        ppl$recovery.day[ppl$state=="I"]<-round(runif(length(ppl$recovery.day[ppl$state=="I"]),infectious.period.min,infectious.period.max))
        
        
        ## summary data frame
        summit<-data.frame(t=1:simdur,S=0,E=0,I=0,R=0,new.infs=0,tests=0,retests=0,missed=0,caught=0,quarantines=0,
                           false_quarantines=0,p_expose=NA,
                           qS=0,qE=0,qI=0,qR=0,positive.tests=0,negative.tests=0,
                           tS=0,tE=0,tI=0,tR=0,site=NA,test.strategy=NA,prop.testing=NA)
        
        for(t in 1:simdur){## loop of days
          ## fill summary row for day's starting conditions
          summit$S[t]<-length(ppl$person[ppl$state=="S"])
          summit$E[t]<-length(ppl$person[ppl$state=="E"])
          summit$I[t]<-length(ppl$person[ppl$state=="I"])
          summit$R[t]<-length(ppl$person[ppl$state=="R"])
          summit$caught[t]<-length(ppl$person[ppl$ever.true.positive==1 & ppl$state=="R"])
          summit$missed[t]<-length(ppl$person[ppl$state=="R"])-length(ppl$person[ppl$ever.true.positive==1 & ppl$state=="R"])
          summit$tests[t]<-length(ppl$person[ppl$new.test==1])+length(ppl$person[ppl$new.retest==1])
          summit$positive.tests[t]<-length(ppl$person[ppl$test.result==1])
          summit$negative.tests[t]<-summit$tests[t]-summit$positive.tests[t]
          summit$retests[t]<-length(ppl$person[ppl$new.retest==1])
          summit$quarantines[t]<-length(ppl$person[ppl$quarantine==1])
          summit$qS[t]<-length(ppl$person[ppl$quarantine==1 & ppl$state=="S"])
          summit$qE[t]<-length(ppl$person[ppl$quarantine==1 & ppl$state=="E"])
          summit$qI[t]<-length(ppl$person[ppl$quarantine==1 & ppl$state=="I"])
          summit$qR[t]<-length(ppl$person[ppl$quarantine==1 & ppl$state=="R"])
          
          summit$tS[t]<-length(ppl$person[ppl$tested>0 & ppl$state=="S"])
          summit$tE[t]<-length(ppl$person[ppl$tested>0 & ppl$state=="E"])
          summit$tI[t]<-length(ppl$person[ppl$tested>0 & ppl$state=="I"])
          summit$tR[t]<-length(ppl$person[ppl$tested>0 & ppl$state=="R"])
          
          summit$false_quarantines[t]<-length(ppl$person[ppl$quarantine==1 & ppl$state!="I"])
          summit$prop.testing[t]<-daily.testing.proportion
          
          ## reset new.test to zero if tests were returned yesterday
          ppl$new.test[ppl$test.returned==1]<-0
          ## then reset test.returned to zero after tests have been returned and results determined from yesterday
          ppl$test.returned<-0
          ## then reset new.retest
          ppl$new.retest<-0
          ## then reset test.result
          ppl$test.result<-2 ## 2 just means neither 1 nor 0 (to remove problems with using NAs)
          
          ## people who are due to develop symptoms today do so
          ppl$symptom[ppl$symptom.onset.day==t]<-1
          
          ## and some people will have non-covid symptoms (Influenza-Like Illness)
          ppl$ILI<-rbinom(length(ppl$person),1,ILI.incidence) ## this is just random ILI that lasts only 1 day. Not the most realistic, but adds a level of realism to the reflect symptomatic PCR strategy
          
          ## calculate p.expose as the density-dependent force of infection \/\/\/
          p.expose<-((beta*length(ppl$person[ppl$state=="I" & ppl$waiting==0 & ppl$quarantine!=1]))/pop.size)+
            ((beta*q*length(ppl$person[ppl$state=="I" & ppl$waiting==1 & ppl$quarantine!=1]))/pop.size)
          
          summit$p_expose[t]<-mean(p.expose) ## record the mean
          
          ###### testing process
          ## First, we test some new people, who could have already been tested
          ppl$new.test<-rbinom(length(ppl$person),1,p.test)
          ppl$test.day[ppl$new.test==1]<-t
          ppl$tested<-ppl$tested+ppl$new.test
          ## and some others who need retesting
          ppl$new.retest[ppl$needs.retest==1 & ppl$retest.day==t]<-1
          ##clear that those who have a new retest today are in need of retest
          ppl$needs.retest[ppl$new.retest==1]<-0
          ppl$retested[ppl$new.retest==1]<-1+(ppl$retested[ppl$new.retest==1])
          ppl$ever.retested[ppl$new.retest==1]<-1
          
          ## waiting column will be 1 if awaiting test (or retest) results
          ppl$waiting[ppl$new.test==1 | ppl$new.retest==1]<-1
          
          
          ## then, we return test results from today (if lag=0) and previous days (if lag>0) for people who have been tested and awaiting results
          ppl$test.returned[ppl$waiting==1] <- rbinom(length(ppl$person[ppl$waiting==1]),1,theta)#,p.return) ## can replace with p.return if you want. ## if 1, test results have been returned
          
          ## but in this special case of test-switch reflexing, we overwrite test.returned for those retesting
          if(test.type=="retest negative symptomatic PCR"){ 
            ppl$test.returned[ppl$waiting==1 & ppl$new.retest==1] <- rbinom(length(ppl$person[ppl$waiting==1 & ppl$new.retest==1]),1,theta.pcr)## if 1, test results have been returned
          }
          
          
          ## then, we find out the outcome of their test. test.result==1 means positive test, ==0 means negative test. ==2 means nothing (just using to prevent problems due to NAs)
          if(test.type!="retest negative symptomatic PCR"){ ## for all test strategies other than type-switching reflexing
            for(p in 1:length(ppl$person)){ ##loop of people 
              if(ppl$test.returned[p]==1 ){ ## if person just had their test returned
                if(ppl$state[p]=="S"){
                  ppl$test.result[p]<-rbinom(1,1,phi.S) ## test rarely gives false positives
                }
                
                if(ppl$state[p]=="E"){
                  ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (we consider E not infectious, therefore a positive test result is false)
                }
                
                if(ppl$state[p]=="I"){
                  ppl$test.result[p]<-rbinom(1,1,phi.e) ## test rarely gives false negatives
                  if(ppl$test.result[p]==1){ ## if true positive, then record that this infection was caught
                    ppl$ever.true.positive[p]<-1
                  }
                }
                
                if(ppl$state[p]=="R"){
                  if(include.rs==FALSE){
                    ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (again, we consider R like we do E)
                  }
                  if(include.rs==TRUE){
                    if(ppl$recovery.day[p]<=(t-r.time)){
                      ppl$test.result[p]<-rbinom(1,1,phi.e) ## test rarely gives false positives (again, we consider R like we do E)
                    }
                    if(ppl$recovery.day[p]>(t-r.time)){
                      ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (again, we consider R like we do E)
                    }
                  }## end include Rs
                }
              }## end test returned contingency
            }## end loop of people
          }## end if contingency to specify all testing strategies other than type-switching reflex
          
          ###### if using strategy "retest negative symptomatic PCR", then we are going to do something extra
          if(test.type=="retest negative symptomatic PCR"){
            # print("alert") ## for debugging
            for(p in 1:length(ppl$person)){ ##loop of people 
              if(ppl$new.retest[p]==1){ ## if being retested
                if(ppl$test.returned[p]==1 ){ ## if person just had their test returned
                  if(ppl$state[p]=="S"){
                    ppl$test.result[p]<-rbinom(1,1,phi.S) ## test rarely gives false positives
                  }
                  
                  if(ppl$state[p]=="E"){
                    ppl$test.result[p]<-rbinom(1,1,1-phi.p.pcr) ## test rarely gives false positives (we consider E not infectious, therefore a positive test result is false)
                  }
                  
                  if(ppl$state[p]=="I"){
                    ppl$test.result[p]<-rbinom(1,1,phi.e.pcr) ## test rarely gives false negatives
                    if(ppl$test.result[p]==1){ ## if true positive, then record that this infection was caught
                      ppl$ever.true.positive[p]<-1
                    }
                  }
                  
                  if(ppl$state[p]=="R"){
                    if(include.rs==FALSE){
                      ppl$test.result[p]<-rbinom(1,1,1-phi.p.pcr) ## test rarely gives false positives (again, we consider R like we do E)
                    }
                    if(include.rs==TRUE){
                      if(ppl$recovery.day[p]<=(t-r.time)){
                        ppl$test.result[p]<-rbinom(1,1,phi.e.pcr) ## test rarely gives false positives (again, we consider R like we do E)
                      }
                      if(ppl$recovery.day[p]>(t-r.time)){
                        ppl$test.result[p]<-rbinom(1,1,1-phi.p.pcr) ## test rarely gives false positives (again, we consider R like we do E)
                      }
                    } ## end include Rs in sensitivity/specificity
                  }
                }## end test returned contingency
              }## end new retest contingency
              
              if(ppl$new.retest[p]!=1){ ## if not being retested then do as normal::\/
                if(ppl$test.returned[p]==1 ){ ## if person just had their test returned
                  if(ppl$state[p]=="S"){
                    ppl$test.result[p]<-rbinom(1,1,phi.S) ## test rarely gives false positives
                  }
                  
                  if(ppl$state[p]=="E"){
                    ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (we consider E not infectious, therefore a positive test result is false)
                  }
                  
                  if(ppl$state[p]=="I"){
                    ppl$test.result[p]<-rbinom(1,1,phi.e) ## test rarely gives false negatives
                    if(ppl$test.result[p]==1){ ## if true positive, then record that this infection was caught
                      ppl$ever.true.positive[p]<-1
                    }
                  }
                  
                  if(ppl$state[p]=="R"){
                    if(include.rs==FALSE){
                      ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (again, we consider R like we do E)
                    }
                    if(include.rs==TRUE){
                      if(ppl$recovery.day[p]<=(t-r.time)){
                        ppl$test.result[p]<-rbinom(1,1,phi.e) ## test rarely gives false positives (again, we consider R like we do E)
                      }
                      if(ppl$recovery.day[p]>(t-r.time)){
                        ppl$test.result[p]<-rbinom(1,1,1-phi.p) ## test rarely gives false positives (again, we consider R like we do E)
                      }
                    }
                    
                  }
                }## end test returned contingency
              } ## end if not being retested
            }## end loop of people
          }## end special case for test.type==retest negative symptomatic PCR
          
          
          ## people who test positive go into quarantine immediately
          ppl$quarantine[ppl$test.returned==1 & ppl$test.result==1]<-1
          ppl$quarantine[ppl$test.returned==1 & ppl$test.result==0]<-0 ## but people who get negative results back while in quarantine come out of quarantine

          ## and people who have test results back, tested negative but are symptomatic get set up to be retested after delay
          if(test.type=="retest negative symptomatic"){
            for(p in 1:length(ppl$person[p])){
              if(ppl$new.test[p]==1){
                if(ppl$test.returned[p]==1){
                  if(ppl$test.result[p]==0){
                    if(ppl$symptom[p]==1 | ppl$ILI[p]==1){
                      #if(ppl$done.testing[p]!=1){
                      ppl$needs.retest[p]<-1
                      ppl$retest.day[p]<-t+retest.delay
                      #}
                    }
                  }
                }
              }
            }
          } ## end retest negative symptomatic stuff
          
          ## same as above, but with pcr
          if(test.type=="retest negative symptomatic PCR"){
            for(p in 1:length(ppl$person[p])){
              if(ppl$new.test[p]==1){
                if(ppl$test.returned[p]==1){
                  if(ppl$test.result[p]==0){
                    if(ppl$symptom[p]==1 | ppl$ILI[p]==1){
                      #if(ppl$done.testing[p]!=1){
                      ppl$needs.retest[p]<-1
                      ppl$retest.day[p]<-t+retest.delay
                      #}
                    }
                  }
                }
              }
            }
          } ## end retest negative symptomatic stuff
          
          ## and people who have test results back (from an initial test only! not from a retest.), tested negative get set up to be retested after delay
          if(test.type=="retest all negatives"){
            for(p in 1:length(ppl$person)){
              if(ppl$new.test[p]==1){
                if(ppl$test.returned[p]==1){
                  if(ppl$test.result[p]==0){
                    #if(ppl$symptom[p]==1){ ## only need to remove this line relative to other reflex strategy
                    #if(ppl$done.testing[p]!=1){
                    ppl$needs.retest[p]<-1
                    ppl$retest.day[p]<-t+retest.delay
                    #}
                    #}
                  }
                }
              }
            }
          } ## end retest all negatives
          
          
          
          
          ## and people who have results back, and will not be retested are no longer waiting.
          ppl$waiting[ppl$test.returned==1 & ppl$needs.retest==0]<-0
          ##and people who are still waiting are recorded as having waited for another day.
          ppl$time.waiting[ppl$waiting==1]<-ppl$time.waiting[ppl$waiting==1]+1
          
          
          ## next, we remove people from quarantine who have lived it out 
          for(p in 1:length(ppl$person)){ ## instead of a rate, we use this individual-based approrach so that we can account for testing delay
            if(ppl$quarantine[p]==1){
              if(ppl$test.day[p]>0){
                if(ppl$test.day[p]==t-(1/omega)){
                  ppl$quarantine[p]<-0
                  ppl$test.day[p]<-0
                }
              }## end catch errors
              if(ppl$retest.day[p]>0){## catch errors due to NAs
                if(ppl$retest.day[p]==t-(1/omega)){
                  ppl$quarantine[p]<-0
                  ppl$retest.day[p]<-0
                }
              }## end catch errors
            }
          }
          
          
          ######### the virus spreads
          for(p in 1:length(ppl$person)){
            if(ppl$quarantine[p]==0){ ## quarantine is not leaky
              if(ppl$state[p]=="S"){ ## if a person is susceptible
                if(rbinom(1,1,p.expose[p])==1){
                  ppl$state[p]<-"E" ## expose with above prob
                }
              }
            } ## end no quarantine. everything else happens whether or not a person is quarantined
            if(ppl$state[p]=="E"){ ## if a person is exposed
              if(rbinom(1,1,p.infect)){
                ppl$state[p]<-"I" ## infect with above prob
                ppl$ever.infected[p]<-1
                ppl$recovery.day[p]<-round(t+runif(1,infectious.period.min,infectious.period.max))
                ppl$symptom.onset.day[p]<-ifelse(rbinom(1,1,p.symptomatic),round(t+incubation[p]-(1/sigma)),0) ## and become symptomatic with this prob
                summit$new.infs[t]<-1+summit$new.infs[t]
              }
            }
            if(ppl$state[p]=="I"){ ## if a person is infected
              ## instead of rate-based recovery.....
              # if(rbinom(1,1,p.recover)){
              #   ppl$state[p]<-"R" ## recover with above prob
              # }
              ## we'll use an individual-based recovery
              if(ppl$recovery.day[p]==t){
                ppl$state[p]<-"R"
              }
            }
            
            
          } ## end loop of people
          
          
          
        } ## end loop of days
        ## concatenate the summit dfs for each simulation
        summit$site<-site
        summit$test.strategy<-test.type
        sum4plots<-rbind(summit,sum4plots)
        
        ## uncomment following to look at results
        # par(mfrow=c(2,1))
        # plot(x=0,y=0,type="n",ylim=c(0,.5*pop.size),xlim=c(0,simdur),
        #      main=paste0(test.type,"  ",site))
        # #plot(summit$S,col="black",ylim=c(0,pop.size))
        # lines(summit$E,col="pink")
        # lines(summit$I,col="red",lwd=2)
        # #lines(summit$R,col="green")
        # lines(summit$tests,col="purple")
        # lines(summit$quarantines,col="yellow")
        # #lines(summit$false_quarantines,col="brown")
        
        # lines(summit$tS,col="black",lty=2)
        # lines(summit$tE,col="pink",lty=2)
        # lines(summit$tI,col="red",lty=2)
        # lines(summit$tR,col="green",lty=2)
        
        # # look at performance of testing
        # plot(summit$missed,type="l",col="red",main=paste0(test.type,"  ",site))
        # lines(summit$caught,col="green")
        
        ## look at distribution of number of tests
        # plot(density(ppl$retested))
        # plot(density(ppl$tested))
        # par(mfrow=c(1,1))
        
        ########## reporting
        
        megs$test.strategy[i]<-test.type
        megs$location[i]<-site
        megs$gamma[i]<-mean(gamma)
        megs$R0[i]<-mean(R0)
        megs$prop.testing[i]<-daily.testing.proportion
        megs$infections[i]<-sum(ppl$ever.infected)+n.seed.events ## total infections
        megs$missed.infections[i]<-max(summit$missed) ## total missed infections
      ########## HARDCODE FLAG: hospitalization/ICU/fatality rates
        megs$hospitalizations[i]<-megs$infections[i]*ifelse(site=="NH",0.25,0.039) ## hospitalization rate
        megs$ICU[i]<-megs$hospitalizations[i]*ifelse(site=="NH",0.353,0.238) ## ICU entrance rate (of hospitalized cases)
        megs$deaths[i]<-megs$infections[i]*ifelse(site=="NH",0.054,0.0001) ## infection fatality rate
      ########## END hardcode flag
        megs$days.quarantine[i]<-sum(summit$quarantines)
        megs$days.false.quarantine[i]<-sum(summit$false_quarantines)
        megs$n.reflexes[i]<-sum(ppl$retested) ## total number of retests (people can be subject to multiple retests)
        megs$n.tests[i]<-sum(summit$tests) ## total # tests, including retests (people can be tested multiple times, and multiple retests)
        megs$positive.tests[i]<-sum(summit$positive.tests)
        megs$negative.tests[i]<-sum(summit$negative.tests)
        megs$time.waiting[i]<-sum(ppl$time.waiting)
        megs$caught.infs[i]<-max(summit$caught)
        megs$infectious_period[i]<-mean(1/gamma)
        if(test.type=="retest negative symptomatic PCR"){
          megs$n.alt.tests[i]<-sum(ppl$retested) ## if retest strategy uses a different type test than normal, then track that and keep it here.
          megs$n.reflexes[i]<-0 ## because no retests of same type were given
        }
      } ## end full loop of test strategies
      
      megs$test.strategy<-as.factor(megs$test.strategy)
      meg.sum<-data.frame(test.strategy=unique(megs$test.strategy),location=site,infections=NA,
                          missed.infections=NA,hospitalizations=NA,ICU=NA,deaths=NA,days.quarantine=NA,
                          days.false.quarantine=NA,time.waiting=NA,n.reflexes=NA,n.tests=NA,positive.tests=NA,
                          negative.tests=NA,n.alt.tests=NA)
      se <- function(x) sd(x)/sqrt(length(x)) ## function to compute the standard error
      for(i in 1:length(meg.sum$test.strategy)){ ## if you want a summary table:
        meg.sum$infections[i]<-mean(megs$infections[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$missed.infections[i]<-mean(megs$missed.infections[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$hospitalizations[i]<-mean(megs$hospitalizations[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$ICU[i]<-mean(megs$ICU[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$deaths[i]<-mean(megs$deaths[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$days.quarantine[i]<-mean(megs$days.quarantine[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$days.false.quarantine[i]<-mean(megs$days.false.quarantine[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$time.waiting[i]<-mean(megs$time.waiting[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$n.reflexes[i]<-mean(megs$n.reflexes[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$n.tests[i]<-mean(megs$n.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$positive.tests[i]<-mean(megs$positive.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$negative.tests[i]<-mean(megs$negative.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$n.alt.tests[i]<-mean(megs$n.alt.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        
        #s.e.
        meg.sum$se.infections[i]<-se(megs$infections[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.missed.infections[i]<-se(megs$missed.infections[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.hospitalizations[i]<-se(megs$hospitalizations[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.ICU[i]<-se(megs$ICU[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.deaths[i]<-se(megs$deaths[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.days.quarantine[i]<-se(megs$days.quarantine[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.days.false.quarantine[i]<-se(megs$days.false.quarantine[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.time.waiting[i]<-se(megs$time.waiting[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.n.reflexes[i]<-se(megs$n.reflexes[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.n.tests[i]<-se(megs$n.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.positive.tests[i]<-se(megs$positive.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.negative.tests[i]<-se(megs$negative.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        meg.sum$se.n.alt.tests[i]<-se(megs$n.alt.tests[megs$test.strategy==meg.sum$test.strategy[i]])
        
      }
      ### save df of raw data
      prev.wd<-paste0(getwd(),"/")
      setwd(path)
      fwrite(megs,file=paste0(path,"sims_",Sys.Date(),"_",site,"_",daily.testing.proportion,".csv"))
      
      ### save dfs for both sites:
      fwrite(meg.sum,file=paste0(path,"sims_sum_",Sys.Date(),"_",site,"_",daily.testing.proportion,".csv"))
      
      ### save df of all dynamics for making epi curves
      fwrite(sum4plots,file=paste0(path,"sum4plots_",Sys.Date(),"_",site,"_",daily.testing.proportion,".csv"))
      setwd(prev.wd)
    }## end FULL FULL LOOP OF SITES
    
  }## end FULL FULL FULL LOOP OF TESTING PROPORTIONS
  
}## end function


