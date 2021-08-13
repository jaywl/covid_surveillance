###### Supplemental code to reproduce figures in 
###### "Comparison of antigen- and RT-PCR-based testing strategies for detection
###### of SARS-CoV-2 in two high-exposure settings"
###### Authors: Jay Love, Megan T. Wimmer, Damon J.A. Toth, Arthi Chandran, 
###### Dilip Makhija, Charles K. Cooper, Matthew H. Samore, Lindsay T. Keegan"

temp_output_dir<-"/Users/temp_output/"## change to your directory of choice for output of function (if run)
temp_data_dir<-"/Users/plos_data/" ## change to directory of data you downloaded from repo
temp_plot_dir<-"/Users/plots/" ## change to a directory in which to save plots

#### to generate the data yourself (long runtime on normal computers):
setwd(temp_output_dir) 
cvtest_model(nstrats=4, ## number of testing strategies (4, for the manuscript)
           nsims=1000, ## number of simulations to run
           nsites=2, ## number of sites (e.g., nursing home and residence hall)
           ntesting_proportions=5, ## number of surveillance daily testing proportions to simulate (separately)
           testprops=c(0.0,0.01,0.02,0.05,0.10), ## proportions to simulate
           quarantine_days=14, ## days of quarantine following a positive test
           pcr.sensitivity=1.0, ## positive percent agreement wrt viral culture
           pcr.specificity=0.955, ## negative percent agreement wrt viral culture
           antigen.sensitivity=0.964, ## positive percent agreement wrt viral culture
           antigen.specificity=0.987, ## negative percent agreement wrt viral culture
           theta.pcr=(1/3), ## 1/days to return pcr test
           theta.antigen=1, ## =1 is same day return of antigen tests
           prob_symptomatic=0.6, ## proportion of infections that are symptomatic
           ILI_incidence=0.02, ## standing level of non-covid ILI
           waiting.reduction=0.50, ## reduction in mixing due to quarantine
           include.rs=FALSE, ## set to true if recovered class positive tests should be "true" positives, not "false" positives (for exploring e.g., Fig S2) 
           r.time=54, ## if include.rs==TRUE, post-recovery days before which a positive test is considered a true positive.
           path=paste0(getwd(),"/"))## path to save output files

#### To use the data we generated:
setwd(temp_data_dir) ## change to directory you downloaded
fls<-list.files(pattern="sum4plots.")
# length(fls) ##checks
# fls
s1<-read.csv(fls[1])
s2<-read.csv(fls[2])
s3<-read.csv(fls[3])
s4<-read.csv(fls[4])
s5<-read.csv(fls[5])
s6<-read.csv(fls[6])
s7<-read.csv(fls[7])
s8<-read.csv(fls[8])
s9<-read.csv(fls[9])
s10<-read.csv(fls[10])
s11<-rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)

## making figure 2 (and supplemental epidemic curve plots)
tps<-c(0.01,0.02,0.05,0.10) ## proportions testing (same as testprops in function parameters)
stes<-c("NH","dorm")

## save individual files
for(q in 1:length(stes)){
  upperlim<-6
  if(q==2){ ## adjust y limit for dorms
    upperlim<-150
  }
  for(i in 1:length(tps)){
    setwd(temp_plot_dir)#### set to a directory to save plots
    ggplot()+geom_smooth(data=s11[s11$test.strategy=="none"&s11$site==stes[q],],aes(t,I),color="brown",method="gam")+
      geom_smooth(data=s11[with(s11,test.strategy!="none"&site==stes[q]&prop.testing==tps[i]),],aes(t,I,colour=test.strategy),alpha=.1,method="gam")+
      theme_minimal(base_size=12,base_family="Helvetica")+coord_cartesian(xlim=c(0,120),ylim=c(0,upperlim))+
      ggtitle(paste0(tps[i]))+theme(legend.position="none")
    ggsave(filename=paste0("epicurve_",stes[q],"_",tps[i],".eps"),height=4,width=4,dpi=300,units="in")
    ggsave(filename=paste0("epicurve_",stes[q],"_",tps[i],".png"),height=4,width=4,dpi=300,units="in")
    
  }
}

## Define multiplot function (very useful, adapted from www.cookbook-R.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## figure 2 itself:
{q<-1
  upperlim<-6
  if(q==2){ ## for only dorms (y limit problems)
    upperlim<-150
  }
  i<-2
  p1<-ggplot()+geom_smooth(data=s11[s11$test.strategy=="none"&s11$site==stes[q],],aes(t,I),color="brown",method="gam")+
    geom_smooth(data=s11[with(s11,test.strategy!="none"&site==stes[q]&prop.testing==tps[i]),],aes(t,I,colour=test.strategy),alpha=.1,method="gam")+
    theme_minimal(base_size=12,base_family="Helvetica")+coord_cartesian(xlim=c(0,120),ylim=c(0,upperlim))+
    ggtitle(paste0(tps[i]))+theme(legend.position="none")
  
  i<-4
  p2<-ggplot()+geom_smooth(data=s11[s11$test.strategy=="none"&s11$site==stes[q],],aes(t,I),color="brown",method="gam")+
    geom_smooth(data=s11[with(s11,test.strategy!="none"&site==stes[q]&prop.testing==tps[i]),],aes(t,I,colour=test.strategy),alpha=.1,method="gam")+
    theme_minimal(base_size=12,base_family="Helvetica")+coord_cartesian(xlim=c(0,120),ylim=c(0,upperlim))+
    ggtitle(paste0(tps[i]))+theme(legend.position="none")
  
  q<-2
  upperlim<-6
  if(q==2){ ## for only dorms (y limit problems)
    upperlim<-150
  }
  i<-2
  p3<-ggplot()+geom_smooth(data=s11[s11$test.strategy=="none"&s11$site==stes[q],],aes(t,I),color="brown",method="gam")+
    geom_smooth(data=s11[with(s11,test.strategy!="none"&site==stes[q]&prop.testing==tps[i]),],aes(t,I,colour=test.strategy),alpha=.1,method="gam")+
    theme_minimal(base_size=12,base_family="Helvetica")+coord_cartesian(xlim=c(0,120),ylim=c(0,upperlim))+
    ggtitle(paste0(tps[i]))+theme(legend.position="none")
  
  i<-4
  p4<-ggplot()+geom_smooth(data=s11[s11$test.strategy=="none"&s11$site==stes[q],],aes(t,I),color="brown",method="gam")+
    geom_smooth(data=s11[with(s11,test.strategy!="none"&site==stes[q]&prop.testing==tps[i]),],aes(t,I,colour=test.strategy),alpha=.1,method="gam")+
    theme_minimal(base_size=12,base_family="Helvetica")+coord_cartesian(xlim=c(0,120),ylim=c(0,upperlim))+
    ggtitle(paste0(tps[i]))+theme(legend.position="none")
  multiplot(p1,p2,p3,p4,cols=2)
}

## making figure 3
setwd(temp_data_dir)## change to directory you downloaded
fls1<-list.files(pattern="sims_2021")
length(fls1)
fls1
m1<-read.csv(fls1[1])
m2<-read.csv(fls1[2])
m3<-read.csv(fls1[3])
m4<-read.csv(fls1[4])
m5<-read.csv(fls1[5])
m6<-read.csv(fls1[6])
m7<-read.csv(fls1[7])
m8<-read.csv(fls1[8])
m9<-read.csv(fls1[9])
m10<-read.csv(fls1[10])
m11<-rbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)

mss.d<-m11[m11$location=="dorm",]
mss.nh<-m11[m11$location=="NH",]
## to get infections averted
baseline.infections.d<-mean(mss.d$infections[with(mss.d,test.strategy=="none")])
baseline.infections.nh<-mean(mss.nh$infections[with(mss.nh,test.strategy=="none")])
## to get hospitalizations averted
baseline.hospitalizations.d<-mean(mss.d$hospitalizations[with(mss.d,test.strategy=="none")])
baseline.hospitalizations.nh<-mean(mss.nh$hospitalizations[with(mss.nh,test.strategy=="none")])
## to get deaths averted
baseline.deaths.d<-mean(mss.d$deaths[with(mss.d,test.strategy=="none")])

##figure 3
require(gplots) ## helpful plotting package. Thanks to Warnes, Bolker, Bonebakker, et al. https://github.com/talgalili/gplots
setwd(temp_plot_dir)#### set to a directory to save plots
{par(mfrow=c(1,2))
  # png(filename=paste0("percent_inf_averted",".png"),height=4,width=8,res=300,units="in")
  setEPS()
  # cairo.ps(paste0("percent_inf_averted",".eps"),height=4,width=8,fallback.res=300)
  postscript(paste0("percent_inf_averted",".eps"),height=7,width=7)
  plotmeans(100*(baseline.infections.nh-infections)/baseline.infections.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.01,],
            connect=FALSE,barcol="gray",pch=16,cex=0.8,ci.label=TRUE,digits=1,ylab="% Infections Averted",xlab="Testing Strategy",
            n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
to PCR"),ylim=c(0,70),main="Nursing Home",gap=0)
  plotmeans(100*(baseline.infections.nh-infections)/baseline.infections.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.02,],
            col="blue",add=TRUE,connect=FALSE,barcol="light blue",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  plotmeans(100*(baseline.infections.nh-infections)/baseline.infections.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.05,],
            col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  plotmeans(100*(baseline.infections.nh-infections)/baseline.infections.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.10,],
            col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  
  plotmeans(100*(baseline.infections.d-infections)/baseline.infections.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.01,],
            connect=FALSE,barcol="gray",pch=16,cex=0.8,ci.label=TRUE,digits=1,ylab="% Infections Averted",xlab="Testing Strategy",
            n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
to PCR"),ylim=c(0,70),main="Residence Hall",las=1,gap=0)
  plotmeans(100*(baseline.infections.d-infections)/baseline.infections.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.02,],
            col="blue",add=TRUE,connect=FALSE,barcol="light blue",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  plotmeans(100*(baseline.infections.d-infections)/baseline.infections.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.05,],
            col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  plotmeans(100*(baseline.infections.d-infections)/baseline.infections.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.10,],
            col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=0.8,ci.label=TRUE,digits=1,n.label=FALSE,xaxt="n",gap=0)
  
  par(mfrow=c(1,1))
  dev.off()
}




##Figure S1

{par(mfrow=c(1,2))
  ## hosps and deaths averted are proportions of infections, so plots look the same. But if curious:
#   plotmeans(100*(baseline.hospitalizations.d-hospitalizations)/baseline.hospitalizations.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.01,],
#             connect=FALSE,barcol="gray",pch=16,cex=2.5,ylab="% Hospitalizations Averted",xlab="Testing Strategy",
#             n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
# to PCR"),ylim=c(0,100),main="Residence Hall")
#   plotmeans(100*(baseline.hospitalizations.d-hospitalizations)/baseline.hospitalizations.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.02,],
#             col="blue",add=TRUE,connect=FALSE,barcol="lightblue",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#   plotmeans(100*(baseline.hospitalizations.d-hospitalizations)/baseline.hospitalizations.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.05,],
#             col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#   plotmeans(100*(baseline.hospitalizations.d-hospitalizations)/baseline.hospitalizations.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.10,],
#             col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#
#   plotmeans(100*(baseline.hospitalizations.nh-hospitalizations)/baseline.hospitalizations.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.01,],
#             connect=FALSE,barcol="gray",pch=16,cex=2.5,ylab="% Hospitalizations Averted",xlab="Testing Strategy",
#             n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
# to PCR"),ylim=c(0,100),main="Nursing Home")
#   plotmeans(100*(baseline.hospitalizations.nh-hospitalizations)/baseline.hospitalizations.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.02,],
#             col="blue",add=TRUE,connect=FALSE,barcol="lightblue",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#   plotmeans(100*(baseline.hospitalizations.nh-hospitalizations)/baseline.hospitalizations.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.05,],
#             col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#   plotmeans(100*(baseline.hospitalizations.nh-hospitalizations)/baseline.hospitalizations.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.10,],
#             col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
#
  ##deaths
  plotmeans(100*(baseline.deaths.d-deaths)/baseline.deaths.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.01,],
            connect=FALSE,barcol="gray",pch=16,cex=2.5,ylab="% Deaths Averted",xlab="Testing Strategy",main="Residence Hall",
            n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
to PCR"),ylim=c(0,100))
  plotmeans(100*(baseline.deaths.d-deaths)/baseline.deaths.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.02,],
            col="blue",add=TRUE,connect=FALSE,barcol="lightblue",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
  plotmeans(100*(baseline.deaths.d-deaths)/baseline.deaths.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.05,],
            col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
  plotmeans(100*(baseline.deaths.d-deaths)/baseline.deaths.d~test.strategy,mss.d[mss.d$test.strategy!="none" & mss.d$prop.testing==0.10,],
            col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=2.5,n.label=FALSE,xaxt="n")

  plotmeans(100*(baseline.deaths.nh-deaths)/baseline.deaths.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.01,],
            connect=FALSE,barcol="gray",pch=16,cex=2.5,ylab="% Deaths Averted",xlab="Testing Strategy",main="Nursing Home",
            n.label=FALSE,legends=c("antigen","PCR","Reflex to Antigen","Reflex Symp.
to PCR"),ylim=c(0,100))
  plotmeans(100*(baseline.deaths.nh-deaths)/baseline.deaths.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.02,],
            col="blue",add=TRUE,connect=FALSE,barcol="lightblue",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
  plotmeans(100*(baseline.deaths.nh-deaths)/baseline.deaths.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.05,],
            col="green",add=TRUE,connect=FALSE,barcol="light green",pch=16,cex=2.5,n.label=FALSE,xaxt="n")
  plotmeans(100*(baseline.deaths.nh-deaths)/baseline.deaths.nh~test.strategy,mss.nh[mss.nh$test.strategy!="none" & mss.nh$prop.testing==0.10,],
            col="red",add=TRUE,connect=FALSE,barcol="pink",pch=16,cex=2.5,n.label=FALSE,xaxt="n")


  par(mfrow=c(1,1))
}


### per test reduction in infections ('C' panels in figures S2-S8)

## dorm
mss.d$percentreduction[mss.d$test.strategy!="none"]<-(100*(baseline.infections.d-mss.d$infections[mss.d$test.strategy!="none"])/baseline.infections.d)
mss.d$pertestreduction[mss.d$test.strategy!="none"]<-(100*(baseline.infections.d-mss.d$infections[mss.d$test.strategy!="none"])/baseline.infections.d)/(mss.d$n.tests[mss.d$test.strategy!="none"])
ptp1<-ggplot()+
  geom_boxplot(data=mss.d[mss.d$test.strategy!="none",],
               aes(x=as.factor(prop.testing),y=pertestreduction,fill=test.strategy),
               position=position_dodge(),outlier.shape=1,outlier.alpha=.1)+
  coord_cartesian(ylim=c(-0.003,0.004))+
  theme_minimal(base_size=12,base_family="Helvetica")+
  xlab("Daily Testing Proportion")+
  ylab("Per-test Percent Infections Averted")+
  ggtitle("Residence Hall")+
  theme(legend.position="none")
## Nursing home
mss.nh$percentreduction[mss.nh$test.strategy!="none"]<-(100*(baseline.infections.nh-mss.nh$infections[mss.nh$test.strategy!="none"])/baseline.infections.nh)
mss.nh$pertestreduction[mss.nh$test.strategy!="none"]<-(100*(baseline.infections.nh-mss.nh$infections[mss.nh$test.strategy!="none"])/baseline.infections.nh)/(mss.nh$n.tests[mss.nh$test.strategy!="none"])
ptp2<-ggplot()+
  geom_boxplot(data=mss.nh[mss.nh$test.strategy!="none",],
               aes(x=as.factor(prop.testing),y=pertestreduction,fill=test.strategy),
               position=position_dodge(),outlier.shape=1,outlier.alpha=.1)+
  coord_cartesian(ylim=c(-0.45,0.6))+
  theme_minimal(base_size=12,base_family="Helvetica")+
  xlab("Daily Testing Proportion")+
  ylab("Per-test Percent Infections Averted")+
  ggtitle("Nursing Home")+
  theme(legend.position="none")
setwd(temp_plot_dir) ##change to directory to save plots
# png(filename=paste0("perTest_inf_averted",".png"),height=4,width=8,res=300,units="in")
setEPS()
postscript(paste0("perTest_inf_averted",".eps"),height=4,width=8)
multiplot(ptp2,ptp1,cols=2)
dev.off()


