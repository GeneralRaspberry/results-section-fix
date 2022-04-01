
library("purrr")

test1<-function(x){
  f<-max(pluck(x,"time"))
  #stti<-sample(1:g)
  #de<-c(f,stti)
}

#maxtimes<-map(s,test1)


datasurexperiment<-data#%>%filter(theta==140)
s<-split(datasurexperiment,datasurexperiment$sim)

funcsur<-function(i=NULL){

n<-seq(15,120,length.out = 7)
fre<-seq(15,120,length.out = 7)
datalist<-list()
indicator<-1
for (f in n){
  for (l in fre){
test4<-function(x,y){
  for(g in y){
    r5<-pluck(x)
    d<-r5[sample(nrow(r5),size=f,replace=FALSE),]
    print(d)
    print(match(g,y))
print(paste0("fre=",f,"n=",l))
    
    if(g>min(d$time)){
      m<-sum(d <= g)
      q<-mean(r5$time <= g)
      t<-min(d$time)
      d<-match(g,y)
      mylist<-list(q,m,t=g,d=d,theta=head(r5$theta,1),beta=head(r5$beta,1),sim=head(r5$sim,1),truesim=indicator,frequency=l,samplesize=f)
      return(mylist)
    }
  }
}




stti<-sample(1:l,length(s),replace=TRUE)
samp.time <- lapply(stti, function(x) seq(from = x, to = max(datasurexperiment$time), by = l))


dftest4<-map2(s,samp.time,test4)
print("....................simulation set completed.............................")
dftestlistdocall<-data.frame(do.call(rbind,dftest4))
datalist[[indicator]]<-dftestlistdocall
      indicator<-indicator+1
  }
}

dfcompleteheatmap<-do.call(rbind,datalist)
}
cl <- makeCluster(mc <- getOption("cl.cores", 20))
clusterCall(cl,function() library("dplyr"))
clusterCall(cl,function() library("purrr"))
clusterExport(cl=cl, varlist=c("s","datasurexperiment"),envir = environment())
par_r1<-parLapply(1,fun=funcsur,cl=cl)
stopCluster(cl)
dfcompleteheatmap <- do.call("rbind", par_r1)


dfcompleteheatmap$q<-as.numeric(dfcompleteheatmap$V1)
sum.q<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(q),list(q_mean = mean))

sum.q$anq<-((rdatamean$r_mean)*(as.numeric(sum.q$frequency)/as.numeric(sum.q$samplesize)))
sum.q$absdif<-abs(sum.q$anq-sum.q$q_mean)
sum.q$reldif<-sum.q$absdif/sum.q$q_mean
dfcompleteheatmap$t1<-as.numeric(dfcompleteheatmap$t)
dfcompleteheatmap$d1<-as.numeric(dfcompleteheatmap$d)
time314<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(t1),list(t_mean = mean))
steps314<-dfcompleteheatmap%>%group_by(frequency,samplesize)%>%summarise_at(vars(d1),list(steps_mean = mean))

####added correctional factor
attach(sum.q)
sum.q$frequency<-as.numeric(sum.q$frequency)
sum.q$samplesize<-as.numeric(sum.q$samplesize)
maxlist<-apply(sum.q[c(1:2)],1,function(x) max(x))
minlist<-apply(sum.q[c(1:2)],1,function(x) min(x))

rtn<-rdatamean$r_mean*sum.q$samplesize/sum.q$frequency

correctionfactor5<-(sum.q$frequency/sum.q$samplesize)*log((rtn),rdatamean$r_mean)

sum.q$predictioncorrection5<-rtn*correctionfactor5

sum.q$abscorrect5<-abs(sum.q$q_mean-sum.q$predictioncorrection5)

sum.q$relcorrect5<-sum.q$abscorrect5/sum.q$q_mean
###############working correction
sum.q$testout2<-log(sum.q$anq,rdatamean$r_mean)
sum.q$testout2r<-log(rdatamean$r_mean,sum.q$anq)
sum.q$testout2.1<-sum.q$anq+((sum.q$anq*sum.q$testout2^2)/2*sum.q$testout2r)
sum.q$correction8abs<-abs(sum.q$testout2.1-sum.q$q_mean)
sum(sum.q$correction8abs)
sum.q$relcorrect8<-sum.q$correction8abs/sum.q$q_mean
################newtest
test1<-sum.q$anq+rdatamean$r_mean*(sum.q$frequency/sum.q$samplesize)
absnewtest<-abs(sum.q$q_mean-test1)

#########no good new test
cor<-1.618*((1+rdatamean$r_mean)^(sum.q$frequency/sum.q$samplesize))
cor<-1.6180339887*exp(sum.q$anq)
sum.q$abscor<-abs(sum.q$q_mean-(sum.q$anq*cor))
sum.q$relcor<-sum.q$abscor/sum.q$q_mean

sum.q$newcor<-(rdatamean$r_mean*sum.q$frequency/sum.q$samplesize)*(1.61803398)
sum.q$abscornew<-abs(sum.q$q_mean-sum.q$newcor)
sum.q$relcornew<-sum.q$abscornew/sum.q$q_mean
sum(sum.q$absdif)
sum(sum.q$abscor)
sum(sum.q$relcor)
sum(sum.q$reldif)
sum(sum.q$relcornew)
##################################################
sum.q$store<-vector(length=length(sum.q$frequency))
  for(i in 1:length(sum.q$frequency)){
  if(sum.q$frequency[i]<=sum.q$samplesize[i]){
    sum.q$store[i]<-(rdatamean$r_mean*sum.q$frequency[i])/(sum.q$samplesize[i])*1.61803395
  }else{
    sum.q$store[i]<-(rdatamean$r_mean*sum.q$frequency[i])/(sum.q$samplesize[i])*1.61803395
  }
  }

sum.q$grabs<-abs(sum.q$q_mean-sum.q$store)
sum.q$grrel<-sum.q$grabs/sum.q$q_mean
sum(sum.q$absdif)
sum(sum.q$reldif)
sum(sum.q$grabs)
sum(sum.q$grrel)
#################################
sum.q$newcorfinal<-(1.61803395)^(sum.q$frequency/sum.q$samplesize)
sum.q$grabs1<-abs(sum.q$q_mean-(sum.q$anq*exp(sum.q$q_mean*1.61803395)))
sum.q$grrel1<-sum.q$grabs1/sum.q$q_mean
sum(sum.q$grabs1)
sum(sum.q$grrel1)
################################lunchtestthatworks quite nicely
sum.q$lunchtest<-sum.q$anq*140*rdatamean$r_mean/(1+(sum.q$frequency/sum.q$samplesize))


sum.q$lunchtest2<-sum.q$anq*(20*rdatamean$r_mean/2+(rdatamean$r_mean^2)/20)


sum.q$lunchtest1<-abs(sum.q$q_mean-(sum.q$lunchtest2))
sum.q$grrel1<-sum.q$lunchtest1/sum.q$q_mean
###############################
sum.q$correction8<-sum.q$anq+(log(rdatamean$r_mean,sum.q$anq)*sum.q$frequency/sum.q$samplesize)
sum.q$correction8abs<-abs(sum.q$correction8-sum.q$q_mean)
sum.q$correction6<-sum.q$anq+rdatamean$r_mean+rdatamean$r_mean^(sum.q$frequency+sum.q$samplesize/sum.q$samplesize)
sum.q$correction6abs<-sum.q$frequency*sum.q$absdif
sum(sum.q$correction6abs)
sum.q$correction6rel<-sum.q$correction6abs/sum.q$q_mean
save(dfcompleteheatmap,file="dfheatmap4testbeta50theta20randomnew.rda")
save(sum.q,file="sum.qbeta50theta20newdatacompleterandomnew.rda")
save(rdatamean,file="rdatasum.qbeta50theta20randomnew.Rda")

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=absdif)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Absolute Difference",
                      limits=c(0,.25))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(absdif, 3)))

ggsave(file="absdiffixedbeta150theta20.pdf")



ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=anq)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="rule of thumb prediction",
                      limits=c(0,.25))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(anq, 3)))

ggsave(file="anqvaluesbeta150theta20.pdf")


ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=q_mean)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Simulated detection",
                      limits=c(0,.25))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(q_mean, 3)))

ggsave(file="simulatedvaluesbeta150theta20.pdf")



ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=reldif)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(reldif, 3)))

ggsave(file="reldifbeta150theta20.pdf")

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=grrel)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(grrel, 3)))

ggsave(file="reldifgoldenratiobeta150theta20.pdf")




ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=abscornew)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Simulated detection",
                      limits=c(0,.025))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(abscornew, 3)))

ggsave(file="heatmapqmeanbeta150theta140.pdf")


ggsave(file="heatmapabsbeta50theta20.pdf")




ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=relcornew)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(relcornew, 3)))

ggsave(file="heatmaprelbeta50theta20withcorrect.pdf")

ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=grrel1)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="Relative Difference",
                      limits=c(0,1))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1) + 
  geom_text(aes(label = round(grrel1, 3)))


ggplot(sum.q,aes(x=as.numeric(frequency),y=as.numeric(samplesize),fill=lunchtest1)) +
  geom_tile() + 
  scale_fill_gradient(low="white",
                      high="darkred",
                      name="lunchtest",
                      limits=c(0,.25))+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),aspect.ratio = 1)+ 
  geom_text(aes(label = round(lunchtest1, 3)))