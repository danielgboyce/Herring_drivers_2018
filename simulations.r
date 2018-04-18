library(reshape2)
library(gplots)
library(forecast)
library(cluster)
library(vegan)
library(ggplot2)
library(hybridHclust)
library(raster)
library(fields)
library(gridExtra)
library(colorRamps)
library(mapdata)
library(scales)
library(MASS)
library(mgcv)
library(maps)
library(plyr)
library(plotrix)
library(lubridate)


datadir<-'/scratch/dboyce/chl_phenology/data'
codedir<-'/scratch/dboyce/chl_phenology/code'
figsdir<-'/scratch/dboyce/chl_phenology/figures'

setwd(codedir)
source('helper.functions.r')

data<-data.frame(day=seq(1,365,1))

#e1<-rnorm(365,mean=0,sd=.5)
#e2<-rnorm(365,mean=0,sd=1)             
e1<-rnorm(365,mean=5,sd=.18)
e2<-rnorm(365,mean=5,sd=.54)             
cos1<-cos((2*3.14*1*data$day)/365.25)*-1
#sin1<-sin((2*3.14*1*data$day)/365.25)
cos2<-cos((2*3.14*2*data$day)/365.25)*-1
#sin2<-sin((2*3.14*2*data$day)/365.25)
#PHENOLOGY WITH NOISE ADDED; 3 PATTERNS: 1) COS: SINGLE BLOOM IN MIDDLE; 2) SIN: SINGLE BLOOM AT START; 3) COMBINATION WITH 2 BLOOMS
data$p1<-rescale(cos1,newrange=c(-2,2))+e1+3
data$p2<-rescale(cos2,newrange=c(-2,2))+e1+3
data$p3<-rescale(cos1,newrange=c(-2,2))+e1+3
data$p4<-rescale(cos2,newrange=c(-2,2))+e2+3
#data$p5<-rescale(cos1,newrange=c(-2,2))+e2
#data$p6<-rescale(sin1,newrange=c(-2,2))+e2
#PHENOLOGY WITH NO NOISE
data$p1mn<-rescale(cos1,newrange=c(-2,2))+3
data$p2mn<-rescale(cos2,newrange=c(-2,2))+3
data$p3mn<-rescale(cos1,newrange=c(-2,2))+3
data$p4mn<-rescale(cos2,newrange=c(-2,2))+3

data$date<-strptime(gsub(' ','',paste(2000,'-',data$day)),"%Y-%j")
data$month<-month(data$date)

plot(data,pch=16)

setwd(datadir)
q<-read.csv('phenology.averages.matrix.tbin10.csv',header=TRUE)
#q<-q[,sample(ncol(q),400,replace=FALSE)]


##SIMULATION TO DETERMINE THE MINIMUM SAMPLING NEEDED TO RECOVER PHENOLOGY TRENDS
##GENERATE SAMPLES RANGING FROM 10 (~1 PER MONTH) TO 100, WITH 1000 REPLICATES EACH
#UNDER ABOVE CRITERIA, SHOULD TAKE ~3 HOURS

ssizes<-seq(4,100,1)
#ssizes<-c(8)
l<-list()
system.time(for(i in 1:length(ssizes)){
print(ssizes[i])

#REPEAT FUNCTION SPECIFIED NUMBER OF TIMES
    dt<-replicate(500, {
        #GENERATE RANDOM SAMPLE OF DATA
        d<-data[sample(nrow(data),ssizes[i],replace=FALSE),]
        nmonths<-length(unique(d$month))
        
###########################
            mfun<-function(d2){
            nm<-gsub('p','f',names(d2)[2])
            names(d2)[1:2]<-c('x','y')            
                                        #FIT MODEL - THE WORKHORSE
            dum<-data.frame(x=sort(unique(d2$x)),
                            sump=tapply(d2$y,d2$x,mean))
            dum$ct1<-dum$x*dum$sump
            ct<-sum(dum$ct1)/sum(dum$sump)
            dout<-data.frame(nmonth=length(unique(d2$month)),
                  n=dim(d2)[1],
                  ct=ct,
                  nm=nm)
            return(dout)
            }
            datout<-rbind(mfun(subset(d,select=c('day','p1','month'))),
                       mfun(subset(d,select=c('day','p3','month'))))
            
            datout
           },simplify=FALSE)

            l[[i]]<-data.frame(do.call('rbind.fill',dt))
            rm(dt)
    })                            
out<-data.frame(do.call('rbind',l))#78 HOURS

plot(out$nmonth,out$ct)
plot(out$n,out$ct)






########################################################################
#SAME BUT ON FALL PEAK SIMULATION
data2<-subset(data,day>184)
ssizes<-seq(4,100,1)
l<-list()
system.time(for(i in 1:length(ssizes)){
print(ssizes[i])

#REPEAT FUNCTION SPECIFIED NUMBER OF TIMES
    dt<-replicate(500, {
        #GENERATE RANDOM SAMPLE OF DATA
        d<-data2[sample(nrow(data2),ssizes[i],replace=FALSE),]
        nmonths<-length(unique(d$month))
        
            mfun<-function(d2){
            nm<-gsub('p','f',names(d2)[2])
            names(d2)[1:2]<-c('x','y')            
                                        #FIT MODEL - THE WORKHORSE
            dum<-data.frame(x=sort(unique(d2$x)),
                            sump=tapply(d2$y,d2$x,mean))
            dum$ct1<-dum$x*dum$sump
            ct<-sum(dum$ct1)/sum(dum$sump)
            dout<-data.frame(nmonth=length(unique(d2$month)),
                  n=dim(d2)[1],
                  ct=ct,
                  nm=nm)
            return(dout)
            }
            datout<-rbind(mfun(subset(d,select=c('day','p2','month'))),
                          mfun(subset(d,select=c('day','p4','month'))))
            
            datout
           },simplify=FALSE)

            l[[i]]<-data.frame(do.call('rbind.fill',dt))
            rm(dt)
    })                            
out2<-data.frame(do.call('rbind',l))#78 HOURS



dt<-replicate(50, {
        d<-data[sample(nrow(data),366,replace=TRUE),]
        nmonths<-length(unique(d$month))
        
            mfun<-function(d2){
            nm<-gsub('p','f',names(d2)[2])
            names(d2)[1:2]<-c('x','y')            
                                        #FIT MODEL - THE WORKHORSE
            dum<-data.frame(x=sort(unique(d2$x)),
                            sump=tapply(d2$y,d2$x,mean))
            dum$ct1<-dum$x*dum$sump
            ct<-sum(dum$ct1)/sum(dum$sump)
            dout<-data.frame(nmonth=length(unique(d2$month)),
                  n=dim(d2)[1],
                  ct=ct,
                  nm=nm)
            return(dout)
            }
            datout<-rbind(mfun(subset(d,select=c('day','p1','month'))),
                          mfun(subset(d,select=c('day','p3','month'))))
            datout
},simplify=FALSE)

outt<-data.frame(do.call('rbind.fill',dt))


par(mfrow=c(2,2))



plot(out$nmonth,out$ct)
mn<-mean(outt$ct)
sdd<-sd(outt$ct)
rect(min(out$nmonth),mn-(1.96*sdd),max(out$nmonth),mn+(1.96*sdd),col=alpha('lightgray',.3),border=NA)

dum<-data.frame(month=sort(unique(out$nmonth)),
                mn=tapply(out$ct,out$nmonth,mean),
                sdd=tapply(out$ct,out$nmonth,sd))
dum$upr<-dum$mn+(1.96*dum$sdd)
dum$lwr<-dum$mn-(1.96*dum$sdd)
f<-function(d){   lines(c(d$month,d$month),c(d$upr,d$lwr),col='red',lwd=2)}
z<-dlply(dum,.(month),.fun=f)
points(dum$month,dum$mn,pch=15,col='red')

plot(out2$nmonth,out2$ct,ylim=c(0,365))
plot(out2$n,out2$ct,ylim=c(0,365))

plot(out$nmonth,out$ct,ylim=c(0,365))
plot(out$n,out$ct,ylim=c(0,365))
rect(min(out$n),mn-(1.96*sdd),max(out$n),mn+(1.96*sdd),col=alpha('lightgray',.3),border=NA)


      
out$ct2<-log10(out$ct)
mod<-lm(ct2~n,data=out)
pdat<-data.frame(xx=seq(min(out$n),max(out$n),length.out=500))
#pdat<-data.frame(xx=sort(unique(out$n)))
#p<-data.frame(predict(mod,newdata=data.frame(n=pdat$xx),interval='prediction',level=0.95))
p<-predict(mod,newdata=data.frame(n=pdat$xx),se.fit=TRUE,type='response')
pdat$p<-10^p$fit
pdat$se<-10^p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)

lines(pdat$xx,pdat$p,col='red')
lines(pdat$xx,pdat$upr,col='red',lty=2)
lines(pdat$xx,pdat$lwr,col='red',lty=2)



mod<-gam(ct~nmonth,data=out,family='quasipoisson'(link='log'))
pdat<-data.frame(xx=seq(min(out$nmonth),max(out$nmonth),length.out=500))
p<-predict(mod,newdata=data.frame(nmonth=pdat$xx),se.fit=TRUE,type='response')
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)

par(mfrow=c(2,2))
lines(pdat$xx,pdat$p,col='red')
lines(pdat$xx,pdat$upr,col='red',lty=2)
lines(pdat$xx,pdat$lwr,col='red',lty=2)






md<-lm(ct2lg~n,data=out)
mds<-summary(md)
k.strt<-mds$coef[2,1]
lwr<-c(k.strt-0.01)
upr<-c(k.strt+0.01)
a.strt<-mds$coef[1,1]
lwra<-c(a.strt-0.01)
upra<-c(a.strt+0.01)

modnl<-nls(ct2 ~ a*exp(-k * n), data=out, start=list(k=k.strt,a=a.strt),control=nls.control(maxiter=500),algorithm='port')
pdat<-data.frame(xx=seq(min(out$n),max(out$n),length.out=500))
p<-predict(modnl,newdata=data.frame(n=pdat$xx),se.fit=TRUE)
pdat$p<-p
plot(out$n,out$ct2,col=alpha('blue',.1),pch=16)
lines(pdat$xx,pdat$p,col='red')





plot(out$n,out$ct2,col=alpha('blue',.1),pch=16)
mod<-glm(ct2~n,data=out,family='quasipoisson'(link='log'))
pdat<-data.frame(xx=seq(min(out$n),max(out$n),length.out=10))
p<-predict(mod,newdata=data.frame(n=pdat$xx),se.fit=TRUE,type='response')
p<-predict(mod,newdata=data.frame(n=pdat$xx),interval='prediction',level=0.95,type='response')
#p<-data.frame(predict(mod,newdata=data.frame(n=pdat$xx),interval='prediction',level=0.95))
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)
lines(pdat$xx,pdat$p,col='red')
lines(pdat$xx,pdat$lwr,col='red',lty=2)
lines(pdat$xx,pdat$upr,col='red',lty=2)



hist(out$ct2)
hist(log(out$ct2+1))

mean(out$ct)
183*2

setwd(datadir)
#write.csv(out,'simulation1.csv',row.names=FALSE)
write.csv(out,'simulation2.csv',row.names=FALSE)
out<-read.csv('simulation2.csv',header=TRUE)

plot(out2$nmonth,out2$ct)










#35 seconds for 1 iteration (50/2)
n<-77*50*25
(n)/60/60



#ADDS 'TRUE' CLUSTER MEMBERSHIPS
clfun<-function(d){
    nm<-names(d)[2]
    names(d)[1:2]<-c('x','y')
    q2<-q
    q2[,868]<-d$y#REPLACE COLUMN WITH CURRENT PHENOLOGY
    names(q2)[868]<-'y'
    
    #RE-RUNS FUZZY CLUSTERING WITH SUBSTITUTION
    k<-6#NUMBER OF CLUSTERS
    dmat<-1-cor(q2,use='pairwise.complete.obs')
    dst2<-as.dist(dmat)
    ff2<-fanny(dst2,k,maxit=500,diss=T)
    memb2<-data.frame(prb=ff2$membership[868,])

    #OUTPUTS RELEVANT INFO
    out<-data.frame(clus=ff2$clustering[868],
                    prb=max(memb2$prb))
    names(out)<-gsub(' ','',paste(nm,names(out)))
    rownames(out)<-NULL    
    return(out)
}
 clout<-cbind(clfun(subset(data,select=c('day','p1mn'))),
              clfun(subset(data,select=c('day','p2mn'))),
              clfun(subset(data,select=c('day','p3mn'))))

out$p1truecl<-clout$p1mnclus
out$p4truecl<-clout$p1mnclus
out$p2truecl<-clout$p2mnclus
out$p5truecl<-clout$p2mnclus
out$p3truecl<-clout$p3mnclus
out$p6truecl<-clout$p3mnclus
out$p1trueprb<-clout$p1mnprb
out$p4trueprb<-clout$p1mnprb
out$p2trueprb<-clout$p2mnprb
out$p5trueprb<-clout$p2mnprb
out$p3trueprb<-clout$p3mnprb
out$p6trueprb<-clout$p3mnprb

out$p1cat<-ifelse(out$cl1==out$p1truecl,1,0)
out$p2cat<-ifelse(out$cl2==out$p2truecl,1,0)
out$p3cat<-ifelse(out$cl3==out$p3truecl,1,0)
out$p4cat<-ifelse(out$cl4==out$p4truecl,1,0)
out$p5cat<-ifelse(out$cl5==out$p5truecl,1,0)
out$p6cat<-ifelse(out$cl6==out$p6truecl,1,0)

out$p1prb<-abs(out$p1trueprb-out$prb1)
out$p2prb<-abs(out$p2trueprb-out$prb2)
out$p3prb<-abs(out$p3trueprb-out$prb3)
out$p4prb<-abs(out$p4trueprb-out$prb4)
out$p5prb<-abs(out$p5trueprb-out$prb5)
out$p6prb<-abs(out$p6trueprb-out$prb6)

out$a1true<-max(data$p1mn,na.rm=TRUE)-min(data$p1mn,na.rm=TRUE)
out$a2true<-max(data$p2mn,na.rm=TRUE)-min(data$p2mn,na.rm=TRUE)
out$a3true<-max(data$p3mn,na.rm=TRUE)-min(data$p3mn,na.rm=TRUE)
out$p1a<-abs(out$a1true-out$a1)
out$p2a<-abs(out$a2true-out$a2)
out$p3a<-abs(out$a3true-out$a3)
out$p4a<-abs(out$a1true-out$a4)
out$p5a<-abs(out$a2true-out$a5)
out$p6a<-abs(out$a3true-out$a6)


out$b1true<-subset(data,p1mn==max(data$p1mn,na.rm=TRUE))$day
out$b2true<-subset(data,p2mn==max(data$p2mn,na.rm=TRUE))$day
out$b3true<-subset(data,p3mn==max(data$p3mn,na.rm=TRUE))$day
out$p1b<-abs(out$b1true-out$b1)
out$p2b<-abs(out$b2true-out$b2)
out$p3b<-abs(out$b3true-out$b3)
out$p4b<-abs(out$b1true-out$b4)
out$p5b<-abs(out$b2true-out$b5)
out$p6b<-abs(out$b3true-out$b6)


mod1<-gam(p1cat~contig + nmonths,data=out,family=quasibinomial,method='ML')
mod2<-gam(p2cat~contig + nmonths,data=out,family=quasibinomial,method='ML')
mod3<-gam(p3cat~contig + nmonths,data=out,family=quasibinomial,method='ML')
mod4<-gam(p4cat~contig + nmonths,data=out,family=quasibinomial,method='ML')
mod5<-gam(p5cat~contig + nmonths,data=out,family=quasibinomial,method='ML')
mod6<-gam(p6cat~contig + nmonths,data=out,family=quasibinomial,method='ML')

mod1<-gam(p1cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')
mod2<-gam(p2cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')
mod3<-gam(p3cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')
mod4<-gam(p4cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')
mod5<-gam(p5cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')
mod6<-gam(p6cat~s(contig,nmonths,k=10),data=out,family=quasibinomial,method='ML')

f<-function(d,ttl){
names(d)[1:3]<-c('x','y','y2')
plot(d$x,d$y,pch=16,col=alpha('darkred',.5),xlab='Day',ylab='Simulated phytoplankton',axes=FALSE)
abline(h=0,lty=2,lwd=.25)
lines(d$x,d$y2,col='blue3',lwd=2)
axis(1,c(1,365),labels=TRUE,tick=TRUE)
axis(2,seq(-3,3,1),las=1)
mtext(ttl,side=3,line=-1,cex=.75)
}
setwd(figsdir)
pdf('simulation_patterns.pdf',height=5,width=9)
par(mfrow=c(2,3),mar=c(4,4,1,1),mgp=c(3,1,.75))
f(subset(data,select=c('day','p1','p1mn')),'Phenotype 1; SD=0.18')
f(subset(data,select=c('day','p2','p2mn')),'Phenotype 2; SD=0.18')
f(subset(data,select=c('day','p3','p3mn')),'Phenotype 3; SD=0.18')
f(subset(data,select=c('day','p4','p1mn')),'Phenotype 1; SD=0.54')
f(subset(data,select=c('day','p5','p2mn')),'Phenotype 2; SD=0.54')
f(subset(data,select=c('day','p6','p3mn')),'Phenotype 3; SD=0.54')
dev.off()
d<-(subset(out,select=c('nmonths','contig','p1cat')))


kk<-10

f<-function(d,ttl){
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')    
d$y<-as.factor(d$y)
mod<-gam(y~s(x1,x2,k=kk),data=d,family=quasibinomial('logit'))
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
cls<-colorRampPalette(c('magenta','blue3','green','yellow','red3'))
xx1<-seq(min(d$x1),max(d$x1),length.out=100)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=100)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-predict(mod,newdata=prc2,type='response')
dat<-acast(prc2,x1~x2,value.var="p")
x<-dat
xmin=0; xmax=1;
x[x<xmin]=xmin; x[x>xmax]=xmax;
collist<-c('magenta','blue3','green','yellow','red3')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
image(x=sort(unique(xx1)),y=sort(unique(xx2)),x,col=ColorRamp_ex,xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(xmin,xmax),col=ColorRamp,legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
points(9,75,col='gold',pch=16,cex=2)
points(9,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=9,x2=75),type='response')
text(9.2,75,round(p,digits=2),col='black',adj=0,cex=1.5)
points(8,75,col='gold',pch=16,cex=2)
points(8,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=8,x2=75),type='response')
text(8.2,80,round(p,digits=2),col='black',adj=0,cex=1.5)

}
pdf('simulation_cluster_categorization.pdf',height=7,width=12)
par(mfrow=c(2,3),mar=c(4,4,2,4))
f(subset(out,select=c('nmonths','contig','p1cat')),'Cluster category [P1; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p2cat')),'Cluster category [P2; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p3cat')),'Cluster category [P3; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p4cat')),'Cluster category [P1; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p5cat')),'Cluster category [P2; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p6cat')),'Cluster category [P3; SD=0.54]')

           
f<-function(d,ttl){
cls<-(colorRampPalette(c('magenta','blue3','green','yellow','red3')))
nm<-names(d)[3]
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')
#d$y<-d$y^.25
mod<-gam(y~s(x1,x2,k=kk),data=d,method='ML')
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
xx1<-seq(min(d$x1),max(d$x1),length.out=150)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=150)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-predict(mod,newdata=prc2,type='response')
dat<-acast(prc2,x1~x2,value.var="p")
print(summary(prc2$p))

x<-dat
xmin=0.0; xmax=0.3;
x[x<xmin]=xmin; x[x>xmax]=xmax;
collist<-c('red3','yellow','green','blue3','magenta')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
image(x=sort(unique(xx1)),y=sort(unique(xx2)),x,col=ColorRamp_ex,xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(xmin,xmax),col=ColorRamp,legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
points(9,75,col='gold',pch=16,cex=2)
points(9,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=9,x2=75),type='response')
text(9.2,75,round(p,digits=2),col='black',adj=0,cex=1.5)
points(8,75,col='gold',pch=16,cex=2)
points(8,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=8,x2=75),type='response')
text(8.2,80,round(p,digits=2),col='black',adj=0,cex=1.5)
}
par(mfrow=c(2,3))
f(subset(out,select=c('nmonths','contig','p1prb')),'Cluster probability [P1; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p2prb')),'Cluster probability [P2; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p3prb')),'Cluster probability [P3; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p4prb')),'Cluster probability [P1; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p5prb')),'Cluster probability [P2; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p6prb')),'Cluster probability [P3; SD=0.54]')


f<-function(d,ttl){
cls<-colorRampPalette(c('magenta','blue3','green','yellow','red3'))
nm<-names(d)[3]
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')
mod<-gam(y~s(x1,x2,k=kk),data=d,method='ML')
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
xx1<-seq(min(d$x1),max(d$x1),length.out=150)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=150)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-predict(mod,newdata=prc2,type='response')
#prc2$p<-ifelse(prc2$p==min(prc2$p),0.5,prc2$p)
dat<-acast(prc2,x1~x2,value.var="p")
print(summary(prc2$p))

x<-dat
xmin=0.5; xmax=1;
x[x<xmin]=xmin; x[x>xmax]=xmax;
collist<-c('magenta','blue3','green','yellow','red3')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
image(x=sort(unique(xx1)),y=sort(unique(xx2)),x,col=ColorRamp_ex,xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(xmin,xmax),col=ColorRamp,legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
points(9,75,col='gold',pch=16,cex=2)
points(9,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=9,x2=75),type='response')
text(9.2,75,round(p,digits=2),col='black',adj=0,cex=1.5)
points(8,75,col='gold',pch=16,cex=2)
points(8,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=8,x2=75),type='response')
text(8.2,80,round(p,digits=2),col='black',adj=0,cex=1.5)
}
par(mfrow=c(2,3))
f(subset(out,select=c('nmonths','contig','r1')),'Correlation [P1; SD=0.18]')
f(subset(out,select=c('nmonths','contig','r2')),'Correlation [P2; SD=0.18]')
f(subset(out,select=c('nmonths','contig','r3')),'Correlation [P3; SD=0.18]')
f(subset(out,select=c('nmonths','contig','r4')),'Correlation [P1; SD=0.54]')
f(subset(out,select=c('nmonths','contig','r5')),'Correlation [P2; SD=0.54]')
f(subset(out,select=c('nmonths','contig','r6')),'Correlation [P3; SD=0.54]')



f<-function(d,ttl,xmn,xmx){
cls<-(colorRampPalette(c('magenta','blue3','green','yellow','red3')))
nm<-names(d)[3]
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')
d$y<-log10(d$y+.01)
mod<-gam(y~s(x1,x2,k=kk),data=d,method='ML')
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
xx1<-seq(min(d$x1),max(d$x1),length.out=150)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=150)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-10^(predict(mod,newdata=prc2,type='response'))-.01
dat<-acast(prc2,x1~x2,value.var="p")
print(summary(prc2$p))

x<-dat
xmin=xmn; xmax=xmx
x[x<xmin]=xmin; x[x>xmax]=xmax;
#collist<-c('magenta','blue3','green','yellow','red3')
collist<-c('red3','yellow','green','blue3','magenta')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
image(x=sort(unique(xx1)),y=sort(unique(xx2)),x,col=ColorRamp_ex,xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(xmin,xmax),col=ColorRamp,legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
points(9,75,col='gold',pch=16,cex=2)
points(9,75,pch=1,cex=2)
p<-10^(predict(mod,newdata=data.frame(x1=9,x2=75),type='response'))-.01
text(9.2,75,round(p,digits=2),col='black',adj=0,cex=1.5)
points(8,75,col='gold',pch=16,cex=2)
points(8,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=8,x2=75),type='response')
text(8.2,80,round(p,digits=2),col='black',adj=0,cex=1.5)
}
par(mfrow=c(2,3))
f(subset(out,select=c('nmonths','contig','p1a')),'Amplitude [P1; SD=0.18]',0.04,1.3)
f(subset(out,select=c('nmonths','contig','p2a')),'Amplitude [P2; SD=0.18]',0.04,1.3)
f(subset(out,select=c('nmonths','contig','p3a')),'Amplitude [P3; SD=0.18]',0.04,1.3)
f(subset(out,select=c('nmonths','contig','p4a')),'Amplitude [P1; SD=0.54]',0.04,1.3)
f(subset(out,select=c('nmonths','contig','p5a')),'Amplitude [P2; SD=0.54]',0.04,1.3)
f(subset(out,select=c('nmonths','contig','p6a')),'Amplitude [P3; SD=0.54]',0.04,1.3)
dev.off()



ptfun<-function(d,xlab,ylab){
    names(d)[1:2]<-c('x','y')
    plot(d$x,d$y,pch=16,col=alpha('darkred',.2),las=1,cex=1.5,xlab=xlab,ylab=ylab)
}
setwd(figsdir)
pdf('simulation_sample_sizes.pdf',height=10,width=5)
par(mfrow=c(3,1),mar=c(4,4,1,1))
ptfun(subset(out,select=c('n','nmonths')),'N','Number of Months')
ptfun(subset(out,select=c('n','contig')),'N','Longest gap in series')
ptfun(subset(out,select=c('nmonths','contig')),'Number of months','Longest gap in series')
dev.off()






d<-subset(out,select=c('nmonths','contig','p1cat'))

d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')
d$x2<-round(d$x2/10,digits=0)
            
f2<-function(a){
    return(data.frame(y=mean(a$y)))
}
dat<-ddply(d,.(x1,x2),.fun=f2)

n<-21
brks<-seq((min(dat$y,na.rm=TRUE)-.01),max(dat$y,na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq((min(dat$y,na.rm=TRUE)-.01),max(dat$y,na.rm=TRUE)+.01,length.out=n),digits=3)
dat$ycat<-cut(dat$y,breaks=brks)
lbls<-sort(unique(dat$ycat))
lbls2<-sort(unique(cut(dat$y,breaks=brks2)))
cls<-colorRampPalette(c('blue3','green','yellow','red3'))

ggplot()+
geom_tile(data=dat, aes(x=x1, y=x2,fill=ycat),color='white',size=.1)+
scale_fill_manual(breaks=as.character(lbls),values=cls(length(lbls)),labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),plot.background = element_rect(fill = 'green', colour = 'red'))+
coord_equal()

+
scale_x_discrete(expand=c(0,0),breaks=seq(5,13,1),labels=seq(5,13,1))+
scale_y_discrete(expand=c(0,0),breaks=seq(0,30,5),labels=seq(0,30,5))+
coord_cartesian(ylim=c(0,30),xlim=c(5,13))+
    xlab('')+
    ylab('')



f<-function(d,ttl){
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')    
d$y<-as.factor(d$y)
mod<-gam(y~s(x1,x2,k=8),data=d,family=quasibinomial('logit'))
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
cls<-colorRampPalette(c('magenta','blue3','green','yellow','red3'))
xx1<-seq(min(d$x1),max(d$x1),length.out=100)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=100)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-predict(mod,newdata=prc2,type='response')
dat<-acast(prc2,x1~x2,value.var="p")
x<-dat
xmin=0; xmax=1;
x[x<xmin]=xmin; x[x>xmax]=xmax;
collist<-c('magenta','blue3','green','yellow','red3')
ColorRamp<-colorRampPalette(collist)(10000)
ColorLevels<-seq(from=xmin, to=xmax, length=10000)
ColorRamp_ex <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
image(x=sort(unique(xx1)),y=sort(unique(xx2)),x,col=ColorRamp_ex,xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(xmin,xmax),col=ColorRamp,legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
points(9,75,col='gold',pch=16,cex=2)
points(9,75,pch=1,cex=2)
p<-predict(mod,newdata=data.frame(x1=9,x2=75),type='response')
text(9.2,75,round(p,digits=2),col='black',adj=0,cex=1.5)
}
pdf('simulation_cluster_categorization.pdf',height=7,width=12)
par(mfrow=c(2,3),mar=c(4,4,2,4))
f(subset(out,select=c('nmonths','contig','p1cat')),'Cluster category [P1; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p2cat')),'Cluster category [P2; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p3cat')),'Cluster category [P3; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p4cat')),'Cluster category [P1; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p5cat')),'Cluster category [P2; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p6cat')),'Cluster category [P3; SD=0.54]')











###
####
f<-function(d,ttl){
d<-na.omit(d)    
names(d)[1:3]<-c('x1','x2','y')    
d$y<-as.factor(d$y)
mod<-gam(y~s(x1,x2,k=10),data=d,family=quasibinomial('logit'))
s<-summary(mod)
r2<-round(s$dev.expl,digits=2)
cls<-colorRampPalette(c('magenta','blue3','green','yellow','red3'))
xx1<-seq(min(d$x1),max(d$x1),length.out=100)#XAXIS
xx2<-seq(min(d$x2),max(d$x2),length.out=100)#YAXIS
prc2<-expand.grid(x1=xx1,x2=xx2)
prc2$p<-predict(mod,newdata=prc2,type='response')
dat<-acast(prc2,x1~x2,value.var="p")
image(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,col=cls(100),xlab='Number of months',ylab='Number of missing days',las=1,cex.axis=.8,main=ttl,useRaster=TRUE)
contour(x=sort(unique(xx1)),y=sort(unique(xx2)),dat,add=TRUE,lwd=.75,col='gray40',labcex=.4,nlevels=6)
image.plot(legend.only=TRUE,zlim=c(min(prc2$p),max(prc2$p)),col=cls(100),legend.lab="",horizontal=FALSE,smallplot=c(.88,.90,.215,.85),cex=.25)
legend('topleft',legend=gsub(' ','',paste('r2 = ',r2)),fill=TRUE,bg='white')
}
pdf('simulation_cluster_categorization.pdf',height=7,width=12)
par(mfrow=c(2,3),mar=c(4,4,2,4))
f(subset(out,select=c('nmonths','contig','p1cat')),'Cluster category [P1; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p2cat')),'Cluster category [P2; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p3cat')),'Cluster category [P3; SD=0.18]')
f(subset(out,select=c('nmonths','contig','p4cat')),'Cluster category [P1; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p5cat')),'Cluster category [P2; SD=0.54]')
f(subset(out,select=c('nmonths','contig','p6cat')),'Cluster category [P3; SD=0.54]')
dev.off()
