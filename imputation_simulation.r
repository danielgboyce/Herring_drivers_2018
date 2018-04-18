library(randomForest)
library(rcompanion)
library(randomForest)
library(DataCombine)
library(lattice)
library(tidyr)
library(lavaan)
library(akima)
library(corrplot)
library(stringr)
library(RColorBrewer)
library(segmented)
library(splines)
library(strucchange)
library(data.table)
library(psych)
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
library(fossil)
library(bnstruct)
library(mice)
library(VIM)
library(progress)

datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-'N://cluster_2017//scratch//spera//data//finaldat_v2'
datadirout<-'N://cluster_2017//scratch//spera//data//dataoutput'
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
figsdir<-'C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures'





setwd(datadir)
#load('SPERA_andata.RData')
load('SPERA_andata_new.RData')
data$sst.t12<-ifelse(data$sst.t12==Inf,NA,data$sst.t12)
nms<-names(data)
nms<-nms[grepl('.se',nms)==FALSE]
nms<-nms[grepl('.tplus',nms)==FALSE]
nms<-nms[grepl('.tmin',nms)==FALSE]
dat<-data[,names(data) %in% nms]
dat<-subset(dat,year>=1969 & year<=2011)
dat<-dat[,!(names(dat)%in%c('her.durng.rv','her.dumx.rv','herjuv.durng.rv','herjuv.dumxrv','herjuv.dvm2.rv','herjuv.dur2.rv','herlrv.sp.r','herjuv.dep.rv','herjuv.sp.n','her.sp.n','herlrv.sp.n','her.szdiv.rv','her.dur2.rv','herjuv.sprng','herjuv.spcv','herjuv.spvar','herjuv.spnug','herlrv.sprng','herlrv.spcv','herlrv.spvar','herlrg.spnug','herlrv.georng','herlrv.spnug'))]


#Z-STANDARDIZE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

mth<-'rf'
tempdata<-mice(dat[,-1],m=5,maxit=50,meth=mth,seed=500,print=FALSE,diagnostics=TRUE)
#GET COMPLETE IMPUTED DATASET
impdat<-complete(tempdata,1)
impdat$year<-dat$year

#GET AVERAGE OF ALL 5 IMPUTATIONS
id1<-as.matrix(complete(tempdata,1))
id2<-as.matrix(complete(tempdata,2))
id3<-as.matrix(complete(tempdata,3))
id4<-as.matrix(complete(tempdata,4))
id5<-as.matrix(complete(tempdata,5))
impdatav<-(id1+id2+id3+id4+id5)/5
impdatav<-data.frame(impdatav)
impdatav$year<-dat$year

save(impdatav,file='impdatav.RData')


####################################################

#SAME BUT NORMALIZE THE DATA BEFORE TRANSFORMING
setwd(datadir)
#load('SPERA_andata.RData')
load('SPERA_andata_new.RData')

nms<-names(data)
nms<-nms[grepl('.se',nms)==FALSE]
nms<-nms[grepl('.tplus',nms)==FALSE]
nms<-nms[grepl('.tmin',nms)==FALSE]
dat<-data[,names(data) %in% nms]
dat<-subset(dat,year>=1969 & year<=2011)
dat<-dat[,!(names(dat)%in%c('her.durng.rv','her.dumx.rv','herjuv.durng.rv','herjuv.dumxrv','herjuv.dvm2.rv','herjuv.dur2.rv','herlrv.sp.r','herjuv.dep.rv','herjuv.sp.n','her.sp.n','herlrv.sp.n','her.szdiv.rv','her.dur2.rv','herjuv.sprng','herjuv.spcv','herjuv.spvar','herjuv.spnug','herlrv.sprng','herlrv.spcv','herlrv.spvar','herlrg.spnug','herlrv.georng','herlrv.spnug'))]

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
f<-function(x){transformTukey(x,plotit=FALSE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

#Z-STANDARDIZE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

mth<-'rf'
tempdata<-mice(dat[,-1],m=5,maxit=50,meth=mth,seed=500,print=FALSE)
#GET COMPLETE IMPUTED DATASET
impdat<-complete(tempdata,1)
impdat$year<-dat$year

#GET AVERAGE OF ALL 5 IMPUTATIONS
id1<-as.matrix(complete(tempdata,1))
id2<-as.matrix(complete(tempdata,2))
id3<-as.matrix(complete(tempdata,3))
id4<-as.matrix(complete(tempdata,4))
id5<-as.matrix(complete(tempdata,5))
impdatav<-(id1+id2+id3+id4+id5)/5
impdatav<-data.frame(impdatav)
impdatav$year<-dat$year

impdatav.std<-impdatav
save(impdatav.std,file='impdatav.std.RData')






######################################################
## CODE TO LAG SOME VARS
rspnms<-c('her.ssb','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.spnug','her.spcv','her.spvar','her.prod','her.metai.rv','her.state','her.ssbc')
for(i in 1:length(rspnms)){
nm<-rspnms[i]
dat<-slide(dat,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
dat<-slide(dat,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
dat<-slide(dat,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
dat<-slide(dat,Var=nm,slideBy=4,NewVar=gsub(' ','',paste(nm,'.t4')))
dat<-slide(dat,Var=nm,slideBy=5,NewVar=gsub(' ','',paste(nm,'.t5')))
dat<-slide(dat,Var=nm,slideBy=6,NewVar=gsub(' ','',paste(nm,'.t6')))
}

nms<-names(data)
clmnms<-nms[grepl('asl',nms)==TRUE |
         grepl('gs',nms)==TRUE |
         grepl('ss\\.dist',nms)==TRUE |
         grepl('sst',nms)==TRUE |
         grepl('wind',nms)==TRUE |
         grepl('wnd',nms)==TRUE |
         grepl('temp',nms)==TRUE |
         grepl('nut',nms)==TRUE |
         grepl('t50',nms)==TRUE |
         grepl('nit',nms)==TRUE |
         grepl('phos',nms)==TRUE |
         grepl('sil',nms)==TRUE |
         grepl('had\\.pi',nms)==TRUE]
clmnms<-clmnms[grepl('\\.se',clmnms)==FALSE]
clmnms<-clmnms[grepl('tmin',clmnms)==FALSE]
for(i in 1:length(clmnms)){
nm<-clmnms[i]
dat<-slide(dat,Var=nm,slideBy=-1,NewVar=gsub(' ','',paste(nm,'.t-1')))
dat<-slide(dat,Var=nm,slideBy=-2,NewVar=gsub(' ','',paste(nm,'.t-2')))
}

nms<-names(data)
rspnms<-nms[grepl('land',nms)==TRUE |
         grepl('bfish',nms)==TRUE |
         grepl('expr',nms)==TRUE]
rspnms<-rspnms[grepl('\\.se',rspnms)==FALSE]
for(i in 1:length(rspnms)){
nm<-rspnms[i]
dat<-slide(dat,Var=nm,slideBy=-1,NewVar=gsub(' ','',paste(nm,'.t-1')))
dat<-slide(dat,Var=nm,slideBy=-2,NewVar=gsub(' ','',paste(nm,'.t-2')))
}






setwd(datadir)
load('SPERA_andata.RData')
load('SPERA_andata_new.RData')

nms<-names(data)
nms<-nms[grepl('.se',nms)==FALSE]
nms<-nms[grepl('.tplus',nms)==FALSE]
nms<-nms[grepl('.tmin',nms)==FALSE]
nms<-nms[grepl('.sabs',nms)==FALSE]
nms<-nms[grepl('state',nms)==FALSE]
dat<-data[,names(data) %in% nms]
dat<-dat[,!(names(dat)%in%c('her.durng.rv','her.dumx.rv','sst.stalt','her.szpe.rv','herjuv.durng.rv'))]


#GET PROPORTION OF EACH COLUMN THAT IS MISSING
dat<-subset(dat,year>=1970 & year<=2006)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
dm<-data.frame(var=names(dat),
               pmis=apply(dat,2,pMiss))
dm<-dm[order(dm$pmis),]
#apply(dat,1,pMiss)
dm2<-subset(dm,pmis<3)
dat2<-dat[,names(dat) %in% as.character(dm2$var)]
dat2<-na.omit(dat2)

#EXAMINE DISTRIBUTION OF INDIVIDUAL VARIABLES STANDARDIZED TO Z
#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
b<-data.frame(cbind(subset(dat2,select='year'),apply(dat2[,2:dim(dat2)[2]],2,f)))
b<-b %>% gather(var, y, her.spnug:sst.bbay)
histogram( ~ y | var,data=b)

#image(md.pattern(q))
aggr_plot<-aggr(q, col=c('dodgerblue4','firebrick3'), numbers=TRUE, sortVars=TRUE, labels=names(q), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))



################################################

################################################
#setwd(figsdir)
#pdf('imputed_scatterplot.pdf',height=10,width=10)

cr<-cor(dat2,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
names(cr.t)[3]<-'r.real'


#########################################################

##############  IMPUTATION SIMULATION REPEATED 100 TIMES

#########################################################
impsim<-function(){

smp<-c(.25,.2,.15,.1)
lll<-list()
zzz<-list()
for(g in 1:length(smp)){

#GETS RANDOM SAMPLE (DATAFRAME) ACCORDING TO SPECIFIED MISSINNESS [G]
nms<-names(dat2)[2:dim(dat2)[2]]
l<-list()
for(i in 1:length(nms)){
n<-dim(dat2)[1]-ceiling(38*smp[g])#SAMPLE SIZE TO SAMPLE
    d<-subset(dat2,select=c('year',paste(nms[i])))
d<-d[sample(nrow(d),n,replace=FALSE),]
l[[i]]<-d
}
q<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), l)#COMBINE

#FUNCTION GETS TRUE VALUES FROM ORIGINAL DATASET THAT WERE REMOVED AND THEN IMPUTED; NON IMPUTED VALUES ARE NA'S
nms<-names(dat2)[2:dim(dat2)[2]]
l<-list()
for(i in 1:length(nms)){
d<-subset(dat2,select=c('year',paste(nms[i])))#TRUE DATA
d2<-subset(q,select=c('year',paste(nms[i])))#SAMPLED DATA
names(d2)[2]<-'y'
dout<-merge(d,d2,by=c('year'),all=TRUE)
dout<-dout[is.na(dout[,3])==TRUE,]
l[[i]]<-dout[,1:2]
}
qreal<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), l)#COMBINE

#DIFFERENT IMPUTATION METHODS TO TEST
#mth<-c('norm.predict','cart','pmm','sample','rf')
mth<-c('cart','pmm','sample','rf')
ll<-list()
zz<-list()
for(h in 1:length(mth)){

print(smp[g])
print(mth[h])
tempdata<-mice(q[,-1],m=5,maxit=50,meth=mth[h],seed=500,print=FALSE)
#GET COMPLETE IMPUTED DATASET
impdat<-complete(tempdata,1)
impdat$year<-q$year

#GET AVERAGE OF ALL 5 IMPUTATIONS
id1<-as.matrix(complete(tempdata,1))
id2<-as.matrix(complete(tempdata,2))
id3<-as.matrix(complete(tempdata,3))
id4<-as.matrix(complete(tempdata,4))
id5<-as.matrix(complete(tempdata,5))
impdatav<-(id1+id2+id3+id4+id5)/5
impdatav<-data.frame(impdatav)
impdatav$year<-q$year

#AVERAGE CORRELATION BETWEEN ALL IMPUTED
cr<-cor(impdat,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
crimp = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
crimp=crimp[crimp$Var1 != crimp$Var2,]
crimp=crimp[paste(crimp$Var1,crimp$Var2,sep="_") %in% combinations,]
names(crimp)[3]<-'r.imp1'

#AVERAGE CORRELATION BETWEEN AVERAGE OF ALL IMPUTED
cr<-cor(impdatav,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
crimpav = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
crimpav=crimpav[crimpav$Var1 != crimpav$Var2,]
crimpav=crimpav[paste(crimpav$Var1,crimpav$Var2,sep="_") %in% combinations,]
names(crimpav)[3]<-'r.imp5'

aa<-merge(cr.t,crimp,by=c('Var1','Var2'),all=TRUE)
aa<-merge(aa,crimpav,by=c('Var1','Var2'),all=TRUE)

nms<-names(dat2)[2:dim(qreal)[2]]
l<-list()
for(i in 1:length(nms)){
d<-subset(qreal,select=c('year',paste(nms[i])))#TRUE VALUES THAT WERE LATER IMPUTED
d2<-subset(impdat,select=c('year',paste(nms[i])))#IMPUTED VALUES
d3<-subset(dat2,select=c('year',paste(nms[i])))#ORIGINAL DATA
d4<-subset(impdatav,select=c('year',paste(nms[i])))#AVERAGE IMPUTED
names(d2)[2]<-'imputed'
dout<-merge(d,d2,by=c('year'),all=TRUE)
dout<-na.omit(dout)
r<-round(cor(dout[,2],dout[,3]),digits=2)#CORRELAION BETWEEN ONLY 'TRUE' VALUES AND THOSE THAT WERE IMPUTED

#plot(dout[,2],dout[,3],pch=16,las=1)
#legend('topright',legend=paste('r=',r),bty='n',title=tempdata$meth[1])
#legend('topleft',legend=paste('nmiss=',smp[g]),bty='n')
#abline(lm(dout[,3]~dout[,2]))

names(d4)[2]<-'imputed'
dout<-merge(d,d4,by=c('year'),all=TRUE)
dout<-na.omit(dout)
r2<-round(cor(dout[,2],dout[,3]),digits=2)#CORRELATION BETWEEN ONLY 'TRUE' VALUES AND AVERAGE OF TOP 5 IMPUTATIONS

#points(dout[,2],dout[,3],pch=16,las=1,col='red3')
#legend('topright',legend=paste('r=',r),bty='n',title=tempdata$meth[1])
#abline(lm(dout[,3]~dout[,2]),col='red')

l[[i]]<-data.frame(var=nms[i],
                   r=r,
                   r2=r2,
                   sddata=round(sd(d3[,2],na.rm=TRUE),digits=2),
                   meth=tempdata$meth[1],
                   smp=smp[g])
}
ll[[h]]<-data.frame(do.call('rbind',l))
zz[[h]]<-aa
}
lll[[g]]<-data.frame(do.call('rbind',ll))#RETURNS R AMONG IMPUTED
zzz[[g]]<-data.frame(do.call('rbind',zz))#RETURNS R AMONG ALL
}

list(data.frame(do.call('rbind',lll)),data.frame(do.call('rbind',zzz)))
}


#R.IMP1 IS AVERAGE OF ALL INDIVIDUAL IMPUTATIONS
#R.IMP5 IS AVERAGE OF AVERAGE OF ALL IMPUTED
#R IS CORRELATION BETWEEN TRUE VALUES AND IMPUTED
#R2 IS CORRELATION BETWEEN TRUE VALUES AND AVERAGE OF TOP 5 IMPUTATIONS
library(pbapply)
dout00<-pbreplicate(2,impsim(),simplify=FALSE)
dout00<-pbreplicate(100,impsim(),simplify=FALSE)
#dout00<-replicate(2,impsim(),simplify=FALSE)

#GETS CORRELATION BETWEEN IMPUTED AND REAL
f<-function(d){ return(data.frame(d[1]))}
dout<-ldply(dout00,.fun=f)
#GETS AVERAGE CORRELATION BETWEEN 1) REAL DATA; IMPUTED DATA(1) WITH REAL DATA; 3) IMPUTED DATA(5) WITH REAL DATA
f<-function(d){ return(data.frame(d[2]))}
dout2<-ldply(dout00,.fun=f)


#save(dout,file='imputation.simdata.RData')
setwd(figsdir)
load('imputation.simdata.RData')
dout00<-dout

f<-function(d){
    return(data.frame(r=mean(d$r,na.rm=TRUE),
                      r.sd=sd(d$r,na.rm=TRUE),
                      r2=mean(d$r2,na.rm=TRUE),
                      r2.sd=sd(d$r2,na.rm=TRUE)))
}
ot<-ddply(dout,.(meth,smp),.fun=f)
ot$id<-seq(1,dim(ot)[1],1)
dm<-data.frame(meth=sort(unique(ot$meth)),
    cl=c('firebrick3','dodgerblue3','green3','orange','magenta'))
ot<-merge(ot,dm,by=c('meth'),all=TRUE)


setwd(figsdir)
pdf('imputation_simulation_compare.pdf',height=10,width=7)
par(mfrow=c(2,1),mar=c(4,4,1,1))

plot(ot$r,ot$r2,pch=16,xlim=c(-.1,.4),ylim=c(-.1,.4),cex=ot$smp*10,las=1,col=alpha(as.character(ot$cl),.7),xlab='Average r between imputed and real values (top imp)',ylab='Average r between imputed and real values (top 5 imps)')
points(ot$r,ot$r2,pch=1,cex=ot$smp*10,col=alpha('black',.5),lwd=.02)
abline(a=0,b=1,lty=2)
legend('topleft',legend=dm$meth,col=as.character(dm$cl),pch=15,bty='n')
f<-function(d){
    lines(c(d$r,d$r),c(d$r2-(1*d$r2.sd),d$r2+(1*d$r2.sd)))
    lines(c(d$r+(1*d$r.sd),d$r-(1*d$r.sd)),c(d$r2,d$r2))
}
#z<-dlply(ot,.(id),.fun=f)

d<-data.frame(sddata=sort(unique(dout$sddata)),
              r=tapply(dout$r,dout$sddata,mean),
              r2=tapply(dout$r2,dout$sddata,mean))

d$lsddata<-log10(d$sddata)
plot(d$lsddata,d$r2,las=1,xlab='Average log10 variance of timeseries',ylab='Average r between imputed and real data',pch=16,xlim=c(-1.75,2),ylim=c(-.05,.5))
mod<-gam(r2~s(lsddata,k=10),data=subset(d,lsddata<3),gamma=1.4)
pdat<-data.frame(lsddata=seq(min(d$lsddata),max(d$lsddata),length.out=100))
pdat$p<-predict(mod,newdata=pdat)
lines(pdat$lsddata,pdat$p)
s<-summary(mod)
legend('topright',legend=paste('r2=',round(s$r.sq,digits=2)),bty='n')
dev.off()



library(lattice)
library(tidyr)
library(lavaan)
library(akima)
library(corrplot)
library(stringr)
library(RColorBrewer)
library(segmented)
library(splines)
library(strucchange)
library(data.table)
library(psych)
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
library(fossil)
library(bnstruct)
library(mice)
library(VIM)

datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-'N://cluster_2017//scratch//spera//data//finaldat_v2'
datadirout<-'N://cluster_2017//scratch//spera//data//dataoutput'
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
figsdir<-'C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures'



setwd(datadir)
load('SPERA_andata.RData')
load('SPERA_andata_new.RData')

nms<-names(data)
nms<-nms[grepl('.se',nms)==FALSE]
nms<-nms[grepl('.tplus',nms)==FALSE]
nms<-nms[grepl('.tmin',nms)==FALSE]
nms<-nms[grepl('.sabs',nms)==FALSE]
nms<-nms[grepl('state',nms)==FALSE]
dat<-data[,names(data) %in% nms]
dat<-dat[,!(names(dat)%in%c('her.durng.rv','her.dumx.rv','sst.stalt','her.szpe.rv','herjuv.durng.rv'))]


#GET PROPORTION OF EACH COLUMN THAT IS MISSING
dat<-subset(dat,year>=1970 & year<=2006)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
dm<-data.frame(var=names(dat),
               pmis=apply(dat,2,pMiss))
dm<-dm[order(dm$pmis),]
#apply(dat,1,pMiss)
dm2<-subset(dm,pmis<3)
dat2<-dat[,names(dat) %in% as.character(dm2$var)]
dat2<-na.omit(dat2)

#EXAMINE DISTRIBUTION OF INDIVIDUAL VARIABLES STANDARDIZED TO Z
#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
b<-data.frame(cbind(subset(dat2,select='year'),apply(dat2[,2:dim(dat2)[2]],2,f)))
b<-b %>% gather(var, y, her.spnug:sst.bbay)
histogram( ~ y | var,data=b)

#image(md.pattern(q))
aggr_plot<-aggr(q, col=c('dodgerblue4','firebrick3'), numbers=TRUE, sortVars=TRUE, labels=names(q), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))



################################################

################################################
#setwd(figsdir)
#pdf('imputed_scatterplot.pdf',height=10,width=10)

cr<-cor(dat2,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
names(cr.t)[3]<-'r.real'


#########################################################

##############  IMPUTATION SIMULATION REPEATED 100 TIMES

#########################################################
impsim<-function(){

smp<-c(.25,.2,.15,.1)
lll<-list()
zzz<-list()
for(g in 1:length(smp)){

#GETS RANDOM SAMPLE (DATAFRAME) ACCORDING TO SPECIFIED MISSINNESS [G]
nms<-names(dat2)[2:dim(dat2)[2]]
l<-list()
for(i in 1:length(nms)){
n<-dim(dat2)[1]-ceiling(38*smp[g])#SAMPLE SIZE TO SAMPLE
    d<-subset(dat2,select=c('year',paste(nms[i])))
d<-d[sample(nrow(d),n,replace=FALSE),]
l[[i]]<-d
}
q<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), l)#COMBINE

#FUNCTION GETS TRUE VALUES FROM ORIGINAL DATASET THAT WERE REMOVED AND THEN IMPUTED; NON IMPUTED VALUES ARE NA'S
nms<-names(dat2)[2:dim(dat2)[2]]
l<-list()
for(i in 1:length(nms)){
d<-subset(dat2,select=c('year',paste(nms[i])))#TRUE DATA
d2<-subset(q,select=c('year',paste(nms[i])))#SAMPLED DATA
names(d2)[2]<-'y'
dout<-merge(d,d2,by=c('year'),all=TRUE)
dout<-dout[is.na(dout[,3])==TRUE,]
l[[i]]<-dout[,1:2]
}
qreal<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), l)#COMBINE

#DIFFERENT IMPUTATION METHODS TO TEST
#mth<-c('norm.predict','cart','pmm','sample','rf')
mth<-c('cart','pmm','sample','rf')
ll<-list()
zz<-list()
for(h in 1:length(mth)){

print(smp[g])
print(mth[h])
tempdata<-mice(q[,-1],m=5,maxit=50,meth=mth[h],seed=500,print=FALSE)
#GET COMPLETE IMPUTED DATASET
impdat<-complete(tempdata,1)
impdat$year<-q$year

#GET AVERAGE OF ALL 5 IMPUTATIONS
id1<-as.matrix(complete(tempdata,1))
id2<-as.matrix(complete(tempdata,2))
id3<-as.matrix(complete(tempdata,3))
id4<-as.matrix(complete(tempdata,4))
id5<-as.matrix(complete(tempdata,5))
impdatav<-(id1+id2+id3+id4+id5)/5
impdatav<-data.frame(impdatav)
impdatav$year<-q$year

#AVERAGE CORRELATION BETWEEN ALL IMPUTED
cr<-cor(impdat,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
crimp = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
crimp=crimp[crimp$Var1 != crimp$Var2,]
crimp=crimp[paste(crimp$Var1,crimp$Var2,sep="_") %in% combinations,]
names(crimp)[3]<-'r.imp1'

#AVERAGE CORRELATION BETWEEN AVERAGE OF ALL IMPUTED
cr<-cor(impdatav,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
crimpav = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
crimpav=crimpav[crimpav$Var1 != crimpav$Var2,]
crimpav=crimpav[paste(crimpav$Var1,crimpav$Var2,sep="_") %in% combinations,]
names(crimpav)[3]<-'r.imp5'

aa<-merge(cr.t,crimp,by=c('Var1','Var2'),all=TRUE)
aa<-merge(aa,crimpav,by=c('Var1','Var2'),all=TRUE)

nms<-names(dat2)[2:dim(qreal)[2]]
l<-list()
for(i in 1:length(nms)){
d<-subset(qreal,select=c('year',paste(nms[i])))#TRUE VALUES THAT WERE LATER IMPUTED
d2<-subset(impdat,select=c('year',paste(nms[i])))#IMPUTED VALUES
d3<-subset(dat2,select=c('year',paste(nms[i])))#ORIGINAL DATA
d4<-subset(impdatav,select=c('year',paste(nms[i])))#AVERAGE IMPUTED
names(d2)[2]<-'imputed'
dout<-merge(d,d2,by=c('year'),all=TRUE)
dout<-na.omit(dout)
r<-round(cor(dout[,2],dout[,3]),digits=2)#CORRELAION BETWEEN ONLY 'TRUE' VALUES AND THOSE THAT WERE IMPUTED

#plot(dout[,2],dout[,3],pch=16,las=1)
#legend('topright',legend=paste('r=',r),bty='n',title=tempdata$meth[1])
#legend('topleft',legend=paste('nmiss=',smp[g]),bty='n')
#abline(lm(dout[,3]~dout[,2]))

names(d4)[2]<-'imputed'
dout<-merge(d,d4,by=c('year'),all=TRUE)
dout<-na.omit(dout)
r2<-round(cor(dout[,2],dout[,3]),digits=2)#CORRELATION BETWEEN ONLY 'TRUE' VALUES AND AVERAGE OF TOP 5 IMPUTATIONS

#points(dout[,2],dout[,3],pch=16,las=1,col='red3')
#legend('topright',legend=paste('r=',r),bty='n',title=tempdata$meth[1])
#abline(lm(dout[,3]~dout[,2]),col='red')

l[[i]]<-data.frame(var=nms[i],
                   r=r,
                   r2=r2,
                   sddata=round(sd(d3[,2],na.rm=TRUE),digits=2),
                   meth=tempdata$meth[1],
                   smp=smp[g])
}
ll[[h]]<-data.frame(do.call('rbind',l))
zz[[h]]<-aa
}
lll[[g]]<-data.frame(do.call('rbind',ll))#RETURNS R AMONG IMPUTED
zzz[[g]]<-data.frame(do.call('rbind',zz))#RETURNS R AMONG ALL
}

list(data.frame(do.call('rbind',lll)),data.frame(do.call('rbind',zzz)))
}


#R.IMP1 IS AVERAGE OF ALL INDIVIDUAL IMPUTATIONS
#R.IMP5 IS AVERAGE OF AVERAGE OF ALL IMPUTED
#R IS CORRELATION BETWEEN TRUE VALUES AND IMPUTED
#R2 IS CORRELATION BETWEEN TRUE VALUES AND AVERAGE OF TOP 5 IMPUTATIONS
library(pbapply)
dout00<-pbreplicate(2,impsim(),simplify=FALSE)
dout00<-pbreplicate(100,impsim(),simplify=FALSE)
#dout00<-replicate(2,impsim(),simplify=FALSE)

#GETS CORRELATION BETWEEN IMPUTED AND REAL
f<-function(d){ return(data.frame(d[1]))}
dout<-ldply(dout00,.fun=f)
#GETS AVERAGE CORRELATION BETWEEN 1) REAL DATA; IMPUTED DATA(1) WITH REAL DATA; 3) IMPUTED DATA(5) WITH REAL DATA
f<-function(d){ return(data.frame(d[2]))}
dout2<-ldply(dout00,.fun=f)


dout<-data.frame(do.call('rbind',dout00))
#save(dout,file='imputation.simdata.RData')
setwd(figsdir)
load('imputation.simdata.RData')
dout00<-dout

f<-function(d){
    return(data.frame(r=mean(d$r,na.rm=TRUE),
                      r.sd=sd(d$r,na.rm=TRUE),
                      r2=mean(d$r2,na.rm=TRUE),
                      r2.sd=sd(d$r2,na.rm=TRUE)))
}
ot<-ddply(dout,.(meth,smp),.fun=f)
ot$id<-seq(1,dim(ot)[1],1)
dm<-data.frame(meth=sort(unique(ot$meth)),
    cl=c('firebrick3','dodgerblue3','green3','orange','magenta'))
ot<-merge(ot,dm,by=c('meth'),all=TRUE)


setwd(figsdir)
pdf('imputation_simulation_compare.pdf',height=10,width=7)
par(mfrow=c(2,1),mar=c(4,4,1,1))

plot(ot$r,ot$r2,pch=16,xlim=c(-.1,.4),ylim=c(-.1,.4),cex=ot$smp*10,las=1,col=alpha(as.character(ot$cl),.7),xlab='Average r between imputed and real values (top imp)',ylab='Average r between imputed and real values (top 5 imps)')
points(ot$r,ot$r2,pch=1,cex=ot$smp*10,col=alpha('black',.5),lwd=.02)
abline(a=0,b=1,lty=2)
legend('topleft',legend=dm$meth,col=as.character(dm$cl),pch=15,bty='n')
f<-function(d){
    lines(c(d$r,d$r),c(d$r2-(1*d$r2.sd),d$r2+(1*d$r2.sd)))
    lines(c(d$r+(1*d$r.sd),d$r-(1*d$r.sd)),c(d$r2,d$r2))
}
#z<-dlply(ot,.(id),.fun=f)

d<-data.frame(sddata=sort(unique(dout$sddata)),
              r=tapply(dout$r,dout$sddata,mean),
              r2=tapply(dout$r2,dout$sddata,mean))

d$lsddata<-log10(d$sddata)
plot(d$lsddata,d$r2,las=1,xlab='Average log10 variance of timeseries',ylab='Average r between imputed and real data',pch=16,xlim=c(-1.75,2),ylim=c(-.05,.5))
mod<-gam(r2~s(lsddata,k=10),data=subset(d,lsddata<3),gamma=1.4)
pdat<-data.frame(lsddata=seq(min(d$lsddata),max(d$lsddata),length.out=100))
pdat$p<-predict(mod,newdata=pdat)
lines(pdat$lsddata,pdat$p)
s<-summary(mod)
legend('topright',legend=paste('r2=',round(s$r.sq,digits=2)),bty='n')
dev.off()



