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
library(rcompanion)
library(moments)
library(DataCombine)

datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-'N://cluster_2017//scratch//spera//data//finaldat_v2'
datadirout<-'N://cluster_2017//scratch//spera//data//dataoutput'
figsdir<-'C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures'
figsdir<-'C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures'



############################################

###READ IN RAW DATA TO DETERMINE MISSINGNESS
setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE &
         grepl('\\.her',nms)==FALSE &
         grepl('\\.tplus',nms)==FALSE &
         grepl('\\.tmin',nms)==FALSE &
         grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-data[,!(names(data) %in% c('her.ssb'))]

dfs<-df %>% gather(var,y,-year)
dfs$y<-ifelse(dfs$y==Inf,NA,dfs$y)
dfs$y<-ifelse(is.na(dfs$y)==FALSE,1,NA)


#GET PROPORTION OF EACH COLUMN THAT IS MISSING
#dat<-subset(df,year>=1971 & year<=2005)
dat<-subset(df,year>=1975 & year<=2005)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
dm<-data.frame(var=names(dat),
               pmis=apply(dat,2,pMiss))
dm<-dm[order(dm$pmis),]
#apply(dat,1,pMiss)
dm2<-subset(dm,pmis<=40)



#####################################################
###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','her.dvm2.rv','herlrv.dep.bof')
setwd(datadir)
#load("impdatav.std.RData")
load("impdatav.RData")
data<-impdatav
nms<-nms[grepl('\\.se',nms)==FALSE &
         grepl('\\.her',nms)==FALSE &
         grepl('\\.tplus',nms)==FALSE &
         grepl('\\.tmin',nms)==FALSE &
         grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c(omt))]

#RETAINS INSTANCES WHERE MISSINGNESS<30%
df<-df[,names(df) %in% as.character(dm2$var)]

df<-df[,!(names(df)%in%c("her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv",'her.eggprod','her.tbio','herjuv.len.rv','her.georng','herjuv.dumx.rv',"her.cf.cat","her.dep.rv","her.dvm.rv","her.dvm2.rv","her.eggprod","her.georng","her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv","herlrv.dep.bof","her.ajrat.rv", "her.ajrat.rv.t1", "her.ajrat.rv.t2", "her.ajrat.rv.t3","her.cf.rv", "her.cf.rv.t1","her.cf.rv.t2", "her.cf.rv.t3","her.fmass.rv","her.fmass.rv.t1","her.fmass.rv.t2", "her.fmass.rv.t3","her.metai.rv","her.metai.rv.t1","her.metai.rv.t2", "her.metai.rv.t3","her.prod","her.prod.t1", "her.prod.t2","her.prod.t3","her.szpe.rv", "her.szpe.rv.t1", "her.szpe.rv.t2","her.szpe.rv.t3","her.waa","her.waa.t1","her.waa.t2","her.waa.t3", "herjuv.fmass.rv", "herjuv.fmass.rv.t1","herjuv.fmass.rv.t2","herjuv.fmass.rv.t3","herjuv.prey.bof","her.age.cat","her.age5.cat","her.agediv.cat", "her.agepe.cat.x","her.agepe.cat.y","herjuv.metai.rv",'sst.t12','sst.fmax'))]

rspnms<-c('her.ssbc','her.state','her.rec1','her.expr','her.land.sprich','her.land','sst.lur','zp.mn.sabs','nut.state','sst.stalt','t50.ct')
for(i in 1:length(rspnms)){
nm<-rspnms[i]
df<-slide(df,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
df<-slide(df,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
df<-slide(df,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
df<-slide(df,Var=nm,slideBy=4,NewVar=gsub(' ','',paste(nm,'.t4')))
df<-slide(df,Var=nm,slideBy=5,NewVar=gsub(' ','',paste(nm,'.t5')))
df<-slide(df,Var=nm,slideBy=6,NewVar=gsub(' ','',paste(nm,'.t6')))
}


#DATA TO USE FOR NETWORK
mat<-subset(df,year>=1975 & year<=2005)

f<-function(x){transformTukey(x,plotit=FALSE)}
mat<-data.frame(cbind(subset(mat,select='year'),apply(mat[,2:dim(mat)[2]],2,f)))

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
mats<-mat %>% gather(var,y,-her.ssbc)
mats<-mat %>% gather(var,y,-her.ssbc.t3)
mats<-mat %>% gather(var,y,-her.ssbc.t4)


cf<-function(d){
return(data.frame(ra=abs(cor(d$her.ssbc,d$y,use='pairwise.complete.obs')),
                r=cor(d$her.ssbc,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]


mats<-mat %>% gather(var,y,-herlrv.surv)
cf<-function(d){
return(data.frame(ra=abs(cor(d$herlrv.surv,d$y,use='pairwise.complete.obs')),
                  r=cor(d$herlrv.surv,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]

mats<-mat %>% gather(var,y,-herlrv.len)
cf<-function(d){
return(data.frame(ra=abs(cor(d$herlrv.len,d$y,use='pairwise.complete.obs')),
                  r=cor(d$herlrv.len,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]

md<-stepAIC(lm(herlrv.len~herlrv.mn.bof+had.pi+asl+herlrv.mn.bof+her.len.rv +zp.mn.sabs + sst.lur,data=mat))
AIC(md1,md2)
s<-summary(md)

rmat<-cor(mat)

md<-lm(herlrv.surv~herlrv.mn.bof+herlrv.len,data=mat)
plot(mat$herlrv.len,log10(mat$herlrv.surv+2))
md<-lm(her.ssbc.t4~her.expr.t3,data=mat)
plot(mat$her.expr.t3,mat$her.ssbc.t4,pch=15)

mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
her.ssbc.t4 ~ had.pi +her.rec1.t1 + her.expr.t3
her.rec1.t3 ~ herlrv.len + herlrv.mn.bof
herlrv.len ~ had.pi + zp.mn.sabs + her.len.rv + asl + nut.state
herlrv.mn.bof ~ had.pi
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
#herlrv.len ~~ had.pi
'
smod<-sem(mod,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE)

summary(smod,standardized=TRUE,fit=TRUE,rsquare=TRUE)
AIC(smod)

smod<-sem(mod,sample.cov=d,sample.nobs=dim(a)[1],std.lv=FALSE,warn=FALSE)
semPaths(smod)
semPaths(smod,layout='tree2')
library(semPlot)
mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
her.ssbc.t6 ~ her.expr.t5 + her.prod.t3 + her.rec1.t3 + temp.state.t6 + strt.ct.t.1
her.rec1.t3 ~ herlrv.len + herlrv.mn.bof + bfish.state.t.1
herlrv.mn.bof ~ had.pi + strt.t.1 + wnd.ct.t.2 + wind.amp + temp.state
herlrv.len ~ had.pi + strt.t.1 + wnd.ct.t.2 + wind.amp + temp.state
herlrv.surv ~ had.pi + sst.fsd.t.1 + ss.dist.t.2 + nut.state + strt + temp.state
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
#herlrv.len ~~ had.pi
'


modindices(smod)#high values indicate variable shoudl not be in model; 40 would be large
vartable(smod)#variances
inspect(smod, 'theta')#TEST IF RESIDUALS ON OBSERVED VARS ARE CORRELATED - MIGHT BE NEGATIVELY SO!

pv<-fitMeasures(smod,'pvalue')
#TABLE OF PARAMETERS
df<-data.frame(parameterEstimates(mod))
#df<-data.frame(standardizedSolution(mod))

################################################

################################################






semPaths(smod)
semPaths(smod,layout='tree2')
mod<-sem(smod,data=a,std.lv=FALSE,warn=FALSE)
#modindices(mod)#high values indicate variable shoudl not be in model; 40 would be large
#vartable(mod)#variances
#inspect(mod, 'theta')#TEST IF RESIDUALS ON OBSERVED VARS ARE CORRELATED - MIGHT BE NEGATIVELY SO!

#plot(subset(d,select=c('pnit','psstpi','pwndp','pcldfi','pirrad','pmm','pzoop')),pch=16)

pv<-fitMeasures(mod,'pvalue')
#TABLE OF PARAMETERS
df<-data.frame(parameterEstimates(mod))
#df<-data.frame(standardizedSolution(mod))
df<-subset(df,select=c('lhs','rhs','est','se','pvalue'))
names(df)<-c('response','predictor','estimate','std.error','pvalue')
df$modpv<-pv




a<-subset(data,year>1965 & is.na(her.ssb)==FALSE,select=c('her.ssb','her.expr','sst.state'))
d<-cov(a,use='na.or.complete',method='pearson')
dd<-cor(a,use='na.or.complete',method='pearson')
mn<-colMeans(a,na.rm=TRUE)

library(semPlot)

mod<-'
#REGRESSIONS - OBSERVATION EQUATIONS
her.ssb ~ her.expr + sst.state
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
her.ssb ~~ her.ssb
'

mod<-sem(mod,sample.cov=d,sample.nobs=dim(a)[1],std.lv=FALSE,warn=FALSE)
semPaths(mod)

a<-subset(data,herlrv.mn.bof<200)
plot(a$herlrv.len,log10(a$herlrv.mn.bof),pch=15)
cor(a$herlrv.len,a$herlrv.mn.bof,use='pairwise.complete.obs')
plot(data$herlrv.len,data$herlrv.mn.bof,pch=15)
plot(data$herlrv.len,log10(log10(data$herlrv.mn.bof)),pch=15)
cor(data$herlrv.len,data$herlrv.mn.bof,use='pairwise.complete.obs')
