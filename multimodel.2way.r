library(parallel)
library(MuMIn)
library(moments)
library(DataCombine)
library(rcompanion)
library(pbapply)
library(MASS)
library(parallel)
library(MuMIn)
library(permute)
library(rcompanion)
library(tidyr)
library(moments)
#install.packages('ggplot2',repos = c("http://rstudio.org/_packages",
#                           "http://cran.rstudio.com"))
#source('http://bioconductor.org/biocLite.R')
#biocLite('Rgraphviz')

datadirout<-'N://cluster_2017//scratch//spera//data//dataoutput'
datadir1<-'N://cluster_2017//scratch//spera//data//stagingdat'
datadir<-"N:/cluster_2017/scratch/spera/data/finaldat_v2"
figsdir<-"C://Users//sailfish//Documents//aalldocuments//literature//research//active//SPERA//Figures"
figsdir<-"C://Users//copepod//Documents//aalldocuments//literature//research//active//SPERA//Figures"


############################################

###READ IN RAW DATA TO DETERMINE MISSINGNESS
setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-data[,!(names(data) %in% c('her.ssb'))]
#df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv'))]

#CREATE LAGGED VERSIONS OF VARIABLES, AS IN IMPUTED
rspnms<-c('her.ssbc','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.prod','her.metai.rv','herlrv.mn.bof')
#rspnms<-c('her.ssbc','her.len.rv','herlrv.len','her.rec1','her.prod','her.state')
for(i in 1:length(rspnms)){
nm<-rspnms[i]
df<-slide(df,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
df<-slide(df,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
df<-slide(df,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
}

clmnms<-c("asl","gs.dist","nao", "nut.state","ss.ct","ss.dist","sst.fsd","sst.tmax","strt","strt.ct","t50.ct","temp.state","wind.amp",  "wind.state","wind.tmax", "wnd.ct","wnd.strs323fall",'nut.ct','nut.state')
for(i in 1:length(clmnms)){
nm<-clmnms[i]
df<-slide(df,Var=nm,slideBy=-1,NewVar=gsub(' ','',paste(nm,'.t-1')))
df<-slide(df,Var=nm,slideBy=-2,NewVar=gsub(' ','',paste(nm,'.t-2')))
}

rspnms<-c("her.expr","her.land","her.land.pct1","her.land.spdiv","her.land.sprich",'bfish.state')
for(i in 1:length(rspnms)){
nm<-rspnms[i]
df<-slide(df,Var=nm,slideBy=-1,NewVar=gsub(' ','',paste(nm,'.t-1')))
df<-slide(df,Var=nm,slideBy=-2,NewVar=gsub(' ','',paste(nm,'.t-2')))
}

names(df)<-gsub('t-1','t\\.1',names(df))
names(df)<-gsub('t-2','t\\.2',names(df))
names(df)<-gsub('t-3','t\\.3',names(df))
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
dm2<-subset(dm,pmis<=35)



#####################################################
#####################################################
###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','her.dvm2.rv','herlrv.dep.bof')
setwd(datadir)
#load("impdatav.std.RData")
load("impdatav.RData")
#load('SPERA_andata_new.RData')
data<-impdatav
#names(data)<-gsub('t\\.1','t-1',names(data))
#names(data)<-gsub('t\\.2','t-2',names(data))
#names(data)<-gsub('t\\.3','t-3',names(data))
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c(omt))]

#RETAINS INSTANCES WHERE MISSINGNESS<30%
df<-df[,names(df) %in% as.character(dm2$var)]

df<-df[,!(names(df)%in%c("her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv",'her.eggprod','her.tbio','herjuv.len.rv','her.georng','herjuv.dumx.rv',"her.cf.cat","her.dep.rv","her.dvm.rv","her.dvm2.rv","her.eggprod","her.georng","her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv","herlrv.dep.bof","her.ajrat.rv", "her.ajrat.rv.t1", "her.ajrat.rv.t2", "her.ajrat.rv.t3","her.cf.rv", "her.cf.rv.t1","her.cf.rv.t2", "her.cf.rv.t3","her.fmass.rv","her.fmass.rv.t1","her.fmass.rv.t2", "her.fmass.rv.t3","her.metai.rv","her.metai.rv.t1","her.metai.rv.t2", "her.metai.rv.t3","her.prod","her.prod.t1", "her.prod.t2","her.prod.t3","her.szpe.rv", "her.szpe.rv.t1", "her.szpe.rv.t2","her.szpe.rv.t3","her.waa","her.waa.t1","her.waa.t2","her.waa.t3", "herjuv.fmass.rv", "herjuv.fmass.rv.t1","herjuv.fmass.rv.t2","herjuv.fmass.rv.t3","herjuv.prey.bof","her.age.cat","her.age5.cat","her.agediv.cat", "her.agepe.cat.x","her.agepe.cat.y","herjuv.metai.rv"))]

nms<-names(df)
rspnms<-nms[grepl('land',nms)==TRUE |
         grepl('bfish',nms)==TRUE |
         grepl('expr',nms)==TRUE]
rspnms<-rspnms[grepl('\\.se',rspnms)==FALSE]
for(i in 1:length(rspnms)){
nm<-rspnms[i]
df<-slide(df,Var=nm,slideBy=-1,NewVar=gsub(' ','',paste(nm,'.t.1')))
}


df<-df[,!(names(df) %in% c('sst.t12','sst.fmax','sst.fmin','sst.bbay','sst.grg','sst.hal','wind.fmin','wind.min','wind.max','her.land.spdiv','her.land.spdiv.t.1','her.land.spdiv.t.2','lrv.mn.bof','zp.rich.sabs','zp.div.sabs.state','ph.rich.sabs','lrv.rich.bof','phos','sil','sst.min','temp.state','wind.state','sst.amp','her.expr','her.land.pct1','her.land','her.land.sprich','bfish.state'))]


#NORMALIZE SO THAT EFFECTS ARE APPROX LINEAR
f<-function(x){transformTukey(x,plotit=FALSE)}
#df<-data.frame(cbind(subset(df,select='year'),apply(df[,!(names(df) %in% c('year'))],2,f)))

###########################################################

###########################################################
#clusterType<-if(length(find.package('snow',quiet=TRUE))) 'SOCK' else 'PSOCK'
#clust<-try(makeCluster(getOption('cl.cores',8),type=clusterType))
#clusterExport(clust,'df2')
nrep<-2000
df2<-df
df2<-subset(df2,year>=1975 & year<=2005)
sts<-c('her.rec1','her.state','herlrv.mn.bof','herlrv.len','herlrv.surv','year','her.len.rv')
df2<-df2[,!(names(df2) %in% sts)]


mavgfun<-function(){
#GETS SAMPLE OF VARIABLES TO EXAMINE
nms<-names(df2[,!(names(df2) %in% c('year','her.ssbc'))])
nms<-sample(nms,5,replace=FALSE)
d<-subset(df2,select=c('her.ssbc',nms))
options(na.action='na.fail')
#mod<-lm(her.ssbc~.,data=d)#ONE-WAY INTERRACTIONS
#mod<-lm(her.ssbc~.^n,data=d)#ALL INTERACTIONS
mod<-lm(her.ssbc~.*.,data=d)#2-WAY INTERRACTIONS
modl<-dredge(mod,extra=c('R^2','adjR^2'),rank='AICc',trace=FALSE)
ma<-model.avg(modl)
imp<-data.frame(ma$importance)
imp$var<-rownames(imp)
ttab<-data.frame(coefficients(ma))
ttab$var<-rownames(ttab)
ttab<-ttab[-1,]
out<-merge(ttab,imp,by=c('var'),all=TRUE)
out$resp<-'her.ssbc'
return(out)
}
system.time(dout0<-pbreplicate(nrep,mavgfun(),simplify=FALSE))
dssb<-data.frame(do.call('rbind',dout0))





df2<-df
df2<-subset(df2,year>=1975 & year<=2005)
sts<-c('her.rec1','her.ssbc','herlrv.mn.bof','herlrv.len','herlrv.surv','year','her.len.rv')
df2<-df2[,!(names(df2) %in% sts)]

mavgfun<-function(){
#GETS SAMPLE OF VARIABLES TO EXAMINE
nms<-names(df2[,!(names(df2) %in% c('year','her.state'))])
nms<-sample(nms,5,replace=FALSE)
d<-subset(df2,select=c('her.state',nms))
options(na.action='na.fail')
mod<-lm(her.state~.*.,data=d)#2-WAY INTERRACTIONS
#mod<-lm(her.state~.,data=d)
modl<-dredge(mod,extra=c('R^2','adjR^2'),rank='AICc',trace=FALSE)
ma<-model.avg(modl)
imp<-data.frame(ma$importance)
imp$var<-rownames(imp)
ttab<-data.frame(coefficients(ma))
ttab$var<-rownames(ttab)
ttab<-ttab[-1,]
out<-merge(ttab,imp,by=c('var'),all=TRUE)
out$resp<-'her.state'
return(out)
}

system.time(dout0<-pbreplicate(nrep,mavgfun(),simplify=FALSE))
dstate<-data.frame(do.call('rbind',dout0))




df2<-df
df2<-subset(df2,year>=1975 & year<=2005)
sts<-c('her.state','her.ssbc','herlrv.mn.bof','herlrv.len','herlrv.surv','year','her.len.rv')
df2<-df2[,!(names(df2) %in% sts)]


mavgfun<-function(){
#GETS SAMPLE OF VARIABLES TO EXAMINE
nms<-names(df2[,!(names(df2) %in% c('year','her.rec1'))])
nms<-sample(nms,5,replace=FALSE)
d<-subset(df2,select=c('her.rec1',nms))
options(na.action='na.fail')
mod<-lm(her.rec1~.*.,data=d)#2-WAY INTERRACTIONS
#mod<-lm(her.rec1~.,data=d)
modl<-dredge(mod,extra=c('R^2','adjR^2'),rank='AICc',trace=FALSE)
ma<-model.avg(modl)
imp<-data.frame(ma$importance)
imp$var<-rownames(imp)
ttab<-data.frame(coefficients(ma))
ttab$var<-rownames(ttab)
ttab<-ttab[-1,]
out<-merge(ttab,imp,by=c('var'),all=TRUE)
out$resp<-'her.rec1'
return(out)
}
system.time(dout0<-pbreplicate(nrep,mavgfun(),simplify=FALSE))
drec<-data.frame(do.call('rbind',dout0))




df2<-df
df2<-subset(df2,year>=1975 & year<=2005)
sts<-c('her.state','her.ssbc','her.rec1','herlrv.len','herlrv.surv','year','her.len.rv')
df2<-df2[,!(names(df2) %in% sts)]

mavgfun<-function(){
#GETS SAMPLE OF VARIABLES TO EXAMINE
nms<-names(df2[,!(names(df2) %in% c('year','herlrv.mn.bof'))])
nms<-sample(nms,5,replace=FALSE)
d<-subset(df2,select=c('herlrv.mn.bof',nms))
options(na.action='na.fail')
mod<-lm(herlrv.mn.bof~.*.,data=d)#2-WAY INTERRACTIONS
#mod<-lm(herlrv.mn.bof~.,data=d)
modl<-dredge(mod,extra=c('R^2','adjR^2'),rank='AICc',trace=FALSE)
ma<-model.avg(modl)
imp<-data.frame(ma$importance)
imp$var<-rownames(imp)
ttab<-data.frame(coefficients(ma))
ttab$var<-rownames(ttab)
ttab<-ttab[-1,]
out<-merge(ttab,imp,by=c('var'),all=TRUE)
out$resp<-'herlrv.mn.bof'
return(out)
}
system.time(dout0<-pbreplicate(nrep,mavgfun(),simplify=FALSE))
dlrv<-data.frame(do.call('rbind',dout0))





df2<-df
df2<-subset(df2,year>=1975 & year<=2005)
sts<-c('her.state','her.ssbc','her.rec1','herlrv.mn.bof','herlrv.surv','year','her.len.rv')
df2<-df2[,!(names(df2) %in% sts)]

mavgfun<-function(){
#GETS SAMPLE OF VARIABLES TO EXAMINE
nms<-names(df2[,!(names(df2) %in% c('year','herlrv.len'))])
nms<-sample(nms,5,replace=FALSE)
d<-subset(df2,select=c('herlrv.len',nms))
options(na.action='na.fail')
mod<-lm(herlrv.len~.*.,data=d)#2-WAY INTERRACTIONS
#mod<-lm(herlrv.len~.,data=d)
modl<-dredge(mod,extra=c('R^2','adjR^2'),rank='AICc',trace=FALSE)
ma<-model.avg(modl)
imp<-data.frame(ma$importance)
imp$var<-rownames(imp)
ttab<-data.frame(coefficients(ma))
ttab$var<-rownames(ttab)
ttab<-ttab[-1,]
out<-merge(ttab,imp,by=c('var'),all=TRUE)
out$resp<-'herlrv.len'
return(out)
}
system.time(dout0<-pbreplicate(nrep,mavgfun(),simplify=FALSE))
dlrvlen<-data.frame(do.call('rbind',dout0))








library(plyr)
mmdat2<-rbind.fill(dssb,drec,dlrv,dlrvlen,dstate)
mmdat2$id<-seq(1,dim(mmdat2)[1],1)

efun<-function(d){
aa<-data.frame(strsplit(d$var,':'))
d$var1<-as.character(aa[1,])
if(dim(aa)[1]>1){
d$var2<-as.character(aa[2,])
} else {
d$var2<-NA
}
return(d)
}
mmdat2<-ddply(mmdat2,.(id),.fun=efun,.progress='text')

setwd(datadir)
#save(mmdat2,file='mmdat2.RData')



a1<-subset(mmdat2,select=c('var','var1','resp','ma.importance','coefficients.ma.'))
names(a1)[1:2]<-c('var2w','var')

a2<-subset(mmdat2,select=c('var','var2','resp','ma.importance','coefficients.ma.'))
names(a2)[1:2]<-c('var2w','var')
mmdat3<-rbind(a1,a2)

ff<-function(d){
    return(data.frame(var=unique(d$var),
                      imp=mean(d$ma.importance,na.rm=TRUE)))
}
ot<-ddply(subset(mmdat3,is.na(var)==FALSE),.(resp,var),.fun=ff,.progress='text')

f<-function(d){
    return(d[order(d$imp),])
}
ot<-ddply(ot,.(resp),.fun=f)

subset(ot,resp=='her.ssbc')
subset(ot,resp=='her.rec1')
subset(ot,resp=='herlrv.mn.bof')
subset(ot,resp=='herlrv.len')




ff<-function(d){
    return(data.frame(var2w=unique(d$var2w),
                      imp=mean(d$ma.importance,na.rm=TRUE)))
}
ot<-ddply(subset(mmdat3,is.na(var2w)==FALSE),.(resp,var2w),.fun=ff,.progress='text')

f<-function(d){    return(d[order(d$imp),]) }
ot<-ddply(ot,.(resp),.fun=f)

subset(ot,resp=='her.ssbc')
subset(ot,resp=='her.rec1')
subset(ot,resp=='herlrv.mn.bof')
subset(ot,resp=='herlrv.len')


f<-function(d){
d<-subset(d,var!='her.len.rv')
nms<-unique(d$var)
nms<-nms[grepl('sst\\.',nms)==TRUE | grepl('temp\\.',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE &
         grepl('sst\\.sha',nms)==FALSE &
         grepl('\\.fsd',nms)==FALSE &
         grepl('\\.tmax',nms)==FALSE &
         grepl('\\.amp',nms)==FALSE &
         grepl('dur12',nms)==FALSE]
d1<-subset(d,var %in% nms)
o1<-data.frame(grp='Temperature',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('\\.ct',nms)==TRUE |
         grepl('\\.dur12',nms)==TRUE |
         grepl('\\.fsd',nms)==TRUE |
         grepl('\\.tmax',nms)==TRUE|
         grepl('\\.amp',nms)==TRUE]
d1<-subset(d,var %in% nms)
o2<-data.frame(grp='Phenology',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('chl\\.',nms)==TRUE |
         grepl('ph\\.',nms)==TRUE |
         grepl('phyt\\.',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE]
d1<-subset(d,var %in% nms)
o3<-data.frame(grp='Phytoplankton',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('had\\.pi',nms)==TRUE]
d1<-subset(d,var %in% nms)
o4<-data.frame(grp='Predation',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('jf',nms)==TRUE]
d1<-subset(d,var %in% nms)
o5<-data.frame(grp='Competition',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('land',nms)==TRUE | grepl('expr',nms)==TRUE| grepl('bfish',nms)==TRUE]
d1<-subset(d,var %in% nms)
o6<-data.frame(grp='Exploitation',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))


nms<-unique(d$var)
nms<-nms[grepl('asl',nms)==TRUE | grepl('wind',nms)==TRUE| grepl('wnd',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE &
         grepl('\\.tmax',nms)==FALSE]
d1<-subset(d,var %in% nms)
o7<-data.frame(grp='Currents',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('zp\\.',nms)==TRUE |
         grepl('lrv\\.rich',nms)==TRUE|
         grepl('lrv\\.mn',nms)==TRUE|
         grepl('preytot',nms)==TRUE|
         grepl('lrv\\.pe',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE]
d1<-subset(d,var %in% nms)
o8<-data.frame(grp='Zooplankton',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('strt',nms)==TRUE |
         grepl('sil',nms)==TRUE|
         grepl('phos',nms)==TRUE|
         grepl('nit',nms)==TRUE|
         grepl('ph\\.mn',nms)==TRUE|
         grepl('nut\\.',nms)==TRUE|
         grepl('phyt\\.mn',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE]
d1<-subset(d,var %in% nms)
o9<-data.frame(grp='Mixing/production',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

nms<-unique(d$var)
nms<-nms[grepl('nao',nms)==TRUE |
         grepl('\\.dist',nms)==TRUE]
nms<-nms[grepl('\\.ct',nms)==FALSE]
d1<-subset(d,var %in% nms)
o10<-data.frame(grp='Climate',
               imp=mean(d1$ma.importance),
               b=mean(d1$coefficients.ma.))

dout<-rbind(o1,o2,o3,o4,o5,o6,o7,o8,o9,o10)
return(dout)
}
idat<-ddply(mmdat2,.(resp),.fun=f)




library(fmsb)
library(scales)
rfun<-function(d){
rsp<-unique(d$resp)
d0<-spread(data=d,key=grp, value=imp)
d1<-data.frame(grp=d$grp,
               imp=.5)
d1<-spread(data=d1,key=grp, value=imp)
d2<-data.frame(grp=d$grp,
               imp=0)
d2<-spread(data=d2,key=grp, value=imp)
dd<-rbind.fill(d1,d2,d0)
dd<-subset(dd,select=c('Temperature','Mixing/production','Climate','Currents','Phenology','Predation','Competition','Zooplankton','Phytoplankton','Exploitation'))
radarchart(dd,maxmin=TRUE,pfcol=alpha('dodgerblue',.5),pcol='dodgerblue',title=rsp)
}
par(mfrow=c(2,3))
z<-dlply(subset(idat,select=c('grp','imp','resp')),.(resp),.fun=rfun)



