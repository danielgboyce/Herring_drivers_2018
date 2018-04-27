library(tidyr)
library(lavaan)
library(scales)
library(MASS)
library(plyr)
library(plotrix)
library(lubridate)
library(bnstruct)
library(rcompanion)
library(moments)
library(DataCombine)
library(semPlot)

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
dm2<-subset(dm,pmis<=35)



#####################################################
###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','her.dvm2.rv','herlrv.dep.bof')
setwd(datadir)
#load("impdatav.std.RData")
#data<-impdatav.std
load("impdatav.RData")
data<-impdatav
nms<-names(data)
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

#rspnms<-c('her.ssbc','her.state','her.rec1','her.expr','her.land.sprich','her.land','sst.lur','nut.state','sst.stalt','t50.ct','herlrv.surv')
rspnms<-names(df)
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

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
mat<-data.frame(cbind(subset(mat,select='year'),apply(mat[,!(names(mat) %in% c('year'))],2,f)))

f<-function(x){transformTukey(x,plotit=FALSE)}
mat2<-data.frame(cbind(subset(mat,select='year'),apply(mat[,2:dim(mat)[2]],2,f)))


#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
mats<-mat %>% gather(var,y,-her.ssbc)
mats<-mat %>% gather(var,y,-her.ssbc.t3)

mats<-mat %>% gather(var,y,-her.ssbc.t4)
cf<-function(d){
return(data.frame(ra=abs(cor(d$her.ssbc.t4,d$y,use='pairwise.complete.obs')),
          r=cor(d$her.ssbc.t4,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra,decreasing=FALSE),]


mats<-mat %>% gather(var,y,-herlrv.surv.t1)
cf<-function(d){
return(data.frame(ra=abs(cor(d$herlrv.surv.t1,d$y,use='pairwise.complete.obs')),
          r=cor(d$herlrv.surv.t1,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]

mats<-mat %>% gather(var,y,-herlrv.len.t1)
cf<-function(d){
return(data.frame(ra=abs(cor(d$herlrv.len.t1,d$y,use='pairwise.complete.obs')),
                  r=cor(d$herlrv.len.t1,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]


mats<-mat %>% gather(var,y,-herlrv.mn.bof.t1)
cf<-function(d){
return(data.frame(ra=abs(cor(d$herlrv.mn.bof.t1,d$y,use='pairwise.complete.obs')),
          r=cor(d$herlrv.mn.bof.t1,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]


mats<-mat %>% gather(var,y,-her.rec1.t4)
cf<-function(d){
return(data.frame(ra=abs(cor(d$her.rec1.t4,d$y,use='pairwise.complete.obs')),
          r=cor(d$her.rec1.t4,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]


md<-stepAIC(lm(herlrv.len~herlrv.mn.bof+had.pi+asl+herlrv.mn.bof+her.len.rv +zp.mn.sabs + sst.lur,data=mat))
AIC(md1,md2)
s<-summary(md)
a<-subset(mat,select=c('her.ssbc.t4','herlrv.len','herlrv.len.t1','herlrv.len.t2','herlrv.len.t3','herlrv.len.t4'))
round(cor(a,use='pairwise.complete.obs',method='spearman'),digits=2)


a<-subset(mat,select=c('her.ssbc.t4','had.pi','her.expr.t3','herlrv.len','sst.stalt.t4','sst.sha.t4','herlrv.surv','herlrv.mn.bof','her.rec1','her.rec1.t1','her.rec1.t2','her.rec1.t3'))
#a<-subset(mat,select=c('herlrv.len.t1','had.pi.t1','herlrv.mn.bof.t1','gs.dist.t1','nao.t1','wind.amp.t1','sst.amp.t1','sst.lur.t1'))
md<-stepAIC(lm(her.ssbc.t4 ~.,data=a))

md1<-lm(her.ssbc.t4 ~ her.expr.t3+ herlrv.len + herlrv.mn.bof + sst.stalt.t4 + had.pi + herlrv.surv,data=a)
#md2<-lm(her.ssbc.t4 ~ her.expr.t3+ herlrv.len + herlrv.mn.bof + sst.stalt.t4 + had.pi + herlrv.surv + her.rec1.t1,data=a)
md3<-lm(her.ssbc.t4 ~ her.expr.t3+ herlrv.len + herlrv.mn.bof + sst.stalt.t4 + had.pi + herlrv.surv + her.rec1.t2,data=a)
md4<-lm(her.ssbc.t4 ~ her.expr.t3+ herlrv.len + herlrv.mn.bof + sst.stalt.t4 + had.pi + herlrv.surv + her.rec1.t3,data=a)
AIC(md1,md2,md3,md4)

md<-stepAIC(lm(herlrv.len.t1 ~.*.,data=a))
md<-stepAIC(lm(her.ssbc.t4 ~ had.pi + her.expr.t3 + sst.amp + t50.ct.t4 + herlrv.len + sst.stalt.t4 + sst.sha.t4 + herlrv.surv + herlrv.mn.bof,data=mat))
md<-lm(her.ssbc.t4 ~ had.pi + her.expr.t3 + sst.amp + t50.ct.t4 + herlrv.len + herlrv.mn.bof + sst.stalt.t4+sst.sha.t4+herlrv.surv,data=mat)
summary(md)
summary(md2)
AIC(md,md2)

library(Hmisc)
rmat<-cor(mat)
rmat<-cov(mat)
rmat2<-cor(mat2)

#PREYTOT IS HIGHLY CORRELATED WITH LARVAL SERIES - MAY BE DUE TO COMMON FACTOR DRIVIGN ALL LARVAE
mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
her.ssbc.t5 ~ had.pi.t1 + her.expr.t4 + sst.amp.t1 + t50.ct.t5 + herlrv.len.t1 + sst.stalt.t5 + sst.sha.t5 + herlrv.surv.t1 + herlrv.mn.bof.t1
her.rec1.t3 ~ herlrv.len + herlrv.mn.bof
herlrv.len.t1 ~ had.pi.t1  + herlrv.mn.bof.t1 + gs.dist.t1
herlrv.mn.bof.t1 ~ had.pi.t1 + her.ssbc.t1 + wind.amp.t1 + sst.lur.t1 + herlrv.len.t1
herlrv.surv.t1 ~ had.pi + herlrv.mn.bof + herlrv.len + sst.amp + sst.fsd + t50.ct + sst.stalt + gs.dist + her.ssbc
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
#herlrv.len ~~ had.pi
'
smod<-sem(mod,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE)
smod2<-sem(mod,sample.cov=rmat2,sample.nobs=dim(mat2)[1],std.lv=FALSE,warn=FALSE)

pe<-parameterEstimates(smod)
ge<-as.matrix(smod,matrix.type='edgelist')

summary(smod,standardized=TRUE,fit=TRUE,rsquare=TRUE)
                 Estimate
    her.ssbc.t5       0.937
    her.rec1.t3       0.226
    herlrv.len.t1     0.440
    hrlrv.mn.bf.t1    0.655
    herlrv.surv.t1    0.724

summary(smod2,standardized=TRUE,fit=TRUE,rsquare=TRUE)
                   Estimate
    her.ssbc.t5       0.931
    her.rec1.t3       0.211
    herlrv.len.t1     0.419
    hrlrv.mn.bf.t1    0.053
    herlrv.surv.t1    0.655
AIC(smod,smod2)

#fit simple reg
fit<-smod
fit.params = fit %>% parameterestimates()
fit.betas = fit.params %>% dplyr::filter(op == "~") %>% `[`("est") %>% unlist

#adj mats
beta.mat = fit@Model@GLIST$beta
log.mat = (fit@Model@GLIST$beta != 0)
adj.mat = log.mat %>% as.integer %>% matrix(nrow = nrow(fit@Model@GLIST$beta), byrow = T)
rownames(adj.mat) = colnames(adj.mat) = fit@Data@ov.names[[1]]
betas = beta.mat[log.mat]
#to igraph
ga<-graph_from_adjacency_matrix(adj.mat)
plot.igraph(ga,edge.label = round(betas,digits=2))





##############################################

#SETS COLOURS FOR GROUPS
fit.params = smod %>% parameterestimates()
grps<-list(Resp=c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),Herring=c('herlrv.mn.bof','herlrv.len','her.ssbc.t1','her.ssbc'),Predation=c('had.pi','had.pi.t1'),Exploitation=c('her.expr.t4'),Environment=c('sst.amp.t1','t50.ct.t5','sst.stalt.t5','sst.sha.t5','gs.dist.t1','wind.amp.t1','sst.lur.t1','sst.amp','sst.fsd','t50.ct','sst.stalt','gs.dist'))

#GETS EDGE LABELS
fit.params$sig<-ifelse(fit.params$pvalue<0.05,'*','')
fit.params$sig<-ifelse(fit.params$pvalue<0.01,'**',fit.params$sig)
fit.params$sig<-ifelse(fit.params$pvalue<0.001,'***',fit.params$sig)
fit.params$sig<-ifelse(is.na(fit.params$pvalue)==TRUE,'',fit.params$sig)
fit.params$lab<-gsub(' ','',paste(round(fit.params$est,digits=2),fit.params$sig))

#WIDTH 0F EDGES
#fit.params$esz<-abs(log10(abs(fit.params$est)))*30
fit.params$esz<-rescale(log10(abs(fit.params$est)),newrange=c(.2,10))

#COLOURS FOR GROUPS
cls<-c(Resp='white',Herring='gray90',Predation='mediumpurple2',Exploitation='hotpink',Environment='skyblue')
#cls<-c(Resp='gray10',Herring='gray50',Predation='mediumpurple4',Exploitation='maroon4',Environment='dodgerblue4')

#NODE NAMES
nms<-data.frame(nms1=smod@Data@ov.names[[1]])
nms$id<-seq(1,dim(nms)[1],1)
nms$var<-nms$nms1
nms$var<-gsub('\\.t1','',nms$var)
nms$var<-gsub('\\.t3','',nms$var)
nms$var<-gsub('\\.t4','',nms$var)
nms$var<-gsub('\\.t5','',nms$var)

nms$lag<-ifelse(grepl('\\.t1',nms$nms1)==TRUE,'\n(t+1)','nolag')
nms$lag<-ifelse(grepl('\\.t2',nms$nms1)==TRUE,'\n(t+2)',nms$lag)
nms$lag<-ifelse(grepl('\\.t3',nms$nms1)==TRUE,'\n(t+3)',nms$lag)
nms$lag<-ifelse(grepl('\\.t4',nms$nms1)==TRUE,'\n(t+4)',nms$lag)
nms$lag<-ifelse(grepl('\\.t5',nms$nms1)==TRUE,'\n(t+5)',nms$lag)

setwd(figsdir)
dm<-read.csv('dm.csv',header=TRUE)
dm<-subset(dm,select=c('var','lbl'))
nms<-merge(nms,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
nms<-nms[order(nms$id),]

nms$nodenms<-gsub('nolag','',paste(nms$lbl,nms$lag))
nms$nodenms<-gsub('Her ','',paste(nms$nodenms))
nms$nodenms<-capitalize(nms$nodenms)

nms$bw<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),5,.5)#WIDTH OF NODE OUTLINES
nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),'red3','gray20')#COLOR OF NODE OUTLINES
nms$bsz<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),1.45,1.1)#SIZE OF NODES

#SUBSET ONLY DIRECTED EFFECTS (NOT VARIANCES/COVARIANCES)
fp<-subset(fit.params,op=='~')

setwd(figsdir)
pdf('SEM_untransf.pdf',height=8,width=8)
pdf('SEM_untransf.cov.pdf',height=8,width=8)
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.5,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,edgeLabels=fp$lab,rainbowStart=1,esize=fp$esz,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,color=cls,nodeLabels=nms$nodenms,mar=c(1.5,1.5,1.5,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30')
dev.off()



















mats<-mat %>% gather(var,y,-her.rec1.t3)
cf<-function(d){
return(data.frame(ra=abs(cor(d$her.rec1.t3,d$y,use='pairwise.complete.obs')),
          r=cor(d$her.rec1.t3,d$y,use='pairwise.complete.obs')))
}
ot<-ddply(mats,.(var),.fun=cf)
ot<-ot[order(ot$ra),]




rmat<-cor(mat)

#PREYTOT IS HIGHLY CORRELATED WITH LARVAL SERIES - MAY BE DUE TO COMMON FACTOR DRIVIGN ALL LARVAE
mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
#her.ssbc.t6 ~ her.expr.t5 + sst.amp.t2 + her.rec1.t3 + herlrv.surv.t3
her.ssbc.t6 ~ her.expr.t5 + her.rec1.t3 + sst.amp.t2
her.rec1.t3 ~ herlrv.surv.t3 + sst.amp.t2 + sst.max.t2
herlrv.len.t2 ~ had.pi.t2  + herlrv.mn.bof.t2 + gs.dist.t2
herlrv.mn.bof.t2 ~ had.pi.t2 + her.ssbc.t2 + wind.amp.t2 + sst.lur.t2 + herlrv.len.t2
herlrv.surv.t3 ~ had.pi.t2 + herlrv.mn.bof.t2 + herlrv.len.t2 + sst.amp.t2 + sst.fsd.t2 + t50.ct.t2 + sst.stalt.t2 + gs.dist.t2 + her.ssbc.t2
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
#herlrv.len ~~ had.pi
'

mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
her.ssbc.t6 ~ her.expr.t5 + her.rec1.t3 + sst.amp.t2
her.rec1.t3 ~ herlrv.surv.t3 + sst.amp.t2 + sst.max.t2
herlrv.len.t2 ~ had.pi.t2  + herlrv.mn.bof.t2 + gs.dist.t2
herlrv.mn.bof.t2 ~ had.pi.t2 + her.ssbc.t2 + wind.amp.t2 + sst.lur.t2 + herlrv.len.t2
herlrv.surv.t3 ~ had.pi.t2 + herlrv.mn.bof.t2 + herlrv.len.t2 + sst.amp.t2 + sst.fsd.t2 + t50.ct.t2 + sst.stalt.t2 + gs.dist.t2 + her.ssbc.t2
#ESTIMATE COVARIANCE BETWEEN RESIDUALS
#sst.stalt.t2 ~~ sst.amp.t2 + sst.lur.t2 + sst.amp.t2
#her.ssbc.t6 ~~ her.ssbc.t2
#herlrv.surv.t3~~0*herlrv.mn.bof.t2 + 0*had.pi.t2 + 0*wind.amp.t2
#t50.ct.t2~~0*sst.lur.t2
'

smod<-sem(mod,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE)
smod.rnd<-sem(mod,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE,fixed.x=FALSE)
#smod2<-sem(mod2,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE)

summary(smod.rnd,standardized=TRUE,fit=TRUE,rsquare=TRUE)
summary(smod,standardized=TRUE,fit=TRUE,rsquare=TRUE)
    her.ssbc.t6       0.877
    her.rec1.t3       0.659
    herlrv.len.t2     0.547
    hrlrv.mn.bf.t2    0.637
    herlrv.surv.t3    0.696

summary(smod2,standardized=TRUE,fit=TRUE,rsquare=TRUE)
AIC(smod,smod2)

a<-subset(mat,select=c('her.rec1','her.rec1.t1','her.rec1.t2','her.rec1.t3','her.rec1.t4','her.rec1.t5','her.ssbc.t6'))
round(cor(a),digits=2)

#fit simple reg
fit<-smod
fit.params = fit %>% parameterestimates()
fit.betas = fit.params %>% dplyr::filter(op == "~") %>% `[`("est") %>% unlist

#adj mats
beta.mat = fit@Model@GLIST$beta
log.mat = (fit@Model@GLIST$beta != 0)
adj.mat = log.mat %>% as.integer %>% matrix(nrow = nrow(fit@Model@GLIST$beta), byrow = T)
rownames(adj.mat) = colnames(adj.mat) = fit@Data@ov.names[[1]]
betas = beta.mat[log.mat]
#to igraph
ga<-graph_from_adjacency_matrix(adj.mat)
plot.igraph(ga,edge.label = round(betas,digits=2))





##############################################

#SETS COLOURS FOR GROUPS
fit.params = smod %>% parameterestimates()
grps<-list(Resp=c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),Herring=c('her.ssbc.t2'),Predation=c('had.pi.t2'),Exploitation=c('her.expr.t5'),Environment=c('sst.amp.t2','t50.ct.t6','sst.stalt.t6','sst.max.t2','gs.dist.t2','wind.amp.t2','sst.lur.t2','sst.fsd.t2','t50.ct.t2','sst.stalt.t2'))

#GETS EDGE LABELS
fit.params$sig<-ifelse(fit.params$pvalue<0.05,'*','')
fit.params$sig<-ifelse(fit.params$pvalue<0.01,'**',fit.params$sig)
fit.params$sig<-ifelse(fit.params$pvalue<0.001,'***',fit.params$sig)
fit.params$sig<-ifelse(is.na(fit.params$pvalue)==TRUE,'',fit.params$sig)
fit.params$lab<-gsub(' ','',paste(round(fit.params$est,digits=2),fit.params$sig))

#WIDTH 0F EDGES
#fit.params$esz<-abs(log10(abs(fit.params$est)))*30
fit.params$esz<-rescale(log10(abs(fit.params$est)),newrange=c(.2,12))

#COLOURS FOR GROUPS
cls<-c(Resp='white',Herring='gray90',Predation='mediumpurple2',Exploitation='hotpink',Environment='lightskyblue2')
#cls<-c(Resp='gray10',Herring='gray50',Predation='mediumpurple4',Exploitation='maroon4',Environment='dodgerblue4')

#NODE NAMES
nms<-data.frame(nms1=smod@Data@ov.names[[1]])
nms$id<-seq(1,dim(nms)[1],1)
nms$var<-nms$nms1
nms$var<-gsub('\\.t1','',nms$var)
nms$var<-gsub('\\.t2','',nms$var)
nms$var<-gsub('\\.t3','',nms$var)
nms$var<-gsub('\\.t4','',nms$var)
nms$var<-gsub('\\.t5','',nms$var)
nms$var<-gsub('\\.t6','',nms$var)

nms$lag<-ifelse(grepl('\\.t1',nms$nms1)==TRUE,'\n(t+1)','nolag')
nms$lag<-ifelse(grepl('\\.t2',nms$nms1)==TRUE,'\n(t+2)',nms$lag)
nms$lag<-ifelse(grepl('\\.t3',nms$nms1)==TRUE,'\n(t+3)',nms$lag)
nms$lag<-ifelse(grepl('\\.t4',nms$nms1)==TRUE,'\n(t+4)',nms$lag)
nms$lag<-ifelse(grepl('\\.t5',nms$nms1)==TRUE,'\n(t+5)',nms$lag)
nms$lag<-ifelse(grepl('\\.t6',nms$nms1)==TRUE,'\n(t+6)',nms$lag)

setwd(figsdir)
dm<-read.csv('dm.csv',header=TRUE)
dm<-subset(dm,select=c('var','lbl'))
nms<-merge(nms,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
nms<-nms[order(nms$id),]

nms$nodenms<-gsub('nolag','',paste(nms$lbl,nms$lag))
nms$nodenms<-gsub('Her ','',paste(nms$nodenms))
nms$nodenms<-capitalize(nms$nodenms)

nms$bw<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),7,.5)#WIDTH OF NODE OUTLINES
nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),'red3','gray100')
nms$bsz<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),1.6,1.1)

#SUBSET ONLY DIRECTED EFFECTS (NOT VARIANCES/COVARIANCES)
fp<-subset(fit.params,op=='~')

setwd(figsdir)
pdf('SEM_untransf_mod2.pdf',height=8,width=8)
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),'red3','gray100')
semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.6,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2.5,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,edgeLabels=fp$lab,rainbowStart=1,esize=fp$esz,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,color=cls,nodeLabels=nms$nodenms,mar=c(1.5,1.5,2,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30',residuals=FALSE,threshholds=TRUE)


nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),'black','gray100')
semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.6,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2.5,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,edgeLabels=fp$lab,rainbowStart=1,esize=fp$esz,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,color=cls,nodeLabels=nms$nodenms,mar=c(1.5,1.5,2,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30',residuals=FALSE,threshholds=TRUE)

nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t6','her.rec1.t3','herlrv.len.t2','herlrv.mn.bof.t2','herlrv.surv.t3'),'gold','gray100')
semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.6,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2.5,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,edgeLabels=fp$lab,rainbowStart=1,esize=fp$esz,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,color=cls,nodeLabels=nms$nodenms,mar=c(1.5,1.5,2,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30',residuals=FALSE,threshholds=TRUE)
dev.off()

#TRIED TO PLOT INDIVIDUAL MODELS OF THE SEM ON SEPARATE PAGES BUT DOESN'T WORK
grps2<-list(mod1=c('her.ssbc.t6','her.expr.t5','sst.amp.t2','her.rec1.t3'),
           mod2=c('her.rec1.t3','sst.fmax.t2','herlrv.surv.t3','sst.amp.t2'),
           mod3=c('herlrv.len.t2','gs.dist.t2','had.pi.t2','herlrv.mn.bof.t2'),
           mod4=c('herlrv.mn.bof.t2','her.ssbc.t2','had.pi.t2','sst.lur.t2','wind.amp.t2'),
           mod5=c('herlrv.surv.t3','t50.ct.t2','sst.stalt.t2','sst.fsd.t2','her.ssbc.t2','herlrv.mn.bof.t2','herlrv.len.t2','had.pi.t2'))

grps2<-list(mod1=c('her.ssbc.t6','her.expr.t5','her.rec1.t3','sst.amp.t2'),
mod2=c('her.rec1.t3','herlrv.surv.t3','sst.amp.t2','sst.max.t2'),
mod3=c('herlrv.len.t2','had.pi.t2','herlrv.mn.bof.t2','gs.dist.t2'),
mod4=c('herlrv.mn.bof.t2','had.pi.t2','her.ssbc.t2','wind.amp.t2','sst.lur.t2','herlrv.len.t2'),
mod5=c('nherlrv.surv.t3','had.pi.t2','herlrv.mn.bof.t2','herlrv.len.t2','sst.amp.t2','sst.fsd.t2','t50.ct.t2','sst.stalt.t2','gs.dist.t2','her.ssbc.t2'))

semPaths(smod,groups=grps2,layout='spring',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',edge.label.cex=.6,borders=TRUE,include=c('mod1'),allVars=FALSE,combineGroups=FALSE)


semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.6,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=TRUE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2.5,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,rainbowStart=1,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,mar=c(1.5,1.5,2,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30',residuals=FALSE,threshholds=TRUE,title=TRUE,include='Herring')

pdat0<-subset(mat,select=c(as.character(nms$nms1)))
pdat0<-data.frame(y=colMeans(pdat0))
pdat0$var<-rownames(pdat0)
pdat0<-spread(data=pdat0,key=var, value=y)
nm<-'sst.amp.t2'
pdat0<-pdat0[,!(names(pdat0) %in% c(nm))]
pdat<-data.frame(sst.amp.t2=seq(min(mat$sst.amp.t2,na.rm=TRUE),max(mat$sst.amp.t2,na.rm=TRUE),length.out=50))
pdat<-cbind(pdat,pdat0)

p<-predict(smod,newdata=pdat)



setwd(datadir)
load('SPERA_andata_new.RData')
andata<-data
nms<-names(andata)
nms<-nms[grepl('\\.se',nms)==FALSE &
         grepl('\\.her',nms)==FALSE &
         grepl('\\.tplus',nms)==FALSE &
         grepl('\\.tmin',nms)==FALSE &
         grepl('sp\\.n',nms)==FALSE]
andata<-andata[,names(andata) %in% c('year',nms)]
andata<-andata[,!(names(andata) %in% c('her.ssb'))]

andata<-andata[,!(names(andata)%in%c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','her.dvm2.rv','herlrv.dep.bof',"her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv",'her.eggprod','her.tbio','herjuv.len.rv','her.georng','herjuv.dumx.rv',"her.cf.cat","her.dep.rv","her.dvm.rv","her.dvm2.rv","her.eggprod","her.georng","her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv","herlrv.dep.bof","her.ajrat.rv", "her.ajrat.rv.t1", "her.ajrat.rv.t2", "her.ajrat.rv.t3","her.cf.rv", "her.cf.rv.t1","her.cf.rv.t2", "her.cf.rv.t3","her.fmass.rv","her.fmass.rv.t1","her.fmass.rv.t2", "her.fmass.rv.t3","her.metai.rv","her.metai.rv.t1","her.metai.rv.t2", "her.metai.rv.t3","her.prod","her.prod.t1", "her.prod.t2","her.prod.t3","her.szpe.rv", "her.szpe.rv.t1", "her.szpe.rv.t2","her.szpe.rv.t3","her.waa","her.waa.t1","her.waa.t2","her.waa.t3", "herjuv.fmass.rv", "herjuv.fmass.rv.t1","herjuv.fmass.rv.t2","herjuv.fmass.rv.t3","herjuv.prey.bof","her.age.cat","her.age5.cat","her.agediv.cat", "her.agepe.cat.x","her.agepe.cat.y","herjuv.metai.rv",'sst.t12','sst.fmax','zp.cpr.ct','wnd.ct'))]

andatas<-andata %>% gather(var,y,-year)
andatas$y<-ifelse(andatas$y==Inf,NA,andatas$y)
andatas$y<-ifelse(is.na(andatas$y)==FALSE,1,NA)

rspnms<-names(andata)
for(i in 1:length(rspnms)){
nm<-rspnms[i]
andata<-slide(andata,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
andata<-slide(andata,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
andata<-slide(andata,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
andata<-slide(andata,Var=nm,slideBy=4,NewVar=gsub(' ','',paste(nm,'.t4')))
andata<-slide(andata,Var=nm,slideBy=5,NewVar=gsub(' ','',paste(nm,'.t5')))
andata<-slide(andata,Var=nm,slideBy=6,NewVar=gsub(' ','',paste(nm,'.t6')))
}

andata<-subset(andata,year>=1965)

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
andata<-data.frame(cbind(subset(andata,select='year'),apply(andata[,!(names(andata) %in% c('year'))],2,f)))

p<-lavPredict(smod,newdata=andata)
p<-lavPredict(smod,newdata=andata,type='ov',method='regression',label=TRUE,se.fit=TRUE)

plot








########################################################################


############################################################################
###READ IN RAW DATA TO DETERMINE MISSINGNESS
setwd(datadir)
load('SPERA_andata_new.RData')
andata<-data
nms<-names(andata)
nms<-nms[grepl('\\.se',nms)==FALSE &
         grepl('\\.her',nms)==FALSE &
         grepl('\\.tplus',nms)==FALSE &
         grepl('\\.tmin',nms)==FALSE &
         grepl('sp\\.n',nms)==FALSE]
andata<-andata[,names(andata) %in% c('year',nms)]
andata<-andata[,!(names(andata) %in% c('her.ssb'))]


andata<-andata[,!(names(andata)%in%c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','her.dvm2.rv','herlrv.dep.bof',"her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv",'her.eggprod','her.tbio','herjuv.len.rv','her.georng','herjuv.dumx.rv',"her.cf.cat","her.dep.rv","her.dvm.rv","her.dvm2.rv","her.eggprod","her.georng","her.spcv","her.spnug","her.sprng","her.spvar","her.tbio","her.totno.rv","her.totwgt.rv","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.len.rv","herjuv.totno.rv","herjuv.totwgt.rv","herlrv.dep.bof","her.ajrat.rv", "her.ajrat.rv.t1", "her.ajrat.rv.t2", "her.ajrat.rv.t3","her.cf.rv", "her.cf.rv.t1","her.cf.rv.t2", "her.cf.rv.t3","her.fmass.rv","her.fmass.rv.t1","her.fmass.rv.t2", "her.fmass.rv.t3","her.metai.rv","her.metai.rv.t1","her.metai.rv.t2", "her.metai.rv.t3","her.prod","her.prod.t1", "her.prod.t2","her.prod.t3","her.szpe.rv", "her.szpe.rv.t1", "her.szpe.rv.t2","her.szpe.rv.t3","her.waa","her.waa.t1","her.waa.t2","her.waa.t3", "herjuv.fmass.rv", "herjuv.fmass.rv.t1","herjuv.fmass.rv.t2","herjuv.fmass.rv.t3","herjuv.prey.bof","her.age.cat","her.age5.cat","her.agediv.cat", "her.agepe.cat.x","her.agepe.cat.y","herjuv.metai.rv",'sst.t12','sst.fmax','zp.cpr.ct','wnd.ct'))]

andatas<-andata %>% gather(var,y,-year)
andatas$y<-ifelse(andatas$y==Inf,NA,andatas$y)
andatas$y<-ifelse(is.na(andatas$y)==FALSE,1,NA)


rspnms<-names(andata)
for(i in 1:length(rspnms)){
nm<-rspnms[i]
andata<-slide(andata,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
andata<-slide(andata,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
andata<-slide(andata,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
andata<-slide(andata,Var=nm,slideBy=4,NewVar=gsub(' ','',paste(nm,'.t4')))
andata<-slide(andata,Var=nm,slideBy=5,NewVar=gsub(' ','',paste(nm,'.t5')))
#andata<-slide(andata,Var=nm,slideBy=6,NewVar=gsub(' ','',paste(nm,'.t6')))
}


#DATA TO USE FOR NETWORK
#mat<-subset(andata,year>=1975 & year<=2005)
mat<-subset(andata,year>=1965)

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
mat<-data.frame(cbind(subset(mat,select='year'),apply(mat[,!(names(mat) %in% c('year'))],2,f)))

f<-function(x){transformTukey(x,plotit=FALSE)}
#mat2<-data.frame(cbind(subset(mat,select='year'),apply(mat[,2:dim(mat)[2]],2,f)))


#GET PROPORTION OF EACH COLUMN THAT IS MISSING
#dat<-subset(df,year>=1971 & year<=2005)
#dat<-subset(andata,year>=1975 & year<=2005)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
dmm<-data.frame(var=names(mat),
               pmis=apply(mat,2,pMiss))
dmm<-dmm[order(dmm$pmis),]
#apply(dat,1,pMiss)
dmm2<-subset(dmm,pmis<=75)

mat<-mat[,names(mat) %in% dmm2$var]


library(Hmisc)
rmat<-cor(mat,use='pairwise.complete.obs')

#PREYTOT IS HIGHLY CORRELATED WITH LARVAL SERIES - MAY BE DUE TO COMMON FACTOR DRIVIGN ALL LARVAE
mod<-'
#MEASUREMENT MODELS - OBSERVATION EQUATIONS
#rec1 =~ her.rec1
##REGRESSION EQUATIONS
her.ssbc.t5 ~ had.pi.t1 + her.expr.t4 + sst.amp.t1 + t50.ct.t5 + herlrv.len.t1 + sst.stalt.t5 + sst.sha.t5 + herlrv.surv.t1 + herlrv.mn.bof.t1
her.rec1.t3 ~ herlrv.len + herlrv.mn.bof
herlrv.len.t1 ~ had.pi.t1  + herlrv.mn.bof.t1 + gs.dist.t1
herlrv.mn.bof.t1 ~ had.pi.t1 + her.ssbc.t1 + wind.amp.t1 + sst.lur.t1 + herlrv.len.t1
herlrv.surv.t1 ~ had.pi + herlrv.mn.bof + herlrv.len + sst.amp + sst.fsd + t50.ct + sst.stalt + gs.dist + her.ssbc
#ESTIMATE VARIANCES OF SINGLE OBSERVED MODELS
#herlrv.len ~~ had.pi
'
smod<-sem(mod,sample.cov=rmat,sample.nobs=dim(mat)[1],std.lv=FALSE,warn=FALSE)



pe<-parameterEstimates(smod)
ge<-as.matrix(smod,matrix.type='edgelist')

summary(smod,standardized=TRUE,fit=TRUE,rsquare=TRUE)
                 Estimate
    her.ssbc.t5       0.937
    her.rec1.t3       0.226
    herlrv.len.t1     0.440
    hrlrv.mn.bf.t1    0.655
    herlrv.surv.t1    0.724

summary(smod2,standardized=TRUE,fit=TRUE,rsquare=TRUE)
                   Estimate
    her.ssbc.t5       0.931
    her.rec1.t3       0.211
    herlrv.len.t1     0.419
    hrlrv.mn.bf.t1    0.053
    herlrv.surv.t1    0.655
AIC(smod,smod2)

#fit simple reg
fit<-smod
fit.params = fit %>% parameterestimates()
fit.betas = fit.params %>% dplyr::filter(op == "~") %>% `[`("est") %>% unlist

#adj mats
beta.mat = fit@Model@GLIST$beta
log.mat = (fit@Model@GLIST$beta != 0)
adj.mat = log.mat %>% as.integer %>% matrix(nrow = nrow(fit@Model@GLIST$beta), byrow = T)
rownames(adj.mat) = colnames(adj.mat) = fit@Data@ov.names[[1]]
betas = beta.mat[log.mat]
#to igraph
ga<-graph_from_adjacency_matrix(adj.mat)
plot.igraph(ga,edge.label = round(betas,digits=2))





##############################################

#SETS COLOURS FOR GROUPS
fit.params = smod %>% parameterestimates()
grps<-list(Resp=c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),Herring=c('herlrv.mn.bof','herlrv.len','her.ssbc.t1','her.ssbc'),Predation=c('had.pi','had.pi.t1'),Exploitation=c('her.expr.t4'),Environment=c('sst.amp.t1','t50.ct.t5','sst.stalt.t5','sst.sha.t5','gs.dist.t1','wind.amp.t1','sst.lur.t1','sst.amp','sst.fsd','t50.ct','sst.stalt','gs.dist'))

#GETS EDGE LABELS
fit.params$sig<-ifelse(fit.params$pvalue<0.05,'*','')
fit.params$sig<-ifelse(fit.params$pvalue<0.01,'**',fit.params$sig)
fit.params$sig<-ifelse(fit.params$pvalue<0.001,'***',fit.params$sig)
fit.params$sig<-ifelse(is.na(fit.params$pvalue)==TRUE,'',fit.params$sig)
fit.params$lab<-gsub(' ','',paste(round(fit.params$est,digits=2),fit.params$sig))

#WIDTH 0F EDGES
#fit.params$esz<-abs(log10(abs(fit.params$est)))*30
fit.params$esz<-rescale(log10(abs(fit.params$est)),newrange=c(.2,10))

#COLOURS FOR GROUPS
cls<-c(Resp='white',Herring='gray90',Predation='mediumpurple2',Exploitation='hotpink',Environment='skyblue')
#cls<-c(Resp='gray10',Herring='gray50',Predation='mediumpurple4',Exploitation='maroon4',Environment='dodgerblue4')

#NODE NAMES
nms<-data.frame(nms1=smod@Data@ov.names[[1]])
nms$id<-seq(1,dim(nms)[1],1)
nms$var<-nms$nms1
nms$var<-gsub('\\.t1','',nms$var)
nms$var<-gsub('\\.t3','',nms$var)
nms$var<-gsub('\\.t4','',nms$var)
nms$var<-gsub('\\.t5','',nms$var)

nms$lag<-ifelse(grepl('\\.t1',nms$nms1)==TRUE,'\n(t+1)','nolag')
nms$lag<-ifelse(grepl('\\.t2',nms$nms1)==TRUE,'\n(t+2)',nms$lag)
nms$lag<-ifelse(grepl('\\.t3',nms$nms1)==TRUE,'\n(t+3)',nms$lag)
nms$lag<-ifelse(grepl('\\.t4',nms$nms1)==TRUE,'\n(t+4)',nms$lag)
nms$lag<-ifelse(grepl('\\.t5',nms$nms1)==TRUE,'\n(t+5)',nms$lag)

setwd(figsdir)
dm<-read.csv('dm.csv',header=TRUE)
dm<-subset(dm,select=c('var','lbl'))
nms<-merge(nms,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
nms<-nms[order(nms$id),]

nms$nodenms<-gsub('nolag','',paste(nms$lbl,nms$lag))
nms$nodenms<-gsub('Her ','',paste(nms$nodenms))
nms$nodenms<-capitalize(nms$nodenms)

nms$bw<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),5,.5)#WIDTH OF NODE OUTLINES
nms$bcl<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),'red3','gray20')#COLOR OF NODE OUTLINES
nms$bsz<-ifelse(nms$nms1 %in% c('her.ssbc.t5','her.rec1.t3','herlrv.len.t1','herlrv.mn.bof.t1','herlrv.surv.t1'),1.45,1.1)#SIZE OF NODES

#SUBSET ONLY DIRECTED EFFECTS (NOT VARIANCES/COVARIANCES)
fp<-subset(fit.params,op=='~')

setwd(figsdir)
pdf('SEM_untransf.pdf',height=8,width=8)
par(oma=c(0,0,0,0),mar=c(0,0,0,0))
semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='OpenMx',what='std',rotation=2,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.5,cardinal=TRUE,equalizeManifests=TRUE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,border.color=nms$bcl,border.width=nms$bw,fade=TRUE,asize=2,curvePivotShape=.8,groups=grps,posCol=c("forestgreen","firebrick3"),legend=FALSE,edgeLabels=fp$lab,rainbowStart=1,esize=fp$esz,edge.label.bg=FALSE,edge.label.position=.5,vTrans=255,color=cls,nodeLabels=nms$nodenms,mar=c(1.5,1.5,1.5,1.5),node.width=nms$bsz,node.height=nms$bsz,label.color='black',label.cex=.55,label.scale=FALSE,edge.label.color='gray30')
dev.off()



semPaths(smod,curvePivot=TRUE,intercepts=FALSE,sizeMan=8,layout='spring',style='lisrel',what='std',rotation=4,nCharNodes=0,nCharEdges=0,shapeMan='ellipse',reorder=TRUE,pastel=FALSE,edge.label.cex=.5,cardinal=TRUE,equalizeManifests=FALSE,intAtSide=TRUE,exoCov=FALSE,exoVar=FALSE,curve=TRUE,curveAdjacent='->',borders=TRUE,mar=c(1,1,1,1),fade=FALSE)











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
