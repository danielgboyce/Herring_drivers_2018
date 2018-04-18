#library(Rgraphviz)
#library(bnstruct)
library(rcompanion)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(maptools)
library(akima)
library(segmented)
library(splines)
library(strucchange)
library(data.table)
library(reshape2)
library(gplots)
library(cluster)
library(vegan)
library(ggplot2)
library(forecast)
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
library(rgdal)
library(car)
library(factoextra)
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



setwd('C:/Users/copepod/Documents/aalldocuments/literature/research/active/ESS_trophic_control/data')
setwd('C:/Users/sailfish/Documents/aalldocuments/literature/research/active/ESS_trophic_control/data')
plg<-readShapePoly('polygons_ecnasap.shp')#COMMAND TO READ BACK IN
plg<-subset(plg,region=='NS')


mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N://data/shapefiles/naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position="right",
                         plot.title = element_text(size=16)))



setwd(datadirout)
load('pdat.int.RData')#ANNUAL AVERAGES, NO SES
load('pdat.fac.RData')#INTERPOLATED 1980-2005
load('pdat.lin.RData')#MODEL PREDICTED
load('mpdat.RData');#ALL MODEL OUTPUT
load('bpdat.RData')#CONTAINS ALL BREAKPOINTS, IF PRESENT
load('data.RData')#ANNUAL ESTIMATES AND SES AS LONG FORM
load('data2.RData')#ANNUAL ESTIMATES AND SES AS STACKED

#ST ANDREWS SERIES RUNS 1988 ONWARD SO REMOVE
nms<-names(data)
nms<-nms[grepl('\\.sabs',nms)==FALSE]
data<-data[,names(data) %in% nms]
pdat.lin<-pdat.lin[,names(pdat.lin) %in% nms]
pdat.int<-pdat.int[,names(pdat.int) %in% nms]
pdat.fac<-pdat.fac[,names(pdat.fac) %in% nms]
data2<-subset(data2,var %in% nms)
bpdat<-subset(bpdat,var %in% nms)
mpdat<-subset(mpdat,var %in% nms)
mpdat0<-mpdat





#####################################################

#####################################################
setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv','her.durng.rv','her.dumx.rv','herjuv.durng.rv','herjuv.dumxrv','herjuv.dvm2.rv','herjuv.dur2.rv','herlrv.sp.r','herjuv.sp.n','her.sp.n','herlrv.sp.n','her.dur2.rv','herjuv.sprng','herjuv.spcv','herjuv.spvar','herjuv.spnug','herlrv.sprng','herlrv.spcv','herlrv.spvar','herlrv.spnug','herlrv.georng','herlrv.spnug','had.dep.rv','had.len.rv','her.cf.rv','herlrv.dep.bof'))]

#df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv','her.durng.rv','her.dumx.rv','herjuv.durng.rv','herjuv.dumxrv','herjuv.dvm2.rv','herjuv.dur2.rv','herlrv.sp.r','herjuv.dep.rv','herjuv.sp.n','her.sp.n','herlrv.sp.n','her.szdiv.rv','her.dur2.rv','herjuv.sprng','herjuv.spcv','herjuv.spvar','herjuv.spnug','herlrv.sprng','herlrv.spcv','herlrv.spvar','herlrg.spnug','herlrv.georng','herlrv.spnug','had.dep.rv','had.len.rv','her.cf.rv'))]

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
#f<-function(x){transformTukey(x,plotit=FALSE)}
#df<-data.frame(cbind(subset(df,select='year'),apply(df[,2:dim(df)[2]],2,f)))

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
df<-data.frame(cbind(subset(df,select='year'),apply(df[,2:dim(df)[2]],2,f)))

dfs<-df %>% gather(var,y,-year)
dfs<-subset(dfs,y!=Inf & is.na(y)==FALSE)

dum<-data.frame(var=sort(unique(dfs$var)))
dum$ct<-ifelse(grepl('lrv\\.',dum$var)==TRUE,'Larvae','Climate')
dum$ct<-ifelse(grepl('herlrv\\.',dum$var)==TRUE,'Herring larvae',dum$ct)
dum$ct<-ifelse(grepl('her\\.',dum$var)==TRUE,'Herring',dum$ct)
dum$ct<-ifelse(grepl('herjuv\\.',dum$var)==TRUE,'Juvenile herring',dum$ct)
dum$ct<-ifelse(grepl('ph\\.',dum$var)==TRUE,'Phytoplankton',dum$ct)
dum$ct<-ifelse(grepl('chl\\.',dum$var)==TRUE,'Phytoplankton',dum$ct)
dum$ct<-ifelse(grepl('phyt\\.',dum$var)==TRUE,'Phytoplankton',dum$ct)
dum$ct<-ifelse(grepl('zp\\.',dum$var)==TRUE,'Zooplankton',dum$ct)
dum$ct<-ifelse(grepl('land',dum$var)==TRUE,'Anthropogenic',dum$ct)
dum$ct<-ifelse(grepl('expr',dum$var)==TRUE,'Anthropogenic',dum$ct)
dum$ct<-ifelse(grepl('bfish\\.',dum$var)==TRUE,'Anthropogenic',dum$ct)
dum$ct<-ifelse(grepl('had\\.',dum$var)==TRUE,'Predation',dum$ct)
dum$ct<-ifelse(grepl('cop\\.',dum$var)==TRUE,'Zooplankton',dum$ct)

dm<-unique(subset(dfs,select=c('var')))
setwd(figsdir)
#write.csv(dm,'dm.csv')
dm<-read.csv('dm.csv',header=TRUE)

#ADDS LABELS
dm$ct2<-ifelse(dm$var %in% c('sst.mn','sst.amp','sst.bbay','sst.dur12','sst.fmax','sst.fmin','sst.fsd','sst.grg','sst.hal','sst.lur','sst.max','sst.min','sst.prin','sst.sha','sst.stalt','sst.tmax','t50','t50.ct','temp.state'),'Temperature',NA)
dm$ct2<-ifelse(dm$var %in% c('wind.state','wind.amp','wind.fmax','wind.fmin','wind.max','wind.min','wind.tmax','wnd.ct','wnd.mn','wnd.pctg10fall','wnd.strs323fall'),'Wind',dm$ct2)
dm$ct2<-ifelse(dm$var %in% c('ph.div.cpr','ph.div.sabs','ph.div.state','ph.mn.cpr','ph.mn.sabs','ph.mn.state','ph.pe.cpr','ph.pe.sabs','ph.pe.state','ph.rich.cpr','ph.rich.sabs','ph.to.sabs','chl.ct','phyt.mn.insitu','phyt.mn.me','phyt.mn.mo','phyt.mn.sw'),'Phytoplankton',dm$ct2)
dm$ct2<-ifelse(dm$var %in% c('asl','gs.dist','nao','ss.ct','ss.dist','strt','strt.ct'),'Mixing',dm$ct2)
dm$ct2<-ifelse(dm$var %in% c('nit','nut.ct','nut.state','phos','sil'),'Mixing',dm$ct2)
dm$ct2<-ifelse(grepl('had\\.',dm$var)==TRUE,'Predation',dm$ct2)
dm$ct2<-ifelse(grepl('her\\.',dm$var)==TRUE,'Herring',dm$ct2)
dm$ct2<-ifelse(grepl('herlrv\\.',dm$var)==TRUE,'Herring larvae',dm$ct2)
dm$ct2<-ifelse(grepl('herjuv\\.',dm$var)==TRUE,'Herring juvenile',dm$ct2)
dm$ct2<-ifelse(dm$var %in% c('bfish.state','her.land','her.land.pct1','her.land.spdiv','her.land.sprich','her.expr'),'Anthropogenic',dm$ct2)
dm$ct2<-ifelse(dm$var %in% c('zp.cop.cpr','zp.cpr.ct','zp.div.cpr','zp.div.cpr.state','zp.div.sabs','zp.div.sabs.state','zp.mn.cpr','zp.mn.cpr.ct','zp.mn.sabs','zp.pe.cpr','cp.cop.cpr.ct','zp.pe.sabs','zp.rich.cpr','zp.rich.sabs','zp.to.sabs','lrv.rich.bof','lrv.pe.bof','lrv.mn.bof','her.preytot.bof','her.jf.bof'),'Zooplankton-larvae',dm$ct2)


#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-lm(y~year,data=d)
    s<-summary(mod)
y1<-subset(d,year==min(d$year))$y
y2<-subset(d,year==max(d$year))$y
return(data.frame(b=s$coef[2,1],
                pv=s$coef[2,4],
                year1=min(d$year),
                year2=max(d$year),
                y1=y1,
                y2=y2))
}
bdat<-ddply(dfs,.(var),.fun=f,.progress='text')
bdat<-merge(bdat,dum,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat$span<-bdat$year2-bdat$year1
bdat<-subset(bdat,span>20)

#COLOR SCHEME
bdat$cb<-cut(bdat$b,breaks=seq(-.12,.12,length.out=21))
rdat<-data.frame(x=seq(-max(abs(bdat$b)),max(abs(bdat$b)),length.out=10000))
rdat$cb<-cut(rdat$x,breaks=seq(-.12,.12,length.out=21))
rdat<-unique(subset(rdat,select=c('cb')))
#rdat$cl<-c(blue2red(20)
cl<-colorRampPalette(c(blue2red(19),'darkred'))
rdat$cl<-cl(20)
bdat<-merge(bdat,rdat,by=c('cb'),all.x=TRUE,all.y=FALSE)

bdat<-merge(bdat,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat<-bdat[order(bdat$ct2),]

aa<-unique(subset(bdat,select=c('ct','ct2')))
aa<-aa[order(aa$ct),]

dum2<-data.frame(ct2=c('Herring','Herring larvae','Herring juvenile','Temperature','Mixing','Wind','Anthropogenic','Predation','Zooplankton-larvae','Phytoplankton'),
                 id1=seq(length(unique(bdat$ct2)),1,-1))
bdat<-merge(bdat,dum2,by=c('ct2'),all.x=TRUE,all.y=FALSE)

dum2<-bdat[order(bdat$id1, bdat$year2, bdat$year1),]
dum2<-unique(subset(dum2,select=c('id1','span','var')))
dum2$id<-seq(1,dim(dum2)[1],1)

bdat<-merge(bdat,dum2,by=c('id1','span','var'),all.x=TRUE,all.y=FALSE)
ydat<-unique(subset(bdat,select=c('id','lbl')))
ydat$lbl<-gsub('Juvenile','Juv.',ydat$lbl)



ff2<-function(d2){
cx1<-1
cx2<-1.5
#lines(c(d2$year1,d2$year2),c(d2$id,d2$id),col=as.character(d2$cl),lwd=2)
rect(xleft=d2$year1,ybottom=d2$id-.15,xright=d2$year2,ytop=d2$id+.15,col=as.character(d2$cl),border=NA)

points(d2$year1,d2$id,pch=16,cex=cx2,col=ifelse(d2$pv<.05,'gray25','white'))
points(d2$year1,d2$id,pch=1,lwd=.1,cex=cx2)

points(d2$year2,d2$id,pch=16,cex=cx2,col=ifelse(d2$pv<.05,'gray25','white'))
points(d2$year2,d2$id,pch=1,lwd=.1,cex=cx2)
}
setwd(figsdir)
pdf('indices_timetrends.pdf',height=24,width=10)
par(mar=c(4,.5,.5,14))
plot(0,0,ylim=c(min(bdat$id)+4,max(bdat$id)-4),las=1,xlim=c(1910,2017),axes=FALSE)
abline(v=1984)
axis(1,at=seq(1920,2020,10))
abline(v=seq(1920,2010,10),col=alpha('black',.2),lty=2)
axis(4,at=ydat$id,labels=ydat$lbl,las=1)
zz<-dlply(bdat,.(id),.fun=ff2)
dev.off()




#ABSOLUTE CHANGES
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-lm(y~year,data=d)
    s<-summary(mod)
y1<-subset(d,year==min(d$year))$y
y2<-subset(d,year==max(d$year))$y
return(data.frame(b=s$coef[2,1],
                pv=s$coef[2,4],
                year1=min(d$year),
                year2=max(d$year),
                y1=y1,
                y2=y2))
}
bdat<-ddply(dfs,.(var),.fun=f,.progress='text')
bdat<-merge(bdat,dum,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat$span<-bdat$year2-bdat$year1
bdat$b<-abs(bdat$b)*bdat$span
bdat<-subset(bdat,span>20)

#COLOR SCHEME
bdat$cb<-cut(bdat$b,breaks=seq(0,3.1,length.out=21))
rdat<-data.frame(x=seq(0.001,max(bdat$b),length.out=10000))
rdat$cb<-cut(rdat$x,breaks=seq(0,3.1,length.out=21))
rdat<-unique(subset(rdat,select=c('cb')))
#rdat$cl<-c(blue2red(20)
cl<-colorRampPalette(c(blue2red(19),'darkred'))
rdat$cl<-cl(20)
bdat<-merge(bdat,rdat,by=c('cb'),all.x=TRUE,all.y=FALSE)

bdat<-merge(bdat,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat<-bdat[order(bdat$ct2),]

aa<-unique(subset(bdat,select=c('ct','ct2')))
aa<-aa[order(aa$ct),]

dum2<-data.frame(ct2=c('Herring','Herring larvae','Herring juvenile','Temperature','Mixing','Wind','Anthropogenic','Predation','Zooplankton-larvae','Phytoplankton'),
                 id1=seq(length(unique(bdat$ct2)),1,-1))
bdat<-merge(bdat,dum2,by=c('ct2'),all.x=TRUE,all.y=FALSE)

dum2<-bdat[order(bdat$id1, bdat$year2, bdat$year1),]
dum2<-unique(subset(dum2,select=c('id1','span','var')))
dum2$id<-seq(1,dim(dum2)[1],1)

bdat<-merge(bdat,dum2,by=c('id1','span','var'),all.x=TRUE,all.y=FALSE)
ydat<-unique(subset(bdat,select=c('id','lbl')))
ydat$lbl<-gsub('Juvenile','Juv.',ydat$lbl)

setwd(figsdir)
pdf('indices_timetrends_abs.pdf',height=24,width=10)
par(mar=c(4,.5,.5,14))
plot(0,0,ylim=c(min(bdat$id)+4,max(bdat$id)-4),las=1,xlim=c(1910,2017),axes=FALSE)
abline(v=1984)
axis(1,at=seq(1920,2020,10))
abline(v=seq(1920,2010,10),col=alpha('black',.2),lty=2)
axis(4,at=ydat$id,labels=ydat$lbl,las=1)
zz<-dlply(bdat,.(id),.fun=ff2)
dev.off()



#RANKED TIME TRENDS IN ORDER OF LARGEST-SMALLEST
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-lm(y~year,data=d)
    s<-summary(mod)
y1<-subset(d,year==min(d$year))$y
y2<-subset(d,year==max(d$year))$y
return(data.frame(b=s$coef[2,1],
                pv=s$coef[2,4],
                year1=min(d$year),
                year2=max(d$year),
                y1=y1,
                y2=y2))
}
bdat<-ddply(dfs,.(var),.fun=f,.progress='text')
bdat<-merge(bdat,dum,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat$span<-bdat$year2-bdat$year1
bdat$b<-abs(bdat$b)*bdat$span
bdat<-subset(bdat,span>20)

#COLOR SCHEME
bdat$cb<-cut(bdat$b,breaks=seq(0,3.1,length.out=21))
rdat<-data.frame(x=seq(0.001,max(bdat$b),length.out=10000))
rdat$cb<-cut(rdat$x,breaks=seq(0,3.1,length.out=21))
rdat<-unique(subset(rdat,select=c('cb')))
#rdat$cl<-c(blue2red(20)
cl<-colorRampPalette(c(blue2red(19),'darkred'))
rdat$cl<-cl(20)
bdat<-merge(bdat,rdat,by=c('cb'),all.x=TRUE,all.y=FALSE)

bdat<-merge(bdat,dm,by=c('var'),all.x=TRUE,all.y=FALSE)
bdat<-bdat[order(bdat$ct2),]

aa<-unique(subset(bdat,select=c('ct','ct2')))
aa<-aa[order(aa$ct),]

dum2<-data.frame(ct2=c('Herring','Herring larvae','Herring juvenile','Temperature','Mixing','Wind','Anthropogenic','Predation','Zooplankton-larvae','Phytoplankton'),
                 id1=seq(length(unique(bdat$ct2)),1,-1))
bdat<-merge(bdat,dum2,by=c('ct2'),all.x=TRUE,all.y=FALSE)

dum2<-bdat[order(bdat$b, bdat$year2),]
dum2<-unique(subset(dum2,select=c('id1','span','var')))
dum2$id<-seq(1,dim(dum2)[1],1)

bdat<-merge(bdat,dum2,by=c('id1','span','var'),all.x=TRUE,all.y=FALSE)
ydat<-unique(subset(bdat,select=c('id','lbl')))
ydat$lbl<-gsub('Juvenile','Juv.',ydat$lbl)

setwd(figsdir)
pdf('indices_timetrends_ranked.pdf',height=24,width=10)
par(mar=c(4,.5,.5,14))
plot(0,0,ylim=c(min(bdat$id)+4,max(bdat$id)-4),las=1,xlim=c(1910,2017),axes=FALSE)
abline(v=1984)
axis(1,at=seq(1920,2020,10))
abline(v=seq(1920,2010,10),col=alpha('black',.2),lty=2)
axis(4,at=ydat$id,labels=ydat$lbl,las=1)
zz<-dlply(bdat,.(id),.fun=ff2)
dev.off()







as.character(subset(aa,is.na(c1)==TRUE,select=c('var')))

mmdat<-bdat
mmdat$c1<-ifelse(mmdat$var %in% c('bfish.state.t.1','bfish.state.t.2','bfish.state','her.expr','her.expr.t.1','her.expr.t.2','her.land.pct1','her.land.pct1.t.1','her.land.pct1.t.2','her.land.spdiv','her.land.spdiv.pct1','her.land.spdiv.t.1','her.land.spdiv.t.2','her.land.sprich.pct1','her.land.sprich.t.1','her.land.sprich.t.2','her.land.sprich','her.land','her.land.t.1','her.land.t.2'),'Exploitation',NA)
mmdat$c1<-ifelse(mmdat$var %in% c('had.pi'),'Predation',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('her.jf.bof'),'Competition',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('sst.bbay','sst.had','sst.grg','sst.lur','sst.max','sst.min','sst.mn','sst.prin','sst.fmin','sst.stalt','sst.hal','temp.state','t50'),'Temperature',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('sil','phos','nit','nut.state'),'Nutrients',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('ph.mn.sabs','ph.mn.state','phyt.mn.insitu','ph.div.sabs','ph.div.state','ph.rich.sabs','ph.pe.sabs','ph.pe.state','ph.div.cpr','ph.mn.cpr','ph.pe.cpr','ph.pe.cpr','ph.rich.cpr','ph.to.sabs','phyt.mn.me','phyt.mn.mo','phyt.mn.sw','ph.div.cpr'),'Phytoplankton',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('lrv.pe.bof','zp.pe.sabs','zp.rich.sabs','zp.div.sabs','zp.div.sabs.state','lrv.rich.bof','her.preytot.bof','lrv.mn.bof','zp.mn.sabs','zp.cop.cpr','zp.div.cpr','zp.div.cpr.state','zp.mn.cpr','zp.pe.cpr','zp.rich.cpr','zp.to.sabs', "herjuv.prey.bof"),'Zooplankton',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('chl.ct','ss.ct','gs.ct','strt.ct','t50.ct','wnd.ct','sst.tmax','sst.fsd','sst.sha','sst.amp','sst.dur12','wind.amp','wind.tmax','cp.cop.cpr.ct','zp.cpr.ct','zp.mn.cpr.ct','nut.ct','sst.fmax'),'Phenology',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c('asl','gs.dist','ss.dist','wind.fmax','wind.fmin','wind.max','wind.min','wnd.mn','wnd.pctg10fall','wnd.strs323fall','wind.state','nao','strt'),'Transport',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c("her.age.cat", "her.age5.cat", "her.agediv.cat", "her.agepe.cat.x", "her.agepe.cat.y", "her.ajrat.rv", "her.cf.cat", "her.dep.rv", "her.dvm2.rv", "her.eggprod", "her.fmass.rv", "her.georng", "her.len.rv", "her.metai.rv", "her.prod", "her.rec1", "her.spcv", "her.spnug", "her.sprng", "her.spvar", "her.ssb", "her.ssbc", "her.state", "her.szpe.rv", "her.tbio", "her.totno.rv", "her.totwgt.rv", "her.waa"),'Adult herring',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c("herjuv.dumx.rv", "herjuv.dvm.rv", "herjuv.fmass.rv", "herjuv.len.rv", "herjuv.metai.rv", "herjuv.totno.rv", "herjuv.totwgt.rv"),'Juvenile herring',mmdat$c1)
mmdat$c1<-ifelse(mmdat$var %in% c("herlrv.dep.bof", "herlrv.len", "herlrv.mn.bof", "herlrv.surv"),'Larval herring',mmdat$c1)

mmdat<-merge(mmdat,dm,by=c('var'),all.x=TRUE,all.y=FALSE)



##################################################

# LOLIPOP FIGURES SHOWING TOTAL CHANGE IN INDICES OVER TIME
d<-mmdat
d$span<-d$year2-d$year1
d<-subset(d,span>=20)
d<-subset(d,!(var %in% c('her.dvm2.rv','herlrv.dep.bof','her.agepe.cat.x','her.agepe.cat.y','her.state','nut.state','temp.state','bfish.state','ph.div.state','ph.pe.state','wind.state','zp.div.cpr.state','zp.div.sabs.state','her.ssb')))
dmm<-data.frame(c1=c('Adult herring','Juvenile herring','Larval herring','Zooplankton','Phytoplankton','Phenology','Exploitation','Temperature','Nutrients','Transport','Competition','Predation'),
               ord=seq(length(unique(d$c1)),1,-1),
               cl1=brewer.pal(n=length(unique(d$c1)),'Paired'),
               cl=rainbow(n=length(unique(d$c1))))
d<-merge(d,dmm,by=c('c1'),all.x=TRUE,all.y=FALSE)
library(RColorBrewer)

d<-d[order(d$ord,d$b),]
d$id<-seq(1,dim(d)[1],1)
d2<-subset(d,pv<0.05)
d2$b<-ifelse(d2$b>0,d2$b+.004,d2$b-.004)

pl1<-ggplot(d, aes(x=id, y=b)) +
geom_point(size=4,alpha=0.7,aes(x=id,y=b),color=d$cl1) +
geom_segment(aes(x=id, xend=id, y=0, yend=b)) +
geom_point(data=d2,size=5,alpha=0.9,aes(x=id,y=b),shape='*',color='gray20') +
#  theme_light() +
  coord_flip() +
  theme(panel.grid.major.y = element_blank())+
    scale_x_continuous(expand=c(0,0),breaks=seq(1,dim(d)[1],1),labels=d$lbl)+
    expand_limits(x=c(0,dim(d)[1]+1),y=c(-.07,.11))+
    xlab('')+
    ylab('Slope')



d<-mmdat
d$span<-d$year2-d$year1
d<-subset(d,span>=20)
d<-subset(d,!(var %in% c('her.dvm2.rv','herlrv.dep.bof','her.agepe.cat.x','her.agepe.cat.y','her.state','nut.state','temp.state','bfish.state','ph.div.state','ph.pe.state','wind.state','zp.div.cpr.state','zp.div.sabs.state','her.ssb')))
d$b<-d$b*d$span
dmm<-data.frame(c1=c('Adult herring','Juvenile herring','Larval herring','Zooplankton','Phytoplankton','Phenology','Exploitation','Temperature','Nutrients','Transport','Competition','Predation'),
               ord=seq(length(unique(d$c1)),1,-1),
               cl1=brewer.pal(n=length(unique(d$c1)),'Paired'),
               cl=rainbow(n=length(unique(d$c1))))
d<-merge(d,dmm,by=c('c1'),all.x=TRUE,all.y=FALSE)

d<-d[order(d$ord,d$b),]
d$id<-seq(1,dim(d)[1],1)
d2<-subset(d,pv<0.05)
d2$b<-ifelse(d2$b>0,d2$b+.1,d2$b-.1)

pl2<-ggplot(d, aes(x=id, y=b)) +
geom_point(size=4,alpha=0.9,aes(x=id,y=b),color=d$cl1) +
geom_segment(aes(x=id, xend=id, y=0, yend=b)) +
geom_point(data=d2,size=5,alpha=0.9,aes(x=id,y=b),shape='*',color='gray20') +
  coord_flip() +
  theme(panel.grid.major.y = element_blank())+
  scale_x_continuous(expand=c(0,0),breaks=seq(1,dim(d)[1],1),labels=d$lbl)+
    expand_limits(x=c(0,dim(d)[1]+1),y=c(-3,3))+
    xlab('')+
    ylab('Total change')



d<-mmdat
d$span<-d$year2-d$year1
d<-subset(d,span>=20)
d<-subset(d,!(var %in% c('her.dvm2.rv','herlrv.dep.bof','her.agepe.cat.x','her.agepe.cat.y','her.state','nut.state','temp.state','bfish.state','ph.div.state','ph.pe.state','wind.state','zp.div.cpr.state','zp.div.sabs.state','her.ssb')))
d$b<-abs(d$b*d$span)
dmm<-data.frame(c1=c('Adult herring','Juvenile herring','Larval herring','Zooplankton','Phytoplankton','Phenology','Exploitation','Temperature','Nutrients','Transport','Competition','Predation'),
               ord=seq(length(unique(d$c1)),1,-1),
               cl1=brewer.pal(n=length(unique(d$c1)),'Paired'),
               cl=rainbow(n=length(unique(d$c1))))
d<-merge(d,dmm,by=c('c1'),all.x=TRUE,all.y=FALSE)

d<-d[order(d$ord,d$b),]
d$id<-seq(1,dim(d)[1],1)
d2<-subset(d,pv<0.05)
d2$b<-d2$b+.1

pl3<-ggplot(d, aes(x=id, y=b)) +
geom_point(size=4,alpha=0.9,aes(x=id,y=b),color=d$cl1) +
geom_segment(aes(x=id, xend=id, y=0, yend=b)) +
geom_point(data=d2,size=5,alpha=0.9,aes(x=id,y=b),shape='*',color='gray20') +
  coord_flip() +
  theme(panel.grid.major.y = element_blank())+
  scale_x_continuous(expand=c(0,0),breaks=seq(1,dim(d)[1],1),labels=d$lbl)+
    expand_limits(x=c(0,dim(d)[1]+1),y=c(0,3))+
    xlab('')+
    ylab('Total change')



d<-mmdat
d$span<-d$year2-d$year1
d<-subset(d,span>=20)
d<-subset(d,!(var %in% c('her.dvm2.rv','herlrv.dep.bof','her.agepe.cat.x','her.agepe.cat.y','her.state','nut.state','temp.state','bfish.state','ph.div.state','ph.pe.state','wind.state','zp.div.cpr.state','zp.div.sabs.state','her.age.cat','her.age5.cat','her.agediv.cat','her.ssb')))
d$b<-d$b*d$span
d<-subset(d,c1 %in% c('Adult herring','Juvenile herring','Larval herring'))
dmm<-data.frame(c1=c('Adult herring','Juvenile herring','Larval herring'),
               ord=seq(length(unique(d$c1)),1,-1),
               cl1=brewer.pal(n=length(unique(d$c1)),'Paired'),
               cl=rainbow(n=length(unique(d$c1))))
d<-merge(d,dmm,by=c('c1'),all.x=TRUE,all.y=FALSE)

d<-d[order(d$ord,d$b),]
d$id<-seq(1,dim(d)[1],1)
d2<-subset(d,pv<0.05)
d2$b<-ifelse(d2$b>0,d2$b+.1,d2$b-.1)

pl4<-ggplot(d, aes(x=id, y=b)) +
geom_point(size=4,alpha=0.9,aes(x=id,y=b),color=d$cl) +
geom_segment(aes(x=id, xend=id, y=0, yend=b)) +
geom_point(data=d2,size=5,alpha=0.9,aes(x=id,y=b),shape='*',color='gray20') +
  coord_flip() +
  theme(panel.grid.major.y = element_blank())+
  scale_x_continuous(expand=c(0,0),breaks=seq(1,dim(d)[1],1),labels=d$lbl)+
    expand_limits(x=c(0,dim(d)[1]+1),y=c(-3,3))+
    xlab('')+
    ylab('Total change')


setwd(figsdir)
pdf('indices_timetrends_lolipop.pdf',height=14,width=10)
par(mar=c(4,.5,.5,14))
pl1
pl2
pl3
grid.arrange(pl4,pl4,ncol=1)
dev.off()


5*1000
5000/60
2^30

plot(0,0,col='#00FFFF',pch=15)


#################################################################

#REDUNDANCY ANALYSIS OPERATES ON INTERPOLATED DATA - NEEDS NON-MISSING DATA
setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv'))]

#PREDICTORS
preds<-df[,!(names(df) %in% c('sst.t12','herjuv.spcv','herjuv.spnug','herjuv.spvar','herjuv.sprng','cp.cop.cpr.ct','zp.mn.cpr.ct','zp.cpr.ct',"her.cf.rv", "her.dep.rv","her.dumx.rv", "her.dur2.rv",  "her.durng.rv", "her.dvm.rv","her.dvm2.rv", "her.fmass.rv", "her.georng", "her.len.rv","her.metai.rv", "her.prod","her.rec1","her.spcv",  "her.spnug","her.sprng", "her.spvar", "her.ssb","her.szdiv.rv","her.szpe.rv","her.tbio","her.totno.rv","her.totwgt.rv","her.waa",'herjuv.totno.rv','herjuv.totwgt.rv','herjuv.metai','herjuv.dur2.rv','herjuv.len.rv','herlrv.len','herlrv.georng','herjuv.fmass.rv','her.cf.cat','herlrv.sprng','herlrv.spvar','herlrv.surv','herlrv.mn.bof','herjuv.durng.rv','herjuv.dvm2.rv','herlrv.spnug','herjuv.metai.rv','herjuv.dvm.rv','herlrv.dep.bof','herjuv.dumx.rv','herlrv.spcv','herjuv.dep.rv','zp.cop.cpr','her.ajrat.rv','her.eggprod'))]


#RESPONSE
resp<-df[,(names(df) %in% c('her.ssb','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.spnug','her.spcv','her.spvar','her.prod','her.metai.rv'))]
resp$her.spnug<-resp$her.spnug*-1
resp$her.spcv<-resp$her.spcv*-1
resp$her.spvar<-resp$her.spvar*-1


preds<-nms[grepl('her\\.',nms)==FALSE]
preds<-nms[grepl('\\.tplus',nms)==FALSE]
preds<-nms[grepl('\\.tmin',nms)==FALSE]
preds<-preds[!(preds %in% c('year','her.fmass.rv','her.metai.rv','her.sprng','her.spnug','her.spcv','her.spvar','her.waa','her.len.rv','her.szdiv.rv','her.szpe.rv','her.ssb'))]
rownames(resp)<-df$year
rownames(preds)<-df$year



#names(prds)<-gsub(' ','_',names(prds))
srda<-rda(resp~.,data=preds,scale=TRUE)
plot(srda,xlim=c(-1,2))





##########################################################

##########################################################

setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv'))]

#HERRING INDICES TO REMOVE
df<-df[,!(names(df) %in% c('sst.t12','herjuv.spcv','herjuv.spnug','herjuv.spvar','herjuv.sprng','cp.cop.cpr.ct','zp.mn.cpr.ct','zp.cpr.ct',"her.cf.rv", "her.dep.rv","her.dumx.rv", "her.dur2.rv",  "her.durng.rv", "her.dvm.rv","her.dvm2.rv", "her.fmass.rv", "her.georng", "her.len.rv","her.metai.rv", "her.prod","her.rec1","her.spcv",  "her.spnug","her.sprng", "her.spvar", "her.ssb","her.szdiv.rv","her.szpe.rv","her.tbio","her.totno.rv","her.totwgt.rv","her.waa",'herjuv.totno.rv','herjuv.totwgt.rv','herjuv.metai','herjuv.dur2.rv','herjuv.len.rv','herlrv.len','herlrv.georng','herjuv.fmass.rv','her.cf.cat','herlrv.sprng','herlrv.spvar','herlrv.surv','herlrv.mn.bof','herjuv.durng.rv','herjuv.dvm2.rv','herlrv.spnug','herjuv.metai.rv','herjuv.dvm.rv','herlrv.dep.bof','herjuv.dumx.rv','herlrv.spcv','herjuv.dep.rv','zp.cop.cpr','her.ajrat.rv','her.eggprod','had.len.rv','had.dep.rv'))]


#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)

#LOOK AT SKEWNESS OF VARIABLES
sk<-data.frame(var=sort(unique(dfs$var)),
         sk=tapply(dfs$y,dfs$var,function(x) skewness(x,na.rm=TRUE)))
sk<-sk[order(abs(sk$sk),decreasing=TRUE),]

#TRANSFORM VARIABLES TO OPTIMIZE NORMALITY
f<-function(d){
    d$y<-transformTukey(d$y,plotit=FALSE)
    return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)


#LOOK AT SKEWNESS OF VARIABLES
sk<-data.frame(var=sort(unique(dfs$var)),
         sk=tapply(dfs$y,dfs$var,function(x) skewness(x,na.rm=TRUE)))
sk<-sk[order(abs(sk$sk),decreasing=TRUE),]

#CONVERT BACK TO LONG FORM
df<-spread(data=dfs,key=var, value=y)

df2<-data.frame(t(df))
names(df2)<-df2[1,]
df2<-df2[-1,]

#histogram(~y | var,data=dfs)

f<-function(d){
    d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
    return(d)
}
pdats<-ddply(dfs,.(var),.fun=f)

pdatl<-spread(data=pdats,key=var, value=y)

#histogram(~y | var,data=pdats)
#xyplot(y~year|var,data=pdats,pch=15,type=c('p','spline'))
#xyplot(y~year|var,data=pdats,pch=15)

dff<-pdatl[,-1]
dff<-df[,-1]
#dff2<-pdats
dmat<-1-(cor(dff,use='pairwise.complete.obs',method='spearman'))
#dmat<-abs(cor(dff,use='pairwise.complete.obs',method='spearman'))
dst<-as.dist(dmat)


#########################################################
setwd(figsdir)
pdf('pressures_state_hcluster.pdf',height=14,width=10)
par(mar=c(4,4,4,4),oma=c(1,1,1,1))

K<-6
mds<-as.data.frame(cmdscale(dst))
colnames(mds)<-c('Dim.1','Dim.2')
p1<-ggscatter(mds,x='Dim.1',y='Dim.2',label=names(dff),size=2,repel=TRUE,alpha=.5)

#KMEANS CLUSTERING
#dcol<-data.frame(cl=c('gold3','firebrick3','forestgreen','dodgerblue3','purple','orange'),clust=seq(1,k,1))
dcol<-data.frame(cl=c('gold3','firebrick3','magenta3','green3','darkblue','lightblue'),clust=seq(1,K,1))

clust <- kmeans(mds, K)$cluster %>%
as.factor()
mds <- mds %>%
mutate(groups = clust)
# Plot and color by groups
p2<-ggscatter(mds, x = "Dim.1", y = "Dim.2",
          label = names(dff),
          color = "groups",
#          palette = "jco",
          palette =as.character(dcol$cl),
          size = 1,
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
grid.arrange(p1,p2,ncol=1)
#grid.arrange(p1,p2,ncol=2)



#HCLUST PLOTS
par(mfrow=c(2,1),mar=c(4,6,1,6))
hc<-hclust(dst,method='ward.D2')
sub_grp<-cutree(hc,k=K)

hc2<-as.dendrogram(hc)
plot(hc2)
rect.hclust(hc, k = K, border = 2:5)


#PLOT CLUST AS HORIZONTAL
agh<-as.hclust(hc)
labelColors<-c('gold3','firebrick3','magenta3','green3','darkblue','lightblue')
clusMember = cutree(agh, k)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
clusDendro = dendrapply(as.dendrogram(agh), colLab)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19))
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.01,pch=16,nodePar=nodePar)


#GETS NAMES OF VARIABLES IN EACH CLUSTER
cldat<-data.frame(var=names(dff),
               clust=clust)
cldat<-cldat[order(cldat$clust),]
cldat<-merge(cldat,dcol,by=c('clust'),all=FALSE)


#GETS AVERAGE CORREALTION, RATE OF CHANGE AND  PLOTS TS FOR EACH CLUSTER
par(mfrow=c(4,2))
for(i in 1:6){
print(i)
cl<-subset(cldat,clust==i)
    d<-subset(pdatl,select=as.character(cl$var))
    d2<-na.omit(subset(pdats,var %in% as.character(cl$var)))
    cr<-cor(d,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
print(mean(cr.t$Freq))
mod<-gamm(y~year,data=d2,random=list(var=~1))
modf<-gamm(y~as.factor(year),data=d2,random=list(var=~1))
mods<-gamm(y~s(year),data=d2,random=list(var=~1))
pdat<-data.frame(year=sort(unique(d2$year)))
pdatsm<-data.frame(year=seq(min(d2$year),max(d2$year),length.out=1000))
psm<-predict(mods$gam,newdata=pdatsm,se.fit=TRUE)
pdatsm$p<-psm$fit
pdatsm$se<-psm$se.fit
pdatsm$upr<-pdatsm$p+(1.96*pdatsm$se)
pdatsm$lwr<-pdatsm$p-(1.96*pdatsm$se)
p<-predict(modf$gam,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*p$se)
pdat$lwr<-pdat$p-(1.96*p$se)
plot(pdat$year,pdat$p,las=1,pch=15,xlim=c(1900,2018),ylim=c(-2,2),col=unique(as.character(cl$cl)),xaxt='n',xlab='Year',ylab='Mixed model average')
axis(1,seq(1900,2018,5))
abline(h=0,lty=2)
polygon(c(pdatsm$year,pdatsm$year[length(pdatsm$year):1]),c(pdatsm$upr,pdatsm$lwr[length(pdatsm$lwr):1]),col=alpha(as.character(unique(cl$cl)),.3),border=NA)
s<-summary(mod$gam)
r2<-round(s$r.sq,digits=2)
f2<-function(g){    lines(c(g$year,g$year),c(g$upr,g$lwr),col=alpha(unique(as.character(cl$cl)),.3))}
#zz<-dlply(pdat,.(year),.fun=f2)
legend('topright',c(paste('r2=',r2),paste('av r=',round(mean(cr.t$Freq),digits=2))),bty='n')
lines(pdatsm$year,pdatsm$p,col=unique(as.character(cl$cl)),lwd=2)
}
dev.off()









######################################################

######################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")
#library(gRain)
library(bnlearn)
#library(network)
#library(ggnet)
library(igraph)
#data(marks)
#bn.gs = bnlearn::gs(marks)


#########################################

###  FIRST LOOK AT MISSINGNESS OF DATA
setwd(datadir)
load('SPERA_andata_new.RData')
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
#df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv'))]


library(DataCombine)
library(rcompanion)
rspnms<-c('her.ssb','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.spnug','her.spcv','her.spvar','her.prod','her.metai.rv')
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

dfs<-df %>% gather(var,y,-year)
dfs$y<-ifelse(is.na(dfs$y)==FALSE,1,NA)
dfs<-subset(dfs,year>=1960)

p1<-ggplot()+
geom_tile(data=dfs,aes(x=year,y=as.factor(var),fill=y),col='gray')+
#scale_fill_distiller(palette='Spectral')+
scale_fill_gradient(low='blue',high='red',na.value='white')+
#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(1960,2015,5),labels=seq(1960,2015,5))+
    expand_limits(x=c(1960,2015),y=c(0,length(unique(dfs$var))))+
    xlab('Year')+
    ylab('')

y1<-1970
y2<-2010
dfss<-subset(dfs,year>=1970 & year<=2010)
p2<-ggplot()+
geom_tile(data=dfss,aes(x=year,y=as.factor(var),fill=y),col='gray')+
#scale_fill_distiller(palette='Spectral')+
scale_fill_gradient(low='blue',high='red',na.value='white')+
#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(y1,y2,5),labels=seq(y1,y2,5))+
    expand_limits(x=c(y1,y2),y=c(0,length(unique(dfs$var))))+
    xlab('Year')+
    ylab('')
setwd(figsdir)
pdf('indices_availability.pdf',height=24,width=12)
grid.arrange(p1,ncol=1)
grid.arrange(p2,ncol=1)
dev.off()


#GET PROPORTION OF EACH COLUMN THAT IS MISSING
dat<-subset(df,year>=1971 & year<=2005)
pMiss <- function(x){sum(is.na(x))/length(x)*100}
dm<-data.frame(var=names(dat),
               pmis=apply(dat,2,pMiss))
dm<-dm[order(dm$pmis),]
#apply(dat,1,pMiss)
dm2<-subset(dm,pmis<=40)


setwd(datadir)
load("impdatav.std.RData")
#load('SPERA_andata_new.RData')
data<-impdatav.std
names(data)<-gsub('t\\.1','t-1',names(data))
names(data)<-gsub('t\\.2','t-2',names(data))
names(data)<-gsub('t\\.3','t-3',names(data))
nms<-names(data)
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.her',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c('wind.fsd','had.totwgt.rv','had.totno.rv'))]
df<-df[,names(df) %in% as.character(dm2$var)]

#PREDICTORS
preds<-df[,!(names(df) %in% c('sst.t12','herjuv.spcv','herjuv.spnug','herjuv.spvar','herjuv.sprng','cp.cop.cpr.ct','zp.mn.cpr.ct','zp.cpr.ct',"her.cf.rv", "her.dep.rv","her.dumx.rv", "her.dur2.rv",  "her.durng.rv", "her.dvm.rv","her.dvm2.rv", "her.fmass.rv", "her.georng", "her.len.rv","her.metai.rv", "her.prod","her.rec1","her.spcv",  "her.spnug","her.sprng", "her.spvar", "her.ssb","her.szdiv.rv","her.szpe.rv","her.tbio","her.totno.rv","her.totwgt.rv","her.waa",'herjuv.totno.rv','herjuv.totwgt.rv','herjuv.metai','herjuv.dur2.rv','herjuv.len.rv','herlrv.len','herlrv.georng','herjuv.fmass.rv','her.cf.cat','herlrv.sprng','herlrv.spvar','herlrv.surv','herlrv.mn.bof','herjuv.durng.rv','herjuv.dvm2.rv','herlrv.spnug','herjuv.metai.rv','herjuv.dvm.rv','herlrv.dep.bof','herjuv.dumx.rv','herlrv.spcv','herjuv.dep.rv','zp.cop.cpr','her.ajrat.rv','her.eggprod','had.len.rv','had.dep.rv','her.state'))]
#RESPONSE
resp<-df[,(names(df) %in% c('her.ssb','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.spnug','her.spcv','her.spvar','her.prod','her.metai.rv','her.state'))]
resp$her.spnug<-resp$her.spnug*-1
resp$her.spcv<-resp$her.spcv*-1
resp$her.spvar<-resp$her.spvar*-1
df<-cbind(resp,preds)

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)



mat<-subset(df,year>=1970)
mat<-mat[,!(names(mat) %in% c('year'))]
nms<-names(mat)
#dataset <- BNDataset(data = mat, discreteness=rep(FALSE,ncol(mat)), variables=nms, node.sizes=rep(2,dim(mat)[2]))


#CONSTRUCT BLACKLIST: NETWORK CONNECTIONS THAT ARE IMPOSSIBLE
#NOTHING CAN EFFECT EXPLOITATION VARIABLES
rspnms<-c("her.expr","her.land","her.land.pct1","her.land.spdiv","her.land.sprich",'bfish.state')
prdnms<-names(mat[,!(names(mat) %in% c('year',rspnms))])
l<-list()
for(i in 1:length(prdnms)){
print(prdnms[i])
l2<-list()
for(j in 1:length(rspnms)){
l2[[j]]<-data.frame(from=prdnms[i],
                    to=rspnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
bl1<-data.frame(do.call('rbind',l))


#BIOLOGICAL VARIABLES CANNOT EFFECT CLIMATE
clmnms<-c("asl","gs.dist","nao", "nut.state","ss.ct","ss.dist","sst.fsd","sst.tmax","strt","strt.ct","t50.ct","temp.state","wind.amp",  "wind.state","wind.tmax", "wnd.ct","wnd.strs323fall")
bnms<-names(mat[,!(names(mat) %in% clmnms)])
l<-list()
for(i in 1:length(bnms)){
print(bnms[i])
l2<-list()
for(j in 1:length(clmnms)){
l2[[j]]<-data.frame(from=bnms[i],
                    to=clmnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
bl2<-data.frame(do.call('rbind',l))
bl<-unique(rbind(bl1,bl2))

#FIT BAYESIAN NETWORK VIA
bn<-gs(mat,blacklist=bl,undirected=FALSE)
#bn<-hc(mat,blacklist=bl)
#bn2<-iamb(mat,blacklist=bl)
#bn3<-fast.iamb(mat,blacklist=bl)
#bn4<-mmpc(mat,blacklist=bl)
#st<-arc.strength(bn,data=mat,criterion='x2',blacklist=bl,algorithm='gs')

#strength.plot(bn,st)
st<-boot.strength(mat,R=3,algorithm='gs',algorithm.args=list(blacklist=bl),cpdag=FALSE)
std<-data.frame(st)

#ADD STRENGTHS TO ARCS
acd<-data.frame(bn$arcs)
acd<-merge(acd,std,by=c('from','to'),all.x=TRUE,all.y=FALSE)
a<-data.frame(bn$arcs)
acd<-merge(a,acd,by=c('from','to'),all.x=TRUE,all.y=FALSE,sort=FALSE)

g<-graph_from_edgelist(bn$arcs,directed=TRUE)
V(g)$color<-'gold3'

rnms<-names(resp)
for(i in 1:length(rnms)){
nm<-rnms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('red3',1)
} else NULL
}

nms<-c('ph.div.state','ph.mn.cpr','ph.mn.sabs','ph.pe.state','ph.to.sabs')
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('lightgreen',1)
} else NULL
}

nms<-c('zp.div.cpr.state','zp.div.sabs.state','zp.mn.cpr','zp.mn.sabs','zp.pe.cpr','zp.pe.sabs','zp.to.sabs','lrv.mn.bof','lrv.pe.bof','lrv.rich.bof','her.preytot.bof','herjuv.prey.bof','her.jf.bof')
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('lightskyblue2',1)
} else NULL
}

nms<-c('her.expr','her.land','her.land.pct1','her.land.spdiv','her.land.spring','bfish.state')
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('lightpink',1)
} else NULL
}



#plot(g,layout=layout_in_circle)
#plot(g,layout=layout.fruchterman.reingold(g, niter=10000),vertex.size=10)
ly<-layout.fruchterman.reingold(g)
ly<-layout_(g,as_star())
ly<-layout_with_fr
ly<-layout_with_kk(g)
ly<-layout_drl
ly<-layout_lgl(g)
layout = g.layout('kamada_kawai')
layout = g.layout("kk")
ly<-layout_nicely(g)

setwd(figsdir)
pdf('herring_pressures_network.pdf',height=10,width=10)
plot(g,
     layout=ly,
     edge.width=acd$strength*10,  # Edge width, defaults to 1
     edge.arrow.size=.5,           # Arrow size, defaults to 1
     edge.arrow.width=1.5,          # Arrow width, defaults to 1
     edge.curved=.7,
     vertex.label.color=alpha('black',.9),
     vertex.label.cex=.75,           #Font size
     vertex.label.dist=1.25,#Distance between label and the vertex
     vertex.shape="circle",
     vertex.size=7,       #Size of the node (default is 15)
     vertex.frame.color='gray60',     # Node border color
     edge.color=alpha('gray',1),
     xlim=c(-1,1),ylim=c(-1,1),asp=0,lwd=.001)
dev.off()

plot(df$year,df$her.ssb,pch=15)
plot(df$year,df$her.rec1,pch=15)

plot(g,axes=TRUE,rescale=FALSE)

nr<-names(resp)
n2<-names(preds)
graphviz.plot(bn,shape='circle',highlight=list(nodes=names(mat),col=ifelse(names(mat) %in% nr,'tomato','blue')))

graphviz.plot(bn,shape='circle',highlight=list(nodes=names(mat),fill='tomato',col='gray'))

graphviz.plot(bn,shape='circle',highlight=list(nodes=names(mat),fill='tomato',col='gray'))

bn.ia<-iamb(mat)
bnh<-hc(mat)
graphviz.plot(bn.ia,shape='circle',highlight=list(nodes=names(mat),fill='tomato',col='gray'))
graphviz.plot(bnh,shape='circle',highlight=list(nodes=names(mat),fill='tomato',col='gray'))
graphviz.chart(bn)


#show(dataset)
#raw.data(dataset)
nw<-learn.network(dataset)
plot(nw,node.col = rep("lightblue", num.nodes(nw)))
dg<-dag(nw)
wpdg<-wpdag(nw)
plot(nw)

wpdag(nw)
parms<-cpts(nw)
#ggnet2(wpdg,node.color='tomato',edge.color='gray',node.size=10,label=nms,mode='kamadakawai')

ggnet2(dg,node.color='tomato',edge.color='gray',node.size=10,mode='kamadakawai',label=nms)

plot(nw,attrs = list(node = list(font=10,fillcolor = "lightblue"), edge = list(arrowsize=0.5)))

plot(nw,y='dot',cex.sub=2,recipEdges='distinct',graph=list(rankdir='RL'))

nw<-learn.dynamic.network(dataset)
data(child)

dataset <- BNDataset(data = mat,  discreteness = rep(F, ncol(mat)), variables = varNames, node.sizes = rep(3, ncol(mat)))

bn <- learn.network(dataset)

dbn <- learn.dynamic.network(dataset, num.time.steps = numEvents)
plot(dbn)







#########################################################
setwd(figsdir)
pdf('statesandpressures_hcluster.pdf',height=14,width=10)
par(mar=c(4,4,4,4),oma=c(1,1,1,1))

K<-6
mds<-as.data.frame(cmdscale(dst))
colnames(mds)<-c('Dim.1','Dim.2')
p1<-ggscatter(mds,x='Dim.1',y='Dim.2',label=names(dff),size=2,repel=TRUE,alpha=.5)

#KMEANS CLUSTERING
#dcol<-data.frame(cl=c('gold3','firebrick3','forestgreen','dodgerblue3','purple','orange'),clust=seq(1,k,1))
dcol<-data.frame(cl=c('gold3','firebrick3','magenta3','green3','darkblue','darkgreen'),clust=seq(1,K,1))

clust <- kmeans(mds, K)$cluster %>%
as.factor()
mds <- mds %>%
mutate(groups = clust)
# Plot and color by groups
p2<-ggscatter(mds, x = "Dim.1", y = "Dim.2",
          label = names(dff),
          color = "groups",
#          palette = "jco",
          palette =as.character(dcol$cl),
          size = 1,
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
grid.arrange(p1,p2,ncol=1)
#grid.arrange(p1,p2,ncol=2)



#HCLUST PLOTS
par(mfrow=c(2,1),mar=c(4,6,1,6))
hc<-hclust(dst,method='ward.D2')
sub_grp<-cutree(hc,k=K)

hc2<-as.dendrogram(hc)
plot(hc2)
rect.hclust(hc, k = K, border = 2:5)


#PLOT CLUST AS HORIZONTAL
agh<-as.hclust(hc)
labelColors<-c('gold3','firebrick3','magenta3','green3','darkblue','lightblue')
clusMember = cutree(agh, k)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
clusDendro = dendrapply(as.dendrogram(agh), colLab)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19))
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.01,pch=16,nodePar=nodePar)


#GETS NAMES OF VARIABLES IN EACH CLUSTER
cldat<-data.frame(var=names(dff),
               clust=clust)
cldat<-cldat[order(cldat$clust),]
cldat<-merge(cldat,dcol,by=c('clust'),all=FALSE)


#GETS AVERAGE CORREALTION, RATE OF CHANGE AND  PLOTS TS FOR EACH CLUSTER
par(mfrow=c(4,2))
for(i in 1:6){
print(i)
cl<-subset(cldat,clust==i)
    d<-subset(pdatl,select=as.character(cl$var))
    d2<-na.omit(subset(pdats,var %in% as.character(cl$var)))
    cr<-cor(d,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
print(mean(cr.t$Freq))
mod<-gamm(y~year,data=d2,random=list(var=~1))
modf<-gamm(y~as.factor(year),data=d2,random=list(var=~1))
mods<-gamm(y~s(year),data=d2,random=list(var=~1))
pdat<-data.frame(year=sort(unique(d2$year)))
pdatsm<-data.frame(year=seq(min(d2$year),max(d2$year),length.out=1000))
psm<-predict(mods$gam,newdata=pdatsm,se.fit=TRUE)
pdatsm$p<-psm$fit
pdatsm$se<-psm$se.fit
pdatsm$upr<-pdatsm$p+(1.96*pdatsm$se)
pdatsm$lwr<-pdatsm$p-(1.96*pdatsm$se)
p<-predict(modf$gam,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*p$se)
pdat$lwr<-pdat$p-(1.96*p$se)
plot(pdat$year,pdat$p,las=1,pch=15,xlim=c(1900,2018),ylim=c(-2,2),col=unique(as.character(cl$cl)),xaxt='n',xlab='Year',ylab='Mixed model average')
axis(1,seq(1900,2018,5))
abline(h=0,lty=2)
polygon(c(pdatsm$year,pdatsm$year[length(pdatsm$year):1]),c(pdatsm$upr,pdatsm$lwr[length(pdatsm$lwr):1]),col=alpha(as.character(unique(cl$cl)),.3),border=NA)
s<-summary(mod$gam)
r2<-round(s$r.sq,digits=2)
f2<-function(g){    lines(c(g$year,g$year),c(g$upr,g$lwr),col=alpha(unique(as.character(cl$cl)),.3))}
#zz<-dlply(pdat,.(year),.fun=f2)
legend('topright',c(paste('r2=',r2),paste('av r=',round(mean(cr.t$Freq),digits=2))),bty='n')
lines(pdatsm$year,pdatsm$p,col=unique(as.character(cl$cl)),lwd=2)
}
dev.off()



















dat<-cbind(resp,preds)
dat$year<-data$year

#GET PROPORTION OF EACH COLUMN THAT IS MISSING
#dat<-subset(data,year>=1965 & year<=2010)
dat<-subset(dat,year>=1970 & year<=2010)
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





setwd(datadir)
load('impdatav.RData')
datai<-impdatav

setwd(datadir)
load('SPERA_andata_new.RData')
datai<-data
nms<-names(datai)
#nms<-nms[grepl('her',nms)==TRUE]
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
#nms<-nms[grepl('\\.cat',nms)==FALSE]
nms<-nms[grepl('sp\\.n',nms)==FALSE]
df<-datai[,names(datai) %in% c('year',nms)]
df<-subset(df,year>=1970 & year<=2010)
df<-df[,!(names(df) %in% c("her.ajrat.rv", "her.cf.cat", "her.cf.rv", "her.dep.rv", "her.dvm.rv", "her.dvm2.rv","her.fmass.rv", "her.georng","her.len.rv","her.metai.rv", "her.prod","her.rec1", "her.spcv", "her.spnug", "her.sprng","her.spvar","her.ssb", "her.szpe.rv","her.tbio", "her.totno.rv","her.totwgt.rv","her.waa","herjuv.dumx.rv","herjuv.dvm.rv","herjuv.fmass.rv","herjuv.len.rv", "herjuv.totno.rv","herjuv.totwgt.rv","herlrv.dep.bof", "herlrv.len", "herlrv.mn.bof","herlrv.surv",'herjuv.spnug','herjuv.spcv','herjuv.sprng','herjuv.spvar'))]

#IMPUTATION DIDN'T WORK FOR THIS
df<-df[,!(names(df) %in% c('sst.t12'))]

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)

#LOOK AT SKEWNESS OF VARIABLES
sk<-data.frame(var=sort(unique(dfs$var)),
         sk=tapply(dfs$y,dfs$var,function(x) skewness(x,na.rm=TRUE)))
sk<-sk[order(abs(sk$sk),decreasing=TRUE),]

#TRANSFORM VARIABLES TO OPTIMIZE NORMALITY
f<-function(d){
    d$y<-transformTukey(d$y,plotit=FALSE)
    return(d)
}
dfs<-ddply(dfs,.(var),.fun=f)


#LOOK AT SKEWNESS OF VARIABLES
sk<-data.frame(var=sort(unique(dfs$var)),
         sk=tapply(dfs$y,dfs$var,function(x) skewness(x,na.rm=TRUE)))
sk<-sk[order(abs(sk$sk),decreasing=TRUE),]

#CONVERT BACK TO LONG FORM
df<-spread(data=dfs,key=var, value=y)

df2<-data.frame(t(df))
names(df2)<-df2[1,]
df2<-df2[-1,]

#histogram(~y | var,data=dfs)

f<-function(d){
    d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
    return(d)
}
pdats<-ddply(dfs,.(var),.fun=f)

pdatl<-spread(data=pdats,key=var, value=y)

histogram(~y | var,data=pdats)
xyplot(y~year|var,data=pdats,pch=15,type=c('p','spline'))
plot(data$year,data$lrv.mn.bof,pch=15)

dff<-pdatl[,-1]
dff2<-pdats
dmat<-1-(cor(dff,use='pairwise.complete.obs',method='spearman'))
#dmat<-abs(cor(dff,use='pairwise.complete.obs',method='spearman'))
dst<-as.dist(dmat)




setwd(figsdir)
#pdf('herring_state_hcluster_spawnar.pdf',height=14,width=10)
#pdf('herring_state_hcluster.pdf',height=14,width=10)
par(mar=c(4,4,4,4),oma=c(1,1,1,1))
#library(tibble)
mds<-as.data.frame(cmdscale(dst))
colnames(mds)<-c('Dim.1','Dim.2')
p1<-ggscatter(mds,x='Dim.1',y='Dim.2',label=names(dff),size=2,repel=TRUE,alpha=.5)

#KMEANS CLUSTERING
k<-12
#dcol<-data.frame(cl=c('gold3','firebrick3','forestgreen','dodgerblue3','purple','orange'),clust=seq(1,k,1))
dcol<-data.frame(cl=rainbow(k),clust=seq(1,k,1))

clust <- kmeans(mds, k)$cluster %>%
as.factor()
mds <- mds %>%
mutate(groups = clust)
# Plot and color by groups
p2<-ggscatter(mds, x = "Dim.1", y = "Dim.2",
          label = names(dff),
          color = "groups",
#          palette = "jco",
          palette =as.character(dcol$cl),
          size = 1,
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
grid.arrange(p1,p2,ncol=2)



#GETS NAMES OF VARIABLES IN EACH CLUSTER
cldat<-data.frame(var=names(dff),
               clust=clust)
cldat<-cldat[order(cldat$clust),]
cldat<-merge(cldat,dcol,by=c('clust'),all=FALSE)

#HCLUST PLOTS
par(mfrow=c(2,1),mar=c(4,6,1,6))
hc<-hclust(dst,method='ward.D2')
sub_grp<-cutree(hc,k=5)

hc2<-as.dendrogram(hc)
plot(hc2)
rect.hclust(hc, k = 3, border = 2:5)


#PLOT CLUST AS HORIZONTAL
agh<-as.hclust(hc)
labelColors = c('red3','royalblue','forestgreen')
clusMember = cutree(agh, 3)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
clusDendro = dendrapply(as.dendrogram(agh), colLab)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19))
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7,pch=16,nodePar=nodePar)

#GETS AVERAGE CORREALTION, RATE OF CHANGE AND  PLOTS TS FOR EACH CLUSTER
par(mfrow=c(4,2))
for(i in 1:4){
print(i)
cl<-subset(cldat,clust==i)
    d<-subset(pdatl,select=as.character(cl$var))
    d2<-na.omit(subset(pdats,var %in% as.character(cl$var)))
    cr<-cor(d,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
print(mean(cr.t$Freq))
mod<-gamm(y~year,data=d2,random=list(var=~1))
modf<-gamm(y~as.factor(year),data=d2,random=list(var=~1))
mods<-gamm(y~s(year),data=d2,random=list(var=~1))
pdat<-data.frame(year=sort(unique(d2$year)))
pdatsm<-data.frame(year=seq(min(d2$year),max(d2$year),length.out=1000))
psm<-predict(mods$gam,newdata=pdatsm,se.fit=TRUE)
pdatsm$p<-psm$fit
pdatsm$se<-psm$se.fit
pdatsm$upr<-pdatsm$p+(1.96*pdatsm$se)
pdatsm$lwr<-pdatsm$p-(1.96*pdatsm$se)
p<-predict(modf$gam,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*p$se)
pdat$lwr<-pdat$p-(1.96*p$se)
plot(pdat$year,pdat$p,las=1,pch=15,xlim=c(1970,2015),ylim=c(-1.5,1.5),col=unique(as.character(cl$cl)),xaxt='n',xlab='Year',ylab='Mixed model average')
axis(1,seq(1970,2015,5))
abline(h=0,lty=2)
polygon(c(pdatsm$year,pdatsm$year[length(pdatsm$year):1]),c(pdatsm$upr,pdatsm$lwr[length(pdatsm$lwr):1]),col=alpha(as.character(unique(cl$cl)),.3),border=NA)
s<-summary(mod$gam)
r2<-round(s$r.sq,digits=2)
f2<-function(g){    lines(c(g$year,g$year),c(g$upr,g$lwr),col=alpha(unique(as.character(cl$cl)),.3))}
#zz<-dlply(pdat,.(year),.fun=f2)
legend('topright',c(paste('r2=',r2),paste('av r=',round(mean(cr.t$Freq),digits=2))),bty='n')
lines(pdatsm$year,pdatsm$p,col=unique(as.character(cl$cl)),lwd=2)
}
dev.off()




#pdf('herring_state_hcluster_heatmap.pdf',height=8,width=12)
#pdf('herring_state_hcluster_heatmap_spawnar.pdf',height=8,width=12)
#PLOTS OUT 2D CLUSTER OF SELECTED VARIABLES
#ALL DATA
pdatll<-subset(pdatl,year<=2014)
pdatll<-pdatll[order(pdatll$year,decreasing=FALSE),]

pdatl2<-data.frame(t(pdatll))
names(pdatl2)<-pdatl2[1,]
pdatl2<-pdatl2[-1,]

pheatmap(pdatl2,cutree_cols=2,cutree_rows=3,clustering_method='ward.D')

#CLUSTERING BASED ON ABSOLUTE CORRELATION MATRIX
dmat2<-abs(cor(dff,use='pairwise.complete.obs',method='spearman'))
pheatmap(dmat2,cutree_cols=2,cutree_rows=2,clustering_method='ward.D2')

#CLUSTERING BASED ON CORRELATION MATRIX
dmat2<-(cor(dff,use='pairwise.complete.obs',method='spearman'))
pheatmap(dmat2,cutree_cols=3,cutree_rows=3,clustering_method='ward.D2')

#ONLY DATA THAT CLUSTER TOGETHER
a<-subset(cldat,clust %in% c(1,3))
pdatll<-subset(pdatl,year<=2014,select=c('year','her.prod','her.metai.rv',as.character(a$var)))
#pdatll<-pdatll[order(pdatll$year,decreasing=TRUE),]
pdatll<-pdatll[order(pdatll$year,decreasing=FALSE),]
pdatl2<-data.frame(t(pdatll))
names(pdatl2)<-pdatl2[1,]
pdatl2<-pdatl2[-1,]

pheatmap(pdatl2,cutree_cols=2,cutree_rows=5,clustering_method='ward.D')

#TRANSFORM SO ALL INDICES MOVE IN SAME DIRECTION
a<-subset(cldat,clust %in% c(1,3))
pdatll<-subset(pdatl,year<=2014,select=c('year','her.prod','her.metai.rv',as.character(a$var)))
#pdatll<-pdatll[order(pdatll$year,decreasing=TRUE),]
#pdatll$herjuv.totwgt.rv<-pdatll$herjuv.totwgt.rv*-1
#pdatll$herjuv.totno.rv<-pdatll$herjuv.totno.rv*-1
#pdatll$her.totwgt.rv<-pdatll$her.totwgt.rv*-1
#pdatll$her.totno.rv<-pdatll$her.totno.rv*-1
pdatll$her.spcv<-pdatll$her.spcv*-1
pdatll$her.spvar<-pdatll$her.spvar*-1
pdatll$her.spnug<-pdatll$her.spnug*-1

pdatll<-pdatll[order(pdatll$year,decreasing=FALSE),]
pdatl2<-data.frame(t(pdatll))
names(pdatl2)<-pdatl2[1,]
pdatl2<-pdatl2[-1,]

pheatmap(pdatl2,cutree_cols=2,cutree_rows=5,clustering_method='ward.D')

#CLUSTERING BASED ON CORRELATION MATRIX
pdatll<-pdatll[,-1]
dmat2<-(cor(pdatll,use='pairwise.complete.obs',method='spearman'))
pheatmap(dmat2,cutree_cols=3,cutree_rows=3,clustering_method='ward.D2')

dev.off()




##########################################################

#  FITS TIMESERIES MODELS TO NORMALIZED SERIES
dati<-pdatl
dati$her.spcv<-dati$her.spcv*-1
dati$her.spvar<-dati$her.spvar*-1
dati$her.spnug<-dati$her.spnug*-1

df<-dati %>% gather(var,y,her.ajrat.rv:herlrv.surv)


d<-subset(df,var=='herjuv.dvm.rv')

f<-function(d){
print(unique(d$var))
d<-na.omit(d)
modlm<-lm(y~year,data=d)
slm<-summary(modlm)
d$resid<-residuals(modlm)
t<-ts(d$resid)
acpar<-mean(acf(t,plot=F)$acf[2])

mod<-gls(y~year,data=d,method='ML',correlation=corAR1(acpar, form = ~year))
s<-summary(mod)

#modg<-gamm(y~s(year),data=d,correlation=corCAR1(value=acpar, form = ~year))
modg<-gamm(y~s(year,k=7),data=d,correlation=corCAR1(form = ~year),gamma=1)
sg<-summary(modg$gam)
plot(d$year,d$y,xlab='Year',ylab='Z-score',las=1,pch=15)
legend('topright',unique(d$var),bty='n')
pdat<-data.frame(year=seq(min(d$year),max(d$year),1))
pdat$p<-predict(modg$gam,newdata=pdat)
lines(pdat$year,pdat$p)

    return(round(data.frame(b=slm$coef[2,1],
                      pv=slm$coef[2,4],
                      bgls=s$tTable[2,1],
                      pvgls=s$tTable[2,4],
                      pvgam=sg$s.table[1,4]),digits=3))
}
setwd(figsdir)
pdf('herring_indices_timetrends.pdf',height=10,width=10)
par(mfrow=c(4,4),mar=c(4,4,1,1))
ot<-ddply(df,.(var),.fun=f,.progress='text')
dev.off()


ot<-ot[order(ot$pvgam),]
write.csv(ot,file='herring_indices_timetrends.csv',row.names=FALSE)

s<-subset(ot,pv<0.05)
s2<-subset(ot,pvgls<0.05)
s3<-subset(ot,pvgam<0.05)
s3<-subset(ot,pvgam<0.1 | pv<.05)

plot(data$year,data$herjuv.spcv,pch=15)




dff<-pdatl[,-1]
dff2<-pdats
dmat<-1-(cor(dff,use='pairwise.complete.obs',method='spearman'))
dst<-as.dist(dmat)


mds<-as.data.frame(cmdscale(dst))
colnames(mds)<-c('Dim.1','Dim.2')
mds$var<-rownames(mds)
mds<-merge(mds,ot,by=c('var'),all=TRUE)
#mds$cl<-ifelse(mds$pvgam<=.05,'firebrick3','black')
mds$cl<-ifelse(mds$pvgam<=.1,'gold3','black')
mds$cl<-ifelse(mds$pvgam<=.05,'firebrick3',mds$cl)


p2<-ggscatter(mds,x='Dim.1',y='Dim.2',label=names(dff),size=5,repel=TRUE,alpha=1,color=mds$cl,font.label=c(14,'plain','gray50'))

p2<-ggscatter(mds,x='Dim.1',y='Dim.2',label='',size=5,repel=TRUE,alpha=1,color=mds$cl,font.label=c(14,'plain','gray50'))

clust <- kmeans(mds, 4)$cluster %>%
as.factor()
mds <- mds %>%
mutate(groups = clust)
p1<-ggscatter(mds, x = "Dim.1", y = "Dim.2",
          label = '',
          color = "groups",
          palette=rep('black',4),
          size = 1,
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)

setwd(figsdir)
pdf('herring_state_timeseries_cluster.pdf',height=14,width=10)
par(mar=c(4,4,4,4),oma=c(1,1,1,1))
grid.arrange(p2,p1,ncol=1)
dev.off()




















#################################################

######### HERRING STATE DERIVATION
setwd(datadir)
load('SPERA_andata_new.RData')

#LOAD LANDINGS FROM FAO
d<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/TS_FI_CAPTURE.csv',header=TRUE)
spc<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/CL_FI_SPECIES_GROUPS.csv',header=TRUE)
names(d)<-tolower(names(d))
names(spc)<-tolower(names(spc))
names(spc)[1]<-'species'
dher<-subset(d,species=='HER' & country==124)

d<-subset(dher,select=c('year','quantity'))
names(d)[2]<-'her.land2'
data<-merge(data,d,by=c('year'),all.x=TRUE,all.y=FALSE)


#SELECT VARIABLES FOR STATE CALCULATION: ADD PRODUCTION EVEN THOUGH NOT IN CLUSTERS; REMOVE TOTAL BIOMASS BECAUSE CORRELATED WITH SSB
a<-subset(cldat,clust %in% c(1,3))
dat<-subset(data,year>=1965,select=c('year','her.prod','her.metai.rv',as.character(a$var)))
#dat<-dat[,!(names(dat) %in% c('her.totwgt.rv','her.totno.rv','herjuv.totno.rv','herjuv.totwgt.rv','her.georng','her.tbio','her.prod','her.metai.rv'))]
dat<-dat[,!(names(dat) %in% c('her.totwgt.rv','her.totno.rv','herjuv.totno.rv','herjuv.totwgt.rv','her.georng','her.tbio'))]

#abs(round(cor(dat,use='pairwise.complete.obs'),digits=2))
(round(cor(dat,use='pairwise.complete.obs'),digits=2))

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
f<-function(x){transformTukey(x,plotit=FALSE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

f<-function(x){scale(x,center=TRUE,scale=TRUE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

#TRANSFORM TO ENSURE GAUSSIAN DISTRIBUTION
dat$her.spcv<-dat$her.spcv*-1
dat$her.spvar<-dat$her.spvar*-1
dat$her.spnug<-dat$her.spnug*-1

dats<-dat %>% gather(var,y,-year)
#bb<-na.omit(bb)
xyplot(y~year | var,data=dats,pch=15,type=c('p','l'))
histogram(~y | var,data=dats)

mod<-gamm(y~as.factor(year),data=dats,random=list(var=~1))
pdat<-data.frame(year=sort(unique(dats$year)))
p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
pdat$her.state<-p$fit
pdat$her.state.se<-p$se.fit
plot(pdat$year,pdat$her.state,pch=15)
dat<-merge(dat,pdat,by=c('year'),all=TRUE)

cr<-cor(dat,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','her.state.se')) & !(Var1 %in% c('year','her.state.se')))

cr.t<-cr.t[order(cr.t$Freq,decreasing=TRUE),]
cr.t$Var1<-as.character(cr.t$Var1)
cc<-data.frame(var1=sort(unique(cr.t$Var1)),
               r=tapply(cr.t$Freq,cr.t$Var1,mean))

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='her.state' & Var2!='her.state')
mean(cc$Freq);#.47
median(cc$Freq);#.51
#mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.54; MEDIAN=0.56


cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))
n<-length(unique(dats$var))
dum3<-data.frame(var=sort(unique(dats$var)),
                 cls=cls(n+2)[3:(n+2)])
dats<-merge(dats,dum3,by=c('var'),all.x=TRUE,all.y=FALSE)

setwd(figsdir)
pdf('herring_state_derivation_transform.pdf',height=7,width=4)
par(mfcol=c(8,2),mar=c(.5,1,0,0),oma=c(4,4,1,4))
f<-function(d){
d<-na.omit(d)
ylm<-c(floor(min(d$y,na.rm=TRUE)),ceiling(max(d$y,na.rm=TRUE)))
print(ylm)
asp<-aspline(d$year,d$y,xout=seq(min(d$year),max(d$year),length.out=1000))
plot(0,0,ylim=ylm,xlim=c(1965,2015),axes=FALSE)
abline(h=0,col='lightgray',lty=2,lwd=.5)
lines(asp$x,asp$y,col=alpha(as.character(unique(d$cl)),1),lwd=2)
points(d$year,d$y,col=alpha(as.character(unique(d$cl)),.3),cex=1,pch=16)
axis(2,at=ylm,las=1,cex.axis=.5,lwd=.1)
if(unique(d$var) %in% c('her.spcv','herlrv.len')){
axis(1,seq(1965,2015,5),labels=FALSE,cex.axis=.5,lwd=.1)
axis(1,seq(1965,2015,10),cex.axis=.5,lwd=.1)
} else NULL
mod<-lm(y~year,data=d)
pdat<-data.frame(year=seq(min(d$year),max(d$year),length.out=1000))
pdat$p<-predict(mod,newdata=pdat)
#lines(pdat$year,pdat$p,col=alpha('black',.5),lwd=.5)
s<-summary(mod)
b<-round(s$coef[2,1],digits=2)
legend('topright',gsub(' ','',paste(unique(d$var),'=',b)),bty='n',cex=.5)
}
z<-dlply(dats,.(var),.fun=f)



dats$segtrue<-ifelse(dats$var%in% c('herjuv.fmass.rv','her.ajrat.rv'),FALSE,TRUE)
#dats$segtrue<-TRUE

par(mfcol=c(8,2),mar=c(.5,1,0,0),oma=c(4,4,1,4))
f<-function(d){
d<-na.omit(d)
print(unique(d$var))
#CHARACTERIZE TRENDS
    if(length(unique(d$year))<=10){dff<-4
    } else {dff<-5
    }
    bp<-floor(mean(d$year))
    modl<-lm(y~year,data=d)
    m<-breakpoints(y~year,data=d,h=3,breaks=1)
    modst<-lm(y~breakfactor(m,breaks=length(unique(m$breakpoints))),data=d)
    modnl<-lm(y~bs(year,degree=3,df=dff),data=d)
    if(unique(d$segtrue)==TRUE){
        modseg<-segmented(modl,seg.Z = ~year,psi=bp,control=seg.control(it.max=200))
        dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
            dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modst',dt$AIC,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
    modl<-modnl
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
    d3<-data.frame(year=seq(min(d$year),max(d$year),1))
} else if (dt$md=='modseg'){
    modl<-modseg
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
    d3<-data.frame(year=seq(min(d$year),max(d$year),1))
} else if (dt$md=='modst'){
    modl<-modst
    d2<-data.frame(year=sort(unique(d$year)))
    d3<-data.frame(year=sort(unique(d$year)))
} else {
    modl<-modl
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
    d3<-data.frame(year=seq(min(d$year),max(d$year),1))
}


    p<-predict(modl,newdata=d3,se.fit=TRUE,type='response')
    d3$p<-p$fit
    d3$se<-p$se.fit
    d3$upr<-d3$p+(1.96*d3$se)
    d3$lwr<-d3$p-(1.96*d3$se)

    s<-summary(modl)
    r2<-round(s$r.squared,digits=2)
    p<-predict(modl,newdata=d2,se.fit=TRUE,type='response')
    d2$p<-p$fit
    d2$se<-p$se.fit
    d2$upr<-d2$p+(1.96*d2$se)
    d2$lwr<-d2$p-(1.96*d2$se)
    ylm<-c(min(c(d$y,d2$lwr)),max(c(d$y,d2$upr)))

ylm<-c(floor(min(d$y,na.rm=TRUE)),ceiling(max(d$y,na.rm=TRUE)))
print(ylm)
plot(0,0,ylim=ylm,xlim=c(1965,2015),axes=FALSE)
polygon(c(d2$year,d2$year[length(d2$year):1]),c(d2$upr,d2$lwr[length(d2$lwr):1]),col=alpha(as.character(unique(d$cl)),.3),border=NA)
abline(h=0,col='lightgray',lty=2,lwd=.5)
points(d$year,d$y,col=alpha(as.character(unique(d$cl)),.3),cex=1,pch=16)
axis(2,at=ylm,las=1,cex.axis=.5,lwd=.1)
if(unique(d$var) %in% c('her.spcv','herlrv.len')){
axis(1,seq(1965,2015,5),labels=FALSE,cex.axis=.5,lwd=.1)
axis(1,seq(1965,2015,10),cex.axis=.5,lwd=.1)
} else NULL
lines(d2$year,d2$p,col=as.character(unique(d$cl)),lwd=2)
s<-summary(modl)
r2<-round(s$r.squared,digits=2)
legend('topright',gsub(' ','',paste(unique(d$var),'=',r2)),bty='n',cex=.5)
d2$year<-round(d2$year,digits=0)
dout<-d3
dout$chng<-dout$p[1]-dout$p[length(dout$p)]

    #GET BREAKPOINT
    if(dt$md=='modst'){
    dm<-data.frame(m$X)
    dm$id<-seq(1,dim(dm)[1],1)
    dm$bp<-breakfactor(m,breaks=length(unique(m$breakpoints)))
    dout$bpt<-max(subset(dm,bp=='segment1')$year)
    } else if (dt$md=='modseg'){
    dout$bpt<-round(data.frame(modseg$psi)$Est.,digits=0)
    } else {
    dout$bpt<-subset(d2,p==max(d2$p))$year[1]
    }

return(data.frame(dout))
}
z<-ddply(dats,.(var),.fun=f)
z<-z[order(z$chng),]

#IF BREAK IS A MINIMUM, TAKE START OF LINEAR CHANGE INSTEAD
f<-function(d){
pbt<-subset(d,year==unique(d$bpt))
    if(pbt$p<=mean(d$p)){
    bpt<-subset(d,p==max(d$p))$year[1]
    } else { bpt<-unique(d$bpt)
         }

if(unique(d$var)=='her.metai.rv'){
    bpt<-min(d$year,na.rm=TRUE)
} else  bpt<-bpt

d$bpt<-bpt
return(d)
}
z<-ddply(z,.(var),.fun=f)
dev.off()




pdf('herring_state_derivation_heatmap.pdf',height=8,width=12)
p1<-ggplot()+
geom_tile(data=z, aes(x=year, y=as.factor(var),fill=p),col='gray80')+
scale_fill_distiller(palette='Spectral')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))
grid.arrange(p1,ncol=1)

a<-unique(subset(z,select=c('bpt','var')))
a<-a[order(a$bpt,decreasing=TRUE),]
a$id<-seq(1,dim(a)[1],1)

z<-merge(z,a,by=c('var'),all.x=TRUE,all.y=FALSE)
z<-z[order(z$id),]
names(a)[1]<-'year'

ggplot()+
geom_tile(data=z, aes(x=year, y=id,fill=p),col='gray80')+
geom_tile(data=a,aes(x=year,y=id,fill=NULL,alpha=0))+
scale_fill_distiller(palette='Spectral')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.09,.5),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(1965,2015,5),labels=seq(1965,2015,5))+
  scale_y_continuous(expand=c(0,0),breaks=seq(1,16,1),labels=a$var)+
    expand_limits(x=c(1964.5,2016),y=c(0,17))+
    xlab('Year')+
    ylab('Landings [%]')

#TRANSFORM SO ALL INDICES MOVE IN SAME DIRECTION
pdatll<-subset(dat,year<=2014)
pdatll<-pdatll[,!(names(pdatll) %in% c('her.state.se','her.state'))]

pdatll<-pdatll[order(pdatll$year,decreasing=FALSE),]
pdatl2<-data.frame(t(pdatll))
names(pdatl2)<-pdatl2[1,]
pdatl2<-pdatl2[-1,]

pheatmap(pdatl2,cutree_cols=2,cutree_rows=5,clustering_method='ward.D')

#CLUSTERING BASED ON CORRELATION MATRIX
pdatll<-pdatll[,-1]
dmat2<-(cor(pdatll,use='pairwise.complete.obs',method='spearman'))
pheatmap(dmat2,cutree_cols=3,cutree_rows=3,clustering_method='ward.D2')

a$year2<-ifelse(a$var=='her.cf.rv',1969,a$year)
a$year2<-ifelse(a$var=='her.ajrat.rv',1967,a$year2)
a$year2<-ifelse(a$var=='her.szpe.rv',1968,a$year2)
a$year2<-ifelse(a$var=='herlrv.len',1974,a$year2)
a$year2<-ifelse(a$var=='her.metai.rv',1966,a$year2)

a$lbl<-c( 'Herring weight-at-age','Herring spatial variance','Juvenile herring metabolic rate','Herring SSB','Herring productivity','Average mass of herring','Hering spatial covariance','Average length of herring','Average length of herrng larve','Average mass of juvenile herring','Ratio of adult to juvenile herring','Herring condition factor','Herring metabolic rate','Small-scale spatial variance of herring','Size evenness of herring','Recruitment at age 1')
datt<-na.omit(subset(dat,select=c('year','her.ssb')))
pdat<-data.frame(year=seq(min(datt$year),max(datt$year),length.out=1000))
pdat$y<-aspline(datt$year,datt$her.ssb,xout=pdat$year)$y
plot(dat$year,dat$her.ssb,las=1,pch=16,ylim=c(-2,2),cex=2,col=alpha('dodgerblue3',.3),xlim=c(1965,2010),xlab='Year',ylab='Herring SSB')
points(dat$year,dat$her.ssb,pch=1,cex=2,col='dodgerblue3')
points(a$year,rep(-2.15,dim(a)[1]),pch=17,col='firebrick3',cex=2)
text(a$year2,rep(-2,dim(a)[1]),labels=a$lbl,srt=50,adj=0,cex=.8)
lines(pdat$year,pdat$y,col='dodgerblue3',lwd=2)
dev.off()




setwd(figsdir)
pdf('herring_state_plots_transform.pdf',height=8,width=10)
par(mfrow=c(2,2),mar=c(4,4,1,1))
cl0<-'lightskyblue1'
cl1<-'lightskyblue2'
cl2<-'lightskyblue3'
cl3<-'dodgerblue3'
b<-subset(dat,select=c('year','her.state','her.state.se'))
b$upr90<-b$her.state+(1.645*b$her.state.se)
b$lwr90<-b$her.state-(1.645*b$her.state.se)
b$upr95<-b$her.state+(1.96*b$her.state.se)
b$lwr95<-b$her.state-(1.96*b$her.state.se)
b$upr99<-b$her.state+(2.56*b$her.state.se)
b$lwr99<-b$her.state-(2.56*b$her.state.se)
plot(0,0,ylim=c(-1.5,1.75),las=1,xlab='Year',ylab='Herring state',xlim=c(1965,2016))
axis(1,at=seq(1965,2015,5),labels=FALSE)
axis(1,at=seq(1965,2015,10),labels=TRUE)
abline(h=0,col='lightgray')
lines(b$year,b$her.state,pch=15)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha(cl0,.75),border=NA)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr95,b$lwr95[length(b$lwr95):1]),col=alpha(cl1,.75),border=NA)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr90,b$lwr90[length(b$lwr90):1]),col=alpha(cl2,.75),border=NA)
#lines(b$year,b$her.state,ylim=c(-1.5,1.75),xlim=c(1965,2016),col=cl3,lwd=2)
asp<-aspline(b$year,b$her.state,xout=seq(min(b$year),max(b$year),length.out=1000))
lines(asp$x,asp$y,col=cl3,lwd=2)
points(b$year,b$her.state,col=alpha(cl3,.5),pch=16,cex=1)


#CHARACTERIZE TRENDS
bb<-subset(dat,select=c('year','her.state','her.state.se'))
modl<-lm(her.state~year,data=bb,weights=1/bb$her.state.se)
m<-breakpoints(her.state~year,data=bb,h=3,breaks=2)
modst<-lm(her.state~breakfactor(m,breaks=length(unique(m$breakpoints))),data=bb,weights=1/bb$her.state.se)
modnl<-lm(her.state~bs(year,degree=3,df=5),data=bb,weights=1/bb$her.state.se)
bp<-median(bb$year)
modd<-lm(her.state~year,data=bb)
modseg<-segmented(modd,seg.Z = ~year,psi=bp,control=seg.control(it.max=200),weights=1/bb$her.state.se)

dt<-data.frame(AIC(modl,modnl,modseg,modst))
names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
dt$md<-rownames(dt)
dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+2,dt$AIC)
dt1<-subset(dt,AIC==min(dt$AIC))
pdat<-data.frame(year=seq(min(b$year),max(b$year),length.out=100))
pdat2<-data.frame(year=sort(unique(b$year)))
pdat$pseg<-predict(modseg,newdata=pdat)
pdat$pnl<-predict(modnl,newdata=pdat)
pdat$pl<-predict(modl,newdata=pdat)
pdat2$pst<-predict(modst,newdata=pdat2)

plot(0,0,ylim=c(-1.25,1.75),las=1,xlab='Year',ylab='Herring state',xlim=c(1965,2016),xaxt='n')
axis(1,seq(1965,2015,5),labels=FALSE)
axis(1,seq(1965,2015,10),labels=TRUE)
abline(h=0,col='lightgray')
points(b$year,b$her.state,col=alpha('black',.76),pch=15,cex=1)
lines(pdat$year,pdat$pnl,col='black',lty=1)
#lines(pdat$year,pdat$pl,col='black',lty=2)
#lines(pdat$year,pdat$pseg,col='black',lty=2)
#lines(pdat2$year,pdat2$pst,col='black',lty=3)
#lines(pdat$year,pdat$pseg,lty=3,col='red')
rlm<-round(summary(modl)$r.sq,digits=2)
rseg<-round(summary(modseg)$r.sq,digits=2)
rst<-round(summary(modst)$r.sq,digits=2)
rnl<-round(summary(modnl)$r.sq,digits=2)
legend('top',c(paste('Spline (r2=',rnl,')'),paste('Linear (r2=',rlm,')'),paste('Structural (r2=',rst,')')),bty='n',lty=c(1,2,3))

dev.off()


xx<-subset(data,select=c('year','her.expr','her.land','sst.state','sst.mn','had.pi','her.totno.rv','her.totwgt.rv','her.expr','nut.state','ct.state'))
b<-subset(dat,select=c('year','her.state','her.state.se'))
a<-merge(xx,b,by=c('year'),all=FALSE)

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))
n<-length(unique(a$year))
dum3<-data.frame(year=sort(unique(a$year)),
                 cls=cls(n+2)[3:(n+2)])
a<-merge(a,dum3,by=c('year'),all.x=TRUE,all.y=FALSE)

#plot(a,pch=16)
par(mfrow=c(2,2),mar=c(4,6,1,4))
plot(a$sst.mn,a$her.state,pch=16,col=as.character(a$cl),cex=2,las=1,xlab='SST',ylab='Herring state',ylim=c(-1,1),xaxt='n',xlim=c(10.5,14.5))
points(a$sst.mn,a$her.state,pch=1,cex=2,lwd=.1)
axis(1,at=seq(10,15,1))
colorbar.plot(13,.5,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~sst.mn,data=a)
mod2<-lm(her.state~sst.mn +I(sst.mn^2),data=a)
pdat<-data.frame(sst.mn=seq(min(a$sst.mn,na.rm=TRUE),max(a$sst.mn,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$sst.mn,pdat$p,lty=2)


plot(a$her.expr,a$her.state,pch=16,col=alpha(as.character(a$cl),.75),cex=2,las=1,xlab='Exploitation rate',ylab='Herring state')
points(a$her.expr,a$her.state,pch=1,cex=2,lwd=1,col=as.character(a$cl))
colorbar.plot(1,1,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~her.expr,data=a)
mod2<-lm(her.state~her.expr +I(her.expr^2),data=a)
pdat<-data.frame(her.expr=seq(min(a$her.expr,na.rm=TRUE),max(a$her.expr,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$her.expr,pdat$p,lty=2)

a$her.land<-a$her.land/1000
plot(a$her.land,a$her.state,pch=16,col=alpha(as.character(a$cl),.75),cex=2,las=1,xlab='Landings',ylab='Herring state',xlim=c(48,200),xaxt='n')
axis(1,seq(50,200,10),labels=FALSE)
axis(1,seq(50,200,50),labels=TRUE)
points(a$her.land,a$her.state,pch=1,cex=2,lwd=1,col=as.character(a$cl))
colorbar.plot(150000,-1,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~her.land,data=a)
mod2<-lm(her.state~her.land +I(her.land^2),data=a)
pdat<-data.frame(her.land=seq(min(a$her.land,na.rm=TRUE),max(a$her.land,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$her.land,pdat$p,lty=2)

dev.off()






























setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)
rvw<-subset(rvw,month %in% c(6,7,8))
plg<-subset(plg,stratum %in% unique(rvw$strat))
plg$stratum<-as.numeric(as.character(plg$stratum))
plgj<-subset(plg,stratum %in% c('493','494'))
plga<-subset(plg,!(stratum %in% c('493','494')))
plg2<-subset(plga,stratum>=480 & stratum<=495)

setwd(figsdir)
pdf('rv_survey_strata.pdf',height=8,width=8)
par(mar=c(4,4,1,1))
map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-56),ylim=c(41,48))
plot(plga,add=TRUE,col=alpha('dodgerblue3',.1),border='dodgerblue3',lwd=.001,xlim=c(-68,-56),ylim=c(41,48))

map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-56),ylim=c(41,48))
plot(plg2,add=TRUE,col=alpha('firebrick2',.1),border='firebrick3',lwd=.001)

map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-56),ylim=c(41,48))
plot(plgj,add=TRUE,col=alpha('forestgreen',.1),border='forestgreen',lwd=.001)

map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-56),ylim=c(41,48))
plot(plg,add=TRUE,col=alpha('dodgerblue3',.1),border='dodgerblue3',lwd=.001,xlim=c(-68,-56),ylim=c(41,48))
plot(plg2,add=TRUE,col=alpha('firebrick2',.1),border='firebrick3',lwd=.001)

map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-60),ylim=c(43,47),axes=FALSE)
axis(1,seq(-68,-60,.5),labels=FALSE)
axis(1,seq(-68,-60,1),labels=TRUE)
axis(2,seq(43,47,.5),labels=FALSE,las=1)
axis(2,seq(43,47,1),labels=TRUE,las=1)
box()

map('worldHires',col='gray',fill=TRUE,border=NA,xlim=c(-68,-60),ylim=c(43,47),axes=FALSE)
dev.off()




setwd(datadir)
load('SPERA_andata.RData')
nms<-names(data)
nms<-nms[grepl('her',nms)==TRUE]
nms<-nms[grepl('\\.se',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('georng',nms)==FALSE]
nms<-nms[grepl('her.durng.rv',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-subset(df,year>=1965)

dmat<-1-(abs(cor(df,use='pairwise.complete.obs')))+.001
#dmat<-abs(cor(df,use='pairwise.complete.obs'))
dst<-as.dist(dmat)
#dst<-dist(df,method='euclidean')
ord<-metaMDS(dst,trymax=50,method='euclidean',trace=0)
x1<--.5
x2<-.5
plot(0,0,col='white',xlim=c(x1,x2),ylim=c(x1,x2),las=1)
text(ord, display = "sites", cex = 0.8, pch=21,  bg="yellow",adj=-.2)
points(ord, display = "sites", cex =1, pch=16)
points(ord, display = "sites", cex =1,,pch=1,lwd=.1)
cldb<-unique(subset(cldb,select=c('cl','lab')))
legend('topleft',cldb$lab,col=as.character(cldb$cl),pch=16,bty='n',pt.cex=2)
legend('bottomleft',paste(yr1,'-',yr2),bty='n')











setwd(datadir)
load('SPERA_andata.RData')
d<-subset(dher,select=c('year','quantity'))
names(d)[2]<-'her.land2'
data<-merge(data,d,by=c('year'),all.x=TRUE,all.y=FALSE)



mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N:\\cluster_2017\\scratch\\chl_phenology\\data\\naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fcoast.mc<-fortify(coast.mc)

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }




library(tidyr)
###################################################

#MAPS OF SPAWNING LOCATIONS
setwd(datadir)
spawn<-read.csv('spawning_locations_historical.csv',header=TRUE)
spawn$season<-gsub(' ','',as.character(spawn$season))
spawn$lon<-as.numeric(as.character(spawn$lon))
names(spawn)[3]<-'hst'
spawn$hst<-ifelse(is.na(spawn$hst)==TRUE,0,spawn$hst)
spawn<-subset(spawn,select=c('lon','lat','season'))
spawn$cl<-ifelse(spawn$season %in% c('s','es'),'dodgerblue3','red')
spawn$cl<-ifelse(spawn$season=='s/f','orange',spawn$cl)
spawn$cl<-ifelse(spawn$season=='f','firebrick3',spawn$cl)
spawn$lon<-round(spawn$lon,digits=2)
spawn$lat<-round(spawn$lat,digits=2)
spawn<-unique(spawn)
spawn$pch<-ifelse(spawn$season %in% c('s','es'),17,16)
spawn$pch<-ifelse(spawn$season=='s/f',15,spawn$pch)
spawn$pch<-ifelse(spawn$season=='f',16,spawn$pch)
spawn$pch2<-ifelse(spawn$season %in% c('s','es'),2,1)
spawn$pch2<-ifelse(spawn$season=='s/f',0,spawn$pch2)
spawn$pch2<-ifelse(spawn$season=='f',1,spawn$pch2)
spawn$seas2<-ifelse(spawn$season %in% c('s','es'),'spring','fall')
spawn$seas2<-ifelse(spawn$season=='s/f','both',spawn$seas2)
spawn$latc<-cut(spawn$lat,breaks=seq(43,51,.5),labels=seq(43.25,50.75,.5))
spawn$lonc<-cut(spawn$lon,breaks=seq(-67,-52,.5),labels=seq(-66.75,-52.25,.5))

setwd(figsdir)
pdf('spawning_locs.pdf',height=8,width=8)
par(mar=c(4,4,1,1))
map('worldHires',xlim=c(-67,-59.5),ylim=c(43,47),fill=TRUE,col='gray',border=NA)
axis(1,seq(-67,-59,.25),labels=FALSE)
axis(1,seq(-67,-59,1),labels=TRUE)
axis(2,seq(43,47,.25),labels=FALSE)
axis(2,seq(43,47,1),labels=TRUE,las=1)
points(spawn$lon,spawn$lat,col=alpha(as.character(spawn$cl),.3),pch=spawn$pch,cex=1.5)
points(spawn$lon,spawn$lat,col=alpha(as.character(spawn$cl),.9),pch=spawn$pch2,cex=1.5,lwd=.5)
aa<-unique(subset(spawn,select=c('seas2','cl','pch')))
legend('bottomright',legend=aa$seas2,col=as.character(aa$cl),bty='n',pch=aa$pch)
box(lwd=.1)


f<-function(d){    return(data.frame(n=dim(d)[1]))}
ot<-ddply(spawn,.(seas2,latc),.fun=f)

aa<- ot %>% spread(latc, n)
aa[is.na(aa)] <- 0
aaa<-aa[,2:16]
par(mar=c(12,10,12,10))
barplot(as.matrix(aaa),las=1,col=c('orange','firebrick3','dodgerblue3'),names.arg=names(aaa),horiz=TRUE,border=NA,xlab='Number of sites')
aa<-unique(subset(spawn,select=c('seas2','cl','pch')))
legend('topright',legend=aa$seas2,col=as.character(aa$cl),bty='n',pch=15)

f<-function(d){    return(data.frame(n=dim(d)[1]))}
ot<-ddply(spawn,.(seas2,lonc),.fun=f)

aa<- ot %>% spread(lonc, n)
aa[is.na(aa)] <- 0
aaa<-aa[,2:16]
par(mar=c(12,10,12,10))
barplot(as.matrix(aaa),las=1,col=c('orange','firebrick3','dodgerblue3'),names.arg=names(aaa),horiz=FALSE,border=NA,xlab='Number of sites')
aa<-unique(subset(spawn,select=c('seas2','cl','pch')))
legend('topright',legend=aa$seas2,col=as.character(aa$cl),bty='n',pch=15)
dev.off()










setwd(datadir)
load('SPERA_andata_new.RData')
data$her.durng.rv<-ifelse(data$her.durng.rv==-Inf,NA,data$her.durng.rv)
data$herjuv.durng.rv<-ifelse(data$herjuv.durng.rv==-Inf,NA,data$herjuv.durng.rv)

#LOAD LANDINGS FROM FAO
d<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/TS_FI_CAPTURE.csv',header=TRUE)
spc<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/CL_FI_SPECIES_GROUPS.csv',header=TRUE)
names(d)<-tolower(names(d))
names(spc)<-tolower(names(spc))
names(spc)[1]<-'species'
dher<-subset(d,species=='HER' & country==124)

d<-subset(dher,select=c('year','quantity'))
names(d)[2]<-'her.land2'
data<-merge(data,d,by=c('year'),all.x=TRUE,all.y=FALSE)


nms<-names(data)
nms<-nms[grepl('her',nms)==TRUE]
nms<-nms[grepl('\\.se',nms)==FALSE]
dat<-data[,names(data) %in% c('year',nms)]
dat<-subset(dat,year>=1965)
dat<-merge(dat,ds,by=c('year'),all.x=TRUE,all.y=FALSE)

dat<-dat[,!(names(dat) %in% c('her.land2','her.expr','her.land','her.preytot.bof','herjuv.prey.bof','her.jf.bof','her.tbio','her.age5.cat','her.agepe.cat','her.age.cat','her.cf.cat','her.agediv.cat','her.dep.rv','herlrv.dep.bof','herjuv.dep.rv','her.dumx.rv','herlrv.mn.bof','herjuv.sp.n','herlrv.surv','herjuv.spcv','herjuv.spnug','herjuv.spvar','herjuv.sp.n','herjuv.georng','herjuv.sprng','her.totwgt.rv','her.totno.rv','herjuv.totno.rv','herjuv.totwgt.rv','herjuv.dumx.rv','herjuv.durng.rv','herjuv.len.rv','herlrv.georng','herlrv.sprng','herlrv.spnug','her.sprng','herlrv.spcv','her.durng.rv'))]

dat<-subset(data,year>=1965,select=c('year','herlrv.len','her.spnug','her.spcv','her.szdiv.rv','her.szpe.rv','her.len.rv','her.fmass.rv','her.waa','her.spvar','her.cf.rv','her.ssb','her.metai.rv','her.prod','her.rec1'))
#d$her.rec1<-log10(d$her.rec1)
#d$her.land<-log10(d$her.land)
#d$her.land2<-log10(d$her.land2)
dat$her.szdiv.rv<-dat$her.szdiv.rv*-1
dat$her.spcv<-dat$her.spcv*-1
dat$her.spvar<-dat$her.spvar*-1
dat$her.spnug<-dat$her.spnug*-1
dat$her.rec1<-sqrt(dat$her.rec1)
dat$her.prod<-sqrt(dat$her.prod+abs(min(dat$her.prod,na.rm=TRUE))+1)

abs(round(cor(dat,use='pairwise.complete.obs'),digits=2))
(round(cor(dat,use='pairwise.complete.obs'),digits=2))
par(mar=c(4,4,1,1))
plot(dat,pch=16,cex=.75)

plot(dat$year,dat$her.durng.rv,pch=15)
plot(dat$year,dat$her.dep.rv,pch=15)

############################################################


##HERRING STATE DERIVATION USING ONLY DATA FROM SPAWNING AREA
#################################################

######### HERRING STATE DERIVATION
setwd(datadir)
load('SPERA_andata_spawnar.RData')
data<-data.spawnar
data$her.durng.rv<-ifelse(data$her.durng.rv==-Inf,NA,data$her.durng.rv)

#LOAD LANDINGS FROM FAO
d<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/TS_FI_CAPTURE.csv',header=TRUE)
spc<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/CL_FI_SPECIES_GROUPS.csv',header=TRUE)
names(d)<-tolower(names(d))
names(spc)<-tolower(names(spc))
names(spc)[1]<-'species'
dher<-subset(d,species=='HER' & country==124)

d<-subset(dher,select=c('year','quantity'))
names(d)[2]<-'her.land2'
data<-merge(data,d,by=c('year'),all.x=TRUE,all.y=FALSE)


#SAME AS ABOVE BUT EXCLUDES VPA SSB SERIES
#OMIT: HERLRV.SPCV, HERJUV.PREY.BOF, HERLRV.SPVAR, HERJUV.SPRNG, HER.TOTNO.RV, HER.TOTWGT.RV,'her.durng.rv','her.dumx.rv','her.dep.rv',her.georng,her.prod,'herlrv.dep.bof','herlrv.surv','herlrv.mn.bof'
#LANDINGS, EXPLOITATION RATE CORRELATED BUT NOT STATE VARS
#dat<-subset(data,year>=1965,select=c('year','her.cf.rv','her.spvar','her.fmass.rv','her.len.rv','her.metai.rv','her.rec1','her.spcv','her.ssb','her.szpe.rv','her.waa','herlrv.len','her.land','her.land2','her.totwgt.rv','her.state','her.lratio'))
dat<-subset(data,year>=1965,select=c('year','herlrv.len','her.spnug','her.spcv','her.szdiv.rv','her.szpe.rv','her.len.rv','her.fmass.rv','her.waa','her.spvar','her.cf.rv','her.ssb','her.metai.rv','her.prod','her.rec1'))
#d$her.rec1<-log10(d$her.rec1)
#d$her.land<-log10(d$her.land)
#d$her.land2<-log10(d$her.land2)
dat$her.szdiv.rv<-dat$her.szdiv.rv*-1
dat$her.spcv<-dat$her.spcv*-1
dat$her.spvar<-dat$her.spvar*-1
dat$her.spnug<-dat$her.spnug*-1
dat$her.rec1<-sqrt(dat$her.rec1)
dat$her.prod<-sqrt(dat$her.prod+abs(min(dat$her.prod,na.rm=TRUE))+1)

abs(round(cor(dat,use='pairwise.complete.obs'),digits=2))
(round(cor(dat,use='pairwise.complete.obs'),digits=2))
plot(dat,pch=16)


library(tidyr)
#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
dat<-data.frame(cbind(subset(dat,select='year'),apply(dat[,2:dim(dat)[2]],2,f)))

dat<-dat[,!(names(dat) %in% c('her.land2','her.state'))]
bb<-dat %>% gather(var,y,herlrv.len:her.rec1)
#bb<-na.omit(bb)
xyplot(y~year | var,data=bb,pch=15,type=c('p','l'))
histogram(~y | var,data=bb)

mod<-gamm(y~as.factor(year),data=bb,random=list(var=~1))
pdat<-data.frame(year=sort(unique(bb$year)))
p<-predict(mod$gam,newdata=pdat,se.fit=TRUE)
pdat$her.state<-p$fit
pdat$her.state.se<-p$se.fit
plot(pdat$year,pdat$her.state)

dat<-merge(dat,pdat,by=c('year'),all=TRUE)

cr<-cor(dat,use='pairwise.complete.obs',method='spearman')
cr.r<-round(cr,digits=2)
cr.t = as.data.frame(as.table(cr.r))#CORRELATIONS TABLE
combinations=combn(colnames(cr.r),2,FUN=function(x){paste(x,collapse="_")})
cr.t=cr.t[cr.t$Var1 != cr.t$Var2,]
cr.t=cr.t[paste(cr.t$Var1,cr.t$Var2,sep="_") %in% combinations,]
cr.t<-subset(cr.t,!(Var2 %in% c('year','her.state.se')) & !(Var1 %in% c('year','her.state.se')))

cr.t<-cr.t[order(cr.t$Freq,decreasing=TRUE),]
cr.t$Var1<-as.character(cr.t$Var1)
cc<-data.frame(var1=sort(unique(cr.t$Var1)),
               r=tapply(cr.t$Freq,cr.t$Var1,mean))

#AVERAGE CORRELATION OF 'STATE' VARIABLES
cc<-subset(cr.t,Var1!='her.state' & Var2!='her.state')
mean(cc$Freq);#.45
median(cc$Freq);#.48
#mdr<-round(median(cc$Freq),digits=2)
#MEAN=0.54; MEDIAN=0.56


cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))
n<-length(unique(bb$var))
dum3<-data.frame(var=sort(unique(bb$var)),
                 cls=cls(n+2)[3:(n+2)])
bb<-merge(bb,dum3,by=c('var'),all.x=TRUE,all.y=FALSE)

setwd(figsdir)
pdf('herring_state_derivation_spawnar.pdf',height=7,width=4)
par(mfcol=c(7,2),mar=c(.5,1,0,0),oma=c(4,4,1,4))
f<-function(d){
d<-na.omit(d)
ylm<-c(floor(min(d$y,na.rm=TRUE)),ceiling(max(d$y,na.rm=TRUE)))
print(ylm)
asp<-aspline(d$year,d$y,xout=seq(min(d$year),max(d$year),length.out=1000))
plot(0,0,ylim=ylm,xlim=c(1965,2015),axes=FALSE)
abline(h=0,col='lightgray',lty=2,lwd=.5)
lines(asp$x,asp$y,col=alpha(as.character(unique(d$cl)),1),lwd=2)
points(d$year,d$y,col=alpha(as.character(unique(d$cl)),.3),cex=1,pch=16)
axis(2,at=ylm,las=1,cex.axis=.5,lwd=.1)
if(unique(d$var) %in% c('her.spcv','herlrv.len')){
axis(1,seq(1965,2015,5),labels=FALSE,cex.axis=.5,lwd=.1)
axis(1,seq(1965,2015,10),cex.axis=.5,lwd=.1)
} else NULL
mod<-lm(y~year,data=d)
pdat<-data.frame(year=seq(min(d$year),max(d$year),length.out=1000))
pdat$p<-predict(mod,newdata=pdat)
#lines(pdat$year,pdat$p,col=alpha('black',.5),lwd=.5)
s<-summary(mod)
b<-round(s$coef[2,1],digits=2)
legend('topright',gsub(' ','',paste(unique(d$var),'=',b)),bty='n',cex=.5)
}
z<-dlply(bb,.(var),.fun=f)



bb$segtrue<-ifelse(bb$var%in% c('her.len.rv'),FALSE,TRUE)

par(mfcol=c(7,2),mar=c(.5,1,0,0),oma=c(4,4,1,4))
f<-function(d){
d<-na.omit(d)
print(unique(d$var))
#CHARACTERIZE TRENDS
    if(length(unique(d$year))<=10){dff<-4
    } else {dff<-5
    }
    bp<-floor(mean(d$year))
    modl<-lm(y~year,data=d)
    m<-breakpoints(y~year,data=d,h=3,breaks=1)
    modst<-lm(y~breakfactor(m,breaks=length(unique(m$breakpoints))),data=d)
    modnl<-lm(y~bs(year,degree=3,df=dff),data=d)
    if(unique(d$segtrue)==TRUE){
        modseg<-segmented(modl,seg.Z = ~year,psi=bp,control=seg.control(it.max=200))
        dt<-data.frame(AIC(modl,modnl,modseg,modst))
    } else {rm(modseg)
            dt<-data.frame(AIC(modl,modnl,modst))
    }
    names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
    dt$md<-rownames(dt)
    dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+4,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modst',dt$AIC+2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modseg',dt$AIC-2,dt$AIC)
    dt$AIC<-ifelse(dt$md=='modl',dt$AIC-2,dt$AIC)
    dt<-subset(dt,AIC==min(dt$AIC))

if(dt$md=='modnl'){
    modl<-modnl
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
} else if (dt$md=='modseg'){
    modl<-modseg
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
} else if (dt$md=='modst'){
    modl<-modst
    d2<-data.frame(year=sort(unique(d$year)))
} else {
    modl<-modl
    d2<-data.frame(year=seq(min(d$year),max(d$year),length.out=100))
}

    s<-summary(modl)
    r2<-round(s$r.squared,digits=2)
    p<-predict(modl,newdata=d2,se.fit=TRUE,type='response')
    d2$p<-p$fit
    d2$se<-p$se.fit
    d2$upr<-d2$p+(1.96*d2$se)
    d2$lwr<-d2$p-(1.96*d2$se)
    ylm<-c(min(c(d$y,d2$lwr)),max(c(d$y,d2$upr)))

ylm<-c(floor(min(d$y,na.rm=TRUE)),ceiling(max(d$y,na.rm=TRUE)))
print(ylm)
plot(0,0,ylim=ylm,xlim=c(1965,2015),axes=FALSE)
polygon(c(d2$year,d2$year[length(d2$year):1]),c(d2$upr,d2$lwr[length(d2$lwr):1]),col=alpha(as.character(unique(d$cl)),.3),border=NA)
abline(h=0,col='lightgray',lty=2,lwd=.5)
points(d$year,d$y,col=alpha(as.character(unique(d$cl)),.3),cex=1,pch=16)
axis(2,at=ylm,las=1,cex.axis=.5,lwd=.1)
if(unique(d$var) %in% c('her.spcv','herlrv.len')){
axis(1,seq(1965,2015,5),labels=FALSE,cex.axis=.5,lwd=.1)
axis(1,seq(1965,2015,10),cex.axis=.5,lwd=.1)
} else NULL
lines(d2$year,d2$p,col=as.character(unique(d$cl)),lwd=2)
s<-summary(modl)
r2<-round(s$r.squared,digits=2)
legend('topright',gsub(' ','',paste(unique(d$var),'=',r2)),bty='n',cex=.5)
}
z<-dlply(bb,.(var),.fun=f)
dev.off()






library(segmented)
library(splines)
library(strucchange)
setwd(figsdir)
pdf('herring_state_plots_spawnar.pdf',height=8,width=10)
par(mfrow=c(2,2),mar=c(4,4,1,1))
cl0<-'lightskyblue1'
cl1<-'lightskyblue2'
cl2<-'lightskyblue3'
cl3<-'dodgerblue3'
b<-subset(dat,select=c('year','her.state','her.state.se'))
b$upr90<-b$her.state+(1.645*b$her.state.se)
b$lwr90<-b$her.state-(1.645*b$her.state.se)
b$upr95<-b$her.state+(1.96*b$her.state.se)
b$lwr95<-b$her.state-(1.96*b$her.state.se)
b$upr99<-b$her.state+(2.56*b$her.state.se)
b$lwr99<-b$her.state-(2.56*b$her.state.se)
plot(0,0,ylim=c(-1.5,1.75),las=1,xlab='Year',ylab='Herring state',xlim=c(1965,2016))
axis(1,at=seq(1965,2015,5),labels=FALSE)
axis(1,at=seq(1965,2015,10),labels=TRUE)
abline(h=0,col='lightgray')
lines(b$year,b$her.state,pch=15)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr99,b$lwr99[length(b$lwr99):1]),col=alpha(cl0,.75),border=NA)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr95,b$lwr95[length(b$lwr95):1]),col=alpha(cl1,.75),border=NA)
polygon(c(b$year,b$year[length(b$year):1]),c(b$upr90,b$lwr90[length(b$lwr90):1]),col=alpha(cl2,.75),border=NA)
#lines(b$year,b$her.state,ylim=c(-1.5,1.75),xlim=c(1965,2016),col=cl3,lwd=2)
asp<-aspline(b$year,b$her.state,xout=seq(min(b$year),max(b$year),length.out=1000))
lines(asp$x,asp$y,col=cl3,lwd=2)
points(b$year,b$her.state,col=alpha(cl3,.5),pch=16,cex=1)



#CHARACTERIZE TRENDS
bb<-subset(dat,select=c('year','her.state','her.state.se'))
modl<-lm(her.state~year,data=bb,weights=1/bb$her.state.se)
m<-breakpoints(her.state~year,data=bb,h=3,breaks=2)
modst<-lm(her.state~breakfactor(m,breaks=length(unique(m$breakpoints))),data=bb,weights=1/bb$her.state.se)
modnl<-lm(her.state~bs(year,degree=3,df=5),data=bb,weights=1/bb$her.state.se)
bp<-median(bb$year)
modd<-lm(her.state~year,data=bb)
modseg<-segmented(modd,seg.Z = ~year,psi=bp,control=seg.control(it.max=200),weights=1/bb$her.state.se)

dt<-data.frame(AIC(modl,modnl,modseg,modst))
names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
dt$md<-rownames(dt)
dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+2,dt$AIC)
dt1<-subset(dt,AIC==min(dt$AIC))
pdat<-data.frame(year=seq(min(b$year),max(b$year),length.out=100))
pdat2<-data.frame(year=sort(unique(b$year)))
pdat$pseg<-predict(modseg,newdata=pdat)
pdat$pnl<-predict(modnl,newdata=pdat)
pdat$pl<-predict(modl,newdata=pdat)
pdat2$pst<-predict(modst,newdata=pdat2)

plot(0,0,ylim=c(-1.25,1.75),las=1,xlab='Year',ylab='Herring state',xlim=c(1965,2016),xaxt='n')
axis(1,seq(1965,2015,5),labels=FALSE)
axis(1,seq(1965,2015,10),labels=TRUE)
abline(h=0,col='lightgray')
points(b$year,b$her.state,col=alpha('black',.76),pch=15,cex=1)
lines(pdat$year,pdat$pnl,col='black',lty=1)
lines(pdat$year,pdat$pl,col='black',lty=2)
lines(pdat2$year,pdat2$pst,col='black',lty=3)
#lines(pdat$year,pdat$pseg,lty=3,col='red')
rlm<-round(summary(modl)$r.sq,digits=2)
rseg<-round(summary(modseg)$r.sq,digits=2)
rst<-round(summary(modst)$r.sq,digits=2)
rnl<-round(summary(modnl)$r.sq,digits=2)
legend('top',c(paste('Spline (r2=',rnl,')'),paste('Linear (r2=',rlm,')'),paste('Structural (r2=',rst,')')),bty='n',lty=c(1,2,3))




xx<-subset(data,select=c('year','her.expr','her.land','sst.state','sst.mn','had.pi','her.totno.rv','her.totwgt.rv','her.expr','nut.state','ct.state'))
b<-subset(dat,select=c('year','her.state','her.state.se'))
a<-merge(xx,b,by=c('year'),all=FALSE)

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))
n<-length(unique(a$year))
dum3<-data.frame(year=sort(unique(a$year)),
                 cls=cls(n+2)[3:(n+2)])
a<-merge(a,dum3,by=c('year'),all.x=TRUE,all.y=FALSE)

#plot(a,pch=16)
par(mfrow=c(2,2),mar=c(4,6,1,4))
plot(a$sst.mn,a$her.state,pch=16,col=as.character(a$cl),cex=2,las=1,xlab='SST',ylab='Herring state',ylim=c(-1,1),xaxt='n',xlim=c(10.5,14.5))
points(a$sst.mn,a$her.state,pch=1,cex=2,lwd=.1)
axis(1,at=seq(10,15,1))
colorbar.plot(13,.5,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~sst.mn,data=a)
mod2<-lm(her.state~sst.mn +I(sst.mn^2),data=a)
pdat<-data.frame(sst.mn=seq(min(a$sst.mn,na.rm=TRUE),max(a$sst.mn,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$sst.mn,pdat$p,lty=2)


plot(a$her.expr,a$her.state,pch=16,col=alpha(as.character(a$cl),.75),cex=2,las=1,xlab='Exploitation rate',ylab='Herring state')
points(a$her.expr,a$her.state,pch=1,cex=2,lwd=1,col=as.character(a$cl))
colorbar.plot(1,1,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~her.expr,data=a)
mod2<-lm(her.state~her.expr +I(her.expr^2),data=a)
pdat<-data.frame(her.expr=seq(min(a$her.expr,na.rm=TRUE),max(a$her.expr,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$her.expr,pdat$p,lty=2)

a$her.land<-a$her.land/1000
plot(a$her.land,a$her.state,pch=16,col=alpha(as.character(a$cl),.75),cex=2,las=1,xlab='Landings',ylab='Herring state',xlim=c(48,200),xaxt='n')
axis(1,seq(50,200,10),labels=FALSE)
axis(1,seq(50,200,50),labels=TRUE)
points(a$her.land,a$her.state,pch=1,cex=2,lwd=1,col=as.character(a$cl))
colorbar.plot(150000,-1,dum3$year,col=as.character(dum3$cl),horizontal=TRUE,strip.width=.04,strip.length=.4)
mod<-lm(her.state~her.land,data=a)
mod2<-lm(her.state~her.land +I(her.land^2),data=a)
pdat<-data.frame(her.land=seq(min(a$her.land,na.rm=TRUE),max(a$her.land,na.rm=TRUE),length.out=100))
if(AIC(mod)<=AIC(mod2)){pdat$p<-predict(mod,newdata=pdat)
                    } else { pdat$p<-predict(mod2,newdata=pdat)
                         }
lines(pdat$her.land,pdat$p,lty=2)

dev.off()


















setwd(figsdir)
pdf('length_v_abundance.pdf',height=8,width=6)
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(data$her.len.rv,data$her.totwgt.rv,pch=15,las=1,xlab='Herring length',ylab='Herring weight per tow',xlim=c(20,35))
cr<-round(cor(data$her.len.rv,data$her.totwgt.rv,use='pairwise.complete.obs',method='spearman'),digits=2)
legend('topright',legend=paste('r =',cr),bty='n')
plot(data$her.len.rv,data$her.totno.rv,pch=15,las=1,xlab='Herring length',ylab='Herring number per tow',xlim=c(20,35))
cr<-round(cor(data$her.len.rv,data$her.totno.rv,use='pairwise.complete.obs',method='spearman'),digits=2)
legend('topright',legend=paste('r =',cr),bty='n')
dev.off()

plot(data$her.ssb,log10(data$herlrv.mn.bof),pch=15,col='white')
text(data$her.ssb,log10(data$herlrv.mn.bof),labels=data$year)
setwd(datadir)
load('phen2.RData')
load('phen3.RData')
load('stdat.RData')
load('rvw2.RData')
load('ldat.RData')
load("BoF_larval_herring_lengths_spera_spawnar.RData");larvs<-dat


mcrt<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coast<-readOGR('N:\\cluster_2017\\scratch\\chl_phenology\\data\\naturalearthdata_ne_10m_land_poly',layer='ne_10m_land')#works for Alfcoast<-fortify(coast)
coast.mc<-crop(coast,extent(-180,180,-90,90),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
fcoast.mc<-fortify(coast.mc)

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }






#########################################################
library(akima)

#PLOTS HERRINGS SSB FROM ACOUSTIC VERSUS VP
setwd(datadir1)
her<-read.csv('herring_assess_ssb_r_f_w_spera.csv',header=TRUE,na.strings=c('- '))
names(her)<-tolower(names(her))
#names(her)<-gsub('...','.',names(her))
names(her)[56]<-'land2016'

her$lratio<-(her$x4wx.stock.nominal.landings/her$x4wx.stock.tac)*100
plot(her$year,her$lratio,pch=15,type='b')
abline(h=100,lty=2)

#LANDINGS ARE FOR ALL AGES(1-11) AND SSB IS ONLY FOR AGES 4-8: THIS IS WHY EXPLOITATION RATE EXCEEDS 1 IN MANY CASES
her$ssb<-her$ssb/1000
her$landings<-her$landings/1000
her$acoustic<-her$acoustic/1000
her$expr<-her$landings/her$ssb
her$expr.ac<-her$landings/her$acoustic

cl1<-'darkblue'
cl2<-'firebrick3'
cl3<-'black'

setwd(figsdir)
pdf('VPA_vs_acoustic.pdf',height=8,width=11)
par(mfrow=c(2,2),mar=c(4,4,1,1))
#PLOTS VPA AND ACOUSTIC SSB ESTIMATES; CALIBRATES ACOUSTIC
d1<-na.omit(subset(her,select=c('year','ssb')))
d2<-na.omit(subset(her,select=c('year','acoustic')))
#plot(her$year,her$ssb,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB',col=cl1,xaxt='n',ylim=c(0,800))
asp<-aspline(her$year,her$ssb,xout=seq(min(d1$year),max(d1$year),length.out=1000))
asp2<-aspline(her$year,her$acoustic,xout=seq(min(d2$year),max(d2$year),length.out=1000))
plot(asp$x,asp$y,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB',col=cl1,xaxt='n',ylim=c(0,800))

points(her$year,her$ssb,pch=16,xlim=c(1965,2016),col=alpha(cl1,.3))
points(asp2$x,asp2$y,pch=16,col=cl2,type='l')
points(her$year,her$acoustic,pch=16,col=alpha(cl2,.3))
legend('topleft',legend=c('VPA (ages 4-8)','Acoustic'),col=c(cl1,cl2),lwd=2,bty='n')
axis(1,seq(1965,2015,5))
mod1<-lm(ssb~year,data=d1)
pdat<-data.frame(year=seq(min(d1$year),max(d1$year),length.out=100))
pdat$p<-predict(mod1,newdata=pdat)
lines(pdat$year,pdat$p,col=cl1,lty=2)
mod2<-lm(acoustic~year,data=d2)
pdat2<-data.frame(year=seq(min(d2$year),max(d2$year),length.out=100))
pdat2$p<-predict(mod2,newdata=pdat2)
lines(pdat2$year,pdat2$p,col=cl2,lty=2)
abline(h=0,col='gray')


rpr<-round(cor(her$ssb,her$acoustic,use='pairwise.complete.obs'),digits=2)
rknd<-round(cor(her$ssb,her$acoustic,use='pairwise.complete.obs',method='spearman'),digits=2)
plot(her$ssb,her$acoustic,pch=16,las=1,xlab='SSB from VPA (2005 assessment)',ylab='SSB from acoustic additive model',xlim=c(25,140),ylim=c(200,700))
legend('bottomright',legend=c(paste('Pearson r =', rpr),paste('Spearman r =',rknd)),bty='n')
d2<-na.omit(subset(her,select=c('ssb','acoustic')))
mod2<-lm(acoustic~ssb,data=d2)
pdat2<-data.frame(ssb=seq(min(d2$ssb),max(d2$ssb),length.out=100))
pdat2$p<-predict(mod2,newdata=pdat2)
lines(pdat2$ssb,pdat2$p,col=cl3,lty=2)


#PREDICT SSB FROM ACOUSTIC- ADJUST ACOUSTIC DATA DOWNWARD TO MATCH VPA
mod<-lm(ssb~acoustic,data=her)
her$acoustic2<-predict(mod,newdata=data.frame(acoustic=her$acoustic))
d1<-na.omit(subset(her,select=c('year','ssb')))
d2<-na.omit(subset(her,select=c('year','acoustic2')))
#plot(her$year,her$ssb,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB',col=cl1,xaxt='n',ylim=c(0,800))
asp<-aspline(her$year,her$ssb,xout=seq(min(d1$year),max(d1$year),length.out=1000))
asp2<-aspline(her$year,her$acoustic2,xout=seq(min(d2$year),max(d2$year),length.out=1000))
plot(asp$x,asp$y,pch=16,xlim=c(1965,2016),type='l',xlab='Year',ylab='SSB',col=cl1,ylim=c(0,800),xaxt='n',las=1)
points(her$year,her$ssb,pch=16,xlim=c(1965,2016),col=alpha(cl1,.3))
points(asp2$x,asp2$y,pch=16,col=cl2,type='l')
points(her$year,her$acoustic2,pch=16,col=alpha(cl2,.3))
legend('topleft',legend=c('VPA (ages 4-8)','Acoustic calibrated to VPA'),col=c(cl1,cl2),lwd=2,bty='n')
axis(1,seq(1965,2015,5))
abline(h=0,col='gray')

#PREDICT ACOUSTIC FROM SSB- ADJUST SSB DATA UPWARD TO MATCH ACOUSTIC
mod<-lm(acoustic~ssb,data=her)
her$ssb2<-predict(mod,newdata=data.frame(ssb=her$ssb))
d1<-na.omit(subset(her,select=c('year','ssb2')))
d2<-na.omit(subset(her,select=c('year','acoustic')))
asp<-aspline(her$year,her$ssb2,xout=seq(min(d1$year),max(d1$year),length.out=1000))
asp2<-aspline(her$year,her$acoustic,xout=seq(min(d2$year),max(d2$year),length.out=1000))

par(mar=c(4,6,1,1))
plot(asp$x,asp$y,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB',col=cl1,xaxt='n',ylim=c(0,3000))
points(her$year,her$ssb2,pch=16,xlim=c(1965,2016),col=alpha(cl1,.3))
points(asp2$x,asp2$y,pch=16,col=cl2,type='l')
points(her$year,her$acoustic,pch=16,col=alpha(cl2,.3))
legend('topleft',legend=c('VPA (ages 4-8) calibrated to acoustic','Acoustic'),col=c(cl1,cl2),lwd=2,bty='n')
abline(h=0,col='gray')

plot(her$year,her$landings,pch=16,xlim=c(1965,2016),type='l',las=1,xlab='Year',ylab='SSB',col=cl3,xaxt='n',ylim=c(0,4000))
points(her$year,her$landings,pch=16,xlim=c(1965,2016),col=alpha(cl3,.3))
axis(1,seq(1965,2015,5))
abline(h=0)

dev.off()





################################################################

#### PLOTS LANDINGS OVER TIME FOR DIFFERENT FISHERIES
library(ggplot2)
library(tidyr)
library(RColorBrewer)
dt<-data.frame(cbind(her$year,her[,45:56]))

a<- dt %>% gather(stocklb, y, x4w.winter.purse.seine:offshore.scotian.shelf.banks)
names(a)[1]<-'year'
a$y<-a$y/1000

f<-function(d){
    return(data.frame(year=d$year,
                      y=d$y-mean(d$y,na.rm=TRUE)))
}
a<-ddply(a,.(stocklb),.fun=f)


b<-subset(a,stocklb %in% c('x4wx.stock.tac'))
b$num<-as.factor(1)
a<-subset(a,!(stocklb %in% c('x4wx.stock.tac','x4wx.stock.nominal.landings','non.stock.4xs.n.b..weir...shutoff','x4vwx.coastal.nova.scotia','offshore.scotian.shelf.banks','land.2016','x4wx.stock.adjusted.landings.')))

n<-length(unique(a$stocklb))

f<-function(d){
    return(data.frame(y=mean(d$y,na.rm=TRUE),
                      y2=sum(d$y,na.rm=TRUE)))
}
ttl<-ddply(a,.(year),.fun=f)
ttl$stocklb<-NA

dm<-data.frame(stocklb=sort(unique(a$stocklb)))
dm$num<-seq(1,dim(dm)[1],1)
a<-merge(a,dm,by=c('stocklb'),all=FALSE)

cls<-colorRampPalette(c('black','gray','magenta','blue3','green','yellow','orange','red2'))
dum2<-(t(data.frame(cls(n))))
aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
names(dum2)<-as.character(aa$stocklb)
rownames(dum2)<-NULL

#cls<-colorRampPalette(c('black','gray','magenta','blue3','green','yellow','orange','red2'))

cls<-colorRampPalette(c('black','darkgreen','green','lawngreen','palegreen'))
#cls<-colorRampPalette(c('black','darkmagenta','maroon1','orchid1','plum1'))
cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))

aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
dum3<-data.frame(stocklb=aa$stocklb,
                 num=aa$num,
                 cls=cls(n+2)[3:(n+2)])
dum3$stocklb2<-c('4W winter seine','4X summer gillnet','4Xqr summer seine','4Xr NS weir','4Xs fall-winter seine')

num=aa$num,
cls=cls(n+2)[3:(n+2)])
#dum3<-data.frame(stocklb=aa$stocklb,
#                 num=aa$num,
#                 cls=colorRampPalette(brewer.pal(9,'Paired'))(n+2)[3:(n+2)])

a<-a[order(a$num),]
a$stocklb<-factor(a$stocklb)
a$num<-factor(a$num)

p1<-ggplot(a, aes(x=year, fill=num,y=y))+
    geom_bar(data=subset(a,y>0),stat='identity',aes(width=.7,order=num),size=.0001)+
    geom_bar(data=subset(a,y<0),stat='identity',aes(width=.7,order=num),size=.0001)+
    theme(legend.position=c(.8,.8),text = element_text(size=8),axis.text.x = element_text(angle=0, vjust=1), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background = element_blank(), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),axis.line = element_line(color="black", size = .1))+
    scale_fill_manual(values=as.character(dum3$cls),breaks=as.factor(dum3$num),labels=as.character(dum3$stocklb2),name='Stock')+
    scale_x_continuous(expand=c(0,0),breaks=seq(1960,2015,5),labels=seq(1960,2015,5))+
    expand_limits(x=c(1960,2014))+
    xlab('Year')+
    ylab('Purse seine landings')+
    geom_hline(yintercept=0)+
  geom_line(data=b,aes(x=year,y=y),colour='gray30',size=1,alpha=.5)+
  geom_point(data=b,aes(x=year,y=y),colour='gray30',size=5,alpha=.5)






dt<-data.frame(cbind(her$year,her[,45:56]))

a<- dt %>% gather(stocklb, y, x4w.winter.purse.seine:offshore.scotian.shelf.banks)
names(a)[1]<-'year'
a$y<-a$y/1000


b<-subset(a,stocklb %in% c('x4wx.stock.tac'))
b$num<-as.factor(1)
a<-subset(a,!(stocklb %in% c('x4wx.stock.tac','x4wx.stock.nominal.landings','non.stock.4xs.n.b..weir...shutoff','x4vwx.coastal.nova.scotia','offshore.scotian.shelf.banks','land.2016','x4wx.stock.adjusted.landings.')))

n<-length(unique(a$stocklb))

f<-function(d){
    return(data.frame(y=mean(d$y,na.rm=TRUE),
                      y2=sum(d$y,na.rm=TRUE)))
}
ttl<-ddply(a,.(year),.fun=f)
ttl$stocklb<-NA

dm<-data.frame(stocklb=sort(unique(a$stocklb)))
dm$num<-seq(1,dim(dm)[1],1)
a<-merge(a,dm,by=c('stocklb'),all=FALSE)

cls<-colorRampPalette(c('black','gray','magenta','blue3','green','yellow','orange','red2'))
dum2<-(t(data.frame(cls(n))))
aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
names(dum2)<-as.character(aa$stocklb)
rownames(dum2)<-NULL

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))

aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
dum3<-data.frame(stocklb=aa$stocklb,
                 num=aa$num,
                 cls=cls(n+2)[3:(n+2)])
dum3$stocklb2<-c('4W winter seine','4X summer gillnet','4Xqr summer seine','4Xr NS weir','4Xs fall-winter seine')

num=aa$num,
cls=cls(n+2)[3:(n+2)])

a<-a[order(a$num),]
a$stocklb<-factor(a$stocklb)
a$num<-factor(a$num)

p2<-ggplot(a, aes(x=year, fill=num,y=y))+
    geom_bar(data=a,stat='identity',aes(width=.7,order=num),size=.0001)+
    theme(legend.position=c(.8,.8),text = element_text(size=8),axis.text.x = element_text(angle=0, vjust=1), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background = element_blank(), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),axis.line = element_line(color="black", size = .1))+
#    scale_fill_manual(values=dum2,name='Stock')+
    scale_fill_manual(values=as.character(dum3$cls),breaks=as.factor(dum3$num),labels=as.character(dum3$stocklb2),name='Stock')+
    scale_x_continuous(expand=c(0,0),breaks=seq(1960,2015,5),labels=seq(1960,2015,5))+
    expand_limits(x=c(1960,2014))+
    xlab('Year')+
    ylab('Purse seine landings')+
    geom_hline(yintercept=0)+
   geom_line(data=b,aes(x=year,y=y),colour='gray30',size=1,alpha=.5)+
   geom_point(data=b,aes(x=year,y=y),colour='gray30',size=5,alpha=.5)

setwd(figsdir)
pdf('landings_byfishery_barplot.pdf',height=12,width=10)
grid.arrange(p1,p2,ncol=1)
dev.off()
















###LANDIGNS BY AGE CLASS FROM 2016 ASSESS
library(tidyr)
setwd(datadir)
cdat<-read.csv('herring_assess2016_catchatage_pct.csv',header=TRUE,na.strings=c("- "))
dt<-cdat[,1:12]
a<- dt %>% gather(stocklb, y, X1:X11)
names(a)[1]<-'year'
n<-length(unique(a$stocklb))
a$num<-as.numeric(gsub('X','',a$stocklb))

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))

aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
dum3<-data.frame(stocklb=aa$stocklb,
                 num=aa$num,
                 cls=cls(n+2)[3:(n+2)])
dum3<-dum3[order(dum3$num),]

a<-a[order(a$num),]
a$stocklb<-factor(a$stocklb)
a$num<-factor(a$num)
a$num2<-as.factor((as.numeric(as.character(a$num))*-1)+12)

p1<-ggplot(a, aes(x=year, fill=num,y=y))+
geom_bar(data=a,stat='identity',aes(width=.7,order=num2),size=.0001)+
    theme(legend.position=c(.8,.8),text = element_text(size=8),axis.text.x = element_text(angle=0, vjust=1), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background = element_blank(), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),axis.line = element_line(color="black", size = .1))+
#    scale_fill_manual(values=dum2,name='Stock')+
    scale_fill_manual(values=as.character(dum3$cls),breaks=as.factor(dum3$num),labels=as.character(dum3$stocklb),name='Age class')+
    scale_x_continuous(expand=c(0,0),breaks=seq(1965,2015,5),labels=seq(1965,2015,5))+
    scale_y_continuous(expand=c(0,0),breaks=seq(0,100,10),labels=seq(0,100,10))+
    expand_limits(x=c(1964.5,2016),y=c(0,101))+
    xlab('Year')+
    ylab('Landings [%]')+
    geom_hline(yintercept=0)


setwd(datadir)
cdat2<-read.csv('herring_assess2016_catchatage.csv',header=TRUE,na.strings=c("- "))
dt<-cdat2[,1:12]
a<- dt %>% gather(stocklb, y, X1:X11)
names(a)[1]<-'year'

n<-length(unique(a$stocklb))
a$num<-as.numeric(gsub('X','',a$stocklb))

cls<-colorRampPalette(c('black','darkmagenta','hotpink','darkblue','lightskyblue','forestgreen','lawngreen','gold','orange','firebrick1','firebrick4'))

aa<-unique(subset(a,select=c('num','stocklb')))
aa<-aa[order(aa$num),]
dum3<-data.frame(stocklb=aa$stocklb,
                 num=aa$num,
                 cls=cls(n+2)[3:(n+2)])
dum3<-dum3[order(dum3$num),]

a<-a[order(a$num),]
a$stocklb<-factor(a$stocklb)
a$num<-factor(a$num)
a$num2<-as.factor((as.numeric(as.character(a$num))*-1)+12)

p2<-ggplot(a, aes(x=year, fill=num,y=y))+
geom_bar(data=a,stat='identity',aes(width=.7,order=num2),size=.0001)+
    theme(legend.position=c(.8,.8),text = element_text(size=8),axis.text.x = element_text(angle=0, vjust=1), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),panel.background = element_blank(), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6),axis.line = element_line(color="black", size = .1))+
#    scale_fill_manual(values=dum2,name='Stock')+
    scale_fill_manual(values=as.character(dum3$cls),breaks=as.factor(dum3$num),labels=as.character(dum3$stocklb),name='Age class')+
    scale_x_continuous(expand=c(0,0),breaks=seq(1965,2015,5),labels=seq(1965,2015,5))+
    scale_y_continuous(expand=c(0,0),breaks=seq(0,3500,500),labels=seq(0,3500,500))+
    expand_limits(x=c(1964.5,2016),y=c(0,3500))+
    xlab('Year')+
    ylab('Landings')+
    geom_hline(yintercept=0)

setwd(figsdir)
pdf('landings_byageclass_barplot.pdf',height=12,width=10)
grid.arrange(p1,p2,ncol=1)
dev.off()







f<-function(nm,lg,ttl,zt){
setwd(datadir)
print(nm)
d<-read.csv(paste(nm),header=TRUE)
ot<- d %>% gather(age, y, X1:X11)
ot$agen<-as.numeric(gsub('X','',ot$age))

if(zt==TRUE){
ff<-function(d){
    d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
    return(d)
}
ot<-ddply(ot,.(age),.fun=ff)
} else NULL

ot$ly<-log10(ot$y+.01)
dt<-unique(subset(ot,select=c('age','agen')))

if(lg==TRUE){
ott<-subset(ot,select=c('ly','age','year','agen'))
} else {
ott<-subset(ot,select=c('y','age','year','agen'))
}

names(ott)[1]<-'y'
return(ggplot()+
geom_tile(data=ott, aes(x=(agen), y=year,fill=y),col='gray80',size=.0001)+
scale_fill_distiller(palette='Spectral')+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position=c(.1,.2),plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.15, "in"),legend.text=element_text(size=6)) +
scale_x_continuous(expand=c(0,0),breaks=seq(1,11,1),labels=seq(1,11,1),limits=NA)    +
scale_y_continuous(expand=c(0,0),breaks=seq(1970,2015,5),labels=as.character(seq(1970,2015,5)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(1969.5,2015),xlim=c(0.5,11.5))+
    xlab('')+
    ylab('')+
       labs(title = ttl,
       subtitle = "",
       caption = '')
   )
}

p1<-f('herring_assess2006_ssb.csv',TRUE,'SSB at age',FALSE)
p2<-f('herring_assess2006_fatage.csv',FALSE,'F at age',FALSE)
p3<-f('herring_assess2016_weightatage.csv',FALSE,'Weight at age',FALSE)
p4<-f('herring_assess2016_weightatage.csv',FALSE,'Weight at age [Z]',TRUE)


setwd(figsdir)
pdf('herring_assess_trends_byage.pdf',height=8,width=10)
grid.arrange(p1,p2,ncol=2)
grid.arrange(p3,p4,ncol=2)

#EXTRA PLOTS OF WEIGHT AT AGE
setwd(datadir)
dd<-read.csv("herring_assess2016_weightatage.csv",header=TRUE)
ot<- dd %>% gather(age, y, X1:X11)
ot$agen<-as.numeric(gsub('X','',ot$age))

ff<-function(d){
    d$y<-(d$y-mean(d$y,na.rm=TRUE))/sd(d$y,na.rm=TRUE)
    return(d)
}
ot<-ddply(ot,.(age),.fun=ff)

xyplot(y~year | as.factor(agen),data=ot,pch=15,type=c('p','spline'),col='dodgerblue3')

par(mfrow=c(2,2),mar=c(4,4,1,1))
nms<-names(dd)[2:12]
cls<-matlab.like(length(nms))
plot(0,0,col='white',xlim=c(1965,2015),ylim=c(0,.5),las=1,ylab='Weight',xaxt='n')
axis(1,seq(1965,2015,5))
for(i in 1:length(nms)){
    d<-na.omit(d)
    d<-subset(dd,select=c(paste(nms[i]),'year'))
    names(d)[1]<-'y'
    xo<-seq(min(d$year),max(d$year),length.out=1000)
    asp<-aspline(d$year,d$y,xout=xo)
    lines(asp$x,asp$y,col=alpha(cls[i],.75),lwd=2)
}
legend('topright',gsub('X','',nms),col=cls,lwd=2,bty='n',ncol=2)


#ESTIMATES THE FORM OF THE TIME-TRENDS IN WEIGHT AT AGE
dm<-data.frame(year=sort(unique(ot$year)),
y=tapply(ot$y,ot$year,function(x) mean(x,na.rm=TRUE)))

modl<-lm(y~year,data=dm)
m<-breakpoints(y~year,data=dm,h=3,breaks=2)
modst<-lm(y~breakfactor(m,breaks=length(unique(m$breakpoints))),data=dm)
modnl<-lm(y~bs(year,degree=3,df=5),data=dm)
bp<-median(dm$year)
modd<-lm(y~year,data=dm)
modseg<-segmented(modd,seg.Z = ~year,psi=bp,control=seg.control(it.max=200))

dt<-data.frame(AIC(modl,modnl,modseg,modst))
names(dt)<-ifelse(names(dt) %in% c('BIC'),'AIC',names(dt))
dt$md<-rownames(dt)
dt$AIC<-ifelse(dt$md=='modnl',dt$AIC+2,dt$AIC)
dt1<-subset(dt,AIC==min(dt$AIC))
pdat<-data.frame(year=seq(min(dm$year),max(dm$year),length.out=100))
p<-predict(modseg,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
pdat$upr<-pdat$p+(1.96*pdat$se)
pdat$lwr<-pdat$p-(1.96*pdat$se)

s<-summary(modseg)
bp<-modseg$psi[1,2]
bpse<-modseg$psi[1,3]
plot(dm$year,dm$y,pch=16,las=1,xlab='Year',ylab='Average weight [Z-score]',xaxt='n',col='white')
abline(v=bp,col='red',lty=2)
polygon(c(pdat$year,pdat$year[length(pdat$year):1]),c(pdat$upr,pdat$lwr[length(pdat$lwr):1]),col=alpha('gray',.5),border=NA)
points(dm$year,dm$y,pch=16,col=alpha('black',.75))
lines(pdat$year,pdat$p)
axis(1,seq(1965,2015,5))
legend('topright',legend=paste('r2=',round(s$r.squared,digits=2)),bty='n')
dev.off()

library(akima)













###################################################
#HERRING LANDINGS FROM NAFO- PRODUCE WEIRD RESULTS
F<-function(){
setwd(datadir1)
herm<-read.csv('herring_landings_maritimes.csv',header=TRUE)
names(herm)<-tolower(names(herm))
herm2<-subset(herm,division=='4X')
herm<-data.frame(year=sort(unique(herm$year)),
landa=tapply(herm$catch...000.kg.,herm$year,function(x) sum(x,na.rm=TRUE)))

herm2<-data.frame(year=sort(unique(herm2$year)),
land4x=tapply(herm2$catch...000.kg.,herm2$year,function(x) sum(x,na.rm=TRUE)))
a<-merge(herm,herm2,by=c('year'),all=FALSE)
a$pct<-(a$land4x/a$landa)*100
plot(a$year,a$pct,pch=16,las=1,xlab='Year',ylab='Percent of all herring landings comprised of Div. 4X')


land<-read.csv('all_comm_landings_canada.csv',header=TRUE)
names(land)<-tolower(names(land))
land<-data.frame(year=sort(unique(land$year)),
landtot=tapply(land$catch...000.kg.,land$year,function(x) sum(x,na.rm=TRUE)))
a<-merge(a,land,by=c('year'),all=FALSE)
a$pct2<-(a$landa/a$landtot)*100
}




#LOAD LANDINGS FROM FAO
d<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/TS_FI_CAPTURE.csv',header=TRUE)
spc<-read.csv('N:/cluster_2017/scratch/spera/data/stagingdat/fao_fish_landings/CL_FI_SPECIES_GROUPS.csv',header=TRUE)
names(d)<-tolower(names(d))
names(spc)<-tolower(names(spc))
names(spc)[1]<-'species'

dher<-subset(d,species=='HER' & country==124)
dall<-subset(d,country==124 & fishing_area==21)
dall2<-merge(dall,spc,by=c('species'),all.x=TRUE,all.y=FALSE)
unique(dall$name_en)

dall<-data.frame(year=sort(unique(dall$year)),
quantityall=tapply(dall$quantity,dall$year,function(x) sum(x,na.rm=TRUE)))
dd<-merge(dall,dher,by=c('year'),all=FALSE)
dd$pct<-(dd$quantity/dd$quantityall)*100

setwd(figsdir)
pdf('landings_FAO.pdf',height=10,width=8)
par(mfrow=c(2,1),mar=c(4,4,1,4))
cl1<-'black'
cl2<-'firebrick3'
plot(dd$year,dd$pct,pch=16,las=1,xlab='Year',ylab='Herring landings as percentage of total',type='l',col=alpha(cl1,.5),lwd=2)
axis(1,seq(1950,2020,5),labels=FALSE)
axis(1,seq(1950,2020,10),labels=TRUE)
points(dd$year,dd$pct,pch=16,col=alpha(cl1,.5))

plot(dd$year,dd$pct,pch=16,las=1,xlab='Year',ylab='Herring landings as percentage of total',type='l',col=alpha(cl1,.5),lwd=2)
axis(1,seq(1950,2020,5),labels=FALSE)
axis(1,seq(1950,2020,10),labels=TRUE)
points(dd$year,dd$pct,pch=16,col=alpha(cl1,.5))
points(1955,subset(dd,year==1955)$pct,col=cl2,cex=1.5,pch=16)
points(1978,subset(dd,year==1978)$pct,col=cl2,cex=1.5,pch=16)
points(1992,subset(dd,year==1992)$pct,col=cl2,cex=1.5,pch=16)
text(1955,(subset(dd,year==1955)$pct)+1,'Major gear change',srt=50,adj=0,col=cl2,cex=.75)
text(1992,(subset(dd,year==1992)$pct)+1,'Groundfish moratorium',srt=50,adj=0,col=cl2,cex=.75)
text(1978,(subset(dd,year==1978)$pct)+1,'200 mile limit',srt=50,adj=0,col=cl2,cex=.75)


plot(dd$year,dd$quantity/1000,las=1,xlab='Year',ylab='Total herring landings [x1000]',pch=16,type='l',col=alpha(cl1,.5),lwd=2,ylim=c(90,550),xaxt='n')
axis(1,seq(1950,2020,5),labels=FALSE)
axis(1,seq(1950,2020,10),labels=TRUE)
points(dd$year,dd$quantity/1000,pch=16,col=alpha(cl1,.5))
points(1955,subset(dd,year==1955)$quantity/1000,col=cl2,cex=1.5,pch=16)
points(1978,subset(dd,year==1978)$quantity/1000,col=cl2,cex=1.5,pch=16)
points(1992,subset(dd,year==1992)$quantity/1000,col=cl2,cex=1.5,pch=16)
text(1955,(subset(dd,year==1955)$quantity/1000)+15,'Major gear change',srt=50,adj=0,col=cl2,cex=.75)
text(1992,(subset(dd,year==1993)$quantity/1000)+30,'Groundfish moratorium',srt=50,adj=0,col=cl2,cex=.75)
text(1978,(subset(dd,year==1979)$quantity/1000)+70,'200 mile limit',srt=50,adj=0,col=cl2,cex=.75)

cl1<-'darkmagenta'
cl2<-'forestgreen'
cl3<-'firebrick3'
plot(dd$year,dd$quantity/1000,las=1,xlab='Year',ylab='Total herring landings [x1000]',pch=16,type='l',col=alpha(cl2,.5),ylim=c(60,600),xaxt='n',lwd=3)
points(dd$year,dd$quantity/1000,pch=16,col=alpha(cl2,.5),cex=1.5)
par(new=TRUE)
plot(dd$year,dd$pct,pch=16,las=1,type='l',col=alpha(cl1,.5),lwd=3,xaxt='n',yaxt='n',ylim=c(5,37),ylab='',xlab='')
axis(1,seq(1950,2020,5),labels=FALSE)
axis(1,seq(1950,2020,10),labels=TRUE)
axis(4,seq(5,40,5),labels=TRUE,las=1)
points(dd$year,dd$pct,pch=16,col=alpha(cl1,.5),cex=1.5)
pt<-4
points(1955,pt,col=cl3,cex=1.5,pch=17)
points(1978,pt,col=cl3,cex=1.5,pch=17)
points(1992,pt,col=cl3,cex=1.5,pch=17)
pt<-5
text(1955,pt,'Major gear change',srt=0,adj=0,col=cl3,cex=.75)
text(1992,pt,'Groundfish moratorium',srt=0,adj=0,col=cl3,cex=.75)
text(1978,pt,'200 mile limit',srt=0,adj=0,col=cl3,cex=.75)
legend('topright',legend=c('Total landings','Proportion of all landings comprised of herring'),col=c(cl1,cl2),pch=15,bty='n')
dev.off()





herl<-read.csv('herring_landings_spawningarea.csv',header=TRUE)
herl<-herl[,!(names(herl) %in% c('X','X.1','X.2','X.3'))]
herl<- herl %>% gather(year, land, X1985:X2014)
herl$year<-as.numeric(as.character(gsub('X','',herl$year)))
names(herl)<-tolower(names(herl))
a<-subset(herl,stock.areas %in% c('Browns Bank','Gennet,Dry Ledge','German Bank','Lurcher','Long Island','Pollock Point','S.W. Grounds','Seal Island','Trinity'))
dm<-data.frame(year=sort(unique(a$year)),
swns=tapply(a$land,a$year,function(x) sum(x,na.rm=TRUE)))
dm2<-data.frame(year=sort(unique(herl$year)),
tot=tapply(herl$land,herl$year,function(x) sum(x,na.rm=TRUE)))
dma<-merge(dm,dm2,by=c('year'),all=FALSE)
dma$pct<-(dma$swns/dma$tot)*100
plot(dma$year,dma$pct,pch=15)
rowMeans(a$land)








#AVERAGE F PER YEAR
f<-function(d){
    return(data.frame(f=mean(d$fage2,d$fage3,d$fage4,d$fage5,d$fage6,d$fage7,d$fage8,d$fage9,d$fage10,d$fage11,na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)


plot(her$year,her$expr,pch=15)
par(new=TRUE)
plot(her$year,her$f,pch=15,col='red',yaxt='n')
axis(4,at=seq(plot(her$expr,her$f,pch=15)
text(her$expr,her$f,labels=her$year)

#AVERAGE WEIGHT PER YEAR
f<-function(d){
 return(data.frame(wgt=mean(d$wage7,d$wage8,d$wage9,d$wage10,d$wage11,na.rm=TRUE)))
}
d2<-ddply(her,.(year),.fun=f)
her<-merge(her,d2,by=c('year'),all.x=TRUE,all.y=FALSE)

her<-subset(her,select=c('year','r','ssb','s','obs.w','pred.w','f','wgt','expr'))
names(her)<-c('year','her.rec1','her.ssb','herlrv.surv','her.cf.rv','her.cf.cat','her.f','her.waa','her.expr')









setwd(figsdir)
pdatt<-pdat.fac
cr<-cor(pdat.fac[,-1],use='pairwise.complete.obs')
pdf('corrplot.pdf',height=14,width=14)
par(mar=c(4,1,1,1))
corrplot(cr,method='ellipse',type='lower',diag=FALSE,order='AOE',tl.col='black')
dev.off()





################################################################

# 'STOP-LIGHT' PLOT SHOWING MODEL PREDICTED RATES OF CHANGE IN VARIABLES WITH BREAKPOINTS OVERLAID

################################################################
#GET ORDER FOR IMAGE PLOT - SORT FROM RAPID DECLINE TO RAPID INCREASE
p<-pdat.lin
p$year<-as.numeric(rownames(pdat.lin))
vbls<-names(pdat.lin)
ll<-list()
for(i in 1:length(vbls)){
    print(vbls[i])
    d<-na.omit(subset(p,select=c('year',paste(vbls[i]))))
    names(d)[2]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd

    d1<-subset(d,year==min(d$year))
    d2<-subset(d,year==max(d$year))
    a<-data.frame(var=vbls[i],
                  chng=d1$z-d2$z)
    ll[[i]]<-a
}
ord<-data.frame(do.call('rbind',ll))
ord<-ord[order(ord$chng),]


#FORMAT AS MATRIX
vbls<-as.character(ord$var)
ll<-list()
for(i in 1:length(vbls)){
    d<-subset(pdat.lin,select=c(paste(vbls[i])))
    names(d)[1]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd
    dd<-data.frame(year=rownames(pdat.lin),
                   lbl=vbls[i],
                   vbl=i,
                   z=d$z)
    ll[[i]]<-dd
}
mat.lin<-rbind.fill(ll)
mat.lin$year<-as.numeric(as.character(mat.lin$year))


#FORMAT AS MATRIX
bpdat<-bpdat[order(bpdat$bpt,decreasing=TRUE),]
vbls<-as.character(bpdat$var)
ll<-list()
for(i in 1:length(vbls)){
    print(vbls[i])
    d<-subset(pdat.lin,select=c(paste(vbls[i])))
    names(d)[1]<-'y'
    d$mn<-mean(d$y,na.rm=TRUE)
    d$sdd<-sd(d$y,na.rm=TRUE)
    d$z<-(d$y-d$mn)/d$sdd
    dd<-data.frame(year=rownames(pdat.lin),
                       vbl=i,
                       z=d$z,
                   lbl=vbls[i])
    ll[[i]]<-dd
}
mat.lin.bp<-rbind.fill(ll)
mat.lin.bp$year<-as.numeric(as.character(mat.lin.bp$year))


#ADD PLOTTING COORDINATES TO BREAKPOINTS DATA
dm<-unique(subset(mat.lin.bp,select=c('vbl','lbl')))
bpdat<-merge(bpdat,dm,by.x='var',by.y='lbl',all.x=TRUE,all.y=FALSE)
bpdat$bpt<-bpdat$bpt+.5

ggplot <- function(...) { ggplot2::ggplot(...) + theme_bw() }


pltfun<-function(a){
aa<-data.frame(z=seq((-max(abs(a$z),na.rm=TRUE)),max(abs(a$z),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)
n<-21
brks<-seq((min(aa$z,na.rm=TRUE)-.01),max(aa$z,na.rm=TRUE)+.01,length.out=n)
brks2<-round(seq((min(aa$z,na.rm=TRUE)-.01),max(aa$z,na.rm=TRUE)+.01,length.out=n),digits=2)
a$zcat<-cut(a$z,breaks=brks)
a$zcat2<-cut(a$z,breaks=brks2)
lbls<-sort(unique(a$zcat2))
#cls<-colorRampPalette(brewer.pal(name='RdYlBu',8))
aa<-subset(a,is.na(z)==FALSE)
ylm<-seq(min(a$vbl,na.rm=TRUE),max(a$vbl,na.rm=TRUE),1)
lb<-unique(subset(a,select=c('vbl','lbl')))
lb<-na.omit(lb[order(lb$vbl),])

return(
    ggplot()+
    geom_tile(data=a, aes(x=year, y=vbl,fill=zcat),colour='white',size=.001)+
    geom_tile(data=aa, aes(x=year, y=vbl,fill=NA),colour='black',size=.1)+
    scale_fill_manual(breaks=as.character(lbls),values=matlab.like(length(lbls)),labels=lbls,na.value="transparent",guide=guide_legend(title=paste('Z-score'))) +
    scale_alpha(guide = 'none')+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(.2, "in"),legend.text=element_text(size=6))+
    scale_y_continuous(expand=c(0,0),breaks=ylm,labels=lb$lbl,limits=c(0,max(ylm)+1))
)
}

setwd(figsdir)
pdf('linear_trends_image.v2.pdf',height=7,width=9)
pltfun(subset(mat.lin.bp,year>=1965))+ geom_point(aes(x=bpt,y=vbl),data=bpdat,col='black',alpha=1,size=3,shape=23)
dev.off()

setwd(figsdir)
pdf('linear_trends_image2.v2.pdf',height=9,width=8)
pltfun(subset(mat.lin,year>=1965))

pltfun(subset(mat.lin,year>=1965))+ geom_point(aes(x=bpt,y=vbl),data=bpdat,col='black',alpha=1,size=3,shape=23)
dev.off()





################################################################

#  PLOTS RATES OF CHANGE SEPARATED BY SPECIES GROUPS

################################################################
#########################################################

#########################################################
#CONVERT TO TOTAL CHANGE IS Z-UNITS; NEED TO ADD SO THAT ALL POSITIVE
mpdat<-mpdat0
dum<-data.frame(t1=mpdat$t1+300,
                t2=mpdat$t2+300)
mpdat$chng<-dum$t2-dum$t1


brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)

#CATEGORIZE INTO MAJOR GROUPS: PHYTO, ZOOP, HERLARVAE, HERRING, ENVIRONMENTAL
mpdat$grp<-ifelse(regexpr('lrv', mpdat$var)>0 |
                  regexpr('zp', mpdat$var)>0,'d_zoop','a_env')
mpdat$grp<-ifelse(regexpr('her', mpdat$var)>0,'f_her',mpdat$grp)
mpdat$grp<-ifelse(regexpr('herlrv', mpdat$var)>0,'e_herlrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('ph\\.', mpdat$var)>0 |
                  regexpr('chl', mpdat$var)>0,'b_phyt',mpdat$grp)
mpdat$grp<-ifelse(regexpr('prey', mpdat$var)>0,'d_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('had', mpdat$var)>0,'g_had',mpdat$grp)


f<-function(d){d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

setwd(figsdir)
pdf('pedicted_model_changes.v2.pdf',height=10,width=9)
par(mar=c(4,8,1,1))
plot(mpdat$chng,mpdat$mvar,xlim=c(-5,15),col=as.character(mpdat$cl),pch=16,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Standard Deviations]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(10,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=16,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=1,cex=cx,lwd=.01,col='gray')

plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-5,15),col=as.character(mpdat$cl2),pch=16,yaxt='n',ylab='',cex=cx)
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(12,10,zr,col=rev(gray.colors(10)),horizontal=TRUE,strip.width=.04,strip.length=.3)
abline(h=lns+.5,col='gray',lwd=2)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=16,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=1,cex=cx,lwd=.1)



#1965 TO 2015
mpdat<-mpdat0
mpdat$chng<-mpdat$pct.65

brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)
mpdat$pchh<-ifelse(mpdat$pv<=.05,16,17)
mpdat$pchh2<-ifelse(mpdat$pv<=.05,21,2)

#CATEGORIZE INTO MAJOR GROUPS: PHYTO, ZOOP, HERLARVAE, HERRING, ENVIRONMENTAL
mpdat$grp<-ifelse(regexpr('lrv', mpdat$var)>0 |
                  regexpr('zp', mpdat$var)>0,'d_zoop','a_env')
mpdat$grp<-ifelse(regexpr('her', mpdat$var)>0,'f_her',mpdat$grp)
mpdat$grp<-ifelse(regexpr('herlrv', mpdat$var)>0,'e_herlrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('ph\\.', mpdat$var)>0 |
                  regexpr('chl', mpdat$var)>0,'b_phyt',mpdat$grp)
mpdat$grp<-ifelse(regexpr('prey', mpdat$var)>0,'d_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('had', mpdat$var)>0,'g_had',mpdat$grp)

f<-function(d){    d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-400,300),col=as.character(mpdat$cl),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(200,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(200,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-400,300),col=as.character(mpdat$cl2),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(200,5,zr,col=as.character(dm$cl2),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(150,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')




###1980-2005
mpdat<-mpdat0
mpdat$chng<-mpdat$pct.05

brks<-seq(-0.001,1,.1)
mpdat$r2f<-cut(mpdat$r2,breaks=brks)
dm<-data.frame(r2f=sort(unique(mpdat$r2f)),
               cl=matlab.like(length(unique(mpdat$r2f))),
               cl2=rev(gray.colors(length(unique(mpdat$r2f)))),
               alph=seq(0,1,length.out=10))
mpdat<-merge(mpdat,dm,by=c('r2f'),all.x=TRUE,all.y=FALSE)
mpdat$pchh<-ifelse(mpdat$pv<=.05,16,17)
mpdat$pchh2<-ifelse(mpdat$pv<=.05,21,2)

#CATEGORIZE INTO MAJOR GROUPS: PHYTO, ZOOP, HERLARVAE, HERRING, ENVIRONMENTAL
mpdat$grp<-ifelse(regexpr('lrv', mpdat$var)>0 |
                  regexpr('zp', mpdat$var)>0,'d_zoop','a_env')
mpdat$grp<-ifelse(regexpr('her', mpdat$var)>0,'f_her',mpdat$grp)
mpdat$grp<-ifelse(regexpr('herlrv', mpdat$var)>0,'e_herlrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('ph\\.', mpdat$var)>0 |
                  regexpr('chl', mpdat$var)>0,'b_phyt',mpdat$grp)
mpdat$grp<-ifelse(regexpr('prey', mpdat$var)>0,'d_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('had', mpdat$var)>0,'g_had',mpdat$grp)

f<-function(d){    d<-d[order(d$chng),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-100,100),col=as.character(mpdat$cl),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1985-2005 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
colorbar.plot(50,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
legend(50,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')
md<-subset(mpdat,abs(chng)>150)
md$x1<-ifelse(md$chng<150,-110,110)
f<-function(d){
    if(d$x1<0){
        arrows(d$x1+6,d$mvar,d$x1+2,d$mvar,length=.1)
    } else {
        arrows(d$x1-6,d$mvar,d$x1-2,d$mvar,length=.1)
    }
}
zz<-dlply(md,.(mvar),.fun=f)


par(mar=c(4,8,1,1))
plot(mpdat$chng,seq(1,dim(mpdat)[1],1),xlim=c(-100,100),col=as.character(mpdat$cl2),pch=mpdat$pchh,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1985-2005 Change [Percent]')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
colorbar.plot(50,5,zr,col=as.character(dm$cl2),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:4]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl2),pch=mpdat$pchh,cex=cx)
points(mpdat$chng,seq(1,dim(mpdat)[1],1),pch=mpdat$pchh2,cex=cx,lwd=.01,col='gray')
legend(50,10,legend=c('P-value<=0.05','P-value>0.05'),pch=c(16,17),bty='n')
md<-subset(mpdat,abs(chng)>150)
md$x1<-ifelse(md$chng<150,-110,110)
f<-function(d){
    if(d$x1<0){
        arrows(d$x1+6,d$mvar,d$x1+2,d$mvar,length=.1)
    } else {
        arrows(d$x1-6,d$mvar,d$x1-2,d$mvar,length=.1)
    }
}
zz<-dlply(md,.(mvar),.fun=f)

dev.off()






###################################################################

#PLOTS ABSOLUTE RATES OF CHANGE SEPARATED BY SPECIES GROUPS

###################################################################
mpdat<-mpdat0
dum<-data.frame(t1=mpdat$t1+300,
                t2=mpdat$t2+300)
mpdat$chng<-dum$t2-dum$t1
mpdat$chngabs<-abs(mpdat$chng)
mpdat$lchngabs<-log10(abs(mpdat$chng)+1)

#CATEGORIZE INTO MAJOR GROUPS: PHYTO, ZOOP, HERLARVAE, HERRING, ENVIRONMENTAL
mpdat$grp<-ifelse(regexpr('lrv', mpdat$var)>0 |
                    regexpr('zp', mpdat$var)>0,'d_zoop','a_env')
mpdat$grp<-ifelse(regexpr('her', mpdat$var)>0,'g_her',mpdat$grp)
mpdat$grp<-ifelse(regexpr('herlrv', mpdat$var)>0,'f_herlrv',mpdat$grp)
mpdat$grp<-ifelse(regexpr('ph\\.', mpdat$var)>0 |
                    regexpr('chl', mpdat$var)>0,'b_phyt',mpdat$grp)
mpdat$grp<-ifelse(regexpr('prey', mpdat$var)>0,'d_zoop',mpdat$grp)
mpdat$grp<-ifelse(regexpr('had', mpdat$var)>0,'e_had',mpdat$grp)

dm<-data.frame(grp=sort(unique(mpdat$grp)),
               cll=c('firebrick4','forestgreen','magenta4','gold3','dodgerblue4','pink'))
mpdat<-merge(mpdat,dm,by=c('grp'),all.x=TRUE,all.y=FALSE)

f<-function(d){d<-d[order(d$chngabs),]}
mpdat<-ddply(mpdat,.(grp),.fun=f)
mpdat$mvar<-seq(1,dim(mpdat)[1],1)
cx<-1.25

dm<-data.frame(grp=sort(unique(mpdat$grp)),
               cll=c('firebrick4','forestgreen','magenta4','gold3','dodgerblue4','pink'))
dm$lmn<-tapply(mpdat$lchngabs,mpdat$grp, mean)
dm$lmd<-tapply(mpdat$lchngabs,mpdat$grp,median)
dm$x<-(tapply(mpdat$mvar,mpdat$grp,min))
dm$lmn<-tapply(mpdat$lchngabs,mpdat$grp, mean)
dm$lmd<-tapply(mpdat$lchngabs,mpdat$grp,median)
dm$x<-(tapply(mpdat$mvar,mpdat$grp,min))

setwd(figsdir)
pdf('pedicted_absolute_model_changes.pdf',height=11,width=9)
par(mar=c(4,8,1,1))
plot(mpdat$lchngabs,mpdat$mvar,xlim=c(0,1.25),col=alpha(as.character(mpdat$cll),.75),pch=16,yaxt='n',ylab='',cex=cx,ylim=c(2.5,max(mpdat$mvar)-1.5),xlab='1965-2015 Change |[Standard Deviations]|')
axis(2,at=mpdat$mvar,labels=mpdat$var,las=1,cex=.75)
abline(v=0,lty=3)
zr<-dm$r2f
#colorbar.plot(200,5,zr,col=as.character(dm$cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
lns<-tapply(mpdat$mvar,mpdat$grp,max)[1:length(unique(mpdat$grp))]
abline(h=lns+.5,col='gray',lwd=1)
points(mpdat$lchngabs,seq(1,dim(mpdat)[1],1),col=as.character(mpdat$cl),pch=16,cex=cx)
points(mpdat$lchngabs,seq(1,dim(mpdat)[1],1),pch=1,cex=cx,lwd=.01,col='gray')
points(dm$lmd,rep(0,6),col=alpha(as.character(dm$cll),.7),pch=17,cex=2)

barplot(mpdat$lchngabs,horiz=TRUE,col=alpha(as.character(mpdat$cll),.75),border=NA,names.arg=mpdat$var,las=1,xlab='1965-2015 Change [log10|SDs|]',xlim=c(0,1.5),space=.4)
points(dm$lmd,rep(-3,6),col=alpha(as.character(dm$cll),.7),pch=17,cex=2)
legend('bottomright',legend=c('Environmental','Phytoplankton','Zooplankton','Herring larvae','Herring'),col=as.character(dm$cll),pch=15,bty='n')
dev.off()




setwd(datadir)
rvw<-read.csv("herring_weights_RV_survey_spera_allar.csv",header=TRUE)

setwd(figsdir)
pdf('ts_example.pdf',height=6,width=9)
par(mar=c(4,8,1,1))
a<-unique(subset(rvw,select=c('lon','lat')))
plot(a$lon,a$lat,pch=16,col=alpha('firebrick',.4),las=1,xlab='Longitude',ylab='Latitude',cex=.5)
map('world',add=TRUE,col='gray',fill=TRUE,border=NA)

a<-na.omit(subset(data,select=c('her.len.rv','her.len.rv.se','year')))
plot(a$year,a$her.len.rv,pch=15,xlab='Year',ylab='Herring length (cm)',las=1,ylim=c(20,35))
f<-function(d){lines(c(d$year,d$year),c(d$her.len.rv+(1.96*d$her.len.rv.se),d$her.len.rv-(1.96*d$her.len.rv.se)),col='black')}
z<-ddply(a,.(year),.fun=f)

map('world',col='gray20',fill=TRUE,border=NA,xlim=c(-67,-58),ylim=c(40,46))
dev.off()



################################################################

# REDUNDANCY ANALYSIS ON INTERPOLATED DATA

################################################################
################################################################3
#REDUNDANCY ANALYSIS OPERATES ON INTERPOLATED DATA - NEEDS NON-MISSING DATA
nms<-names(pdat.int)
resp<-nms[grepl('her\\.',nms)==TRUE]
resp<-resp[grepl('her\\.f',resp)==FALSE]
resp<-resp[grepl('\\.cat',resp)==FALSE]
resp<-resp[grepl('\\.expr',resp)==FALSE]
resp<-resp[grepl('\\.land',resp)==FALSE]
#resp<-resp[!(resp %in% c('her.agepe.cat','her.age5.cat','her.age.cat','her.f'))]
resp<-c('her.ssb','her.state')

#preds<-nms[!(nms %in% c(resp,'year','her.fmass.rv','her.metai.rv','her.sprng','her.spnug','her.spcv','her.spvar'))]
preds<-nms[grepl('her\\.',nms)==FALSE]
preds<-nms[grepl('\\.tplus',nms)==FALSE]
preds<-nms[grepl('\\.tmin',nms)==FALSE]
preds<-preds[!(preds %in% c('year','her.fmass.rv','her.metai.rv','her.sprng','her.spnug','her.spcv','her.spvar','her.waa','her.len.rv','her.szdiv.rv','her.szpe.rv','her.ssb'))]


pdat.int<-subset(pdat.int,year>=1982)
prds<-pdat.int[,names(pdat.int) %in% preds]
rsp<-pdat.int[,names(pdat.int) %in% resp]
rownames(rsp)<-pdat.int$year
rownames(prds)<-pdat.int$year

#names(prds)<-gsub(' ','_',names(prds))
srda<-rda(rsp~.,data=prds,scale=TRUE)
plot(srda,xlim=c(-1,2))



library(DataCombine)
a<-na.omit(subset(data,select=c('year','her.expr','her.state')))
round(cor(a,use='pairwise.complete.obs'),digits=2)
mod<-lm(her.state~her.expr,data=a)
mod2<-lm(her.state~her.expr+I(her.expr^2),data=a)
mod3<-gam(her.state~s(her.expr,k=10),data=a,gamma=.5)
par(mfrow=c(2,2))
plot(a$her.expr,a$her.state,pch=15)
lines(a$her.expr,predict(mod3))
a$r<-residuals(mod2,type='pearson')
a$r<-residuals(mod3,type='pearson')
plot(a$year,a$r,pch=15)

plot(data$her.expr,data$her.state,pch=15,col='white')
text(data$her.expr,data$her.state,labels=data$year)
cor(data$her.expr,data$her.state,use='pairwise.complete.obs')
#a<-subset(data,!(year %in% c(1965,1966,1967,1968,1970,
RsquareAdj(srda)$r.squared
#biplot(srda,scaling='symmetric',type=c('text','points'))
#ordiellipse(srda,)

setwd(figsdir)
pdf('rdaplot.pdf',height=8,width=9)
xlm <-c(-1.6,1.5)
ylm<-c(-1.25,1.5)
clenv<-'gray10'
clresp<-'springgreen4'
yr<-ifelse(as.numeric(rownames(rsp))<=1988,'red3','green')
cls<-c(rev(heat.colors(dim(rsp)[1])),'darkred')

plot(0,0,xlim=xlm,ylim=ylm,las=1,xlab='RDA1',ylab='RDA2',xaxt='n',yaxt='n',bg='gray')
rect(xlm[1]-1,ylm[1]-1,xlm[2]+1,ylm[2]+1,col='gray80')
points(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=2,col=cls,pch=16)
points(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=2,col='black',pch=1,lwd=.1)
#plot(srda,dislay='sp',add=TRUE)
spe.sc <- scores(srda, choices=1:2, scaling=2, display="sp",xlim=xlm,ylim=ylm)
env.sc <- scores(srda, choices=1:2, scaling=2, display="sites",xlim=xlm,ylim=ylm)
#sc <- scores(srda, choices=1:2, scaling=2, display="bp",xlim=xlm,ylim=ylm)
#arrows(0, 0, sc[, 1], sc[, 2], length=.2, lty=1, col="black",xlim=xlm,ylim=ylm,lwd=.5)
text(srda,display=c('bp'),scaling=2,xlim=xlm,ylim=ylm,cex=.75,col=alpha(clenv,.5))
text(srda,display=c('sp'),scaling=2,xlim=xlm,ylim=ylm,cex=1,col=clresp)
axis(side=1,at=round(seq(xlm[1],xlm[2],.4),digits=1),tick=T,labels=T)
axis(side=2,at=seq(ylm[1],ylm[2],.4),tick=T,labels=T,las=1)
box(lwd=1)
#arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length=.2, lty=1, col="gray40",xlim=xlm,ylim=ylm,lwd=.5)
#arrows(0, 0, env.sc[, 1], env.sc[, 2], length=.2, lty=1, col="red",xlim=xlm,ylim=ylm,lwd=.5)
#text(srda,display=c('sites'),scaling=2,xlim=xlm,ylim=ylm,cex=.6,col=yr)
zr<-sort(as.numeric(rownames(rsp)))
colorbar.plot(1,1.15,zr,col=cls,horizontal=TRUE,strip.width=.04,strip.length=.3)
dev.off()






##############################################################

#  MDA PLOTS

##############################################################
setwd(figsdir)
pdf('mda_plots.pdf',height=12,width=14)
par(mfrow=c(2,2),mar=c(4,4,1,1))
mdfun.cor<-function(yr1,yr2){
df<-subset(pdat.fac,year>=yr1 & year<=yr2)
df<-df[,!(names(df) %in% c('year','Phyt Rich','Phyt Rich Ct','Zoop Rich','Lrv Fish Rich Ct','Lrv Fish Rich'))]
n<-round(dim(df)[1]*0.5,digits=0)
df <- df[,colSums(is.na(df)) <=n]#EXCLUDE TS WHERE >20 MISSING
#colSums(is.na(df))
cldb<-data.frame(var=names(df))
cldb$cl<-ifelse(regexpr('Lrv', cldb$var)>0,'blue3','green3')
cldb$cl<-ifelse(regexpr('Her', cldb$var)>0,'black',cldb$cl)
cldb$cl<-ifelse(regexpr('Zoop', cldb$var)>0,'gold3',cldb$cl)
cldb$cl<-ifelse(regexpr('Dist', cldb$var)>0 |
                regexpr('Nitrate', cldb$var)>0 |
                regexpr('Wind', cldb$var)>0 |
                regexpr('Temp', cldb$var)>0 |
                regexpr('Strat', cldb$var)>0 |
                regexpr('SST', cldb$var)>0 |
                regexpr('Silicate', cldb$var)>0 |
                regexpr('Phosphate', cldb$var)>0 |
                regexpr('NAO', cldb$var)>0
               ,'red3',cldb$cl)
cldb$cl<-as.character(cldb$cl)

cldb$lab<-ifelse(regexpr('Lrv', cldb$var)>0,'Larval Fish','Phytoplankton')
cldb$lab<-ifelse(regexpr('Her', cldb$var)>0,'Herring',cldb$lab)
cldb$lab<-ifelse(regexpr('Zoop', cldb$var)>0,'Zooplankton',cldb$lab)
cldb$lab<-ifelse(regexpr('Dist', cldb$var)>0 |
                regexpr('Nitrate', cldb$var)>0 |
                regexpr('Wind', cldb$var)>0 |
                regexpr('Temp', cldb$var)>0 |
                regexpr('Strat', cldb$var)>0 |
                regexpr('SST', cldb$var)>0 |
                regexpr('Silicate', cldb$var)>0 |
                regexpr('Phosphate', cldb$var)>0 |
                regexpr('NAO', cldb$var)>0
               ,'Environment',cldb$lab)

#dmat<-1-(cor(df,use='pairwise.complete.obs'))
dmat<-1-(abs(cor(df,use='pairwise.complete.obs')))+.001
#dmat<-abs(cor(df,use='pairwise.complete.obs'))
dst<-as.dist(dmat)
#dst<-dist(df,method='euclidean')
ord<-metaMDS(dst,trymax=50,method='euclidean',trace=0)
plot(0,0,col='white',xlim=c(-.5,.6),ylim=c(-.5,.5),las=1)
text(ord, display = "sites", cex = 0.8, pch=21, col=cldb$cl, bg="yellow",adj=-.2)
points(ord, display = "sites", cex =1, col=alpha(cldb$cl,.5),pch=16)
points(ord, display = "sites", cex =1,,pch=1,lwd=.1)
cldb<-unique(subset(cldb,select=c('cl','lab')))
legend('topleft',cldb$lab,col=as.character(cldb$cl),pch=16,bty='n',pt.cex=2)
legend('bottomleft',paste(yr1,'-',yr2),bty='n')
}
par(mfrow=c(2,2))
mdfun.cor(1980,2006)
mdfun.cor(1965,2015)
mdfun.cor(1990,2015)


mdfun.dst<-function(yr1,yr2){
df<-subset(pdat.int,year>=yr1 & year<=yr2)
rownames(df)<-df$year
df<-df[,!(names(df) %in% c('year','Phyt Rich','Phyt Rich Ct','Zoop Rich','Lrv Fish Rich Ct','Lrv Fish Rich'))]

dst<-vegdist(df,method='euclidean')
ord<-metaMDS(dst,trymax=50,trace=0)
#ord<-metaMDS(decostand(df,method='standardize'),trace=0,distance='euclidean')
cl<-rev(heat.colors(dim(df)[1]))
ordiplot(ord,bg='gray',las=1)
rect(-8,-7,7,7,col='gray85')
clps<-ifelse(as.numeric(rownames(df))<=1990,'pre-collapse','post-collapse')
ordihull(ord,groups=clps,draw='lines',label=FALSE,col='black')
text(ord, display = "sites", cex = 1, pch=21, col=as.character(cl),adj=-.2)
points(ord, display = "sites", cex=1.5, col=alpha(cl,.5),pch=16)
points(ord, display = "sites", cex=1.5,,pch=1,lwd=.1)
zr<-as.numeric(rownames(df))
colorbar.plot(4,3.5,zr,col=as.character(cl),horizontal=TRUE,strip.width=.02,strip.length=.2)
legend('bottomleft',paste(yr1,'-',yr2),bty='n')
}
par(mfrow=c(2,2))
mdfun.dst(1982,2005)
mdfun.dst(1990,2015)
dev.off()





df<-subset(lrv,var=='Herring')
f<-function(d){
    names(d)[1]<-'time'
    print(dim(d))
    if(dim(d)[1]>100){
    mod<-gam(stdno~s(lon,lat,k=75),data=d)
    clonep<-seq(min(df$lon),max(df$lon),length.out=500)
    clatep<-seq(min(df$lat),max(df$lat),length.out=500)
    #PREDICTION DATA
    pdat2<-expand.grid(lon=clonep,lat=clatep)
    pdat2$p<-predict(mod,newdata=pdat2,type='response')
    dout<-subset(pdat2,p==max(pdat2$p))[1,]
    dout$n<-dim(d)[1]
    return(dout)
    } else NULL
}
com<-ddply(subset(df,select=c('tbin3','stdno','lon','lat')),.(tbin3),.fun=f,.progress='text')
comy<-ddply(subset(df,select=c('year','stdno','lon','lat')),.(year),.fun=f,.progress='text')

map('world',xlim=c(-70,-64),ylim=c(42,45))
points(comy$lon,comy$lat)
points(com$lon,com$lat)
plot(as.numeric(as.character(com$tbin3)),com$lat)
plot(as.numeric(as.character(com$tbin3)),com$lon)
plot(comy$year,comy$lat)
plot(as.numeric(as.character(com$tbin3)),com$lon)


f<-function(d){
    names(d)[1]<-'time'
    print(dim(d))
    if(dim(d)[1]>100){
    mod<-gam(stdno~s(lon,lat,k=75),data=d)
    clonep<-seq(min(df$lon),max(df$lon),length.out=500)
    clatep<-seq(min(df$lat),max(df$lat),length.out=500)
    #PREDICTION DATA
    pdat2<-expand.grid(lon=clonep,lat=clatep)
    pdat2$p<-predict(mod,newdata=pdat2,type='response')
    dout<-subset(pdat2,p==max(pdat2$p))[1,]
    dout$n<-dim(d)[1]
    return(dout)
    } else NULL
}
cos<-ddply(subset(df,select=c('tbin3','stdno','lon','lat')),.(tbin3),.fun=f,.progress='text')
cosy<-ddply(subset(df,select=c('year','stdno','lon','lat')),.(year),.fun=f,.progress='text')



library(reshape2)
library(zoo)
dat.m <- data.frame(yr=index(dat),value=melt(dat)$value)


dm<-na.omit(subset(her,select=c('year','s','ssb')))
s<-subset(dm,select=c('year','s'))
ssb<-subset(dm,select=c('year','ssb'))
ssb$year<-ssb$year+
s<-as.ts(s$s)
ssb<-as.ts(ssb$ssb)
par(mfrow=c(2,2))
plot(s$year,s$s,type='l')
plot(ss)
lines(lag(ss,k=-5),col='red')
plot(ssb)
lines(lag(ssb,k=-5),col='red')

s<-subset(dm,select=c('year','s'))
sts<-zoo(s$s,order.by=s$year,frequency=1)
plot(sts)
lines(lag(sts,k=5),col='red')

s<-subset(dm,select=c('year','s'))
sts<-zoo(s$s,order.by=s$year,frequency=1)
lt<-list()
for(i in 1:5){
slag<-lag(sts,k=i)
df<-data.frame(year=index(slag),
               surv=as.matrix(slag))
names(df)[2]<-gsub(' ','',paste(names(df)[2],i))
lt[[i]]<-df
}
lgdat<-Reduce(function(x, y) merge(x, y, by=c('year'),all=TRUE), lt)#COMBINE

ssb<-subset(dm,select=c('year','ssb'))
lgdat<-merge(lgdat,ssb,by=c('year'),all.x=TRUE,all.y=FALSE)
lgdat<-subset(lgdat,year>=1965)
round(cor(lgdat,use='pairwise.complete.obs'),digits=2)


ccf(dm$s,dm$ssb,ylim=c(-.5,.5))
abline(v=0,lty=2)
plot(dm$s,dm$ssb)

dst<-vegdist(decostand(pdat.int,method='standardize'),method='euclidean')
ord<-metaMDS(pdat.int,trymax=50,distance='euclidean')
#ord<-isoMDS(dst)
plot(ord)
plot(0,0,col='white',xlim=c(-.5,.5),ylim=c(-.5,.5),las=1)
text(ord, display = "sites", cex = 0.8, pch=21, col=cldb$cl, bg="yellow",adj=0)
#text(ord, display = "sites", cex = 0.8, pch=21, col="blue", bg="yellow",adj=0)
points(ord, display = "sites", cex = 0.8, col=cldb$cl, bg="yellow",pch=15)
text(ord, display = "spec", cex=0.7, col="blue")

stressplot(ord,dst)
decostand(varespec, "standardize")

ord<-isoMDS(dst)
autoplot(ord,label=TRUE)
autoplot(ord,label=TRUE,colour=cldb$cl)
slm<-sammon(dst)
autoplot(slm,label=TRUE)


dd<-subset(pdat.fac,select=c('HerSSB','SST','Wind','Stratification'))
mod<-lm(HerSSB~.^3,data=dd)
mod<-lm(HerSSB~.*.,data=dd)

library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)
library(ggfortify)
pd<-pdat.int
pd<-pd[,!(names(pd) %in% c('Rich1phyto','NAO','SSTduration.10','SSTampitude','Windampitude','CTPtotal','CTSST','PrichnessCT','PPtotalCT','Windduration.10'))]
#era<-ifelse(pd$year<=1987,'pre','collapse')
#era<-ifelse(pd$year>1995,'post',era)
era<-ifelse(pd$year<=1988,'pre','post')
era<-as.factor(era)
pd<-pd[,-1]

lbls<-sort(unique(era))
#c('collapse','pose','pre')
cls<-colorRampPalette(c('green3','magenta4','royalblue'))
cls<-cls(length(lbls))
xlm<-c(-5,8)
ylm<-c(-5,4)
pc<-prcomp(pd,center=TRUE,scale=TRUE)
s<-summary(pc)
plot(seq(1980,2006,1),pc$x[,1],type='l')
plot(seq(1980,2006,1),pc$x[,2],type='l')
plot(seq(1980,2006,1),pc$x[,3],type='l')

setwd(figsdir)
pdf('pcaplot.pdf',height=7,width=7)
ggbiplot(pc, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = FALSE,groups=era, var.axes=TRUE,size=5,shape=15)+
    scale_colour_manual(values=as.character(cls),breaks=lbls,guide=guide_legend(title=paste('Era')))+
    theme(legend.direction = 'horizontal', legend.position = 'top')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='top',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.2, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=xlm,labels=xlm,limits=xlm)+
scale_y_continuous(expand=c(0,0),breaks=ylm,labels=ylm,limits=ylm)+
    coord_equal()+
    coord_cartesian(ylim=ylm,xlim=xlm)
dev.off()


names(pd)<-gsub('>','',names(pd))
rownames(pd)<-seq(1980,2006,1)

pd2<-pd[,!(names(pd) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF'))]
pd2<-pd[,(names(pd) %in% c('SST','WInd','Stratification','Temperature50m'))]
pd2<-pd[,(names(pd) %in% c('ShannonP','TotPhyto','RIchphyto','Pielouphyto','PhytoTO','Richphyto'))]

pd2<-pdat.fac[,!((names(pdat.fac) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF')))]

nms<-c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF')
pd2<-(pdat.fac[,names(pdat.fac) %in% c('HerSSB','HerWeight','HerRecruitment','HerSurvivorship','HerF','LHerring')])[5:54,]



p1<-autoplot(kmeans(pd2,2),data=pd,label=TRUE,label.size=5,,alpha=.1,frame=TRUE)
p2<-autoplot(pam(pd2, 3), frame = TRUE, frame.type = 'norm',label=TRUE)
grid.arrange(p2,p3,ncol=2)
ag<-(agnes(pdat.fac,diss=FALSE,metric='manhattan',stand=TRUE))



#PLOTS ORDIPLOTS
ordfun<-function(pd2,n){
p1<-autoplot(kmeans(pd2,n),data=pd2,label=TRUE,label.size=3,,alpha=.1,frame=TRUE)
p2<-autoplot(pam(pd2, n), frame = TRUE, frame.type = 'norm',label=TRUE)
p3<-autoplot(clara(pd2,n),label=TRUE,alpha=.1,label.size=3,frame=TRUE)
pc2<-prcomp(pd2,center=TRUE,scale=TRUE)
p4<-autoplot(pc2,label=TRUE,alpha=.1,label.size=3,frame=TRUE)

dmat<-1-cor(pd2,use='pairwise.complete.obs')
dst<-as.dist(dmat)
grid.arrange(p1,ncol=1)
grid.arrange(p2,ncol=1)
grid.arrange(p3,ncol=1)
grid.arrange(p4,ncol=1)
}



setwd(figsdir)
pdf('ordiplots_herring.pdf',width=10,height=8)
par(mar=c(8,4,8,1))
psidat<-unique(subset(bpdat,var %in% c('Her F','Her SSB','Her CF obs','Her Weight','L Herring')))
her$ssb2<-her$ssb/1000
plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
cl1<-'gray80'
cl2<-'gray50'
cl3<-'gray10'
mx<-max(her$ssb2,na.rm=TRUE)+100
rect(1960,-100,1988.5,mx,col=alpha(cl1,.75),border=NA)
rect(1988.5,-100,1995,mx,col=alpha(cl2,.75),border=NA)
rect(1995,-100,2010,mx,col=alpha(cl3,.75),border=NA)
cl<-'red3'
rug(psidat$bpt,lwd=2,col=cl)
text(psidat$bpt[1],0,'Herring Condition Factor declines',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[2],0,'F increase',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[3]-.5,0,'SSB collapse',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[4],0,'Herring weights decline',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[5]+.5,0,'Herring larvae decline',srt=40,adj=0,col=cl,cex=.75)
text(1975,750,'Pre-collapse')
text(1991,750,'Collapse')
text(2001,750,'Post-collapse')
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
lines(her$year,her$ssb2,col=cl,lwd=2)
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

her$r2<-her$r/10000
plot(her$year,her$r2,type='l',las=1,xlab='Year',ylab='R [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,r2==max(her$r2,na.rm=TRUE))
mn<-subset(her,year==2004)
lines(c(mx$year-1,mx$year),c(mx$r2,mx$r2))
text(mx$year-1.75,mx$r2,round(mx$r2,digits=0))
text(mn$year+1,mn$r2+10,round(mn$r2,digits=0))
lines(her$year,her$ssb2,type='l',col='red')

par(mar=c(6,4,4,2))
pd2<-pdat.fac[,names(pdat.fac) %in% c('Her SSB','Her Weight','Her Recruitment','Her Survivorship','Her F')][5:54,]
names(pd2)<-gsub(' ','_',names(pd2))
ag<-(agnes(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
agh<-as.hclust(ag)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(agh, 4)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(agh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)

dn<-(diana(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
dnh<-as.hclust(dn)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(dnh, 2)
clusDendro = dendrapply(as.dendrogram(dnh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7,type='triangle')

par(mar=c(8,4,8,2))
plot(clusDendro, main = "",las=1,horiz=FALSE,cex=.7)

#par(mfrow=c(2,2),mar=c(4,4,1,1))
ordfun(na.omit(pd2),3)
dev.off()












par(mar=c(8,4,8,1))
psidat<-unique(subset(bpdat,var %in% c('Her F','Her SSB','Her CF obs','Her Weight','L Herring')))
her$ssb2<-her$ssb/1000
plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

plot(her$year,her$ssb2,type='l',las=1,xlab='Year',ylab='SSB [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
cl1<-'gray80'
cl2<-'gray50'
cl3<-'gray10'
mx<-max(her$ssb2,na.rm=TRUE)+100
rect(1960,-100,1988.5,mx,col=alpha(cl1,.75),border=NA)
rect(1988.5,-100,1995,mx,col=alpha(cl2,.75),border=NA)
rect(1995,-100,2010,mx,col=alpha(cl3,.75),border=NA)
cl<-'red3'
rug(psidat$bpt,lwd=2,col=cl)
text(psidat$bpt[1],0,'Herring Condition Factor declines',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[2],0,'F increase',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[3]-.5,0,'SSB collapse',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[4],0,'Herring weights decline',srt=40,adj=0,col=cl,cex=.75)
text(psidat$bpt[5]+.5,0,'Herring larvae decline',srt=40,adj=0,col=cl,cex=.75)
text(1975,750,'Pre-collapse')
text(1991,750,'Collapse')
text(2001,750,'Post-collapse')
mx<-subset(her,ssb2==max(her$ssb2,na.rm=TRUE))
mn<-subset(her,ssb2==min(her$ssb2,na.rm=TRUE))
lines(c(mx$year-1,mx$year),c(mx$ssb2,mx$ssb2))
lines(her$year,her$ssb2,col=cl,lwd=2)
text(mx$year-1.75,mx$ssb2,round(mx$ssb2,digits=0))
text(mn$year,mn$ssb2+40,round(mn$ssb2,digits=0))

her$r2<-her$r/10000
plot(her$year,her$r2,type='l',las=1,xlab='Year',ylab='R [Tons x 1000]',xaxt='n',xlim=c(1966,2005),ylim=c(0,750),lwd=2,yaxt='n')
axis(1,seq(1965,2010,5))
axis(2,seq(0,750,100),las=1)
mx<-subset(her,r2==max(her$r2,na.rm=TRUE))
mn<-subset(her,year==2004)
lines(c(mx$year-1,mx$year),c(mx$r2,mx$r2))
text(mx$year-1.75,mx$r2,round(mx$r2,digits=0))
text(mn$year+1,mn$r2+10,round(mn$r2,digits=0))
lines(her$year,her$ssb2,type='l',col='red')



library(ggfortify)
#PLOTS ORDIPLOTS
ordfun<-function(pdd,n){
pdd<-na.omit(pdd[,!(names(pdd) %in% c('her.cf.rv'))])
p1<-autoplot(kmeans(pdd,n),data=pdd,label=TRUE,label.size=3,,alpha=.1,frame=TRUE)
p2<-autoplot(pam(pdd, n), frame = TRUE, frame.type = 'norm',label=TRUE)
p3<-autoplot(clara(pdd,n),label=TRUE,alpha=.1,label.size=3,frame=TRUE)
pc2<-prcomp(pdd,center=TRUE,scale=TRUE)
p4<-autoplot(pc2,label=TRUE,alpha=.1,label.size=3,frame=TRUE)

dmat<-1-cor(pdd,use='pairwise.complete.obs')
dst<-as.dist(dmat)
grid.arrange(p1,ncol=1)
grid.arrange(p2,ncol=1)
grid.arrange(p3,ncol=1)
grid.arrange(p4,ncol=1)
}



setwd(figsdir)
pdf('ordiplots_herring_v2.pdf',width=10,height=8)
nms<-names(pdat.fac)
nms<-nms[grepl('her\\.',nms)==TRUE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('\\.prey',nms)==FALSE]
nms<-nms[grepl('\\.totno',nms)==FALSE]
nms<-nms[grepl('\\.totwgt',nms)==FALSE]
#nms<-nms[grepl('\\.cat',nms)==FALSE]
#pd2<-pdat.fac[,names(pdat.fac) %in% c(nms,'year')]
pd2<-pdat.fac[,names(pdat.fac) %in% c(nms)]

par(mar=c(6,4,4,2))
#pd2<-pdat.fac[,names(pdat.fac) %in% c('Her SSB','Her Weight','Her Recruitment','Her Survivorship','Her F')][5:54,]
#names(pd2)<-gsub(' ','_',names(pd2))
ag<-(agnes(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
agh<-as.hclust(ag)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(agh, 4)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(agh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)

dn<-(diana(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
dnh<-as.hclust(dn)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(dnh, 2)
clusDendro = dendrapply(as.dendrogram(dnh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)
#plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7,type='triangle')

par(mar=c(8,4,8,2))
plot(clusDendro, main = "",las=1,horiz=FALSE,cex=.7)

#par(mfrow=c(2,2),mar=c(4,4,1,1))
ordfun(na.omit(pd2),3)
dev.off()





################################################################

#PLOTS OUT CENTRAL TENDANCY TRENDS
nms<-names(pdat.fac)
nms<-nms[grepl('\\.ct',nms)==TRUE]
pd2<-pdat.fac[,names(pdat.fac) %in% c(nms,'wind.tmin','sst.t12','year')]
pd2$sst.t12<-ifelse(pd2$sst.t12==Inf,NA,pd2$sst.t12)
a<-pd2
#STANDARDIZE COLUMNS TO UNIT VARIANCE
f<-function(x){scale(x,center=TRUE,scale=FALSE)}
a<-data.frame(cbind(subset(a,select='year'),apply(a[,2:dim(a)[2]],2,f)))

setwd(figsdir)
pdf('central_tendancy_timetrends.pdf',height=7,width=9)
par(mar=c(4,4,1,1))
plot(0,0,xlim=c(1978,2015),ylim=c(-40,55),xlab='Year',ylab='Day of peak',las=1)
abline(h=0,lty=3)
f<-function(d,cl){
names(d)[2]<-'y'
lines(d$year,d$y,col=alpha(cl,.2),lwd=2)
points(d$year,d$y,pch=16,col=alpha(cl,.4))
}
f(subset(a,select=c('year','sst.t12')),'firebrick4')
f(subset(a,select=c('year','chl.ct')),'forestgreen')
f(subset(a,select=c('year','nut.ct')),'green')
f(subset(a,select=c('year','strt.ct')),'dodgerblue4')
f(subset(a,select=c('year','t50.ct')),'magenta4')
f(subset(a,select=c('year','wind.tmin')),'gold3')

aa<-data.frame(year=sort(unique(a$year)),
               mn=rowMeans(subset(a,select=c('sst.t12','chl.ct','nut.ct','strt.ct','t50.ct')),na.rm=TRUE))
lines(aa$year,aa$mn,col='black',lwd=2)
points(aa$year,aa$mn,col='black',pch=15)
legend('topleft',legend=c('SST','Chl','Nutrients','Stratification','Temperature 50m','Wind'),col=c('firebrick4','forestgreen','green','dodgerblue4','magenta4','gold3'),lwd=3,bty='n')
dev.off()






nms<-names(pdat.fac)
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
df<-pdat.fac[,names(pdat.fac) %in% nms]
df$sst.t12<-ifelse(df$sst.t12==Inf,NA,df$sst.t12)

k<-3
dmat<-abs(cor(df,use='pairwise.complete.obs',method='spearman'))
dst<-as.dist(dmat)
ff<-fanny(dst,k,maxit=10000,diss=T)

dum<-c('red3','forestgreen','gold','cornflowerblue','darkblue')
plot(silhouette(ff),col=dum[1:k],main='')#silhouette plot

dc.pcoa<-cmdscale(dst)
dc.scores<-scores(dc.pcoa,choices=c(1,2))
spefuz.g<-ff$clustering
a<-data.frame(var=as.character(sort(unique(names(df)))),
              clusters=ff$clustering)
aa<-data.frame(ff$membership)
aa$var<-rownames(a)

par(mar=c(1,1,1,8),oma=c(1,1,1,1))
plot(scores(dc.pcoa),asp=1,type='n',xlim=c(-.75,1.1),ylim=c(-1,1.2),las=1,axes=FALSE,xlab='',ylab='')
stars(ff$membership,location=scores(dc.pcoa),draw.segments=T,add=T,scale=F,len=.1,col.segments=alpha(c(dum[1:k]),.75),byt='n',labels=NULL,xlim=c(-1.1,1.4),ylim=c(-1,.5),lwd=.0001,xpd=TRUE,border=NULL,radius=FALSE,col.radius=alpha('white',.1))
for(i in 1:k){ cl<-dum[i]
    gg<-dc.scores[spefuz.g==i,]
    hpts<-chull(gg)
    hpts<-c(hpts,hpts[1])
    lines(gg[hpts,],col=cl,lwd=3,xlim=c(-1.1,1.4),ylim=c(-1,.5))
}








nms<-names(pdat.fac)
nms<-nms[grepl('her',nms)==FALSE]
nms<-nms[grepl('\\.ct',nms)==FALSE]
nms<-nms[grepl('had\\.',nms)==FALSE]
nms<-nms[grepl('\\.tmin',nms)==FALSE]
nms<-nms[grepl('\\.tmax',nms)==FALSE]
nms<-nms[grepl('\\.tplus',nms)==FALSE]
nms<-nms[grepl('lrv',nms)==FALSE]
nms<-nms[grepl('ph\\.',nms)==FALSE]
nms<-nms[grepl('zp\\.',nms)==FALSE]
nms<-nms[grepl('\\.dist',nms)==FALSE]
#pd2<-pdat.fac[,names(pdat.fac) %in% c(nms,'year')]



par(mar=c(6,4,4,2))
#pd2<-pdat.fac[,names(pdat.fac) %in% c('Her SSB','Her Weight','Her Recruitment','Her Survivorship','Her F')][5:54,]
#names(pd2)<-gsub(' ','_',names(pd2))
ag<-(agnes(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
agh<-as.hclust(ag)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(agh, 4)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(agh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)

dn<-(diana(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
dnh<-as.hclust(dn)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(dnh, 2)
clusDendro = dendrapply(as.dendrogram(dnh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)
#plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7,type='triangle')

par(mar=c(8,4,8,2))
plot(clusDendro, main = "",las=1,horiz=FALSE,cex=.7)

#par(mfrow=c(2,2),mar=c(4,4,1,1))
ordfun(na.omit(pd2),3)








nms<-names(pdat.fac)
nms<-nms[grepl('her\\.',nms)==TRUE]
nms<-nms[grepl('prey',nms)==FALSE]
pd2<-na.omit(subset(data,select=c(nms)))

par(mar=c(6,4,4,2))
#pd2<-pdat.fac[,names(pdat.fac) %in% c('Her SSB','Her Weight','Her Recruitment','Her Survivorship','Her F')][5:54,]
#names(pd2)<-gsub(' ','_',names(pd2))
ag<-(agnes(pd2,diss=FALSE,metric='euclidean',stand=TRUE))
agh<-as.hclust(ag)
labelColors = c('red3','royalblue','forestgreen','gold3')
clusMember = cutree(agh, 4)
# function to get color labels
colLab <- function(n) {
    if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
}
# using dendrapply
clusDendro = dendrapply(as.dendrogram(agh), colLab)
plot(clusDendro, main = "",las=1,horiz=TRUE,cex=.7)







####################################################################
#ESTIMATES VARIATION IN LARVAL ABUNDANCE BY DAY FOR EACH YEAR OF SURVEY
f<-function(d){
if(length(unique(d$day))>8){
d$stdno<-d$stdno+.1
print(unique(d$year))
#mod<-gam(stdno~s(day,k=4) + s(lon,lat,k=10),data=d,gamma=.5,family=Gamma('log'))
mod<-gam(stdno~s(day,k=4),data=d,gamma=.5,family=Gamma('log'))
pdat<-data.frame(lon=-67,
                 lat=43.5,
                 day=seq(min(d$day),max(d$day),.2),
                 year=unique(d$year))
 p<-predict(mod,newdata=pdat,se.fit=TRUE,type='response')
 pdat$p<-p$fit
 pdat$p<-(pdat$p-mean(pdat$p))/sd(pdat$p)
 pdat$se<-p$se.fit
 mx<-subset(pdat,p==max(pdat$p))
 pdat$tim<-mx$day
 return(pdat)
} else NULL
}
lphen<-ddply(larv,.(year),.fun=f,.progress='text')

lphen2<-subset(lphen,is.na(p)==FALSE)
dm<-data.frame(year=seq(1975,2010,1),
               day=315)
lphen2<-merge(lphen2,dm,by=c('year','day'),all=TRUE)

cls<-rich.colors(1000)
x2<-sort(unique(na.omit(lphen2$day)))
x1<-sort(unique(na.omit(lphen2$year)))
dm<-expand.grid(day=x2,year=x1)
datplot<-acast(lphen2,day~year,value.var="p",fun.aggregate=mean)

setwd(figsdir)
pdf('larval_herring_phenology.pdf',height=8,width=8)
par(mar=c(6,4,6,1))
image(x=sort(unique(x2)),y=sort(unique(x1)),datplot,col=cls,ylab='Year',xlab='Month',las=1,cex.axis=.8,axes=TRUE,useRaster=FALSE,xaxt='n',xlim=c(220,330))
image.plot(legend.only=TRUE,zlim=c(-2.7,3.7),col=cls,legend.lab="",horizontal=TRUE,smallplot=c(.2,.4,.7,.73),cex=.25)
dm<-data.frame(day=seq(15,345,30),
               month=c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'))
axis(3,at=seq(220,330,5),labels=TRUE)
axis(1,at=dm$day,labels=dm$month,tick=FALSE)
axis(1,at=dm$day-15,labels=FALSE,tick=TRUE)
dev.off()





setwd(figsdir)
pdf('herring_larvae_size_trend_spera.pdf',width=6,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
cl<-'firebrick4'
modf<-gam(length~as.factor(year) + s(lon,lat,k=20) + s(day,bs='cc',k=5),weights=larvs$clen,data=larvs,gamma=1.4)
pdat<-data.frame(year=sort(unique(larvs$year)),
                 lon=median(larvs$lon),
                 lat=median(larvs$lat),
                 day=230)
p<-predict(modf,newdata=pdat,se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit

mods<-gam(length~s(year,k=7) + s(lon,lat,k=20) + s(day,bs='cc',k=5),weights=larvs$clen,data=larvs,gamma=1.4)
pdats<-data.frame(year=seq(min(larvs$year),max(larvs$year),length.out=100),
                 lon=median(larvs$lon),
                 lat=median(larvs$lat),
                 day=230)
p<-predict(mods,newdata=pdats,se.fit=TRUE)
pdats$p<-p$fit
pdats$se<-p$se.fit
pdats$upr<-pdats$p+(1.96*pdats$se)
pdats$lwr<-pdats$p-(1.96*pdats$se)
plot(0,0,ylim=c(5.5,12),las=1,xlab='Year',ylab='Abundance-weighted larvae length',col='white',xlim=c(1975,2010))
polygon(c(pdats$year,pdats$year[length(pdats$year):1]),c(pdats$upr,pdats$lwr[length(pdats$lwr):1]),col=alpha(cl,.2),border=NA)
points(pdat$year,pdat$p,pch=16)
f<-function(d){ lines(c(d$year,d$year),c(d$p-(1.96*d$se),d$p+(1.96*d$se)),col=alpha('black',.2)) }
z<-dlply(pdat,.(year),.fun=f)
lines(pdats$year,pdats$p,col=cl)
dev.off()



##############################################
#PLOTS HERRING ABUNDANCE AND WEIGHT AS A FUNCTION OF AVERAGE TEMPERATURE (SURFACE+BOTTOM)/2 FROM RV SURVEY
setwd(figsdir)
pdf('herring_thermal_range_spera.pdf',width=6,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
cl<-'firebrick4'
plot(stdat$temp,stdat$swgt,pch=16,cex=rescale(stdat$n,newrange=c(.5,2.5)),ylim=c(0,500),las=1,xlab='Temperature',ylab='Total weight',xaxt='n',col=alpha('gray',.65))
#plot(stdat$temp,stdat$swgt,pch=15,ylim=c(0,500),las=1)
mod<-gam(swgt~s(temp,k=10),data=stdat,family=Gamma('log'),weights=stdat$n,gamma=1.4)
pdat<-data.frame(temp=seq(min(stdat$temp),max(stdat$temp),length.out=100))
p<-predict(mod,newdata=data.frame(temp=pdat$temp),type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
lines(pdat$temp,pdat$p,col=cl,lwd=2)
lines(pdat$temp,pdat$p-(1.96*pdat$se),lty=3,col=cl)
lines(pdat$temp,pdat$p+(1.96*pdat$se),lty=3,col=cl)
abline(h=0)
axis(1,at=seq(-2,17,1),labels=F,tick=TRUE)
axis(1,at=seq(0,15,5),labels=TRUE,tick=FALSE)

plot(stdat$temp,stdat$tno,pch=16,cex=rescale(stdat$n,newrange=c(.5,2.5)),ylim=c(0,10000),las=1,xlab='Temperature',ylab='Total number',xaxt='n',col=alpha('gray',.65))
#plot(stdat$temp,stdat$tno,pch=15,ylim=c(0,500),las=1)
mod<-gam(tno~s(temp,k=10),data=stdat,family=Gamma('log'),weights=stdat$n,gamma=1.4)
pdat<-data.frame(temp=seq(min(stdat$temp),max(stdat$temp),length.out=100))
p<-predict(mod,newdata=data.frame(temp=pdat$temp),type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
lines(pdat$temp,pdat$p,col=cl,lwd=2)
lines(pdat$temp,pdat$p-(1.96*pdat$se),lty=3,col=cl)
lines(pdat$temp,pdat$p+(1.96*pdat$se),lty=3,col=cl)
axis(1,at=seq(-2,17,1),labels=F,tick=TRUE)
axis(1,at=seq(0,15,5),labels=TRUE,tick=FALSE)
dev.off()




#USES MORE ROBUST MODEL APPROACH APPLIED TO RAW DATA INSTEAD
setwd(figsdir)
pdf('herring_larvae_thermal_range_spera.pdf',width=6,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,4))
#plot(rvw2$temp,rvw2$tno,pch=15,ylim=c(0,500),las=1)
mod<-gam(totno~s(temp,k=5),data=rvw2,family=nb(link='log'),gamma=1.4)
pdat<-data.frame(temp=seq(min(rvw2$temp,na.rm=TRUE),20,length.out=100))
p<-predict(mod,newdata=data.frame(temp=pdat$temp),type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
plot(rvw2$temp,rvw2$totno,pch=16,cex=1,ylim=c(0,10000),las=1,xlab='Average temperature',ylab='Total number',xaxt='n',col=alpha('gray20',.2),xlim=c(-1,20))
abline(h=0)
par(new=TRUE)
plot(pdat$temp,pdat$p,col=cl,lwd=2,type='l',yaxt='n',ylim=c(0,250),xlab='',ylab='',axes=FALSE)
lines(pdat$temp,pdat$p-(1.96*pdat$se),lty=3,col=cl)
lines(pdat$temp,pdat$p+(1.96*pdat$se),lty=3,col=cl)
a<-subset(pdat,temp>15.5)
lines(a$temp,a$p,col='green',lwd=2,type='l')
lines(a$temp,a$p-(1.96*a$se),lty=3,col='green')
lines(a$temp,a$p+(1.96*a$se),lty=3,col='green')
axis(1,at=seq(-2,20,1),labels=F,tick=TRUE)
axis(1,at=seq(0,20,5),labels=TRUE,tick=FALSE)
axis(4,at=seq(0,250,50),labels=TRUE,las=1)

rvw2$sampwgt2<-rvw2$sampwgt+.01
mod<-gam(sampwgt2~s(temp,k=6),data=subset(rvw2,is.na(sampwgt2)==FALSE),family=Gamma(link='log'),gamma=.7)
pdat<-data.frame(temp=seq(min(rvw2$temp,na.rm=TRUE),20,length.out=100))
p<-predict(mod,newdata=data.frame(temp=pdat$temp),type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
plot(rvw2$temp,rvw2$sampwgt2,pch=16,cex=1,ylim=c(0,100),las=1,xlab='Average temperature',ylab='Total weight',xaxt='n',col=alpha('gray20',.2),xlim=c(-1,20))
abline(h=0)
par(new=TRUE)
plot(pdat$temp,pdat$p,col=cl,lwd=2,type='l',yaxt='n',ylim=c(0,8),xlab='',ylab='',axes=FALSE)
lines(pdat$temp,pdat$p-(1.96*pdat$se),lty=3,col=cl)
lines(pdat$temp,pdat$p+(1.96*pdat$se),lty=3,col=cl)
a<-subset(pdat,temp>15.5)
lines(a$temp,a$p,col='green',lwd=2,type='l')
lines(a$temp,a$p-(1.96*a$se),lty=3,col='green')
lines(a$temp,a$p+(1.96*a$se),lty=3,col='green')
axis(1,at=seq(-2,20,1),labels=F,tick=TRUE)
axis(1,at=seq(0,20,5),labels=TRUE,tick=FALSE)
axis(4,at=seq(0,8,2),labels=TRUE,las=1)

mod<-gam(stdno~s(surftemp,k=5),data=ldat,family=nb(link='log'),gamma=1.4)
pdat<-data.frame(surftemp=seq(min(ldat$surftemp,na.rm=TRUE),20,length.out=100))
p<-predict(mod,newdata=data.frame(surftemp=pdat$surftemp),type='response',se.fit=TRUE)
pdat$p<-p$fit
pdat$se<-p$se.fit
plot(ldat$surftemp,ldat$stdno,pch=16,cex=1,ylim=c(0,2500),las=1,xlab='Surface temperature',ylab='Total number of larvae',xaxt='n',col=alpha('gray20',.2),xlim=c(-1,20))
abline(h=0)
par(new=TRUE)
plot(pdat$surftemp,pdat$p,col=cl,lwd=2,type='l',yaxt='n',ylim=c(0,80),xlab='',ylab='',axes=FALSE)
lines(pdat$surftemp,pdat$p-(1.96*pdat$se),lty=3,col=cl)
lines(pdat$surftemp,pdat$p+(1.96*pdat$se),lty=3,col=cl)
mx<-max(ldat$surftemp)
a<-subset(pdat,surftemp>mx)
lines(a$surftemp,a$p,col='green',lwd=2,type='l')
lines(a$surftemp,a$p-(1.96*a$se),lty=3,col='green')
lines(a$surftemp,a$p+(1.96*a$se),lty=3,col='green')
axis(1,at=seq(-2,20,1),labels=F,tick=TRUE)
axis(1,at=seq(0,20,5),labels=TRUE,tick=FALSE)
axis(4,at=seq(0,80,20),labels=TRUE,las=1)
dev.off()






#################################################
#PLOTS PREDICTED SST PHENOLOGY FOR PRE-POST 1988
setwd(figsdir)
pdf('sst_phenology_pre_post_spera.pdf',width=6,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
cl1<-'dodgerblue4'
cl2<-'firebrick4'
par(mar=c(6,4,1,1))
plot(0,0,xlim=c(5,360),ylim=c(3,16),las=1,col='white',xlab='Day of the year',ylab='Predicted SST',xaxt='n')
axis(1,at=seq(15,350,30),labels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),tick=FALSE,cex.axis=.75)
axis(1,at=seq(1,365,30),labels=FALSE,tick=TRUE)
rect(210,0,290,20,col=alpha('gold',.2),border=NA)
f<-function(d){
lines(d$day,d$p,col=as.character(unique(d$cl)),lwd=2)
mx<-subset(d,p==max(d$p))
points(mx$day,mx$p,pch=16,col=as.character(unique(d$cl)))
}
z<-dlply(phen3,.(tbinpp),.fun=f)
legend('topleft',legend=c('Pre-1988','Post-1988'),col=c(cl2,cl1),lwd=2,bty='n')
a<-subset(phen3,day==210)#DIFFERENCE IS 1.6 DEGREES
dev.off()


#############################################################
#FUNCTION TO RETURN MAP
pltfun<-function(adat,ttl,ct,dg){

nm<-names(adat)[1]
names(adat)[1]<-'y'
adat$y<-ifelse(adat$y>ct,ct,adat$y)
adat$y<-ifelse(adat$y< -ct,-ct,adat$y)
adat$y<-round(adat$y,digits=2)
a<-adat
aa<-data.frame(y=seq(-(max(abs(adat$y),na.rm=TRUE)),max(abs(adat$y),na.rm=TRUE),length.out=100))
a<-rbind.fill(a,aa)

n<-21
mxx<-max(abs(adat$y))
brks<-seq(-mxx-.01,mxx+.01,length.out=n)
brks2<-round(seq(-mxx-.01,mxx+.01,length.out=n),digits=dg)
a$ycat<-cut(a$y,breaks=brks)
lbls<-sort(unique(a$ycat))
lbls2<-sort(unique(cut(a$y,breaks=brks2)))
#cls<-matlab.like(length(lbls))
#cls<-colorRampPalette(c('magenta4','blue3','green','yellow','red3'))
#cls<-(cls(length(lbls)))
cls<-colorRampPalette(c('dodgerblue4','white','firebrick4'))
cls<-cls(length(lbls))


return(
ggplot()+
geom_tile(data=a, aes(x=lon, y=lat,fill=ycat),col='gray80',size=.0001)+
scale_fill_manual(breaks=as.character(lbls),values=cls,labels=lbls2,na.value="transparent",guide=guide_legend(title=paste(ttl)))+
scale_alpha(guide = 'none')+
geom_polygon(aes(long,lat, group=group), fill="grey65", data=coast.mc)+
theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right',plot.background=element_blank(),axis.line = element_line(color = 'black'), legend.key.size =  unit(0.1, "in"),legend.text=element_text(size=6))+
scale_x_continuous(expand=c(0,0),breaks=seq(-70,-62,1),labels=as.character(seq(-70,-62,1)),limits=NA)+
scale_y_continuous(expand=c(0,0),breaks=seq(43,46,1),labels=as.character(seq(43,46,1)),limits=NA)+
coord_equal()+
coord_cartesian(ylim=c(min(adat$lat)-.3,max(adat$lat)+.3),xlim=c(min(adat$lon)-.3,max(adat$lon)+.3))+
    xlab('')+
    ylab('')
)
}
p1<-pltfun(subset(phen2,select=c('delta.mxfall','lon','lat')),'mxfall',2.8,3)
p2<-pltfun(subset(phen2,select=c('delta.sdfall','lon','lat')),'sdfall',.7,3)
p3<-pltfun(subset(phen2,select=c('delta.tim','lon','lat')),'tim max',20,3)
p4<-pltfun(subset(phen2,select=c('delta.amp','lon','lat')),'amp',2.4,3)
p5<-pltfun(subset(phen2,select=c('delta.mn','lon','lat')),'min',.9,3)
p6<-pltfun(subset(phen2,select=c('delta.tim2','lon','lat')),'tim min',15,3)
p7<-pltfun(subset(phen2,select=c('delta.dur20','lon','lat')),'dur14',60,3)


setwd(figsdir)
pdf('sst_phen_maps_spera.pdf',width=5,height=10)
grid.arrange(p1,p2,p3,ncol=1)
grid.arrange(p4,p5,p6,ncol=1)
grid.arrange(p7,p7,p7,ncol=1)
dev.off()


#########################################################
#PLOT RELATIONSHIP BETWEEN LINEAR TIME TREND IN THE TIMING OF SEASONAL MINIMUM VERSUS MAXIMUM
setwd(figsdir)
pdf('sst_phenology_min_max_spera.pdf',width=6,height=8)
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(phen2$beta.tim2,phen2$beta.tim,pch=16,xlim=c(-1,1),ylim=c(-1,1),las=1,xlab='Time trend in timing of min SST',ylab='Time trend in timing of max SST',cex=.1)
f<-function(d){
lines(c(d$beta.tim2-(1.96*d$se.tim2),d$beta.tim2+(1.96*d$se.tim2)),c(d$beta.tim,d$beta.tim),col=alpha('dodgerblue4',.2))
lines(c(d$beta.tim2,d$beta.tim2),c(d$beta.tim-(1.96*d$se.tim),d$beta.tim+(1.96*d$se.tim)),col=alpha('dodgerblue4',.2))
}
z<-dlply(phen2,.(cell),.fun=f)
points(phen2$beta.tim2,phen2$beta.tim,pch=16,col='dodgerblue4',cex=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(a=0,b=1)

#PLOT RELATIONSHIP BETWEEN LINEAR TIME TREND IN SEASONAL MINIMUM VERSUS MAXIMUM
plot(phen2$beta.mn,phen2$beta.mx,pch=16,xlim=c(-.01,.11),ylim=c(-.01,.11),las=1,xlab='Time trend in seasonal min SST',ylab='Time trend in seasonal max SST',cex=.1)
f<-function(d){
lines(c(d$beta.mn-(1.96*d$se.mn),d$beta.mn+(1.96*d$se.mn)),c(d$beta.mx,d$beta.mx),col=alpha('dodgerblue4',.2))
lines(c(d$beta.mn,d$beta.mn),c(d$beta.mx-(1.96*d$se.mx),d$beta.mx+(1.96*d$se.mx)),col=alpha('dodgerblue4',.2))
}
z<-dlply(phen2,.(cell),.fun=f)
points(phen2$beta.mn,phen2$beta.mx,pch=16,col='dodgerblue4',cex=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(a=0,b=1)
dev.off()



