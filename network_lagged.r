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

dt<-subset(df,select=c('her.ssbc','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.prod','her.metai.rv','her.spvar','her.spcv','her.spnug','her.expr','her.land.spdiv','her.land'))
rspnms<-c('her.ssbc')
for(i in 1:length(rspnms)){
nm<-rspnms[i]
dt<-slide(dt,Var=nm,slideBy=1,NewVar=gsub(' ','',paste(nm,'.t1')))
dt<-slide(dt,Var=nm,slideBy=2,NewVar=gsub(' ','',paste(nm,'.t2')))
dt<-slide(dt,Var=nm,slideBy=3,NewVar=gsub(' ','',paste(nm,'.t3')))
dt<-slide(dt,Var=nm,slideBy=4,NewVar=gsub(' ','',paste(nm,'.t4')))
dt<-slide(dt,Var=nm,slideBy=5,NewVar=gsub(' ','',paste(nm,'.t5')))
dt<-slide(dt,Var=nm,slideBy=6,NewVar=gsub(' ','',paste(nm,'.t6')))
}


dts<-dt %>% gather(var,y,-her.ssbc.t6)
dts<-subset(dts,!(var %in% c('her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t2','her.ssbc.t3','her.ssbc.t4','her.ssbc.t5')))
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-gam(her.ssbc.t6~s(y,k=4),data=d)
s<-summary(mod)
return(data.frame(r2=round(s$r.sq,digits=2),
                resp='her.ssbc.t6'))
}
mdat6<-ddply(dts,.(var),.fun=f,.progress='text')


dts<-dt %>% gather(var,y,-her.ssbc.t5)
dts<-subset(dts,!(var %in% c('her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t2','her.ssbc.t3','her.ssbc.t4','her.ssbc.t6')))
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-gam(her.ssbc.t5~s(y,k=4),data=d)
s<-summary(mod)
return(data.frame(r2=round(s$r.sq,digits=2),
                resp='her.ssbc.t5'))
}
mdat5<-ddply(dts,.(var),.fun=f,.progress='text')


dts<-dt %>% gather(var,y,-her.ssbc.t4)
dts<-subset(dts,!(var %in% c('her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t2','her.ssbc.t3','her.ssbc.t5','her.ssbc.t6')))
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-gam(her.ssbc.t4~s(y,k=4),data=d)
s<-summary(mod)
return(data.frame(r2=round(s$r.sq,digits=2),
                resp='her.ssbc.t4'))
}
mdat4<-ddply(dts,.(var),.fun=f,.progress='text')


dts<-dt %>% gather(var,y,-her.ssbc.t3)
dts<-subset(dts,!(var %in% c('her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t2','her.ssbc.t4','her.ssbc.t5','her.ssbc.t6')))
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-gam(her.ssbc.t3~s(y,k=4),data=d)
s<-summary(mod)
return(data.frame(r2=round(s$r.sq,digits=2),
                resp='her.ssbc.t3'))
}
mdat3<-ddply(dts,.(var),.fun=f,.progress='text')

dts<-dt %>% gather(var,y,-her.ssbc.t2)
dts<-subset(dts,!(var %in% c('her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t3','her.ssbc.t4','her.ssbc.t5','her.ssbc.t6')))
#ESTIMATE TIME TRENDS
f<-function(d){
d<-na.omit(d)
nm<-unique(d$var)
print(nm)
mod<-gam(her.ssbc.t2~s(y,k=4),data=d)
s<-summary(mod)
return(data.frame(r2=round(s$r.sq,digits=2),
                resp='her.ssbc.t2'))
}
mdat2<-ddply(dts,.(var),.fun=f,.progress='text')

mdat<-rbind(mdat2,mdat3,mdat4,mdat5,mdat6)
mdat<-mdat[order(mdat$r2,decreasing=TRUE),]
plot(dt$her.expr,dt$her.ssbc.t2)#r2=.66
plot(dt$her.spcv,dt$her.ssbc.t6)#r2=.57
plot(dt$her.fmass.rv,dt$her.ssbc.t6)#r2=.49
plot(dt$her.waa,dt$her.ssbc.t6)#r2=.47



setwd(datadir)
load("impdatav.std.RData")
data<-impdatav.std
f<-function(d,xlb,ylb){
names(d)<-c('x','y')
plot(d$x,d$y,xlab=xlb,ylab=ylb,pch=15,col=alpha('black',.5),cex=.9)
r<-round(cor(d$x,d$y,use='pairwise.complete.obs'),digits=2)
legend('topright',legend=paste('r=',r),bty='n')
}
setwd(figsdir)
pdf('exploitation_scatter.pdf',height=7,width=8)
par(mfrow=c(2,2),mar=c(4,4,0,0),oma=c(4,1,1,1))
f(subset(data,select=c('her.land.sprich.t.1','her.ssb')),'Number of counties reporting catch','Herring SSB')
f(subset(data,select=c('her.expr.t.1','her.ssb')),'Exploitation rate [t-1]','Herring SSB')
f(subset(data,select=c('had.pi','her.ssb')),'predation','Herring SSB')
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
library(moments)
library(DataCombine)
library(rcompanion)
#data(marks)
#bn.gs = bnlearn::gs(marks)


#########################################

#dataset <- BNDataset(data = mat, discreteness=rep(FALSE,ncol(mat)), variables=nms, node.sizes=rep(2,dim(mat)[2]))



#################################################################

#CONSTRUCT BLACKLIST: NETWORK CONNECTIONS THAT ARE IMPOSSIBLE
#NOTHING CAN EFFECT EXPLOITATION VARIABLES
rspnms<-c("her.expr","her.land","her.land.pct1","her.land.spdiv","her.land.sprich",'bfish.state')
rnms<-c(rspnms,
          gsub(' ','',paste(rspnms,'.t-1')),
          gsub(' ','',paste(rspnms,'.t-2')))
prdnms<-names(mat[,!(names(mat) %in% c('year',rspnms))])
l<-list()
for(i in 1:length(prdnms)){
l2<-list()
for(j in 1:length(rspnms)){
l2[[j]]<-data.frame(from=prdnms[i],
                    to=rspnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist0<-data.frame(do.call('rbind',l))


#BIOLOGICAL VARIABLES CANNOT EFFECT CLIMATE
clmnms<-c("asl","gs.dist","nao", "nut.state","ss.ct","ss.dist","sst.fsd","sst.tmax","strt","strt.ct","t50.ct","temp.state","wind.amp",  "wind.state","wind.tmax", "wnd.ct","wnd.strs323fall",'nut.ct','nut.state','sst.t12')
clmnms<-c(clmnms,
          gsub(' ','',paste(clmnms,'.t-1')),
          gsub(' ','',paste(clmnms,'.t-2')))
bnms<-names(mat[,!(names(mat) %in% clmnms)])
l<-list()
for(i in 1:length(bnms)){
l2<-list()
for(j in 1:length(clmnms)){
l2[[j]]<-data.frame(from=bnms[i],
                    to=clmnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist1<-data.frame(do.call('rbind',l))


##############################################
#VARS AT T3 CAN'T EFFECT THOS IN THE PAST
nms<-names(mat)
nms<-nms[grepl('\\.t3',nms)==TRUE]
rnms<-names(mat[,!(names(mat) %in% nms)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist2<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+2 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t2',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist3<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t1',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE],
      nm[grepl('\\.t1',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist4<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T CAN'T EFFECT  THOSE IN PAST
nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
       nm[grepl('\\.t2',nm)==TRUE],
       nm[grepl('\\.t1',nm)==TRUE],
       nm[grepl('\\.t-1',nm)==TRUE],
       nm[grepl('\\.t-2',nm)==TRUE])
nms<-names(mat[,!(names(mat) %in% nm)])

nm<-names(mat)
nm<-c(nm[grepl('\\.t-2',nm)==TRUE],
      nm[grepl('\\.t-1',nm)==TRUE])
rnms<-names(mat[,names(mat) %in% nm])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist5<-data.frame(do.call('rbind',l))


#################################################
#VARS AT T-1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-c(nms[grepl('\\.t-1',nms)==TRUE])

nm<-names(mat)
rnms<-nm[grepl('\\.t-2',nm)==TRUE])
l<-list()
for(i in 1:length(nms)){
print(nms[i])
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist6<-data.frame(do.call('rbind',l))



#################################################
#OMITS LAGGED VERSIONS OF SELF-SAME VARIABLES
nms<-names(mat)
l<-list()
for(i in 1:length(nms)){
vc<-data.frame(y=c(gsub(' ','',paste(nm[i],'.t-2')),
                   gsub(' ','',paste(nm[i],'.t-1')),
                   gsub(' ','',paste(nm[i],'.t1')),
                   gsub(' ','',paste(nm[i],'.t2')),
                   gsub(' ','',paste(nm[i],'.t3')),
                   nm[i]))
l[[i]]<-expand.grid(from=vc$y,to=vc$y)
}
blist7<-data.frame(do.call('rbind',l))


bl<-unique(rbind(blist0,blist1,blist2,blist3,blist4,blist5,blist6,blist7))
bl<-subset(bl,from %in% names(mat) & to %in% names(mat))


########################################################

#FIT BAYESIAN NETWORK WHILE EXCLUSING BLACKLISTED EDGES
bn<-gs(mat,blacklist=bl,undirected=FALSE)
#bn<-hc(mat,blacklist=bl)
#bn<-iamb(mat,blacklist=bl)
#bn<-fast.iamb(mat,blacklist=bl)
#bn<-mmpc(mat,blacklist=bl)
#st<-arc.strength(bn,data=mat,criterion='x2',blacklist=bl,algorithm='gs')

#GET STRENGTH OF ASSOCIATIONS BY BOOTSTRAPPING
#strength.plot(bn,st)
st<-boot.strength(mat,R=5,algorithm='gs',algorithm.args=list(blacklist=bl),cpdag=FALSE)
#st<-boot.strength(mat,R=5,algorithm='iamb',algorithm.args=list(blacklist=bl),cpdag=FALSE)
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

nms<-c('ph.div.state','ph.mn.cpr','ph.mn.sabs','ph.pe.state','ph.to.sabs','ph.mn.state')
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


nm<-c('her.expr','her.land','her.land.pct1','her.land.spdiv','her.land.sprich','bfish.state')
nms<-c(nm,
       gsub(' ','',paste(nm,'.t-1')),
       gsub(' ','',paste(nm,'.t-2')))
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('lightpink',1)
} else NULL
}


#ly<-layout.fruchterman.reingold(g)
#ly<-layout_(g,as_star())
#ly<-layout_with_fr
#ly<-layout_with_kk(g)
#ly<-layout_drl
#ly<-layout_lgl(g)
#layout = g.layout('kamada_kawai')
#layout = g.layout("kk")
ly<-layout_nicely(g)

setwd(figsdir)
pdf('herring_pressures_network_iamb.pdf',height=12,width=12)
plot(g,
     layout=ly,
     edge.width=acd$strength*10,  # Edge width, defaults to 1
     edge.arrow.size=.4,           # Arrow size, defaults to 1
     edge.arrow.width=1.5,          # Arrow width, defaults to 1
     edge.curved=.7,
     vertex.label.color=alpha('black',.9),
     vertex.label.cex=.75,           #Font size
     vertex.label.dist=1.2,#Distance between label and the vertex
     vertex.shape="circle",
     vertex.size=7,       #Size of the node (default is 15)
     vertex.frame.color='gray60',     # Node border color
     edge.color=alpha('black',1),
     xlim=c(-1,1),ylim=c(-1,1),asp=0,lwd=.001)
dev.off()
























######################################################

#SAME AS ABOVE BUT OMIT MOST HERRING INDICES

######################################################
#source("https://bioconductor.org/biocLite.R")
#biocLite("RBGL")
#library(gRain)
library(bnlearn)
library(igraph)
library(DataCombine)
library(rcompanion)


#########################################

###  FIRST LOOK AT MISSINGNESS OF RAW DATA
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
#####################################################
#ADDS NEW CALIBRATED SSB
ddat<-subset(df,select=c('year','her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t3'))

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
f<-function(x){transformTukey(x,plotit=FALSE)}
#ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))

#Z-STANDARDIZE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))


###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv')
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
nms<-nms[grepl('her.ssb',nms)==FALSE]
nms<-nms[grepl('her.state',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c(omt))]

#ADDS CALIBRATED SSB
df<-merge(df,ddat,by=c('year'),all.x=TRUE,all.y=FALSE)

df<-df[,names(df) %in% as.character(dm2$var)]
#df$her.spnug<-df$her.spnug*-1
#df$her.spcv<-df$her.spcv*-1
#df$her.spvar<-df$her.spvar*-1
#PREDICTORS
prd<-c("asl", "bfish.state", "chl.ct",  "cp.cop.cpr.ct", "gs.dist", "had.dep.rv", "had.len.rv","had.pi","had.totno.rv", "had.totwgt.rv", "her.expr","her.jf.bof", "her.land","her.land.pct1","her.land.sprich", "her.preytot.bof", "herjuv.prey.bof", "lrv.mn.bof","lrv.pe.bof", "lrv.rich.bof","nao", "nut.ct", "nut.state", "ph.div.state","ph.mn.cpr","ph.mn.sabs","ph.mn.state", "ph.pe.state","ph.to.sabs","ss.ct", "ss.dist","sst.fsd", "sst.t12", "sst.tmax", "strt","strt.ct", "t50.ct", "temp.state","wind.amp","wind.fsd", "wind.state","wind.tmax", "wnd.ct", "wnd.strs323fall", "zp.cop.cpr", "zp.cpr.ct", "zp.div.cpr.state","zp.div.sabs.state","zp.mn.cpr", "zp.mn.cpr.ct","zp.mn.sabs", "zp.pe.cpr", "zp.pe.sabs","zp.to.sabs",'ph.mn.state')
prdnms<-c(prd,
          gsub(' ','',paste(prd,'.t-1')),
          gsub(' ','',paste(prd,'.t-2')))
preds<-df[,(names(df) %in% c(prdnms))]


#RESPONSE
rnm<-c('her.ssbc','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.prod','her.metai.rv','herlrv.mn.bof')
#rnm<-c('her.ssbc','her.prod','her.state')
rnms<-c(rnm,
          gsub(' ','',paste(rnm,'.t1')),
          gsub(' ','',paste(rnm,'.t2')),
          gsub(' ','',paste(rnm,'.t3')))
resp<-df[,(names(df) %in% c('year',rnms))]
df<-cbind(resp,preds)

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)


#DATA TO USE FOR NETWORK
mat<-subset(df,year>=1975 & year<=2005)
mat<-mat[,!(names(mat) %in% c('year'))]
nms<-names(mat)
#dataset <- BNDataset(data = mat, discreteness=rep(FALSE,ncol(mat)), variables=nms, node.sizes=rep(2,dim(mat)[2]))



#################################################################

#CONSTRUCT BLACKLIST: NETWORK CONNECTIONS THAT ARE IMPOSSIBLE
#NOTHING CAN EFFECT EXPLOITATION VARIABLES
rspnms<-c("her.expr","her.land","her.land.pct1","her.land.spdiv","her.land.sprich",'bfish.state')
rnms<-c(rspnms,
          gsub(' ','',paste(rspnms,'.t-1')),
          gsub(' ','',paste(rspnms,'.t-2')))
prdnms<-names(mat[,!(names(mat) %in% c('year',rspnms))])
l<-list()
for(i in 1:length(prdnms)){
l2<-list()
for(j in 1:length(rspnms)){
l2[[j]]<-data.frame(from=prdnms[i],
                    to=rspnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist0<-data.frame(do.call('rbind',l))


#BIOLOGICAL VARIABLES CANNOT EFFECT CLIMATE
clmnms<-c("asl","gs.dist","nao", "nut.state","ss.ct","ss.dist","sst.fsd","sst.tmax","strt","strt.ct","t50.ct","temp.state","wind.amp",  "wind.state","wind.tmax", "wnd.ct","wnd.strs323fall",'nut.ct','nut.state')
clmnms<-c(clmnms,
          gsub(' ','',paste(clmnms,'.t-1')),
          gsub(' ','',paste(clmnms,'.t-2')))
bnms<-names(mat[,!(names(mat) %in% clmnms)])
l<-list()
for(i in 1:length(bnms)){
l2<-list()
for(j in 1:length(clmnms)){
l2[[j]]<-data.frame(from=bnms[i],
                    to=clmnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist1<-data.frame(do.call('rbind',l))


##############################################
#VARS AT T3 CAN'T EFFECT THOS IN THE PAST
nms<-names(mat)
nms<-nms[grepl('\\.t3',nms)==TRUE]
rnms<-names(mat[,!(names(mat) %in% nms)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist2<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+2 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t2',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist3<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t1',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE],
      nm[grepl('\\.t1',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist4<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T CAN'T EFFECT  THOSE IN PAST
nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
       nm[grepl('\\.t2',nm)==TRUE],
       nm[grepl('\\.t1',nm)==TRUE],
       nm[grepl('\\.t-1',nm)==TRUE],
       nm[grepl('\\.t-2',nm)==TRUE])
nms<-names(mat[,!(names(mat) %in% nm)])

nm<-names(mat)
nm<-c(nm[grepl('\\.t-2',nm)==TRUE],
      nm[grepl('\\.t-1',nm)==TRUE])
rnms<-names(mat[,names(mat) %in% nm])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist5<-data.frame(do.call('rbind',l))


#################################################
#VARS AT T-1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-c(nms[grepl('\\.t-1',nms)==TRUE])

nm<-names(mat)
rnms<-nm[grepl('\\.t-2',nm)==TRUE])
l<-list()
for(i in 1:length(nms)){
print(nms[i])
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist6<-data.frame(do.call('rbind',l))



#################################################
#OMITS LAGGED VERSIONS OF SELF-SAME VARIABLES
nms<-names(mat)
l<-list()
for(i in 1:length(nms)){
vc<-data.frame(y=c(gsub(' ','',paste(nm[i],'.t-2')),
                   gsub(' ','',paste(nm[i],'.t-1')),
                   gsub(' ','',paste(nm[i],'.t1')),
                   gsub(' ','',paste(nm[i],'.t2')),
                   gsub(' ','',paste(nm[i],'.t3')),
                   nm[i]))
l[[i]]<-expand.grid(from=vc$y,to=vc$y)
}
blist7<-data.frame(do.call('rbind',l))


bl<-unique(rbind(blist0,blist1,blist2,blist3,blist4,blist5,blist6,blist7))
bl<-subset(bl,from %in% names(mat) & to %in% names(mat))


########################################################

#FIT BAYESIAN NETWORK WHILE EXCLUSING BLACKLISTED EDGES
bn<-gs(mat,blacklist=bl,undirected=FALSE)
#bn<-hc(mat,blacklist=bl)
#bn<-iamb(mat,blacklist=bl)
#bn3<-fast.iamb(mat,blacklist=bl)
#bn4<-mmpc(mat,blacklist=bl)
#st<-arc.strength(bn,data=mat,criterion='x2',blacklist=bl,algorithm='gs')

#GET STRENGTH OF ASSOCIATIONS BY BOOTSTRAPPING
#strength.plot(bn,st)
st<-boot.strength(mat,R=2,algorithm='gs',algorithm.args=list(blacklist=bl),cpdag=FALSE)
#st<-boot.strength(mat,R=1,algorithm='iamb',algorithm.args=list(blacklist=bl),cpdag=FALSE)
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
V(g)[nm]$color<-alpha('orangered',1)
} else NULL
}

nms<-c('ph.div.state','ph.mn.cpr','ph.mn.sabs','ph.pe.state','ph.to.sabs','ph.mn.state')
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


nm<-c('her.expr','her.land','her.land.pct1','her.land.spdiv','her.land.spring','bfish.state')
nms<-c(nm,
       gsub(' ','',paste(nm,'.t-1')),
       gsub(' ','',paste(nm,'.t-2')))
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('lightpink',1)
} else NULL
}


nm<-c('her.ssbc')
nms<-c(nm,
       gsub(' ','',paste(nm,'.t1')),
       gsub(' ','',paste(nm,'.t2')),
       gsub(' ','',paste(nm,'.t3')))
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('darkred',1)
} else NULL
}



ly<-layout_nicely(g)

setwd(figsdir)
pdf('herring_pressures_network_reduced_iamb.pdf',height=12,width=12)
plot(g,
     layout=ly,
     edge.width=acd$strength*10,  # Edge width, defaults to 1
     edge.arrow.size=.5,           # Arrow size, defaults to 1
     edge.arrow.width=1.5,          # Arrow width, defaults to 1
     edge.curved=.7,
     vertex.label.color=alpha('black',.9),
     vertex.label.cex=.75,           #Font size
     vertex.label.dist=1.2,#Distance between label and the vertex
     vertex.shape="circle",
     vertex.size=7,       #Size of the node (default is 15)
     vertex.frame.color='gray60',     # Node border color
     edge.color=alpha('black',1),
     xlim=c(-1,1),ylim=c(-1,1),asp=0,lwd=.001)
dev.off()






























######################################################

#SAME BUT EXCLUDE EXPLOITATION


#########################################

###  FIRST LOOK AT MISSINGNESS OF RAW DATA
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
#####################################################
#ADDS NEW CALIBRATED SSB
ddat<-subset(df,select=c('year','her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t3'))

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
f<-function(x){transformTukey(x,plotit=FALSE)}
#ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))

#Z-STANDARDIZE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))


###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv','bfish.state','her.expr','her.land','her.land.pct1','her.land.sprich')
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
nms<-nms[grepl('her.ssb',nms)==FALSE]
nms<-nms[grepl('her.state',nms)==FALSE]
nms<-nms[grepl('her.expr',nms)==FALSE]
nms<-nms[grepl('her.land',nms)==FALSE]
nms<-nms[grepl('bfish',nms)==FALSE]
df<-data[,names(data) %in% c('year',nms)]
df<-df[,!(names(df) %in% c(omt))]

#ADDS CALIBRATED SSB
df<-merge(df,ddat,by=c('year'),all.x=TRUE,all.y=FALSE)

df<-df[,names(df) %in% as.character(dm2$var)]
#PREDICTORS
prd<-c("asl", "chl.ct",  "cp.cop.cpr.ct", "gs.dist", "had.pi",'her.jf.bof', "her.preytot.bof", "herjuv.prey.bof", "lrv.mn.bof","lrv.pe.bof", "lrv.rich.bof","nao", "nut.ct", "nut.state", "ph.div.state","ph.mn.cpr","ph.mn.sabs","ph.mn.state", "ph.pe.state","ph.to.sabs","ss.ct", "ss.dist","sst.fsd", "sst.t12", "sst.tmax", "strt","strt.ct", "t50.ct", "temp.state","wind.amp","wind.fsd", "wind.state","wind.tmax", "wnd.ct", "wnd.strs323fall", "zp.cop.cpr", "zp.cpr.ct", "zp.div.cpr.state","zp.div.sabs.state","zp.mn.cpr", "zp.mn.cpr.ct","zp.mn.sabs", "zp.pe.cpr", "zp.pe.sabs","zp.to.sabs",'ph.mn.state')
prdnms<-c(prd,
          gsub(' ','',paste(prd,'.t-1')),
          gsub(' ','',paste(prd,'.t-2')))
preds<-df[,(names(df) %in% c(prdnms))]


#RESPONSE
rnm<-c('her.ssbc','her.len.rv','her.fmass.rv','herlrv.len','her.cf.rv','her.rec1','her.prod','herlrv.mn.bof')
#rnm<-c('her.ssbc','her.prod','her.state')
rnms<-c(rnm,
          gsub(' ','',paste(rnm,'.t1')),
          gsub(' ','',paste(rnm,'.t2')),
          gsub(' ','',paste(rnm,'.t3')))
resp<-df[,(names(df) %in% c('year',rnms))]
df<-cbind(resp,preds)

#CONVERT TO SHORT FORM, EXAMINE NORMALITY, TRANSFORM USING BEST EXPONENTIAL TRANSFORM
dfs<-df %>% gather(var,y,-year)


#DATA TO USE FOR NETWORK
mat<-subset(df,year>=1975 & year<=2005)
mat<-mat[,!(names(mat) %in% c('year'))]
nms<-names(mat)
#dataset <- BNDataset(data = mat, discreteness=rep(FALSE,ncol(mat)), variables=nms, node.sizes=rep(2,dim(mat)[2]))



#################################################################

#CONSTRUCT BLACKLIST: NETWORK CONNECTIONS THAT ARE IMPOSSIBLE
#BIOLOGICAL VARIABLES CANNOT EFFECT CLIMATE
clmnms<-c("asl","gs.dist","nao", "nut.state","ss.ct","ss.dist","sst.fsd","sst.tmax","strt","strt.ct","t50.ct","temp.state","wind.amp",  "wind.state","wind.tmax", "wnd.ct","wnd.strs323fall",'nut.ct','nut.state')
clmnms<-c(clmnms,
          gsub(' ','',paste(clmnms,'.t-1')),
          gsub(' ','',paste(clmnms,'.t-2')))
bnms<-names(mat[,!(names(mat) %in% clmnms)])
l<-list()
for(i in 1:length(bnms)){
l2<-list()
for(j in 1:length(clmnms)){
l2[[j]]<-data.frame(from=bnms[i],
                    to=clmnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist1<-data.frame(do.call('rbind',l))


##############################################
#VARS AT T3 CAN'T EFFECT THOS IN THE PAST
nms<-names(mat)
nms<-nms[grepl('\\.t3',nms)==TRUE]
rnms<-names(mat[,!(names(mat) %in% nms)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist2<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+2 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t2',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist3<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T+1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-nms[grepl('\\.t1',nms)==TRUE]

nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
      nm[grepl('\\.t2',nm)==TRUE],
      nm[grepl('\\.t1',nm)==TRUE])
rnms<-names(mat[,!(names(mat) %in% nm)])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist4<-data.frame(do.call('rbind',l))


#############################################
#VARS AT T CAN'T EFFECT  THOSE IN PAST
nm<-names(mat)
nm<-c(nm[grepl('\\.t3',nm)==TRUE],
       nm[grepl('\\.t2',nm)==TRUE],
       nm[grepl('\\.t1',nm)==TRUE],
       nm[grepl('\\.t-1',nm)==TRUE],
       nm[grepl('\\.t-2',nm)==TRUE])
nms<-names(mat[,!(names(mat) %in% nm)])

nm<-names(mat)
nm<-c(nm[grepl('\\.t-2',nm)==TRUE],
      nm[grepl('\\.t-1',nm)==TRUE])
rnms<-names(mat[,names(mat) %in% nm])
l<-list()
for(i in 1:length(nms)){
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist5<-data.frame(do.call('rbind',l))


#################################################
#VARS AT T-1 CAN'T EFFECT  THOSE IN PAST
nms<-names(mat)
nms<-c(nms[grepl('\\.t-1',nms)==TRUE])

nm<-names(mat)
rnms<-nm[grepl('\\.t-2',nm)==TRUE])
l<-list()
for(i in 1:length(nms)){
print(nms[i])
l2<-list()
for(j in 1:length(rnms)){
l2[[j]]<-data.frame(from=nms[i],
                    to=rnms[j])
}
l[[i]]<-data.frame(do.call('rbind',l2))
}
blist6<-data.frame(do.call('rbind',l))



#################################################
#OMITS LAGGED VERSIONS OF SELF-SAME VARIABLES
nms<-names(mat)
l<-list()
for(i in 1:length(nms)){
vc<-data.frame(y=c(gsub(' ','',paste(nm[i],'.t-2')),
                   gsub(' ','',paste(nm[i],'.t-1')),
                   gsub(' ','',paste(nm[i],'.t1')),
                   gsub(' ','',paste(nm[i],'.t2')),
                   gsub(' ','',paste(nm[i],'.t3')),
                   nm[i]))
l[[i]]<-expand.grid(from=vc$y,to=vc$y)
}
blist7<-data.frame(do.call('rbind',l))


bl<-unique(rbind(blist1,blist2,blist3,blist4,blist5,blist6,blist7))
bl<-subset(bl,from %in% names(mat) & to %in% names(mat))


########################################################

#FIT BAYESIAN NETWORK WHILE EXCLUSING BLACKLISTED EDGES
bn<-gs(mat,blacklist=bl,undirected=FALSE)
#bn<-hc(mat,blacklist=bl)
#bn<-iamb(mat,blacklist=bl)
#bn3<-fast.iamb(mat,blacklist=bl)
#bn4<-mmpc(mat,blacklist=bl)
#st<-arc.strength(bn,data=mat,criterion='x2',blacklist=bl,algorithm='gs')

#GET STRENGTH OF ASSOCIATIONS BY BOOTSTRAPPING
#strength.plot(bn,st)
st<-boot.strength(mat,R=2,algorithm='gs',algorithm.args=list(blacklist=bl),cpdag=FALSE)
#st<-boot.strength(mat,R=1,algorithm='iamb',algorithm.args=list(blacklist=bl),cpdag=FALSE)
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
V(g)[nm]$color<-alpha('orangered',1)
} else NULL
}

nms<-c('ph.div.state','ph.mn.cpr','ph.mn.sabs','ph.pe.state','ph.to.sabs','ph.mn.state')
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

nm<-c('her.ssbc')
nms<-c(nm,
       gsub(' ','',paste(nm,'.t1')),
       gsub(' ','',paste(nm,'.t2')),
       gsub(' ','',paste(nm,'.t3')))
for(i in 1:length(nms)){
nm<-nms[i]
if(nm %in% names(V(g))){
V(g)[nm]$color<-alpha('darkred',1)
} else NULL
}



ly<-layout_nicely(g)

setwd(figsdir)
pdf('herring_pressures_network_reduced_noexpl_gs.pdf',height=12,width=12)
plot(g,
     layout=ly,
     edge.width=acd$strength*10,  # Edge width, defaults to 1
     edge.arrow.size=.5,           # Arrow size, defaults to 1
     edge.arrow.width=1.5,          # Arrow width, defaults to 1
     edge.curved=.7,
     vertex.label.color=alpha('black',.9),
     vertex.label.cex=.75,           #Font size
     vertex.label.dist=1.2,#Distance between label and the vertex
     vertex.shape="circle",
     vertex.size=7,       #Size of the node (default is 15)
     vertex.frame.color='gray60',     # Node border color
     edge.color=alpha('black',1),
     xlim=c(-1,1),ylim=c(-1,1),asp=0,lwd=.001)
dev.off()



install.packages('mgm')
library(mgm)
library(qgraph)
p<-ncol(mat)
fito<-mgm(data=mat,type=rep('g',p),lev=rep(1,p),rule.reg='OR',k=2)
prd<-predict(fito,mat,error.continuous='VarExpl')
qgraph(fito$pairwise$wadj,layout='spring',pie=prd$obj$error$Error,pieColor=rep('#377EB8',p),node.color=fito$edgecolor,labels=colnames(mat))


qgraph(fito$pairwise$wadj,layout='spring')













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
#####################################################
#ADDS NEW CALIBRATED SSB
ddat<-subset(df,select=c('year','her.ssbc','her.ssbc.t1','her.ssbc.t2','her.ssbc.t3'))

#EXPONENTIAL TRANSFORM TO 'NORMALIZE'
#f<-function(x){transformTukey(x,plotit=FALSE)}
#ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))

#Z-STANDARDIZE
f<-function(x){scale(x,center=TRUE,scale=TRUE)}
ddat<-data.frame(cbind(subset(ddat,select='year'),apply(ddat[,2:dim(ddat)[2]],2,f)))


###############################################################
#1 LOAD IMPUTED DATA AND SUBSET BASED ON MISSINGNESS
omt<-c('had.dep.rv','wind.fsd','had.totwgt.rv','had.totno.rv','had.len.rv')
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
df<-df[,!(names(df) %in% c(omt))]

#ADDS CALIBRATED SSB
#df<-merge(df,ddat,by=c('year'),all.x=TRUE,all.y=FALSE)

df<-df[,names(df) %in% as.character(dm2$var)]
#PREDICTORS
prd<-c("asl", "bfish.state", "chl.ct",  "cp.cop.cpr.ct", "gs.dist", "had.dep.rv", "had.len.rv","had.pi","had.totno.rv", "had.totwgt.rv", "her.expr","her.jf.bof", "her.land","her.land.pct1","her.land.sprich", "her.preytot.bof", "herjuv.prey.bof", "lrv.mn.bof","lrv.pe.bof", "lrv.rich.bof","nao", "nut.ct", "nut.state", "ph.div.state","ph.mn.cpr","ph.mn.sabs","ph.mn.state", "ph.pe.state","ph.to.sabs","ss.ct", "ss.dist","sst.fsd", "sst.t12", "sst.tmax", "strt","strt.ct", "t50.ct", "temp.state","wind.amp","wind.fsd", "wind.state","wind.tmax", "wnd.ct", "wnd.strs323fall", "zp.cop.cpr", "zp.cpr.ct", "zp.div.cpr.state","zp.div.sabs.state","zp.mn.cpr", "zp.mn.cpr.ct","zp.mn.sabs", "zp.pe.cpr", "zp.pe.sabs","zp.to.sabs",'ph.mn.state')
prdnms<-c(prd,
          gsub(' ','',paste(prd,'.t-1')),
          gsub(' ','',paste(prd,'.t-2')))
preds<-df[,(names(df) %in% c(prdnms))]

#RESPONSE
rnm<-c('her.ssbc','her.szpe.rv','herjuv.fmass.rv','her.len.rv','her.fmass.rv','her.waa','herlrv.len','herjuv.metai.rv','her.cf.rv','her.rec1','her.ajrat.rv','her.prod','her.metai.rv','herlrv.mn.bof','herlrv.surv','her.state')
rnms<-c(rnm,
          gsub(' ','',paste(rnm,'.t1')),
          gsub(' ','',paste(rnm,'.t2')),
          gsub(' ','',paste(rnm,'.t3')))
resp<-df[,(names(df) %in% c('year',rnms))]
df<-cbind(resp,preds)




dff<-df
names(dff)<-gsub('-','.',names(dff))
dff2<-dff %>% gather(var,y,-her.ssbc.t3)
dff2<-subset(dff2,!(var %in% c('her.ssbc.t2','her.ssbc.t1','her.ssbc','her.state.t2','her.state.t1','her.state','year')))

nms<-c('her.ssbc.t3','her.rec1.t1','herlrv.mn.bof','herlrv.surv','herlrv.len')
l<-list()
for(i in 1:length(nms)){
print(nms[i])
d0<-subset(dff,select=c('year',nms[i]))
d1<-dff %>% gather(var,x,-year)
d1<-merge(d1,d0,by=c('year'),all.x=TRUE,all.y=FALSE)
names(d1)[4]<-c('y')

f<-function(d){
print(unique(d$var))
d<-na.omit(d)
modlm<-lm(y~x,data=d)
slm<-summary(modlm)
d$resid<-residuals(modlm)
t<-ts(d$resid)
acpar<-mean(acf(t,plot=F)$acf[2])
    return(data.frame(b=round(slm$coef[2,1],digits=2),
                      pv=round(slm$coef[2,4],digits=2),
                      r2=round(slm$r.squared,digits=2)))
}
ot<-ddply(d1,.(var),.fun=f,.progress='text')
ot$rvar<-nms[i]
l[[i]]<-data.frame(ot)
}

mdat<-data.frame(do.call('rbind',l))
write.csv(mdat,'C:/Users/sailfish/Downloads/prac.csv')
mdat<-read.csv('C:/Users/sailfish/Downloads/prac.csv',header=TRUE)


a<-subset(mdat,rvar=='her.ssbc.t3')
a<-a[order(a$r2,decreasing=TRUE),]

a<-subset(mdat,rvar=='her.rec1.t1')
a<-a[order(a$r2,decreasing=TRUE),]

a<-subset(mdat,rvar=='herlrv.surv')
a<-a[order(a$r2,decreasing=TRUE),]

s<-subset(ot,pv<0.05)
dat1<-subset(dff,select=c(s$var,'her.ssbc.t3'))
omt<-c('her.expr.t.1','her.expr.t.2','her.fmass.rv.t1','her.fmass.rv.t2','her.fmass.rv.t3','her.land.sprich.t.1','her.len.rv.t1','her.len.rv.t2','her.len.rv.t3','her.prod.t1','her.szpe.rv.t1','her.szpe.rv.t2','her.szpe.rv.t3','her.waa.t1','her.waa.t2','her.waa.t3','herjuv.fmass.rv.t1','herjuv.fmass.rv.t3','herjuv.fmass.rv.t3','t50.ct.t.1','t50.ct.t.2')
dat1<-dat1[,!(names(dat1) %in% omt)]


mod<-stepAIC(lm(her.ssbc.t3~.,data=dat1))


cr<-round(cor(dff,use='pairwise.complete.obs'),digits=2)
crt = as.data.frame(as.table(cr))#CORRELATIONS TABLE
combinations=combn(colnames(cr),2,FUN=function(x){paste(x,collapse="_")})
crt=crt[crt$Var1 != crt$Var2,]
crt=crt[paste(crt$Var1,crt$Var2,sep="_") %in% combinations,]
names(crt)<-tolower(names(crt))
a<-subset(crt,var1=='her.rec1.t3')
a<-subset(crt,var2=='her.rec1.t3')
a<-a[order(abs(a$freq),decreasing=TRUE),]
plot(dff$herlrv.mn.bof,dff$her.ssbc.t3,pch=15)
text(dff$herlrv.mn.bof,dff$her.rec1.t3,labels=dff$year,adj=0)


mod<-gam(her.len.rv.t3~  s(wind.amp,k=4)+ s(t50.ct,k=4)+ s(herlrv.len.t1,k=4)+s(her.state,k=4)+ s('asl.t-1',k=4),data=df)

mod<-gam(her.len.rv.t3~s(wind.amp,k=4)+s(wind.amp,k=4),data=mat)

blst<-subset(bl,to=='her.ssb.t3' & from!='her.ssb.t3')
dmat<-mat[,!(names(mat) %in% blst$from)]

mod<-stepAIC(lm(her.ssb.t3~.,data=dmat))











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




