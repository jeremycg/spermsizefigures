## ----,echo=FALSE,hide=TRUE,message=FALSE,warning=FALSE-------------------
#packages
library(phytools)
library(ape)
library(geiger)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(picante)
library(adephylo)
library(caper)
library(ggtree)
library(psych)
library(knitr)

## ----,echo=FALSE,hide=TRUE,message=FALSE---------------------------------
#tree
tr<-read.tree(text="(C. sp. 1:0.356,(C. plicata:0.3,((C. guadeloupensis:0.223,((C. portoensis:0.158,C. virilis:0.201):0.005,((C. sp. 8:0.169,(C. angaria:0.031,C. castelli:0.037):0.131):0.073,(C. drosophilae:0.027,C. sp. 2:0.03):0.161):0.047):0.033):0.074,((C. kamaaina:0.14,(C. japonica:0.16,((C. imperialis:0.082,C. afra:0.182):0.009,((C. yunquensis:0.086,C. nouraguensis:0.085):0.038,C. macrosperma:0.159):0.057):0.012):0.059):0.041,(C. elegans:0.166,(((C. brenneri:0.129,C. doughertyi:0.111):0.042,(C. wallacei:0.083,C. tropicalis:0.07):0.033):0.013,(C. remanei:0.129,(C. sp. 5:0.127,(C. briggsae:0.047,C. nigoni:0.032):0.088):0.02):0.017):0.046):0.02):0.195):0.07):0.259);")

tr$tip.label<-c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni")

## ----,echo=FALSE,hide=TRUE,message=FALSE---------------------------------
#read in a preprocess data
setwd("C:/Users/jeremy/Desktop/geiger")
spermsizes<-read.table("compileddata.csv",header=T,sep=",")
oocytesize<-read.table("oocyteSizeCaeno.csv",header=T,sep=",")
primaryspermatocyte<-read.table("PrimarySpermatocyteSize.csv",header=T,sep=",")
bodysize<-read.csv("CaenoBodySize.csv")
eggsize<-read.csv("eggvolumesarea.csv")
#treat as ellipses
eggsize$area<-0.5*eggsize$Length.Mean*0.5*eggsize$Width.Mean*pi
#remove ones we don't have phylogenetic data for
eggsize<-eggsize[!eggsize$Species %in% c("latens","Ocheius tipulae"),]
#ugh, bodysize has a row with pasted in headers at line 1507 - removed it in my version of the file
#ugh, anagaria has "PS1010" and "PS1010 " replaced in my file
oocytesize<-oocytesize[oocytesize$oocyte.stage==-1,]
primaryspermatocyte<-primaryspermatocyte[primaryspermatocyte$Sex=="male",]
oocytesizearea<-read.csv("oocyteSizearea.csv")
oocytesizearea<-oocytesizearea[oocytesizearea$oocyte.stage==-1,]

weneednewnames <- list("brenneri"="C. brenneri",
                       "C. brenneri"="C. brenneri",
                       "C. angaria"="C. angaria",
                       "angaria"="C. angaria",
                       "angaria "="C. angaria",
                       "C. briggsae"="C. briggsae",
                       "briggsae"="C. briggsae",
                       "C. drosophilae"="C. drosophilae",
                       "drosophilae"="C. drosophilae",
                       "C. elegans"="C. elegans",
                       "elegans"="C. elegans",
                       "C. japonica"="C. japonica",
                       "japonica"="C. japonica",
                       "C. plicata"="C. plicata",
                       "plicata"="C. plicata",
                       "C. remanei"="C. remanei",
                       "remanei"="C. remanei",
                       "C. sp1"="C. sp. 1",
                       "C. sp 1"="C. sp. 1",
                       "sp1"="C. sp. 1",
                       "C. sp8"="C. sp. 8",
                       "C. sp 8"="C. sp. 8",
                       "sp8"="C. sp. 8",
                       "C. sp10"="C. doughertyi",
                       "C. sp. 10"="C. doughertyi",
                       "sp10"="C. doughertyi",
                       "doughertyi"="C. doughertyi",
                       "C. sp.11"="C. tropicalis",
                       "C. sp11"="C. tropicalis",
                       "C. sp. 11"="C. tropicalis",
                       "sp11"="C. tropicalis",
                       "tropicalis"="C. tropicalis",
                       "C. sp12"="C. castelli",
                       "C. sp. 12"="C. castelli",
                       "sp12"="C. castelli",
                       "castelli"="C. castelli",
                       "C. sp13"="C. virilis",
                       "C. sp. 13"="C. virilis",
                       "sp13"="C. virilis",
                       "virilis"="C. virilis",
                       "C. sp14"="C. imperialis",
                       "C. sp. 14"="C. imperialis",
                       "sp14"="C. imperialis",
                       "imperialis"="C. imperialis",
                       "C. sp15"="C. kamaaina",
                       "C. sp. 15"="C. kamaaina",
                       "sp15"="C. kamaaina",
                       "kamaaina"="C. kamaaina",
                       "C. sp16"="C. wallacei",
                       "C. sp. 16"="C. wallacei",
                       "sp16"="C. wallacei",
                       "wallacei"="C. wallacei",
                       "C. sp17"="C. nouraguensis",
                       "C. sp. 17"="C. nouraguensis",
                       "sp17"="C. nouraguensis",
                       "nouraguensis"="C. nouraguensis",
                       "C. sp 18"="C. macrosperma",
                       "C. sp18"="C. macrosperma",
                       "C. sp. 18"="C. macrosperma",
                       "sp18"="C. macrosperma",
                       "macrosperma"="C. macrosperma",
                       "C. sp19"="C. yunquensis",
                       "C. sp. 19"="C. yunquensis",
                       "sp19"="C. yunquensis",
                       "yunquensis"="C. yunquensis",
                       "C. sp2"="C. sp. 2",
                       "C. sp. 2"="C. sp. 2",
                       "sp2"="C. sp. 2",
                       "C. sp. 20"="C. guadeloupensis",
                       "sp20"="C. guadeloupensis",
                       "guadeloupensis"="C. guadeloupensis",
                       "C. sp5"="C. sinica",
                       "C. sp. 5"="C. sinica",
                       "sp5"="C. sinica",
                       "sinica"="C. sinica",
                       "C. sp6"="C. portoensis",
                       "C. sp. 6"="C. portoensis",
                       "sp6"="C. portoensis",
                       "portoensis"="C. portoensis",
                       "C. sp7"="C. afra",
                       "C. sp. 7"="C. afra",
                       "sp7"="C. afra",
                       "afra"="C. afra",
                       "C. sp9"="C. nigoni",
                       "C. sp. 9"="C. nigoni",
                       "nigoni"="C. nigoni",
                       "sp9"="C. nigoni")

renameworms<-function(df){
  df$Species1<-as.character(df$Species)
  df$Species1<-weneednewnames[df$Species1]
  df$Species<-as.factor(unlist(df$Species1))
  df
}

bodysize<-renameworms(bodysize)
spermsizes<-renameworms(spermsizes)
oocytesize<-renameworms(oocytesize)
primaryspermatocyte<-renameworms(primaryspermatocyte)
oocytesizearea<-renameworms(oocytesizearea)
eggsize<-renameworms(eggsize)

clusterer<-function(df){
  x<-df %>% group_by(Species) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
  rownames(x) <- x$Species
  x
}

clusterer2<-function(df){
  x<-df %>% group_by(Species,Strain) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
  x
}

clusterer3<-function(df){
  x<-df %>% group_by(Species,Strain,individual) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
  x
}
spermsizelog<-spermsizes
spermsizelog$area<-log(spermsizes$area)
spermsize1<-clusterer(spermsizes)
spermsize2<-clusterer2(spermsizes)
spermsize3<-clusterer3(spermsizes)
spermsizelog1<-clusterer(spermsizelog)
spermsizelog2<-clusterer2(spermsizelog)
spermsizelog3<-clusterer3(spermsizelog)

spermvolume<-spermsizes
spermvolume$area<-4/3*pi*(0.5*(0.5*spermvolume$dia1.um+0.5*spermvolume$dia2.um))^3
spermvolume1<-clusterer(spermvolume)
spermvolume2<-clusterer2(spermvolume)
spermvolume3<-clusterer3(spermvolume)

names(oocytesize)<-c("Species","Strain","individual","oocyte.stage","length.um","width.um","area")
oocytesizelog<-oocytesize
oocytesizelog$area<-log(oocytesizelog$area)
oocytesize1<-clusterer(oocytesize)
oocytesize2<-clusterer2(oocytesize)
oocytesize3<-clusterer3(oocytesize)
oocytesizelog1<-clusterer(oocytesizelog)
oocytesizelog2<-clusterer2(oocytesizelog)
oocytesizelog3<-clusterer3(oocytesizelog)

names(primaryspermatocyte)<-c("Species","Strain","Sex","individual","area")
primaryspermatocytelog<-primaryspermatocyte
primaryspermatocytelog$area<-log(primaryspermatocytelog$area)
primaryspermatocyte1<-clusterer(primaryspermatocyte)
primaryspermatocyte2<-clusterer2(primaryspermatocyte)
primaryspermatocyte3<-clusterer3(primaryspermatocyte)
primaryspermatocytelog1<-clusterer(primaryspermatocytelog)
primaryspermatocytelog2<-clusterer2(primaryspermatocytelog)
primaryspermatocytelog3<-clusterer3(primaryspermatocytelog)

names(oocytesizearea)<-c("Species","Strain","individual","oocyte.stage","length.um","width.um","area")
oocytesizearealog<-oocytesizearea
oocytesizearealog$area<-log(oocytesizearealog$area)
oocytesizearea1<-clusterer(oocytesizearea)
oocytesizearea2<-clusterer2(oocytesizearea)
oocytesizearea3<-clusterer3(oocytesizearea)
oocytesizearealog1<-clusterer(oocytesizearealog)
oocytesizearealog2<-clusterer2(oocytesizearealog)
oocytesizearealog3<-clusterer3(oocytesizearealog)

bodysizelength<-bodysize[,1:4]
names(bodysizelength)<-c("Species","Strain","Sex","area")
bodysizelengthlog<-bodysizelength
bodysizelengthlog$area<-log(bodysizelengthlog$area)
bodysizelength1<-bodysizelength %>% group_by(Species,Strain,Sex) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
bodysizelengthlog1<-bodysizelengthlog %>% group_by(Species,Strain,Sex) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))

bodysizearea<-bodysize[,c(1:3,8)]
names(bodysizearea)<-c("Species","Strain","Sex","area")
bodysizearealog<-bodysizearea
bodysizearealog$area<-log(bodysizearealog$area)
bodysizearea1<-bodysizearea %>% group_by(Species,Strain,Sex) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
bodysizearealog1<-bodysizearealog %>% group_by(Species,Strain,Sex) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))

bodysizewidth<-bodysize[,c(1:3,7)]
names(bodysizewidth)<-c("Species","Strain","Sex","area")
bodysizewidthlog<-bodysizewidth
bodysizewidthlog$area<-log(bodysizewidthlog$area)
bodysizewidth1<-bodysizewidth %>% group_by(Species,Strain,Sex) %>%
  summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
bodysizewidthlog1<-bodysizewidthlog %>% group_by(Species,Strain,Sex) %>%
  summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))

eggvolume<-eggsize[,c(1,2,7)]
names(eggvolume)<-c("Species","Strain","area")
eggvolumelog<-eggvolume
eggvolumelog$area<-log(eggvolumelog$area)
#we have only means - assume equal repeats
eggvolume1<-eggvolume %>% group_by(Species) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
eggvolumelog1<-eggvolumelog %>% group_by(Species) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))

eggsize<-eggsize[,c(1,2,8)]
eggsizelog<-eggsize
eggsizelog$area<-log(eggsize$area)
#we have only means - assume equal repeats
eggsize1<-eggsize %>% group_by(Species) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
eggsizelog1<-eggsizelog %>% group_by(Species) %>%
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))

## ----,echo=FALSE,comment=""----------------------------------------------
print("these numbers are a little different - I'm using the data from the newer Farhadifar paper")
#need to have range of variation
print(c("the range of sperm sizes is",max(spermsize1$means)/min(spermsize1$means)))
#need to have % of egg size for different species
ztemp<-merge(eggsize1,spermsize1,by="Species")
ztemp$ratio<-ztemp$means.y/ztemp$means.x
print("NB this is the ratio of area to area, not volume to volume")
print(c("the max ratio of sperm to egg is:",max(ztemp$ratio)))
print(c("c. elegans ratio is:" ,ztemp$ratio[ztemp$Species=="C. elegans"]))
print("volume to volume")
ztemp<-merge(eggvolume1,spermvolume1,by="Species")
ztemp$ratio<-ztemp$means.y/ztemp$means.x
print(c("the max ratio of sperm to egg is:",max(ztemp$ratio)))
print(c("c. elegans ratio is:" ,ztemp$ratio[ztemp$Species=="C. elegans"]))

## ----,echo=FALSE,comment=""----------------------------------------------
print(c("the range of sperm sizes is",max(spermsize1$means)/min(spermsize1$means)))
print(c("C. briggsae sperm",spermsize1$means[spermsize1$Species=="C. briggsae"]))
print(c("C. plicata sperm",spermsize1$means[spermsize1$Species=="C. plicata"]))
#surface analysis
print("surface analysis is written up a little in the methods, might be worth a sentence or two here")

## ----,echo=FALSE,warning=FALSE-------------------------------------------
#figure 1a
#incredibly ugly code, but it works
ff<-read.table("C:/Users/jeremy/Desktop/geiger/compileddata.csv",header=T,sep=",")
ff<-renameworms(ff)
graphspermsize2<-ff %>% group_by(Species,Strain) %>% summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T),sem=mean(area,na.rm=T)/sqrt(n()))
countspec=function(dataframe){
  holding=c()
	for(i in 1:nrow(dataframe)){
		holding=c(holding,sum(dataframe$Species==dataframe[i,]$Species))
	}
	return(holding)
}
graphspermsize2$Species<-ordered(graphspermsize2$Species,levels=c(tr$tip.label))
graphspermsize2<-graphspermsize2[order(graphspermsize2$Species),]
graphspermsize2$Strain<-factor(graphspermsize2$Strain, levels = graphspermsize2$Strain)
graphspermsize2$specnum=countspec(graphspermsize2)
graphspermsize2$cols<-graphspermsize2$Species=="C. plicata"|graphspermsize2$Species=="C. drosophilae"|graphspermsize2$Species=="C. macrosperma"|graphspermsize2$Species=="C. japonica"
findchange<-function(df){
    holding=c()
    for(i in 1:length(df$Species)){
        if(i==1){
            holding=c(holding,0)
        } else if(df$Species[i]==df$Species[i-1]){
            holding=c(holding,0)
        } else {
            holding=c(holding,1)
        }
    }
	return(holding)
}

graphspermsize2$change<-findchange(graphspermsize2)
w<-(1-(graphspermsize2$specnum-1)*0.05)/graphspermsize2$specnum
pos <- 0.5 * (cumsum(w) + cumsum(c(0, w[-length(w)])))

gapsizes<-function(df,pos){
	pos1<-pos
	for(i in 1:length(df$specnum)){
		if(df$specnum[i]==1){
			pos1[i:length(pos1)]<-pos1[i:length(pos1)]+0.1
		} else if(df$change[i]==1){
			pos1[i:length(pos1)]<-pos1[i:length(pos1)]+0.1
		} else {
			pos1[i:length(pos1)]<-pos1[i:length(pos1)]+0.05
		}
	}
	return(pos1)
}

pos<-gapsizes(graphspermsize2,pos)
ggplot(data=graphspermsize2,aes(x = pos, width = w, y = means, fill=cols,ymin=means-sem,ymax=means+sem)) + scale_fill_brewer()+coord_cartesian(ylim=c(0,280))+
    geom_bar(stat = "identity",colour="black") + scale_x_continuous(labels = graphspermsize2$Strain, breaks = pos)+geom_linerange()+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(breaks= seq(from=0,to=280,by=10), 
                  labels = c(0,rep("",4), 50, rep("",4), 100,rep("",4),150,rep("",4),200,rep("",4),250,"","",""))

## ----,echo=FALSE,warning=FALSE-------------------------------------------
#figure 1b
spermsize1graph<-ff %>% group_by(Species) %>% summarise(means=mean(area,na.rm=T))
spermsize2graph<-ff %>% group_by(Species,Strain) %>% summarise(means=mean(area,na.rm=T))
spermsize3graph<-ff %>% group_by(Species,Strain,individual) %>% summarise(means=mean(area,na.rm=T))
fitace <- function(spermsizes,tree){
  holding<-spermsizes$means
  names(holding)<-spermsizes$Species
  return(fastAnc(tree,holding,CI=T))
}
x<-fitace(spermsize1graph,tr)
AncSperm=x$ace
#or AncSperm=apply(x,2,mean) if it's bootstrapped
spermsize1graph$Species <- factor(spermsize1graph$Species, levels = c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni"))
spermsizegraph=spermsize1graph[order(spermsize1graph$Species),]$means
names(spermsizegraph)<-spermsize1graph[order(spermsize1graph$Species),]$Species
plotfancytree<-function(tree,nodesizes,terminalsizes){
  nodesizes<<-nodesizes
  terminalsizes<<-terminalsizes
  print(ggtree(tree)+geom_text(subset=.(isTip),aes(label=label),hjust=-0.19,fontface="italic")+
    geom_point(subset=.(!isTip),color="#ffae1a", size=sqrt(nodesizes/pi))+
    geom_point(subset=.(isTip),color="#7570b3", size=sqrt(terminalsizes/pi))+
    geom_text(subset=.(isTip),aes(label=round(terminalsizes,1)), hjust=1, vjust=-0.4, size=3)+
    geom_text(subset=.(!isTip),aes(label=round(nodesizes,1)), hjust=1.5, vjust=-0.4, size=3)+
    scale_x_continuous(expand = c(0.1, 0.1)))
  rm(list=c("terminalsizes","nodesizes"),pos=".GlobalEnv")
}
plotfancytree(tr,AncSperm,spermsizegraph)

## ----,echo=FALSE,warning=FALSE-------------------------------------------
hist(spermsize1$means,breaks=13,main="",xlab=expression(paste("Sperm Size (", mu, m^{2},")")),ylab="Number of species",col="gray")
minor.tick(ny=5)

## ----,echo=FALSE,warning=FALSE-------------------------------------------
#sperm versus body size in males
pglsdata<-function(df1,df2,tree){
  df<-merge(df1,df2,by="Species")
  rownames(df)<-df$Species
  tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% df$Species)])
  df <- df[match(tree$tip.label,rownames(df)),]
  out<-as.data.frame(summary(gls(means.y ~ means.x, data=df,correlation=corBrownian(1,tree)))$tTable)
  out$p[2]
}
print(c("p-value from pgls of sperm to male body size (volume):",pglsdata(spermsizelog1,bodysizearealog1[bodysizearealog1$Sex=="male",],tr)))
print(c("p-value from pgls of sperm to male body size (length):",pglsdata(spermsizelog1,bodysizelengthlog1[bodysizelengthlog1$Sex=="male",],tr)))
print(c("p-value from pgls of sperm to male body size (width):",pglsdata(spermsizelog1,bodysizewidthlog1[bodysizewidthlog1$Sex=="male",],tr)))
#sperm vs egg variation
print(c("p-value from pgls of sperm to egg size (from Farhadifar):",pglsdata(spermsizelog1,eggsizelog1,tr)))
print(c("p-value from pgls of sperm to oocyte size:",pglsdata(spermsizelog1,oocytesizearealog1,tr)))
#volumetric ratios
ztemp<-merge(spermsize1,eggsize1,by="Species")
ztemp$ratio<-ztemp$means.x/ztemp$means.y
print(c("sperm:egg ratio: Farhadifar elegans",1/ztemp$ratio[ztemp$Species=="C. elegans"]))
print(c("sperm:egg ratio: Farhadifar drosophilae",1/ztemp$ratio[ztemp$Species=="C. drosophilae"]))
ztemp<-merge(spermsize1,oocytesizearea1,by="Species")
ztemp$ratio<-ztemp$means.x/ztemp$means.y
print(c("sperm:oocyte ratio: elegans",1/ztemp$ratio[ztemp$Species=="C. elegans"]))
print(c("sperm:oocyte: drosophilae",1/ztemp$ratio[ztemp$Species=="C. drosophilae"]))

## ----},echo=FALSE,warning=FALSE------------------------------------------
print("I think they did the analysis on sperm transfer as they had very few species pairs - happy to do PIC otherwise")
print(c("p-value from pgls of sperm to male body size (width):",pglsdata(spermsizelog1,bodysizewidthlog1[bodysizewidthlog1$Sex=="male",],tr)))

## ----,echo=FALSE---------------------------------------------------------
plotpic2<-function(df1,df2,tree,title,xlab,ylab){
  df<-merge(df1,df2,by="Species")
  rownames(df)<-df$Species
  tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% df$Species)])
  df <- df[match(tree$tip.label,rownames(df)),]
  rownames(df)<-df$Species
  x1 <- pic(df$means.x, tree)
  y1 <- pic(df$means.y, tree)
  plot(y1~x1,xlab=xlab,ylab=ylab,main=title,pch=19)
  fit1=lm(y1~x1)
  abline(fit1$coefficients,lwd=2,col=8)
  p=pglsdata(df1,df2,tree)
  rp = vector('expression',1)
  rp[2]=substitute(expression(italic(p) == MYVALUE), list(MYVALUE = format(p,dig=3)))[2]
  legend("bottomright",bty="n",legend=rp)
}
plotpic2(eggsizelog1,spermsizelog1,tr,"Sperm vs Eggsize","oocyte","sperm")

## ----,echo=FALSE---------------------------------------------------------
plotpic2(bodysizelengthlog1[bodysizelengthlog1$Sex=="female",],eggsizelog1,tr,"oocyte vs body length","body length","oocyte size")

## ----,echo=FALSE---------------------------------------------------------
ztemp<-merge(eggvolume1,spermvolume1,by="Species")
ztemp$ratio<-ztemp$means.y/ztemp$means.x
ggplot(ztemp,aes(reorder(Species, -ratio), ratio))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_continuous(limits = c(0, 0.05),breaks=seq(from=0,to=0.05,by=0.001),labels = c(0,rep("",4), 0.005, rep("",4), 0.010,rep("",4),0.015,rep("",4),0.020,rep("",4),0.0025,rep("",4),0.030,rep("",4),0.035,rep("",4),0.040,rep("",4),0.045,rep("",4),0.050),expand = c(0, 0))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+xlab("Species")+ylab("Ratio of Sperm volume to Egg volume")

## ----,echo=FALSE---------------------------------------------------------
plotpic<-function(df,tree,title,xlab,ylab){
  tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% df$Species)])
  rownames(df)<-df$Species
  df <- df[match(tree$tip.label,rownames(df)),]
  rownames(df)<-df$Species
  x1 <- pic(df$means, tree)
  y1 <- pic(df$cv, tree)
  plot(y1~x1,xlab=xlab,ylab=ylab,main=title,pch=19)
  fit1=lm(y1~x1)
  abline(fit1$coefficients,lwd=2,col=8)
  p=as.data.frame(summary(gls(means ~ cv, data=df,correlation=corBrownian(1,tree)))$tTable)$p[2]
  rp = vector('expression',1)
  rp[2]=substitute(expression(italic(p) == MYVALUE), list(MYVALUE = format(p,dig=3)))[2]
  legend("bottomright",bty="n",legend=rp)
}
plotpic(spermsizelog1,tr,"Sperm size VS Sperm CV","sperm size","sperm CV")

## ----,echo=FALSE---------------------------------------------------------
print(c("pgls of sperm size vs spermatocyte size",pglsdata(spermsizelog1,primaryspermatocytelog1,tr)))

## ----,echo=FALSE---------------------------------------------------------
ggplot(primaryspermatocyte1,aes(x=Species,y=means,ymin=c(means-cv),ymax=c(means+cv)))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_linerange()

## ----,echo=FALSE---------------------------------------------------------
plotpic2(spermsizelog1,primaryspermatocytelog1,tr,"sperm vs primary spermatocyte","sperm","spermatocyte")

## ----,echo=FALSE---------------------------------------------------------
violinit<-function(f,specieslist){
  g=f[f$Species%in%specieslist,]
  p<-ggplot(g,aes(factor(Strain), area))
  p + facet_wrap(~Species,ncol = 4,scales = "free_x")+geom_violin(aes(fill = Species)) + geom_jitter(position =
    position_jitter(width = .2),alpha=0.5,)+theme(legend.position = "none") +
    ylab( expression(paste("area (", mu, m^{2},")")))+xlab("Strain")+
    theme(panel.margin = unit(0.05, "lines"))+theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_y_continuous(breaks=seq(from=0,to=280,by=10), 
                  labels = c(0,rep("",4), 50, rep("",4), 100,rep("",4),150,rep("",4),200,rep("",4),250,"","","")) 
}

violinit(spermsizes,c("C. brenneri","C. remanei","C. sp. 8","C. macrosperma"))

## ----,echo=FALSE---------------------------------------------------------
clumper<-function(species,dfindv,dfstrain,dfspecies){
  holding=dfindv[dfindv$Species==species,]
  holding=rbind(holding,dfstrain[dfstrain$Species==species,])
  holding=rbind(holding,dfspecies[dfspecies$Species==species,])
  holding$plotter<-factor(holding$plotter,levels=c("intra individual","between individuals","between strains"))
  holding$plotcolour<-as.numeric(factor(holding$Strain))
  return(holding)
}

plotcvbysplits<-function(df,strainlist){
  dfstrain<-df%>% summarise(cv=sd(means)/mean(means),means=mean(means),plotter="between individuals")
  dfstrain$individual<-"dropped"
  dfspecies<-dfstrain %>% summarise(cv=sd(means)/mean(means),means=mean(means),plotter="between strains")
  dfspecies$individual="dropped"
  dfspecies$Strain="dropped"
  dfindv<-df
  dfindv$plotter<-"intra individual"
  x<-clumper(strainlist[1],dfindv,dfstrain,dfspecies)
  if(length(strainlist)>1){
    for(i in 2:length(strainlist)){
      x=rbind(x,clumper(strainlist[i],dfindv,dfstrain,dfspecies))
    }
  }
  p<-ggplot(x,aes(factor(plotter),cv))
  p+geom_jitter(aes(colour = factor(plotcolour),size = 5),position = position_jitter(width = .2))+facet_wrap(~Species,ncol = length(strainlist),scales = "free_x")+
    theme(panel.margin = unit(0.05, "lines"),legend.position = "none")+xlab("")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_y_continuous(limits = c(0, 0.45),breaks=seq(from=0,to=0.45,by=0.01),labels = c(0,rep("",4), 0.5, rep("",4), 1.0,rep("",4),1.5,rep("",4),2.0,rep("",4),2.5,rep("",4),3.0,rep("",4),3.5,rep("",4),4.5,rep("",4),5.5))
}
plotcvbysplits(spermsize3,c("C. brenneri","C. remanei","C. sp. 8","C. macrosperma"))

## ----,echo=FALSE---------------------------------------------------------
ff<-read.table("C:/Users/jeremy/Desktop/geiger/malevhermcomp.csv",header=T,sep=",")
ff<-renameworms(ff)
groupedmeans=ff %>% group_by(Species,sex) %>% summarise(mean=mean(area,na.rm=T),n=n(),sd=sd(area,na.rm=T))
limits <- aes(ymax = mean + 1.96*sd/sqrt(n), ymin=mean - 1.96*sd/sqrt(n))
dodge <- position_dodge(width=0.9)
ggplot(data=groupedmeans, aes(x=Species, y=mean, fill=sex)) + geom_bar(stat="identity", position=position_dodge(), colour="black")+geom_errorbar(limits,position=dodge, width=0.25)+
  scale_fill_manual(values=c("gray","white"))+ylab( expression(paste("area (", mu, m^{2},")")))+xlab("Species")+coord_cartesian(ylim=c(0,25))+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = "none")+scale_y_continuous(limits = c(0, 25),breaks=seq(from=0,to=25,by=1),labels = c(0,rep("",4), 5, rep("",4), 10,rep("",4),15,rep("",4),20,rep("",4),25))

## ----,echo=FALSE---------------------------------------------------------
groupedstrains=ff %>% group_by(Species,Strain) %>%
    do(data.frame(meanM=mean(.[which(.$sex=="M"),]$area,na.rm=T),meanH=mean(.[which(.$sex=="H"),]$area,na.rm=T),nM=sum(.$sex=="M"),nH=sum(.$sex=="H"),
                  sdM=sd(.[which(.$sex=="M"),]$area,na.rm=T),sdH=sd(.[which(.$sex=="H"),]$area,na.rm=T)))
groupedstrains$sem=1.96*groupedstrains$sdM/sqrt(groupedstrains$nM)
groupedstrains$seh=1.96*groupedstrains$sdH/sqrt(groupedstrains$nH)
p=ggplot(data=groupedstrains, aes(x=meanM, y=meanH, colour=Species, ymin = meanH - seh,ymax=meanH + seh,xmin = meanM - sem,xmax=meanM + sem))
p+geom_point(size=5)+coord_cartesian(ylim=c(0,15),xlim=c(0,33))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = "none")+ geom_errorbar(width=0.5)+geom_errorbarh(height=0.5)+
   ylab(expression(paste("Hermaphrodite Sperm area (", mu, m^{2},")")))+xlab(expression(paste("Male Sperm area (", mu, m^{2},")")))+
	 geom_smooth(aes(group=Species), method="lm",size=1.5,colour="black")+scale_y_continuous(limits = c(0, 25),breaks=seq(from=0,to=25,by=1),labels = c(0,rep("",4), 5, rep("",4), 10,rep("",4),15,rep("",4),20,rep("",4),25))+scale_x_continuous(breaks=seq(from=0,to=35,by=1),labels = c(0,rep("",4), 5, rep("",4), 10,rep("",4),15,rep("",4),20,rep("",4),25,rep("",4),30,rep("",4),35))

## ----,echo=FALSE---------------------------------------------------------
bodysize<-read.csv("C:/Users/jeremy/Desktop/geiger/CaenoBodySize.csv")
bodysizelength<-bodysize[,1:4]
bodysizelength<-renameworms(bodysizelength)
names(bodysizelength)<-c("Species","Strain","Sex","area")
ff<-bodysizelength %>% group_by(Species,Strain,Sex) %>%
  summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
tr<-read.tree(text="(C. sp. 1:0.356,(C. plicata:0.3,((C. guadeloupensis:0.223,((C. portoensis:0.158,C. virilis:0.201):0.005,((C. sp. 8:0.169,(C. angaria:0.031,C. castelli:0.037):0.131):0.073,(C. drosophilae:0.027,C. sp. 2:0.03):0.161):0.047):0.033):0.074,((C. kamaaina:0.14,(C. japonica:0.16,((C. imperialis:0.082,C. afra:0.182):0.009,((C. yunquensis:0.086,C. nouraguensis:0.085):0.038,C. macrosperma:0.159):0.057):0.012):0.059):0.041,(C. elegans:0.166,(((C. brenneri:0.129,C. doughertyi:0.111):0.042,(C. wallacei:0.083,C. tropicalis:0.07):0.033):0.013,(C. remanei:0.129,(C. sp. 5:0.127,(C. briggsae:0.047,C. nigoni:0.032):0.088):0.02):0.017):0.046):0.02):0.195):0.07):0.259);")
tr$tip.label<-c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni")

tr1<-drop.tip(tr,tr$tip.label[!(tr$tip.label %in% ff$Species)])
fff<-ff[ff$Sex=="female",]
ffm<-ff[ff$Sex=="male",]
xf<-fitace(fff,tr)
xm<-fitace(ffm,tr)
AncSizef=xf$ace
AncSizem=xm$ace

fff$Species <- factor(fff$Species, levels = c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni"))
fff1=fff[order(fff$Species),]$means
ffm$Species <- factor(ffm$Species, levels = c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni"))
ffm1=ffm[order(ffm$Species),]$means
names(fff1)<-fff[order(fff$Species),]$Species
names(ffm1)<-ffm[order(ffm$Species),]$Species
plotfancytree1<-function(tree,nodesizes,terminalsizes){
  nodesizes<<-nodesizes
  terminalsizes<<-terminalsizes
  print(ggtree(tree)+geom_text(subset=.(isTip),aes(label=label),hjust=-0.19,fontface="italic")+
          geom_point(subset=.(!isTip),color="#ffae1a", size=nodesizes/80,shape=126)+
          geom_point(subset=.(isTip),color="#7570b3", size=terminalsizes/80,shape=126)+
          geom_text(subset=.(isTip),aes(label=round(terminalsizes,1)), hjust=1, vjust=-0.4, size=3)+
          geom_text(subset=.(!isTip),aes(label=round(nodesizes,1)), hjust=1.5, vjust=-0.4, size=3)+
          scale_x_continuous(expand = c(0.1, 0.1)))
  rm(list=c("terminalsizes","nodesizes"),pos=".GlobalEnv")
}


## ----,echo=FALSE---------------------------------------------------------
plotfancytree1(tr1,AncSizef,fff1)

## ----,echo=FALSE---------------------------------------------------------
plotfancytree1(tr1,AncSizem,ffm1)

## ----,echo=FALSE---------------------------------------------------------
eggsizetograph1<-eggsize1
tr1<-drop.tip(tr,tr$tip.label[!(tr$tip.label %in% eggsizetograph1$Species)])
x<-fitace(eggsizetograph1,tr1)
AncSize<-x$ace
eggsizetograph1$Species <- factor(eggsizetograph1$Species, levels = c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni"))
eggsizetograph11=eggsizetograph1[order(eggsizetograph1$Species),]$means
names(eggsizetograph11)<-eggsizetograph1[order(eggsizetograph1$Species),]$Species
plotfancytree2<-function(tree,nodesizes,terminalsizes){
  nodesizes<<-nodesizes
  terminalsizes<<-terminalsizes
  print(ggtree(tree)+geom_text(subset=.(isTip),aes(label=label),hjust=-0.19,fontface="italic")+
    geom_point(subset=.(!isTip),color="#ffae1a", size=sqrt(nodesizes/pi)/10)+
    geom_point(subset=.(isTip),color="#7570b3", size=sqrt(terminalsizes/pi)/10)+
    geom_text(subset=.(isTip),aes(label=round(terminalsizes,1)), hjust=1, vjust=-0.4, size=3)+
    geom_text(subset=.(!isTip),aes(label=round(nodesizes,1)), hjust=1.5, vjust=-0.4, size=3)+
    scale_x_continuous(expand = c(0.1, 0.1)))
  rm(list=c("terminalsizes","nodesizes"),pos=".GlobalEnv")
}
plotfancytree2(tr1,AncSize,eggsizetograph11)

## ----,echo=FALSE,warning=FALSE,message=FALSE-----------------------------
library(surface)

tr<-read.tree(text="(C. sp. 1:0.356,(C. plicata:0.3,((C. guadeloupensis:0.223,((C. portoensis:0.158,C. virilis:0.201):0.005,((C. sp. 8:0.169,(C. angaria:0.031,C. castelli:0.037):0.131):0.073,(C. drosophilae:0.027,C. sp. 2:0.03):0.161):0.047):0.033):0.074,((C. kamaaina:0.14,(C. japonica:0.16,((C. imperialis:0.082,C. afra:0.182):0.009,((C. yunquensis:0.086,C. nouraguensis:0.085):0.038,C. macrosperma:0.159):0.057):0.012):0.059):0.041,(C. elegans:0.166,(((C. brenneri:0.129,C. doughertyi:0.111):0.042,(C. wallacei:0.083,C. tropicalis:0.07):0.033):0.013,(C. remanei:0.129,(C. sp. 5:0.127,(C. briggsae:0.047,C. nigoni:0.032):0.088):0.02):0.017):0.046):0.02):0.195):0.07):0.259);")

tr$tip.label<-c("C. sp. 1","C. plicata","C. guadeloupensis","C. portoensis","C. virilis","C. sp. 8","C. angaria","C. castelli","C. drosophilae","C. sp. 2","C. kamaaina","C. japonica","C. imperialis","C. afra","C. yunquensis","C. nouraguensis","C. macrosperma","C. elegans","C. brenneri","C. doughertyi","C. wallacei","C. tropicalis","C. remanei","C. sinica","C. briggsae","C. nigoni")
tr<-nameNodes(tr)
spermsizes<-read.table("C:/Users/jeremy/Desktop/geiger/compileddata.csv",header=T,sep=",")
spermsizes<-renameworms(spermsizes)
clusterer<-function(df){
  x<-df %>% group_by(Species) %>% 
    summarise(means=mean(area,na.rm=T),cv=sd(area,na.rm=T)/mean(area,na.rm=T))
  rownames(x) <- x$Species
  x
}
surfacesperm<-clusterer(spermsizes)

rownames(surfacesperm)<-surfacesperm$Species
surfacesperm <- surfacesperm[match(tr$tip.label,rownames(surfacesperm)),]
valuesin<-surfacesperm$means
names(valuesin)<-surfacesperm$Species
z<-runSurface(nameNodes(tr),as.data.frame(valuesin))

runonesim<-function(tr,fit,shifts){
  runSurface(tr,surfaceSimulate(tr, type="hansen-fit", hansenfit=fit, 
                  shifts=shifts, sample_optima=TRUE)$dat)
}


getvals<-function(tr,fit,shifts){
  sim<-runonesim(tr,fit,shifts)
  lengthofsim<-length(sim$bwd)
  sim$bwd[[lengthofsim]]$n_regimes
}

distributionofvals<-function(n,tr,fit,shifts){
  rdply(n,getvals(tr,fit,shifts),.progress = "text")
}

surfaceTreePlot(tr,z$bwd[[2]])
surfaceAICPlot(fwd=z$fwd,bwd=z$bwd)


## ----,echo=FALSE---------------------------------------------------------
plotpic2(oocytesizearealog1,spermsizelog1,tr,"Sperm vs Eggsize","oocyte","sperm")

## ----,echo=FALSE---------------------------------------------------------
plotpic2(bodysizelengthlog1[bodysizelengthlog1$Sex=="female",],oocytesizearealog1,tr,"oocyte vs body volume","body length","oocyte size")

## ----,echo=FALSE---------------------------------------------------------
plotpictable<-function(df,tree){
  tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% df$Species)])
  rownames(df)<-df$Species
  df <- df[match(tree$tip.label,rownames(df)),]
  rownames(df)<-df$Species
  x1 <- pic(df$means, tree)
  y1 <- pic(df$cv, tree)
  summary(gls(means ~ cv, data=df,correlation=corBrownian(1,tree)))$tTable
}
kable(plotpictable(spermsizelog1,tr),digits=3,caption="PGLS,all")
spermsizelog1noherm<-spermsizelog1[spermsizelog1$Species=="C. elegans"|spermsizelog1$Species=="C. briggsae"|spermsizelog1$Species=="C. tropicalis",]
kable(plotpictable(spermsizelog1noherm,tr),digits=3,caption="PGLS, no herm")
plotnopictable<-function(df,tree){
  tree<-drop.tip(tree,tree$tip.label[!(tree$tip.label %in% df$Species)])
  rownames(df)<-df$Species
  df <- df[match(tree$tip.label,rownames(df)),]
  rownames(df)<-df$Species
  x1 <- pic(df$means, tree)
  y1 <- pic(df$cv, tree)
  summary(gls(means ~ cv, data=df))$tTable
}
kable(plotnopictable(spermsizelog1,tr),digits=3,caption="GLS,all")
kable(plotnopictable(spermsizelog1noherm,tr),digits=3,caption="GLS, no herm")


