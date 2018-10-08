###############################################################################
###############################################################################
#Analizing Helmintosporiose strain resistance tests from Arvalis field trial
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#load the global dataset
creadat<-read.table("creative_data.txt",header=T,sep="\t")


###############################################################################
#Analysis for the boscalid
###############################################################################

#subsetting the global dataset
bosc.dat<-creadat[creadat$active_substance=="boscalid",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
bosc_rez<-as.character(bosc.dat[bosc.dat$dose=="30" & bosc.dat$perc_croiss>50,
                                "sample_ID"])
REZbos<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#bosc.dat<-bosc.dat[!(bosc.dat$sample_ID %in% bosc_rez),]
bosc.dat<-drop.levels(bosc.dat)
pdf("boscalid.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(bosc.dat$strain_ID))[1]) {
  datatemp<-bosc.dat[bosc.dat$strain_ID==names(table(bosc.dat$strain_ID))[i],]
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZbos<-rbind(REZbos,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30), 
         main=names(table(bosc.dat$strain_ID))[i])
  } else {
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),
         main=names(table(bosc.dat$strain_ID))[i])
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
                      "ED50"=temp[1])
    REZbos<-rbind(REZbos,tempx)
  }
}
dev.off()

REZbos$ED50[REZbos$ED50>30]<-30
plot(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZbos,file="REZbos.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
abline(v=10,col="red",lwd=3)
#export to pdf 4.5 x 9 inches


###############################################################################
#Analysis for the dodine
###############################################################################

#subsetting the global dataset
dodi.dat<-creadat[creadat$active_substance=="dodine",]

#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
dodi_rez<-as.character(dodi.dat[dodi.dat$dose=="30" & dodi.dat$perc_croiss>50,
                                "sample_ID"])
REZdod<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
#we limit the dataset to the sample that reach somehow a IC of 50%
#dodi.dat<-dodi.dat[!(dodi.dat$sample_ID %in% bosc_rez),]
dodi.dat<-drop.levels(dodi.dat)
pdf("dodine.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(dodi.dat$strain_ID))[1]) {
  datatemp<-dodi.dat[dodi.dat$strain_ID==names(table(dodi.dat$strain_ID))[i],]
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(dodi.dat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZbos<-rbind(REZbos,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30), 
         main=names(table(dodi.dat$strain_ID))[i])
  } else {
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),
         main=names(table(dodi.dat$strain_ID))[i])
    temp<-ED(temp.m1,50,type="absolute")
    tempx<-data.frame("strain_ID"=names(table(dodi.dat$strain_ID))[i],
                      "ED50"=temp[1])
    REZbos<-rbind(REZbos,tempx)
  }
}
dev.off()

REZbos$ED50[REZbos$ED50>30]<-30
plot(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="Souches ID",
     ylab="FR",las=1)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)
#export to pdf 10 x 6 inches
write.table(REZbos,file="REZbos.txt",quote=FALSE,sep="\t",row.names=FALSE)

hist(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="FR Classes",
     breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
     las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
abline(v=10,col="red",lwd=3)
#export to pdf 4.5 x 9 inches



###############################################################################
#combined plot
###############################################################################

op<-par(mfrow=c(2,2))

plot(REZbos$ED50[order(REZbos$ED50)],main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)

plot(REZbix$ED50[order(REZbix$ED50)],main="Bixafen")
abline(0.08,0,col="green3",lwd=2)
abline(0.8,0,col="orange3",lwd=2)

plot(REZflo$ED50[order(REZflo$ED50)],main="Fluopyram")
abline(0.44,0,col="green3",lwd=2)
abline(4.4,0,col="orange3",lwd=2)

plot(REZflx$ED50[order(REZflx$ED50)],main="Fluxapyroxade")
abline(0.21,0,col="green3",lwd=2)
abline(2.1,0,col="orange3",lwd=2)

par(op)


###############################################################################
#plot with individual group by pop
###############################################################################

EC50_pop<-read.table("EC50_byPOP.txt",header=TRUE)

op<-par(mfrow=c(2,2),mar=c(2,2.5,3,1))
EC50bosc<-EC50_pop[EC50_pop$SA_ID=="boscalid",]
plot(EC50bosc$ED50,col=as.character(EC50bosc$pop_col),main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)

EC50bixa<-EC50_pop[EC50_pop$SA_ID=="bixafen",]
plot(EC50bixa$ED50,col=as.character(EC50bixa$pop_col),main="Bixafen")
abline(0.08,0,col="green3",lwd=2)
abline(0.8,0,col="orange3",lwd=2)

EC50fluo<-EC50_pop[EC50_pop$SA_ID=="fluopyram",]
plot(EC50fluo$ED50,col=as.character(EC50fluo$pop_col),main="Fluopyram")
abline(0.44,0,col="green3",lwd=2)
abline(4.4,0,col="orange3",lwd=2)

EC50flux<-EC50_pop[EC50_pop$SA_ID=="fluxapyroxade",]
plot(EC50flux$ED50,col=as.character(EC50flux$pop_col),main="Fluxapyroxade")
abline(0.21,0,col="green3",lwd=2)
abline(2.1,0,col="orange3",lwd=2)

par(op)


###############################################################################
#END
###############################################################################
