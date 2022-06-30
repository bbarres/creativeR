##############################################################################/
##############################################################################/
#Code for different figures with IC50
##############################################################################/
##############################################################################/

##loading the data set and the necessary library
source("creative_load_data.R")

#creating a list of correspondence between strain ID and species
speStrai<-unique(creadat[,c(1,3)])

#in order for this code to work, you need first to run the script for 
#IC50 estimation
alterCI50<-read.table("output/Alternaria_mycelial.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)

VenGerCI50<-read.table("output/Venturia_spore.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)

VenMycCI50<-read.table("output/Venturia_mycelial.txt",header=TRUE,
                       sep="\t",stringsAsFactors=TRUE)


##############################################################################/
#barplot to compare the ED50 of the different samples for Alternaria####
##############################################################################/

#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
#adding the species name to the data
alterCI50<-merge(alterCI50,speStrai)
alterCI50$ED50.abs<-as.character(alterCI50$ED50.abs)
alterCI50[alterCI50$ED50.abs==">30","ED50.abs"]<-32
alterCI50$ED50.abs<-as.numeric(as.character(alterCI50$ED50.abs))


#boscalid
byprod<-alterCI50[alterCI50$ActiveSub==levels(alterCI50$ActiveSub)[1],]
activSub<-levels(alterCI50$ActiveSub)[1]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 4 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:10])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
AltFR<-byprod[,c(1:4,24:26)]


#cyprodinil
byprod<-alterCI50[alterCI50$ActiveSub==levels(alterCI50$ActiveSub)[4],]
activSub<-levels(alterCI50$ActiveSub)[4]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 4 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:10])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
AltFR<-rbind(AltFR,byprod[,c(1:4,24:26)])
plot(AltFR[AltFR$ActiveSub=="cyprodinil","factres"])


#difenoconazole
byprod<-alterCI50[alterCI50$ActiveSub==levels(alterCI50$ActiveSub)[5],]
activSub<-levels(alterCI50$ActiveSub)[5]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 10 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:10])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
AltFR<-rbind(AltFR,byprod[,c(1:4,24:26)])
plot(AltFR[AltFR$ActiveSub=="difenoconazole","factres"])


#dodine
byprod<-alterCI50[alterCI50$ActiveSub==levels(alterCI50$ActiveSub)[7],]
activSub<-levels(alterCI50$ActiveSub)[7]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 10 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
AltFR<-rbind(AltFR,byprod[,c(1:4,24:26)])
plot(AltFR[AltFR$ActiveSub=="dodine","factres"])


#fluopyram
byprod<-alterCI50[alterCI50$ActiveSub==levels(alterCI50$ActiveSub)[8],]
activSub<-levels(alterCI50$ActiveSub)[8]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 10 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:10])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
AltFR<-rbind(AltFR,byprod[,c(1:4,24:26)])
plot(AltFR[AltFR$ActiveSub=="fluopyram","factres"])

AltFR<-drop.levels(AltFR,reorder=FALSE)


##############################################################################/
#barplot to compare the ED50 of the different samples for Venturia germin####
##############################################################################/

#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
#adding the species name to the data
VenGerCI50<-merge(VenGerCI50,speStrai)
VenGerCI50$ED50.abs<-as.character(VenGerCI50$ED50.abs)
VenGerCI50[VenGerCI50$ED50.abs==">30","ED50.abs"]<-32
VenGerCI50$ED50.abs<-as.numeric(as.character(VenGerCI50$ED50.abs))


#boscalid
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[1],]
activSub<-levels(VenGerCI50$ActiveSub)[1]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,35),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:15])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-byprod[,c(1:4,24:26)]
plot(VenGFR[VenGFR$ActiveSub=="boscalid","factres"])


#captane
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[2],]
activSub<-levels(VenGerCI50$ActiveSub)[2]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,10),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:15])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="captane","factres"])


#cyprodinil
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[3],]
activSub<-levels(VenGerCI50$ActiveSub)[3]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:15])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="cyprodinil","factres"])


#difénoconazole
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[4],]
activSub<-levels(VenGerCI50$ActiveSub)[4]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="difenoconazole","factres"])


#dithianon
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[5],]
activSub<-levels(VenGerCI50$ActiveSub)[5]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1.5),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:15])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="dithianon","factres"])


#dodine
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[6],]
activSub<-levels(VenGerCI50$ActiveSub)[6]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,0.5),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="dodine","factres"])


#mancozeb
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[7],]
activSub<-levels(VenGerCI50$ActiveSub)[7]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="mancozeb","factres"])


#manebe
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[8],]
activSub<-levels(VenGerCI50$ActiveSub)[8]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="manebe","factres"])


#thirame
byprod<-VenGerCI50[VenGerCI50$ActiveSub==levels(VenGerCI50$ActiveSub)[9],]
activSub<-levels(VenGerCI50$ActiveSub)[9]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,1),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:15])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenGFR<-rbind(VenGFR,byprod[,c(1:4,24:26)])
plot(VenGFR[VenGFR$ActiveSub=="thirame","factres"])

VenGFR<-drop.levels(VenGFR,reorder=FALSE)


##############################################################################/
#barplot to compare the ED50 of the different samples for Venturia mycelial####
##############################################################################/

#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
#adding the species name to the data
VenMycCI50<-merge(VenMycCI50,speStrai)
VenMycCI50$ED50.abs<-as.character(VenMycCI50$ED50.abs)
VenMycCI50[VenMycCI50$ED50.abs==">30","ED50.abs"]<-32
VenMycCI50$ED50.abs<-as.numeric(as.character(VenMycCI50$ED50.abs))


#difenoconazole
byprod<-VenMycCI50[VenMycCI50$ActiveSub==levels(VenMycCI50$ActiveSub)[3],]
activSub<-levels(VenMycCI50$ActiveSub)[3]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,5),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenMFR<-byprod[,c(1:4,24:26)]
plot(VenMFR[VenMFR$ActiveSub=="difenoconazole","factres"])


#trifloxystrobine
byprod<-VenMycCI50[VenMycCI50$ActiveSub==levels(VenMycCI50$ActiveSub)[5],]
activSub<-levels(VenMycCI50$ActiveSub)[5]
#creating categories of CI50
byprod$catCI50<-cut(byprod$ED50.abs,
                    breaks=c(0,4,8,12,16,20,24,28,32),
                    include.lowest=TRUE)
#defining the colors of the points
levels(byprod$catCI50)<-brewer.pal(11,"RdYlGn")[8:1]
op<-par(mfrow=c(2,1))
hist(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     breaks=seq(0,32,by=2),bty="l",freq=FALSE,las=1,
     main=activSub,xlim=c(0,32),
     col=brewer.pal(11,"RdYlGn")[rep(8:1,each=2)],
     xlab="CI50 Class",ylab="Pourcentage")
box(bty="l")
plot(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))]),
     bg=as.character(byprod$catCI50[order(as.numeric(byprod$ED50.abs))]),
     cex=1.5,las=1,ylim=c(0,5),pch=21,
     ylab="CI50",xlab="",
     main="",frame=FALSE)
box(bty="l")
par(op)
#adding RF based on the mean CI50 of the 15 most sensitive individual (
#approximately 10% of the tested strains)
sensit<-mean(as.numeric(byprod$ED50.abs[order(as.numeric(byprod$ED50.abs))])[1:20])
byprod$factres<-byprod$ED50.abs/sensit
byprod<-byprod[order(as.numeric(byprod$factres)),]
VenMFR<-rbind(VenMFR,byprod[,c(1:4,24:26)])
plot(VenMFR[VenMFR$ActiveSub=="trifloxystrobine","factres"])

VenMFR<-drop.levels(VenMFR,reorder=FALSE)


##############################################################################/
#Multi FR plot####
##############################################################################/

cooloor<- brewer.pal(12,"Set3")
options(scipen=10000)
op<-par(mfrow=c(1,3))
bornes<-c(min(AltFR$factres,na.rm=TRUE),max(AltFR$factres,na.rm=TRUE))
i<-1
plot(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],"factres"],log="y",
     pch=as.numeric(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],
                          "species"])+20,
     main="Distribution FR Alternaria",ylim=bornes,cex=1,las=1,
     ylab="FR",bg=cooloor[i])
abline(h=10,col=grey(0.3,0.8),lwd=5,lty=2)
abline(h=100,col=grey(0.3,0.8),lwd=5,lty=2)
i<-2
points(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],"factres"],
       pch=as.numeric(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-3
points(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],"factres"],
       pch=as.numeric(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-4
points(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],"factres"],
       pch=as.numeric(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-5
points(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],"factres"],
       pch=as.numeric(AltFR[AltFR$ActiveSub==levels(AltFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
legend(20,0.5,legend=levels(AltFR$ActiveSub),
       cex=1,pt.cex=1.3,
       y.intersp=1,x.intersp=1.2,
       pch=c(15),
       col=cooloor[c(1,2,3,4,5)],
       bty="n")
#export to .pdf 6 x 6

bornes<-c(min(VenGFR$factres,na.rm=TRUE),max(VenGFR$factres,na.rm=TRUE))
i<-1
plot(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],log="y",
     pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                          "species"])+20,
     main="Distribution FR Venturia Germination",ylim=bornes,cex=1,las=1,
     ylab="FR",bg=cooloor[i])
abline(h=10,col=grey(0.3,0.8),lwd=5,lty=2)
abline(h=100,col=grey(0.3,0.8),lwd=5,lty=2)
i<-2
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-3
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-4
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-5
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                            "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-6
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                             "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-7
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                             "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-8
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                             "species"])+20,
       cex=1,las=1,bg=cooloor[i])
i<-9
points(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenGFR[VenGFR$ActiveSub==levels(VenGFR$ActiveSub)[i],
                             "species"])+20,
       cex=1,las=1,bg=cooloor[i])
legend(30,0.8,legend=levels(VenGFR$ActiveSub),
       cex=1,pt.cex=1.3,
       y.intersp=1,x.intersp=1.2,
       pch=c(15),
       col=cooloor[c(1,2,3,4,5,6,7,8,9)],
       bty="n")
legend(10,0.2,legend=levels(VenGFR$species),cex=1,
       pt.cex=1.3,y.intersp=1,x.intersp=1.2,
       pch=c(21,22,23),bty="n")
#export to .pdf 6 x 6

bornes<-c(min(VenMFR$factres,na.rm=TRUE),max(VenMFR$factres,na.rm=TRUE))
i<-1
plot(VenMFR[VenMFR$ActiveSub==levels(VenMFR$ActiveSub)[i],"factres"],log="y",
     pch=as.numeric(VenMFR[VenMFR$ActiveSub==levels(VenMFR$ActiveSub)[i],
                           "species"])+20,
     main="Distribution FR Venturia Germination",ylim=bornes,cex=1,las=1,
     ylab="FR",bg=cooloor[i])
abline(h=10,col=grey(0.3,0.8),lwd=5,lty=2)
abline(h=100,col=grey(0.3,0.8),lwd=5,lty=2)
i<-2
points(VenMFR[VenMFR$ActiveSub==levels(VenMFR$ActiveSub)[i],"factres"],
       pch=as.numeric(VenMFR[VenMFR$ActiveSub==levels(VenMFR$ActiveSub)[i],
                             "species"])+20,
       cex=1,las=1,bg=cooloor[i])
legend(40,0.8,legend=levels(VenMFR$ActiveSub),
       cex=1,pt.cex=1.3,
       y.intersp=1,x.intersp=1.2,
       pch=c(15),
       col=cooloor[c(1,2)],
       bty="n")
legend(20,0.8,legend=levels(VenMFR$species),cex=1,
       pt.cex=1.3,y.intersp=1,x.intersp=1.2,
       pch=c(21,22,23),bty="n")
#export to .pdf 6 x 6


#same plot but combine on one figure and with log(EC50)
plot(log(temp[order(c(temp$CYPROCONAZOLE)),"CYPROCONAZOLE"]),
     main="log EC50 distribution",bg=cooloor[1],pch=21,cex=1,las=1,
     ylab="log(EC50)",ylim=c(-6,6))
points(log(temp[order(c(temp$DIFENOCONAZOLE)),"DIFENOCONAZOLE"]),
       bg=cooloor[2],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$EPOXICONAZOLE)),"EPOXICONAZOLE"]),
       bg=cooloor[3],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$FLUTRIAFOL)),"FLUTRIAFOL"]),
       bg=cooloor[5],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$MEFENTRIFLUCONAZOLE)),"MEFENTRIFLUCONAZOLE"]),
       bg=cooloor[6],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$METCONAZOLE)),"METCONAZOLE"]),
       bg=cooloor[7],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$PROCHLORAZE)),"PROCHLORAZE"]),
       bg=cooloor[8],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$`PROTHIOCONAZOLE-DESTHIO`)),
                "PROTHIOCONAZOLE-DESTHIO"]),
       bg=cooloor[9],pch=21,cex=1,las=1)
points(log(temp[order(c(temp$TEBUCONAZOLE)),"TEBUCONAZOLE"]),
       bg=cooloor[10],pch=21,cex=1,las=1,)
points(log(temp[order(c(temp$TETRACONAZOLE)),"TETRACONAZOLE"]),
       bg=cooloor[11],pch=21,cex=1,las=1)
legend(60,-1,legend=c("CYPROCONAZOLE","DIFENOCONAZOLE","EPOXICONAZOLE",
                      "FLUTRIAFOL","MEFENTRIFLUCONAZOLE","METCONAZOLE",
                      "PROCHLORAZE","PROTHIOCONAZOLE-DESTHIO","TEBUCONAZOLE",
                      "TETRACONAZOLE"),
       cex=1,pt.cex=1.3,
       y.intersp=0.7,x.intersp=1.2,
       pch=c(15),
       col=cooloor[c(1,2,3,5,6,7,8,9,10,11)],
       bty="n")
#export to .pdf 8 x 7 inches


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 2), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#correlation for Alternaria analyses
temp<-alterCI50[,c(1,2,4)]
#removing SA that didn't work
temp<-temp[temp$ActiveSub!="captane" & 
             temp$ActiveSub!="carbendazim" &
             temp$ActiveSub!="dithianon" & 
             temp$ActiveSub!="trifloxystrobine_sham100",]
temp<-spread(temp,ActiveSub,ED50.abs)
#all combination
pairs(log(temp[,c(2:6)]),las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)


#correlation for Venturia analyses
temp<-VenGerCI50[,c(1,2,4)]
temp<-spread(temp,ActiveSub,ED50.abs)
#all combination
pairs(log(temp[,c(2:11)]),las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

#comparison germination vs mycelial growth in Venturia bioassays
temp<-rbind(VenGerCI50,VenMycCI50)
#for difenoconazole
temp1<-temp[temp$ActiveSub=="difenoconazole",]
temp1<-temp1[,c(1,3,4)]
temp1<-spread(temp1,method,ED50.abs)
#all combination
pairs(log(temp1[,c(2:3)]),las=1,main="Correlation Methods difenoconazole",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#for trifloxystrobine
temp1<-temp[temp$ActiveSub=="trifloxystrobine",]
temp1<-temp1[,c(1,3,4)]
temp1<-spread(temp1,method,ED50.abs)
#all combination
pairs(log(temp1[,c(2:3)]),las=1,main="Correlation Methods trifloxystrobine",
      lower.panel=panel.smooth, upper.panel=panel.cor)


##############################################################################/
#Analyzing the multisensitivity profile of the strains####
##############################################################################/

#Clusterization based on scoring of 10 SA
row.names(temp)<-temp$sample_ID

#PCA for the scoring on 10 SA
truc<-dudi.pca(temp[,-c(1)],
               scannf=FALSE,nf=3)
scatter(truc)
#determining the optimal number of clusters
fviz_nbclust(temp[,c(2:13)],kmeans,method="gap_stat")
clust<-kmeans(temp[,c(2:13)],5)
fviz_cluster(clust,data=temp[,c(2:13)])
plot(truc$li[,c(1,2)],col=brewer.pal(5,"Dark2")[clust$cluster],
     pch=19,cex=2)

#we remove FENTINE HYDROXYDE and TOLNAFTATE because it is of no 
#interest here as well as individuals 39 to 44 that have too many 
#missing values
hclu<-hclust(dist(scale(temp[-c(39:44),c(2:4,6:12)]),
                  method="euclidean"),
             method="ward.D2")
plot(hclu)
fviz_dend(hclu,k=5,cex=0.5,rect=TRUE,
          k_colors=brewer.pal(5,"Dark2"))
#export to pdf 9 x 6 inches


##############################################################################/
#plot of the distribution of IC50 for each active substance####
##############################################################################/

#preparing the data set
temp<-CompRez[,c(1,2,4)]
temp<-spread(temp,Subs_Act,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(2,5))
plot(temp[order(c(temp$CYPROCONAZOLE)),"CYPROCONAZOLE"],
     main="CYPROCONAZOLE IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$DIFENOCONAZOLE)),"DIFENOCONAZOLE"],
     main="DIFENOCONAZOLE IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$EPOXICONAZOLE)),"EPOXICONAZOLE"],
     main="EPOXICONAZOLE IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$FLUTRIAFOL)),"FLUTRIAFOL"],
     main="FLUTRIAFOL IC50",bg=cooloor[5],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$MEFENTRIFLUCONAZOLE)),"MEFENTRIFLUCONAZOLE"],
     main="MEFENTRIFLUCONAZOLE IC50",bg=cooloor[6],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$METCONAZOLE)),"METCONAZOLE"],
     main="METCONAZOLE IC50",bg=cooloor[7],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$PROCHLORAZE)),"PROCHLORAZE"],
     main="PROCHLORAZE IC50",bg=cooloor[8],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$`PROTHIOCONAZOLE-DESTHIO`)),"PROTHIOCONAZOLE-DESTHIO"],
     main="PROTHIOCONAZOLE-DESTHIO IC50",bg=cooloor[9],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$TEBUCONAZOLE)),"TEBUCONAZOLE"],
     main="TEBUCONAZOLE IC50",bg=cooloor[10],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$TETRACONAZOLE)),"TETRACONAZOLE"],
     main="TETRACONAZOLE IC50",bg=cooloor[11],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
par(op)
#export to pdf 14 x 10 inches




##############################################################################/
#END
##############################################################################/