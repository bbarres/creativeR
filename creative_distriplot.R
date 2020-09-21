##############################################################################/
##############################################################################/
#Some plot of the distribution of IC50
##############################################################################/
##############################################################################/

library(RColorBrewer)
library(tidyr)


temp<-ventuMyc.dat[ventuMyc.dat$dose==30 & ventuMyc.dat$active_substance=="boscalid",]
temp<-cbind(as.character(temp$strain_ID),
            as.character(temp$species))[order(as.character(temp$strain_ID)),]

REZ_VeMyBosc<-data.frame(REZ_VeMyBosc,temp)
REZ_VeMyBosc<-REZ_VeMyBosc[order(as.numeric(as.character(REZ_VeMyBosc$ED50))),]
plot(as.numeric(as.character(REZ_VeMyBosc$ED50)),col=REZ_VeMyBosc$X2,pch=19,cex=2, 
     main="Venturia boscalid mycelium")

REZ_VeMyDife<-data.frame(REZ_VeMyDife,temp)
REZ_VeMyDife<-REZ_VeMyDife[order(as.numeric(as.character(REZ_VeMyDife$ED50))),]
plot(as.numeric(as.character(REZ_VeMyDife$ED50)),col=REZ_VeMyDife$X2,pch=19,cex=2, 
     main="Venturia difenoconazole mycelium")

REZ_VeMyDodi<-data.frame(REZ_VeMyDodi,temp)
REZ_VeMyDodi<-REZ_VeMyDodi[order(as.numeric(as.character(REZ_VeMyDodi$ED50))),]
plot(as.numeric(as.character(REZ_VeMyDodi$ED50)),col=REZ_VeMyDodi$X2,pch=19,cex=2, 
     main="Venturia dodine mycelium")

REZ_VeMyTrif<-data.frame(REZ_VeMyTrif,temp)
REZ_VeMyTrif<-REZ_VeMyTrif[order(as.numeric(as.character(REZ_VeMyTrif$ED50))),]
plot(as.numeric(as.character(REZ_VeMyTrif$ED50)),col=REZ_VeMyTrif$X2,pch=19,cex=2, 
     main="Venturia trifloxystrobine mycelium")





temp<-ventuGer.dat[ventuGer.dat$dose==30 & ventuGer.dat$active_substance=="boscalid",]
temp<-cbind(as.character(temp$strain_ID),
            as.character(temp$species))[order(as.character(temp$strain_ID)),]
temp<-temp[-c(4,6),]

REZ_VeGeBosc<-data.frame(REZ_VeGeBosc,temp)
REZ_VeGeBosc<-REZ_VeGeBosc[order(as.numeric(as.character(REZ_VeGeBosc$ED50))),]
plot(as.numeric(as.character(REZ_VeGeBosc$ED50)),col=REZ_VeGeBosc$X2,pch=19,cex=2, 
     main="Venturia boscalid germination")

REZ_VeGeCapt<-data.frame(REZ_VeGeCapt,temp)
REZ_VeGeCapt<-REZ_VeGeCapt[order(as.numeric(as.character(REZ_VeGeCapt$ED50))),]
plot(as.numeric(as.character(REZ_VeGeCapt$ED50)),col=REZ_VeGeCapt$X2,pch=19,cex=2, 
     main="Venturia captane germination")

REZ_VeGeCypr<-data.frame(REZ_VeGeCypr,temp)
REZ_VeGeCypr<-REZ_VeGeCypr[order(as.numeric(as.character(REZ_VeGeCypr$ED50))),]
plot(as.numeric(as.character(REZ_VeGeCypr$ED50)),col=REZ_VeGeCypr$X2,pch=19,cex=2, 
     main="Venturia cyprodinil germination")

REZ_VeGeDith<-data.frame(REZ_VeGeDith,temp)
REZ_VeGeDith<-REZ_VeGeDith[order(as.numeric(as.character(REZ_VeGeDith$ED50))),]
plot(as.numeric(as.character(REZ_VeGeDith$ED50)),col=REZ_VeGeDith$X2,pch=19,cex=2, 
     main="Venturia dithianon germination")

REZ_VeGeDodi<-data.frame(REZ_VeGeDodi,temp)
REZ_VeGeDodi<-REZ_VeGeDodi[order(as.numeric(as.character(REZ_VeGeDodi$ED50))),]
plot(as.numeric(as.character(REZ_VeGeDodi$ED50)),col=REZ_VeGeDodi$X2,pch=19,cex=2, 
     main="Venturia dodine germination")

REZ_VeGeManc<-data.frame(REZ_VeGeManc,temp)
REZ_VeGeManc<-REZ_VeGeManc[order(as.numeric(as.character(REZ_VeGeManc$ED50))),]
plot(as.numeric(as.character(REZ_VeGeManc$ED50)),col=REZ_VeGeManc$X2,pch=19,cex=2, 
     main="Venturia mancozebe germination")

REZ_VeGeMane<-data.frame(REZ_VeGeMane,temp)
REZ_VeGeMane<-REZ_VeGeMane[order(as.numeric(as.character(REZ_VeGeMane$ED50))),]
plot(as.numeric(as.character(REZ_VeGeMane$ED50)),col=REZ_VeGeMane$X2,pch=19,cex=2, 
     main="Venturia manebe germination")

#Alternaria plot
plot(sort(as.numeric(as.character(REZ_AlMyBosc$ED50))),
     main="Alternaria boscalid mycelial growth",
     ylim=c(0,max(as.numeric(as.character(REZ_AlMyBosc$ED50)),na.rm=TRUE))
)

plot(sort(as.numeric(as.character(REZ_AlMyCypr$ED50))),
     main="Alternaria cyprodinil mycelial growth",
     ylim=c(0,max(as.numeric(as.character(REZ_AlMyCypr$ED50)),na.rm=TRUE))
)

plot(sort(as.numeric(as.character(REZ_AlMyDife$ED50))),
     main="Alternaria difenoconazole mycelial growth",
     ylim=c(0,max(as.numeric(as.character(REZ_AlMyDife$ED50)),na.rm=TRUE))
)

plot(sort(as.numeric(as.character(REZ_AlMyDodi$ED50))),
     main="Alternaria dodine mycelial growth",
     ylim=c(0,max(as.numeric(as.character(REZ_AlMyDodi$ED50)),na.rm=TRUE))
)

plot(sort(as.numeric(as.character(REZ_AlMyFluop$ED50))),
     main="Alternaria fluopyram mycelial growth",
     ylim=c(0,max(as.numeric(as.character(REZ_AlMyFluop$ED50)),na.rm=TRUE))
)


##############################################################################/
#barplot to compare the ED50 of the different samples####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
REZ_VeMy$ED50<-as.character(REZ_VeMy$ED50)
REZ_VeMy[REZ_VeMy$ED50==">30","ED50"]<-32
REZ_VeMy$ED50<-as.numeric(as.character(REZ_VeMy$ED50))

pdf(file="output/histo_AllInd_ASA_VeMy.pdf",width=50,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(REZ_VeMy$ED50)),
        ylim=c(0,40),col=cooloor[as.numeric(REZ_VeMy$SubsAct)],
        names.arg=REZ_VeMy$strain_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=32,lty=2)
legend(300,47,levels(REZ_VeMy$SubsAct),fill=cooloor,bty="n")
par(op)
dev.off()

#histogramme by samples
samplelist<-as.character(names(table(REZ_VeMy$strain_ID)))
pdf(file="output/histo_byInd_ASA_VeMy.pdf",width=20,height=14)
op<-par(mfrow=c(5,8))
for (i in (1:length(samplelist))) {
        temp<-merge(as.data.frame(levels(REZ_VeMy$SubsAct)),
                    REZ_VeMy[REZ_VeMy$strain_ID==samplelist[i],],
                    by.x=1,by.y=4,
                    all.x=TRUE)
        barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
                ylim=c(0,32))
}
par(op)
dev.off()
#export to pdf 12 x 16


#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
REZ_VeGe$ED50<-as.character(REZ_VeGe$ED50)
REZ_VeGe[REZ_VeGe$ED50==">30","ED50"]<-32
REZ_VeGe$ED50<-as.numeric(as.character(REZ_VeGe$ED50))

pdf(file="output/histo_AllInd_ASA_VeGe.pdf",width=50,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(REZ_VeGe$ED50)),
        ylim=c(0,7),col=cooloor[as.numeric(REZ_VeGe$SubsAct)],
        names.arg=REZ_VeGe$strain_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=6,lty=2)
legend(300,47,levels(REZ_VeGe$SubsAct),fill=cooloor,bty="n")
par(op)
dev.off()

#histogramme by samples
samplelist<-as.character(names(table(REZ_VeGe$strain_ID)))
pdf(file="output/histo_byInd_ASA_VeGe.pdf",width=20,height=14)
op<-par(mfrow=c(5,8))
for (i in (1:length(samplelist))) {
        temp<-merge(as.data.frame(levels(REZ_VeGe$SubsAct)),
                    REZ_VeGe[REZ_VeGe$strain_ID==samplelist[i],],
                    by.x=1,by.y=4,
                    all.x=TRUE)
        barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
                ylim=c(0,6))
}
par(op)
dev.off()
#export to pdf 12 x 16


#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
REZ_AlMy$ED50<-as.character(REZ_AlMy$ED50)
REZ_AlMy[REZ_AlMy$ED50==">30","ED50"]<-32
REZ_AlMy$ED50<-as.numeric(as.character(REZ_AlMy$ED50))

pdf(file="output/histo_AllInd_ASA_AlMy.pdf",width=50,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(REZ_AlMy$ED50)),
        ylim=c(0,35),col=cooloor[as.numeric(REZ_AlMy$SubsAct)],
        names.arg=REZ_AlMy$strain_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=32,lty=2)
legend(300,47,levels(REZ_AlMy$SubsAct),fill=cooloor,bty="n")
par(op)
dev.off()

#histogramme by samples
samplelist<-as.character(names(table(REZ_AlMy$strain_ID)))
pdf(file="output/histo_byInd_ASA_AlMy.pdf",width=20,height=14)
op<-par(mfrow=c(5,8))
for (i in (1:length(samplelist))) {
        temp<-merge(as.data.frame(levels(REZ_AlMy$SubsAct)),
                    REZ_AlMy[REZ_AlMy$strain_ID==samplelist[i],],
                    by.x=1,by.y=4,
                    all.x=TRUE)
        barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
                ylim=c(0,40))
}
par(op)
dev.off()
#export to pdf 12 x 16


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

#for Venturia mycelium
temp<-REZ_VeMy[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

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

pairs(temp[,c(2:4)],las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[,c(2:4)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 11 x 11 inches


#for Venturia germination
temp<-REZ_VeGe[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

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

pairs(temp[,c(2:10)],las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[,c(2:10)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 11 x 11 inches


#for Venturia germination
temp<-REZ_AlMy[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

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

pairs(temp[,c(2:10)],las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[,c(2:10)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 11 x 11 inches


##############################################################################/
#plot of the distribution of IC50 for each active substance####
##############################################################################/

#preparing the dataset
temp<-REZ_VeMy[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(2,2))
plot(temp[order(c(temp$boscalid)),"boscalid"],
     main="boscalid IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$difenoconazole)),"difenoconazole"],
     main="difenoconazole IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$dodine)),"dodine"],
     main="dodine IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$trifloxystrobine)),"trifloxystrobine"],
     main="trifloxystrobine IC50",bg=cooloor[4],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
par(op)
#export to pdf 8 x 6 inches


#preparing the dataset
temp<-REZ_VeGe[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(3,3))
plot(temp[order(c(temp$boscalid)),"boscalid"],
     main="boscalid IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$captane)),"captane"],
     main="captane IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$cyprodinil)),"cyprodinil"],
     main="cyprodinil IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$difenoconazole)),"difenoconazole"],
     main="difenoconazole IC50",bg=cooloor[4],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$dithianon)),"dithianon"],
     main="dithianon IC50",bg=cooloor[5],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$dodine)),"dodine"],
     main="dodine IC50",bg=cooloor[6],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$mancozebe)),"mancozebe"],
     main="mancozebe IC50",bg=cooloor[7],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$manebe)),"manebe"],
     main="manebe IC50",bg=cooloor[8],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
plot(temp[order(c(temp$thirame)),"thirame"],
     main="thirame IC50",bg=cooloor[9],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,10))
par(op)
#export to pdf 14 x 10 inches


#preparing the dataset
temp<-REZ_AlMy[,c(1,2,4)]
temp<-spread(temp,SubsAct,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(3,3))
plot(temp[order(c(temp$boscalid)),"boscalid"],
     main="boscalid IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$captane)),"captane"],
     main="captane IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$carbendazim)),"carbendazim"],
     main="carbendazim IC50",bg=cooloor[9],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$cyprodinil)),"cyprodinil"],
     main="cyprodinil IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$difenoconazole)),"difenoconazole"],
     main="difenoconazole IC50",bg=cooloor[4],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$dithianon)),"dithianon"],
     main="dithianon IC50",bg=cooloor[5],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$dodine)),"dodine"],
     main="dodine IC50",bg=cooloor[6],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$fluopyram)),"fluopyram"],
     main="fluopyram IC50",bg=cooloor[7],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
plot(temp[order(c(temp$trifloxystrobine)),"trifloxystrobine"],
     main="trifloxystrobine IC50",bg=cooloor[8],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,35))
par(op)
#export to pdf 14 x 10 inches


##############################################################################/
#END
##############################################################################/