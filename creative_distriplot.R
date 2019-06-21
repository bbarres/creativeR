###############################################################################
###############################################################################
#Some plot of the distribution of IC50
###############################################################################
###############################################################################

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




