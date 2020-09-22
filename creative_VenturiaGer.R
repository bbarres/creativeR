##############################################################################/
##############################################################################/
#Analysis for Venturia by mycelial growth bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#subsetting the global dataset
ventuGer.dat<-creadat[creadat$species=="V. inaequalis" | 
                        creadat$species=="V. asperata",]
ventuGer.dat<-ventuGer.dat[ventuGer.dat$test_type=="germination",]
ventuGer.dat$species<-factor(ventuGer.dat$species,
                             levels=rev(levels(ventuGer.dat$species)))
ventuGer.dat<-drop.levels(ventuGer.dat,reorder=FALSE)
#ventuGer.dat<-ventuGer.dat[!is.na(ventuGerc.dat$perc_croiss),]

collist<-c("forestgreen","black")


####boscalild
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="boscalid",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGeboscalid.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeBosc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="boscalid")



####captane
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="captane",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGecaptane.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeCapt<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="captane")



####cyprodinil
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="cyprodinil",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGecyprodinil.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeCypr<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="cyprodinil")



####difenoconazole
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="difenoconazole",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedifenoconazole.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDife<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="difenoconazole")



####dithianon
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="dithianon",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedithianon.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDith<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="dithianon")



####dodine
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="dodine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGedodine.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeDodi<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="dodine")



####mancozebe
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="mancozebe",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGemancozebe.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeManc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="mancozebe")



####manebe
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="manebe",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGemanebe.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeMane<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="manebe")



####thirame
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="thirame",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGethirame.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeThir<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="thirame")



####trifloxystrobine
subdat<-ventuGer.dat[ventuGer.dat$active_substance=="trifloxystrobine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/VeGetrifloxystrobine.pdf",width=12,height=30)
op<-par(mfrow=c(10,4))
for (i in 1: dim(table(subdat$strain_ID))[1]) {
  datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
  couleur<-collist[as.numeric(datatemp$strain_type)]
  typeline<-as.numeric(datatemp$species)[1]
  print(as.character(datatemp$strain_ID[1]))
  if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
    tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                      "ED50"=c("NA"))
    REZsub<-rbind(REZsub,tempx)
    plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         main=names(table(subdat$strain_ID))[i])
  } else { tryCatch({
    temp.m1<-drm(perc_croiss~dose,
                 data=datatemp,
                 fct=LL.4())
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,main=names(table(subdat$strain_ID))[i])
    plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
         lty=typeline,type="confidence",add=TRUE)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    if (!exists("temp.m1")){
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"="ERROR")
    } else {
      temp<-ED(temp.m1,50,type="absolute")
      tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                        "ED50"=as.character(temp[1]))
    }
    REZsub<-rbind(REZsub,tempx)
    rm(temp.m1)
  }
}
par(op)
dev.off()

REZ_VeGeTrif<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Germination",
                         "SubsAct"="trifloxystrobine")


#concatenating the different result files
REZ_VeGe<-rbind(REZ_VeGeBosc,REZ_VeGeCapt,REZ_VeGeCypr,REZ_VeGeDife,
                REZ_VeGeDith,REZ_VeGeDodi,REZ_VeGeManc,REZ_VeGeMane,
                REZ_VeGeThir,REZ_VeGeTrif)
write.table(REZ_VeGe,file="output/Venturia_germination.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/