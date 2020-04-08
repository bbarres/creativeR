##############################################################################/
##############################################################################/
#Analysis for Alternaria sp by mycelial growth bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#subsetting the global dataset
AlterMyc.dat<-creadat[creadat$species=="Alternaria sp.",]
AlterMyc.dat<-AlterMyc.dat[AlterMyc.dat$test_type=="mycelial_growth",]
AlterMyc.dat$species<-factor(AlterMyc.dat$species,
                             levels=rev(levels(AlterMyc.dat$species)))
#removing missing data
AlterMyc.dat<-AlterMyc.dat[!is.na(AlterMyc.dat$perc_croiss),]
AlterMyc.dat<-drop.levels(AlterMyc.dat,reorder=FALSE)
#AlterMyc.dat<-AlterMyc.dat[!is.na(AlterMyc.dat$perc_croiss),]

collist<-c("forestgreen","black")


####boscalild####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="boscalid",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMyboscalid.pdf",width=12,height=30)
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

REZ_AlMyBosc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="boscalid")


####captane####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="captane",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMycaptane.pdf",width=12,height=30)
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

REZ_AlMyCapt<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="captane")



####carbendazim####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="carbendazim",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMycarbendazim.pdf",width=12,height=30)
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

REZ_AlMyCarb<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="carbendazim")



####cyprodinil####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="cyprodinil",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMycyprodinil.pdf",width=12,height=30)
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

REZ_AlMyCypr<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="cyprodinil")



####difenoconazole####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="difenoconazole",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMydifenoconazole.pdf",width=12,height=30)
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

REZ_AlMyDife<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="difenoconazole")



####dithianon####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="dithianon",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMydithianon.pdf",width=12,height=30)
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

REZ_AlMyDith<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="dithianon")


####dodine####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="dodine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMydodine.pdf",width=12,height=30)
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

REZ_AlMyDodi<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="dodine")



####fluopyram####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="fluopyram",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMyfluopyram.pdf",width=12,height=30)
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

REZ_AlMyFluop<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="fluopyram")



####trifloxystrobine####
subdat<-AlterMyc.dat[AlterMyc.dat$active_substance=="trifloxystrobine",]
#some individual never reach an inhibition of 50%, event for the highest 
#tested concentration. 
subdat_rez<-as.character(subdat[subdat$dose=="30" & subdat$perc_croiss>50,
                                "strain_ID"])
subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
subdat_rez<-names(table(subdat_rez))
ifelse(length(subdat_rez)==0,
       REZsub<-data.frame("strain_ID"=as.character(),"ED50"=as.character()),
       REZsub<-data.frame("strain_ID"=subdat_rez,
                          "ED50"=paste(">",max(subdat$dose),sep="")))
#we limit the dataset to the sample that reach somehow a IC of 50%
subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]

subdat<-drop.levels(subdat,reorder=FALSE)
pdf(file="output/AlMytrifloxystrobine.pdf",width=12,height=30)
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

REZ_AlMyTrif<-data.frame(REZsub[sort(REZsub$strain_ID),],
                          "TestType"="Mycelium",
                          "SubsAct"="trifloxystrobine")


#concatenating the different result files
REZ_AlMy<-rbind(REZ_AlMyBosc,REZ_AlMyCapt,REZ_AlMyCarb,REZ_AlMyCypr,
                REZ_AlMyDife,REZ_AlMyDith,REZ_AlMyDodi,REZ_AlMyFluop,
                REZ_AlMyTrif)
write.table(REZ_AlMy,file="output/Alternaria_mycelial.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/