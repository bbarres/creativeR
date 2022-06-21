##############################################################################/
##############################################################################/
#Analysis for Venturia sp by mycelial growth bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#temporary data set loading waiting for the true final data set
VentuMyc.dat<-creadat[creadat$species!="Alternaria sp." & 
                        creadat$test_type=="mycelial_growth",]

#removing missing data
VentuMyc.dat<-VentuMyc.dat[!is.na(VentuMyc.dat$perc_croiss),]
VentuMyc.dat<-drop.levels(VentuMyc.dat,reorder=FALSE)
#VentuMyc.dat<-VentuMyc.dat[!is.na(VentuMyc.dat$perc_croiss),]

collist<-c("forestgreen","black")


##############################################################################/
#Analysis for Venturia sp by mycelial growth bioassay
##############################################################################/

REZ_VentuMy<-data.frame("strain_ID"=as.character(),"ActiveSub"=as.character(),
                        "method"=as.character(),"ED50-abs"=as.numeric(),
                        "ED50-SE"=as.numeric(),"ED50-lower"=as.numeric(),
                        "ED50-upper"=as.numeric(),"b-param"=as.numeric(),
                        "b-SE"=as.numeric(),"b-tval"=as.numeric(),
                        "b-pval"=as.numeric(),"c-param"=as.numeric(),
                        "c-SE"=as.numeric(),"c-tval"=as.numeric(),
                        "c-pval"=as.numeric(),"d-param"=as.numeric(),
                        "d-SE"=as.numeric(),"d-tval"=as.numeric(),
                        "d-pval"=as.numeric(),"e-param"=as.numeric(),
                        "e-SE"=as.numeric(),"e-tval"=as.numeric(),
                        "e-pval"=as.numeric())

for (j in 1: length(levels(VentuMyc.dat$active_substance))) {
  
  tempAS<-levels(VentuMyc.dat$active_substance)[j]
  subdat<-VentuMyc.dat[VentuMyc.dat$active_substance==tempAS,]
  #some individual never reach an inhibition of 50%, event for the highest 
  #tested concentration. 
  subdat_rez<-as.character(subdat[subdat$dose==max(subdat$dose) & 
                                    subdat$perc_croiss>50,
                                  "strain_ID"])
  subdat_rez<-subdat_rez[!(is.na(subdat_rez))]
  subdat_rez<-names(table(subdat_rez))
  ifelse(length(subdat_rez)==0,
         REZsub<-data.frame("strain_ID"=as.character(),
                            "ActiveSub"=as.character(),
                            "method"=as.character(),"ED50-abs"=as.numeric(),
                            "ED50-SE"=as.numeric(),"ED50-lower"=as.numeric(),
                            "ED50-upper"=as.numeric(),"b-param"=as.numeric(),
                            "b-SE"=as.numeric(),"b-tval"=as.numeric(),
                            "b-pval"=as.numeric(),"c-param"=as.numeric(),
                            "c-SE"=as.numeric(),"c-tval"=as.numeric(),
                            "c-pval"=as.numeric(),"d-param"=as.numeric(),
                            "d-SE"=as.numeric(),"d-tval"=as.numeric(),
                            "d-pval"=as.numeric(),"e-param"=as.numeric(),
                            "e-SE"=as.numeric(),"e-tval"=as.numeric(),
                            "e-pval"=as.numeric()),
         REZsub<-data.frame("strain_ID"=subdat_rez,"ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
                            "ED50-abs"=paste(">",max(subdat$dose),sep=""),
                            "ED50-SE"=NA,"ED50-lower"=NA,
                            "ED50-upper"=NA,"b-param"=NA,
                            "b-SE"=NA,"b-tval"=NA,
                            "b-pval"=NA,"c-param"=NA,
                            "c-SE"=NA,"c-tval"=NA,
                            "c-pval"=NA,"d-param"=NA,
                            "d-SE"=NA,"d-tval"=NA,
                            "d-pval"=NA,"e-param"=NA,
                            "e-SE"=NA,"e-tval"=NA,
                            "e-pval"=NA))
  
  #we limit the dataset to the sample that reach somehow a IC of 50%
  subdat<-subdat[!(subdat$strain_ID %in% subdat_rez),]
  subdat<-drop.levels(subdat,reorder=FALSE)
  #the subsequent analyses are unnecessary if all the strains are >maxdose
  if(dim(subdat)[1]==0) {
    REZ_VentuMy<-rbind(REZ_VentuMy,REZsub)
  } else {
    pdf(file=paste("output/VenMy_",tempAS,".pdf",sep=""),
        width=12,height=30)
    op<-par(mfrow=c(10,4))
    for (i in 1: dim(table(subdat$strain_ID))[1]) {
      datatemp<-subdat[subdat$strain_ID==names(table(subdat$strain_ID))[i],]
      couleur<-collist[as.numeric(datatemp$strain_type)]
      typeline<-as.numeric(datatemp$species)[1]
      print(as.character(datatemp$strain_ID[1]))
      if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
        tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                          "ActiveSub"=tempAS,
                          "method"=levels(subdat$test_type),
                          "ED50-abs"=NA,
                          "ED50-SE"=NA,"ED50-lower"=NA,
                          "ED50-upper"=NA,"b-param"=NA,
                          "b-SE"=NA,"b-tval"=NA,
                          "b-pval"=NA,"c-param"=NA,
                          "c-SE"=NA,"c-tval"=NA,
                          "c-pval"=NA,"d-param"=NA,
                          "d-SE"=NA,"d-tval"=NA,
                          "d-pval"=NA,"e-param"=NA,
                          "e-SE"=NA,"e-tval"=NA,
                          "e-pval"=NA)
        REZsub<-rbind(REZsub,tempx)
        plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             main=names(table(subdat$strain_ID))[i])
      } else { tryCatch({
        temp.m1<-drm(perc_croiss~dose,
                     data=datatemp,
                     fct=LL.4())
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             lty=typeline,type="confidence",
             main=names(table(subdat$strain_ID))[i])
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col="red",
             pch=4,type="obs",add=TRUE)
        plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
             lty=typeline,pch=19,add=TRUE)
      },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        if (!exists("temp.m1")){
          tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                            "ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
                            "ED50-abs"="ERROR",
                            "ED50-SE"=NA,"ED50-lower"=NA,
                            "ED50-upper"=NA,"b-param"=NA,
                            "b-SE"=NA,"b-tval"=NA,
                            "b-pval"=NA,"c-param"=NA,
                            "c-SE"=NA,"c-tval"=NA,
                            "c-pval"=NA,"d-param"=NA,
                            "d-SE"=NA,"d-tval"=NA,
                            "d-pval"=NA,"e-param"=NA,
                            "e-SE"=NA,"e-tval"=NA,
                            "e-pval"=NA)
        } else {
          temp<-ED(temp.m1,50,type="absolute",interval="delta")
          tempb<-summary(temp.m1)$coefficients[1,]
          tempc<-summary(temp.m1)$coefficients[2,]
          tempd<-summary(temp.m1)$coefficients[3,]
          tempe<-summary(temp.m1)$coefficients[4,]
          tempx<-data.frame("strain_ID"=names(table(subdat$strain_ID))[i],
                            "ActiveSub"=tempAS,
                            "method"=levels(subdat$test_type),
                            "ED50-abs"=temp[1],
                            "ED50-SE"=temp[2],"ED50-lower"=temp[3],
                            "ED50-upper"=temp[4],"b-param"=tempb[1],
                            "b-SE"=tempb[2],"b-tval"=tempb[3],
                            "b-pval"=tempb[4],"c-param"=tempc[1],
                            "c-SE"=tempc[2],"c-tval"=tempc[3],
                            "c-pval"=tempc[4],"d-param"=tempd[1],
                            "d-SE"=tempd[2],"d-tval"=tempd[3],
                            "d-pval"=tempd[4],"e-param"=tempe[1],
                            "e-SE"=tempe[2],"e-tval"=tempe[3],
                            "e-pval"=tempe[4])
        }
        REZsub<-rbind(REZsub,tempx)
        rm(temp.m1)
      }
      
    }
    par(op)
    dev.off()
    
    REZ_VentuMy<-rbind(REZ_VentuMy,REZsub)
  }
  
}

write.table(REZ_VentuMy,file="output/Venturia_mycelial.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/















##############################################################################/
##############################################################################/
#Analysis for Venturia by mycelial growth bioassay
##############################################################################/
##############################################################################/

source("creative_load_data.R")

#subsetting the global dataset
ventuMyc.dat<-creadat[creadat$species=="V. inaequalis" | 
                        creadat$species=="V. asperata",]
ventuMyc.dat<-ventuMyc.dat[ventuMyc.dat$test_type=="mycelial_growth",]
ventuMyc.dat$species<-factor(ventuMyc.dat$species,
                             levels=rev(levels(ventuMyc.dat$species)))
ventuMyc.dat<-drop.levels(ventuMyc.dat,reorder=FALSE)
#ventuMyc.dat<-ventuMyc.dat[!is.na(ventuMyc.dat$perc_croiss),]

collist<-c("forestgreen","black")


####boscalild
subdat<-ventuMyc.dat[ventuMyc.dat$active_substance=="boscalid",]
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
pdf(file="output/VeMyboscalid.pdf",width=12,height=30)
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

REZ_VeMyBosc<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="boscalid")


####difenoconazole
subdat<-ventuMyc.dat[ventuMyc.dat$active_substance=="difenoconazole",]
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
pdf(file="output/VeMydifenoconazole.pdf",width=12,height=30)
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

REZ_VeMyDife<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="difenoconazole")


####dodine
subdat<-ventuMyc.dat[ventuMyc.dat$active_substance=="dodine",]
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
pdf(file="output/VeMydodine.pdf",width=12,height=30)
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

REZ_VeMyDodi<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="dodine")


####trifloxystrobine
subdat<-ventuMyc.dat[ventuMyc.dat$active_substance=="trifloxystrobine",]
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
pdf(file="output/VeMytrifloxystrobine.pdf",width=12,height=30)
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

REZ_VeMyTrif<-data.frame(REZsub[sort(REZsub$strain_ID),],
                         "TestType"="Mycelium",
                         "SubsAct"="trifloxystrobine")

#concatenating the different result files
REZ_VeMy<-rbind(REZ_VeMyBosc,REZ_VeMyDife,REZ_VeMyDodi,REZ_VeMyTrif)
write.table(REZ_VeMy,file="output/Venturia_mycelial.txt",quote=FALSE,
            col.names=TRUE,row.names=FALSE,sep="\t")


##############################################################################/
#END
##############################################################################/