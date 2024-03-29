##############################################################################/
##############################################################################/
#Loading the packages and dataset needed for analyses
##############################################################################/
##############################################################################/

#loading the libraries
library(drc)
library(plotrix)
library(gdata)
library(RColorBrewer)

#load the global dataset
#creadat<-read.table(file="data/creative_data.txt",header=T,sep="\t")
creadat<-read.table(file="data/creative_data_final.txt",header=T,
                    sep="\t",stringsAsFactors=TRUE)
creadat$species<-factor(creadat$species,levels=rev(levels(creadat$species)))
collist<-c("forestgreen","black")


##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/














#below the code used for the very first round of analyses. To be removed soon

# ###############################################################################
# #Analysis for the boscalid
# ###############################################################################
# 
# #subsetting the global dataset
# bosc.dat<-creadat[creadat$active_substance=="boscalid",]
# bosc.dat<-bosc.dat[!is.na(bosc.dat$perc_croiss),]
# 
# #some individual never reach an inhibition of 50%, event for the highest 
# #tested concentration. 
# bosc_rez<-as.character(bosc.dat[bosc.dat$dose=="30" & bosc.dat$perc_croiss>50,
#                                 "strain_ID"])
# REZbos<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
# #we limit the dataset to the sample that reach somehow a IC of 50%
# #bosc.dat<-bosc.dat[!(bosc.dat$sample_ID %in% bosc_rez),]
# bosc.dat<-drop.levels(bosc.dat,reorder=FALSE)
# pdf(file="output/boscalid.pdf",width=12,height=30)
# op<-par(mfrow=c(10,4))
# for (i in 1: dim(table(bosc.dat$strain_ID))[1]) {
#   datatemp<-bosc.dat[bosc.dat$strain_ID==names(table(bosc.dat$strain_ID))[i],]
#   couleur<-collist[as.numeric(datatemp$strain_type)]
#   typeline<-as.numeric(datatemp$species)[1]
#   if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
#     tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
#                       "ED50"=c("NA"))
#     REZbos<-rbind(REZbos,tempx)
#     plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          main=names(table(bosc.dat$strain_ID))[i])
#   } else {
#     temp.m1<-drm(perc_croiss~dose,
#                  data=datatemp,
#                  fct=LL.4())
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,main=names(table(bosc.dat$strain_ID))[i])
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,type="confidence",add=TRUE)
#     temp<-ED(temp.m1,50,type="absolute")
#     tempx<-data.frame("strain_ID"=names(table(bosc.dat$strain_ID))[i],
#                       "ED50"=temp[1])
#     REZbos<-rbind(REZbos,tempx)
#   }
# }
# par(op)
# dev.off()
# 
# 
# ###############################################################################
# #Analysis for the dodine
# ###############################################################################
# 
# #subsetting the global dataset
# dodi.dat<-creadat[creadat$active_substance=="dodine",]
# dodi.dat<-dodi.dat[!is.na(dodi.dat$perc_croiss),]
# 
# #some individual never reach an inhibition of 50%, event for the highest 
# #tested concentration. 
# dodi_rez<-as.character(dodi.dat[dodi.dat$dose=="30" & dodi.dat$perc_croiss>50,
#                                 "strain_ID"])
# REZdod<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
# #we limit the dataset to the sample that reach somehow a IC of 50%
# #dodi.dat<-dodi.dat[!(dodi.dat$sample_ID %in% dodi_rez),]
# dodi.dat<-drop.levels(dodi.dat,reorder=FALSE)
# pdf("output/dodine.pdf",width=12,height=30)
# op<-par(mfrow=c(10,4))
# for (i in 1: dim(table(dodi.dat$strain_ID))[1]) {
#   datatemp<-dodi.dat[dodi.dat$strain_ID==names(table(dodi.dat$strain_ID))[i],]
#   couleur<-collist[as.numeric(datatemp$strain_type)]
#   typeline<-as.numeric(datatemp$species)[1]
#   if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
#     tempx<-data.frame("strain_ID"=names(table(dodi.dat$strain_ID))[i],
#                       "ED50"=c("NA"))
#     REZdod<-rbind(REZdod,tempx)
#     plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          main=names(table(dodi.dat$strain_ID))[i])
#   } else {
#     temp.m1<-drm(perc_croiss~dose,
#                  data=datatemp,
#                  fct=LL.4())
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,main=names(table(dodi.dat$strain_ID))[i])
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,type="confidence",add=TRUE)
#     temp<-ED(temp.m1,50,type="absolute")
#     tempx<-data.frame("strain_ID"=names(table(dodi.dat$strain_ID))[i],
#                       "ED50"=temp[1])
#     REZdod<-rbind(REZdod,tempx)
#   }
# }
# par(op)
# dev.off()
# 
# 
# ###############################################################################
# #Analysis for the difenoconazole
# ###############################################################################
# 
# #subsetting the global dataset
# dife.dat<-creadat[creadat$active_substance=="difenoconazole",]
# dife.dat<-dife.dat[!is.na(dife.dat$perc_croiss),]
# 
# #some individual never reach an inhibition of 50%, event for the highest 
# #tested concentration. 
# dife_rez<-as.character(dife.dat[dife.dat$dose=="30" & dife.dat$perc_croiss>50,
#                                 "strain_ID"])
# REZdif<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
# #we limit the dataset to the sample that reach somehow a IC of 50%
# #dife.dat<-dife.dat[!(dife.dat$sample_ID %in% dife_rez),]
# dife.dat<-drop.levels(dife.dat,reorder=FALSE)
# pdf("output/difenoconazole.pdf",width=12,height=30)
# op<-par(mfrow=c(10,4))
# for (i in 1: dim(table(dife.dat$strain_ID))[1]) {
#   datatemp<-dife.dat[dife.dat$strain_ID==names(table(dife.dat$strain_ID))[i],]
#   couleur<-collist[as.numeric(datatemp$strain_type)]
#   typeline<-as.numeric(datatemp$species)[1]
#   if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
#     tempx<-data.frame("strain_ID"=names(table(dife.dat$strain_ID))[i],
#                       "ED50"=c("NA"))
#     REZdif<-rbind(REZdif,tempx)
#     plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          main=names(table(dife.dat$strain_ID))[i])
#   } else {
#     temp.m1<-drm(perc_croiss~dose,
#                  data=datatemp,
#                  fct=LL.4())
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,main=names(table(dife.dat$strain_ID))[i])
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,type="confidence",add=TRUE)
#     temp<-ED(temp.m1,50,type="absolute")
#     tempx<-data.frame("strain_ID"=names(table(dife.dat$strain_ID))[i],
#                       "ED50"=temp[1])
#     REZdif<-rbind(REZdif,tempx)
#   }
# }
# par(op)
# dev.off()
# 
# 
# ###############################################################################
# #Analysis for the trifloxystrobine
# ###############################################################################
# 
# #subsetting the global dataset
# trif.dat<-creadat[creadat$active_substance=="trifloxystrobine",]
# trif.dat<-trif.dat[!is.na(trif.dat$perc_croiss),]
# 
# #some individual never reach an inhibition of 50%, event for the highest 
# #tested concentration. 
# trif_rez<-as.character(trif.dat[trif.dat$dose=="30" & trif.dat$perc_croiss>50,
#                                 "strain_ID"])
# REZtri<-data.frame("strain_ID"=as.character(),"ED50"=as.character())
# #we limit the dataset to the sample that reach somehow a IC of 50%
# trif.dat<-trif.dat[!(trif.dat$strain_ID %in% trif_rez),]
# trif.dat<-drop.levels(trif.dat,reorder=FALSE)
# pdf("output/trifloxystrobine.pdf",width=12,height=30)
# op<-par(mfrow=c(10,4))
# for (i in 1: dim(table(trif.dat$strain_ID))[1]) {
#   datatemp<-trif.dat[trif.dat$strain_ID==names(table(trif.dat$strain_ID))[i],]
#   couleur<-collist[as.numeric(datatemp$strain_type)]
#   typeline<-as.numeric(datatemp$species)[1]
#   if (is.na(datatemp[1,"perc_croiss"])==TRUE) {
#     tempx<-data.frame("strain_ID"=names(table(trif.dat$strain_ID))[i],
#                       "ED50"=c("NA"))
#     REZtri<-rbind(REZtri,tempx)
#     plot(0,1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          main=names(table(trif.dat$strain_ID))[i])
#   } else {
#     temp.m1<-drm(perc_croiss~dose,
#                  data=datatemp,
#                  fct=LL.4())
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,main=names(table(trif.dat$strain_ID))[i])
#     plot(temp.m1,ylim=c(-10,120),xlim=c(0,30),col=couleur,
#          lty=typeline,type="confidence",add=TRUE)
#     temp<-ED(temp.m1,50,type="absolute")
#     tempx<-data.frame("strain_ID"=names(table(trif.dat$strain_ID))[i],
#                       "ED50"=temp[1])
#     REZtri<-rbind(REZtri,tempx)
#   }
# }
# par(op)
# dev.off()
# 
# 
# ###############################################################################
# #combined results and export
# ###############################################################################
# 
# REZ<-rbind(REZbos,REZdif,REZdod,REZtri)
# REZ<-cbind(REZ,"active_substance"=c(rep("boscalid",dim(REZbos)[1]),
#                                     rep("difenoconazole",dim(REZdif)[1]),
#                                     rep("dodine",dim(REZdod)[1]),
#                                     rep("trifloxystrobine",dim(REZtri)[1])))
# 
# write.table(REZ,file="output/result_CI50.txt",quote=FALSE,col.names=TRUE, 
#             row.names=FALSE,sep="\t")
# 
# 
# ###############################################################################
# #END
# ###############################################################################
# 
# 
# REZbos$ED50[REZbos$ED50>30]<-30
# plot(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="Souches ID",
#      ylab="FR",las=1)
# abline(0.39/0.39,0,col="green4",lwd=2)
# abline(3.9/0.39,0,col="red",lwd=2)
# #export to pdf 10 x 6 inches
# write.table(REZbos,file="output/REZbos.txt",quote=FALSE,sep="\t",
#             row.names=FALSE)
# 
# hist(REZbos$ED50[order(REZbos$ED50)]/0.39,main="Boscalid",xlab="FR Classes",
#      breaks=c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),
#      las=1,col=heat.colors(8)[8:1],ylim=c(0,450))
# abline(v=10,col="red",lwd=3)
# #export to pdf 4.5 x 9 inches