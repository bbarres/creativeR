##############################################################################/
##############################################################################/
#Code for different figures with IC50
##############################################################################/
##############################################################################/

##loading the dataset and the necessary library
source("recif_load.R")

#in order for this code to work, you need first to run the script for 
#IC50 estimation
alterCI50<-read.table("output/Alternaria_mycelial.txt",header=TRUE,
                      sep="\t",stringsAsFactors=TRUE)


##############################################################################/
#barplot to compare the ED50 of the different samples####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
alterCI50$ED50.abs<-as.character(alterCI50$ED50.abs)
alterCI50[alterCI50$ED50.abs==">30","ED50.abs"]<-32
alterCI50$ED50.abs<-as.numeric(as.character(alterCI50$ED50.abs))

pdf(file="output/histo_AllInd_ASA.pdf",width=60,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(alterCI50$ED50.abs)),
        ylim=c(0,32),
        col=cooloor[as.numeric(as.factor(alterCI50$ActiveSub))],
        names.arg=alterCI50$strain_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=32,lty=2)
legend(300,47,levels(CompRez$Subs_Act),fill=cooloor,bty="n")
par(op)
dev.off()

#histogram by samples
samplelist<-as.character(names(table(CompRez$sample_ID)))
pdf(file="output/histo_byInd_ASA.pdf",width=9,height=20)
op<-par(mfrow=c(9,5))
for (i in (1:length(samplelist))) {
  temp<-merge(as.data.frame(levels(CompRez$Subs_Act)),
              CompRez[CompRez$sample_ID==samplelist[i],],
              by.x=1,by.y=1,
              all.x=TRUE)
  barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
          ylim=c(0,52))
}
par(op)
dev.off()
#export to pdf 12 x 16


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

temp<-CompRez[,c(1,2,4)]
temp<-spread(temp,Subs_Act,ED50)

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

#all 13 active substance
pairs(log(temp[,c(2:14)]),las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 13 x 11 inches
#all 10 DMI active subsatnce
pairs(log(temp[,c(2:4,7:13)]),las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 13 x 11 inches
#only 4 widely investigated DMI
pairs(log(temp[,c(3,8,11,13)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 8 x 6 inches

#just for the difenoconazole and tetraconazole
pairs(log(temp[,c(3,12)]),las=1,main="Correlation between log(ActSubst)",
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
#END
##############################################################################/