#Set directory
setwd("D:/OneDrive - University of Cambridge/IMS_MRL_2014-2018/Collaborations/Sandovici_etal")
#Target Scan files
ts<-list.files(pattern = ".txt")
fl<-list.files(pattern = ".csv")
fl.name<-c("15w_Overexpression","e19_KD")

#conversion function:
library("biomaRt")

#from Human to Mouse
cvrt <- function(x){
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=F)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
#select by score
score <- -0.3
pdf(paste("Sandovici_OverExpression.pdf",sep=""), 
    width = 8.27,
    height = 11.69, 
    paper='special', 
    onefile=T)
par(mfrow=c(5,4))
####overexpression####
rs<-read.csv(file=fl[1],
             header=T)
tgts<-c()
for(i in 7){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="miR-483-5p",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(234,121,14,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(234,121,14,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(234,121,14,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(234,121,14,maxColorValue = 255 ))
tgts<-c()
for(i in 5:6){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="miR-483-3p",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(234,121,14,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(234,121,14,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(234,121,14,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(234,121,14,maxColorValue = 255 ))
tgts<-c()
for(i in 5:length(ts)){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="miR-483",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(234,121,14,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(234,121,14,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(234,121,14,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(234,121,14,maxColorValue = 255 ))
####repression####
rs<-read.csv(file=fl[2],
             header=T)
tgts<-c()
for(i in 5:6){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="miR-483-3p",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(133,125,204,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(133,125,204,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(133,125,204,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(133,125,204,maxColorValue = 255 ))
tgts<-c()
for(i in 7){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="miR-483-5p",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(133,125,204,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(133,125,204,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(133,125,204,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(133,125,204,maxColorValue = 255 ))
tgts<-c()
for(i in 3:4){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="Positively regulated",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(133,125,204,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(133,125,204,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(133,125,204,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(133,125,204,maxColorValue = 255 ))
tgts<-c()
for(i in c(1:2,5:7)){
  x<-read.table(ts[i],
                header=T,
                sep="\t")
  x<-x[which(x$Total.context...score < score),]
  tgts<-c(tgts,x$Ortholog.of.target.gene)
}
target.gene<-cvrt(tgts)
nontg<-rs[which(!rs$Gene.name %in% target.gene),]
tg<-rs[which(rs$Gene.name %in% target.gene & rs$logCPM > 5),]
ks<-ks.test(tg$logFC,nontg$logFC,alternative="two.sided")
ks.p<-round(ks$p.value,digits=3)
#plot
plot(ecdf(nontg$logFC), 
     ylab = "Cumulative Fraction", xlab = "log2(Fold Change)",
     main="Negatively regulated",
     col.01line = rgb(0,0,0,0.5),
     mgp=c(2.5,1,0),
     xlim=c(-2,2),
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     pch = 20, col = rgb(1,1,1), 
     cex = 0.7, cex.lab = 1.2, cex.main = 1,
     font.main = 1,
     yaxt="n"
)
axis(2,seq(0,1,by=0.2),las=1)
abline(v=0,col=rgb(0,0,0,0.5),lty=3)
plot(ecdf(nontg$logFC),
     col.01line = rgb(0,0,0,0), 
     verticals = T, col.ver=rgb(0,0,0,0.5),lty=1,
     add=T, pch = 20, col = rgb(0,0,0), cex=0.5) 
plot(ecdf(tg$logFC), 
     verticals = T, col.ver=rgb(133,125,204,maxColorValue = 255 ),lty=1,
     col.01line = rgb(0,0,0,0), 
     add=T, pch = 20, col = rgb(133,125,204,maxColorValue = 255 ), cex=0.3) 
legend("topleft",c("No target", "miRNA target"), bg=rgb(1,1,1),
       col = c("black", rgb(133,125,204,maxColorValue = 255 )), 
       pch=20,
       cex = 0.7, adj =0)
text(1,0.3,
     labels=paste("P",ks.p,sep="="),
     col=rgb(133,125,204,maxColorValue = 255 ))
dev.off()