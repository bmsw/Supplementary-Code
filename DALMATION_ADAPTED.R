# Link to original code:
# http://www.ufz.de/index.php?en=38440


library("seriation")
# File: Results_Greyscale.txt
dat<-read.delim("Results_Greyscale.txt", row.names=2)
new.names<-gsub(".bmp","",row.names(dat))
row.names(dat)<-new.names

dat.ol<-dat[grep("overlap",row.names(dat)),]
dat.si<-dat[-c(grep("overlap",row.names(dat))),]

out.temp<-data.frame(t(combn(row.names(dat.si), 2)), t(combn(dat.si[,2], 2)))

out<-data.frame(out.temp, dat.ol[,2])
colnames(out)<-c("Scatter.1","Scatter.2","Pixel.Scatter.1","Pixel.Scatter.2","Pixel.overlap")
#rownames(out)<-rownames(dat.ol)

#Jaccard index

jacc.ind<- out[,5]/((out[,3]+ out[,4])-out[,5]) #

#Jaccard distance

jacc.dist<-  jacc.ind #


out.fin<-data.frame(out,Jaccard=jacc.dist) #
size<-nrow(dat.si)
mat<-matrix(nrow=size, ncol=size)
j<-size-1
temp<-0
for (i in 1:j)  {
  mat[,i]<-c(rep(NA,size-length(out.fin[(temp+1):(temp+j),6])), out.fin[(temp+1):(temp+j),6])
  temp<-temp+j
  j<-j-1
}
diag(mat)<-0
colnames(mat)<-rownames(mat)<-rownames(dat.si)

mat.dalmation<-as.dist(mat)

library("gdata")
mat2=as.matrix(mat.dalmation)
mat.dp=mat2
mat.dp[upper.tri(mat.dp, diag = F)]=t(mat.dp)[upper.tri(mat.dp, diag = F)]

# Prepare to rename

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)

# ID's CHIC

id.dp_raw=colnames(mat.dp)
id.dp_raw_regexpr=regexpr("[0-9]{2} dil", colnames(mat.dp), ignore.case = T)
id.dp_clean=regmatches(id.dp_raw,id.dp_raw_regexpr, invert = F)

# Checking order
sum(match(id.dp_clean, id.data_clean)==c(1:31))

# New data
DATA3=cbind("ID"=id.dp_clean, DATA)

# Renamed matrix
colnames(mat.dp)<-rownames(mat.dp)<-DATA3$Lago


# Heatmaps

myheatmap<-function(x){
  myfun <- function(x) hclust(as.dist(x), method = "ward.D")

  heatmap.2(x, hclustfun = myfun,  col=colorRampPalette(colors = c("white","azure","cyan2","blue4"))(100),
            dendrogram = "column", key=T, trace = "none",symm=T,
            lmat=rbind(c(0,3),c(2,1),c(0,4)), lwid=c(0.05,5), lhei = c(1,8,1.1), density.info = "none",
            key.title = NA,labCol=NA, key.xlab = NA, key.par=list(mar=c(2,0,0,8.3), cex=1),
            margins=c(.5,12.6),
            cexRow = 1.7,
            offsetRow = -.2)
}


myheatmap(mat.dp[match(dgge$X1, rownames(mat.dp)),match(dgge$X1, rownames(mat.dp))])

