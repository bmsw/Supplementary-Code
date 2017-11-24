# Link to original code:
# http://www.ufz.de/index.php?en=38441

library("tcltk")
library("seriation")
setwd("/home/brunomsw/Documentos/DOUTORADO/KOCH_PIPELINES/CHIC/Plots")
overlaps<-read.delim("Results_overlaps.txt")     #Reads overlapt results estimated by ImageJ

xor<-read.delim("Results_xor.txt")      #Reads XOR results estimated by ImageJ

overlaps$Label<-gsub("overlap_","",overlaps$Label)    #Removes "overlap_" from file labels

overlaps$Label<-gsub(".bmp","",overlaps$Label)     #Removes ".bmp" from file labels

label.split<-strsplit(as.character(overlaps$Label),"_&_")  #Creates a file with proper label names

nam<-NULL #Creates an empty vector file

for (i in 1:length(label.split)) nam<-rbind(nam, cbind(label.split[[i]][1], label.split[[i]][2])) #creates a file ("nam") with all combinations of label names in it

dats<-data.frame(Label=overlaps$Label, Label.1=nam[,1], Label.2=nam[,2], IntDen=xor$IntDen, Area=overlaps$Area, IntDen.Area=xor$IntDen/overlaps$Area/100)

size<-length(unique(dats$Label.1))+1 #Number of unique data labels is estimated

mat<-matrix(nrow=size, ncol=size)      #Empty matrix ('mat') for disimilariy matrix is created

j=1          #Counter J created and set to 1

for (i in unique(dats$Label.1)) {      #dissimilarity according to formula 1 in MS for each combination is calculated and stored in matrix 'mat'

  mat[,j]<-c(rep(NA,size-length(dats$IntDen.Area[dats$Label.1==i])),dats$IntDen.Area[dats$Label.1==i])

  j<-j+1

  }
diag(mat)<-0           #Diagonals are filled up with '0' values

colnames(mat)<-rownames(mat)<-union(unique(dats$Label.1), unique(dats$Label.2)) #labels names are written into matrix

mat.chic<-as.dist(mat)         #Matrix mat is converted to dissimilarity matrix ('mat.dist') to be used with metaMDS

#Raw matrix
mat.chic2<-as.matrix(mat.chic)
# Prepare to rename

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)

# ID's CHIC

id.chic_raw=colnames(mat.chic2)
id.chic_raw_regexpr=regexpr("[0-9]{2} dil", colnames(mat.chic2), ignore.case = T)
id.chic_clean=regmatches(id.chic_raw,id.chic_raw_regexpr, invert = F)

# Checking order
sum(match(id.chic_clean, id.data_clean)==c(1:31))

# New data
DATA3=cbind("ID"=id.chic_clean, DATA)

# Renamed matrix
colnames(mat.chic2)<-rownames(mat.chic2)<-DATA3$Lago

# Heatmaps
hmap(mat.chic2,method = "HC_ward", col=colorRampPalette(colors = c("white", "blue"))(100), margins=c(8,8), main="CHIC")

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




myheatmap(mat.chic2)

myheatmap(mat.chic2[match(dgge$X1, rownames(mat.chic2)),match(dgge$X1, rownames(mat.chic2))])

