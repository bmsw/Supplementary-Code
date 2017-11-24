# Link to original code:
# http://www.bioconductor.org/packages/devel/bioc/vignettes/flowCyBar/inst/doc/flowCyBar-manual.pdf

# Necessary packages
library("flowCyBar")
library("seriation")
library("vegan")
# Cybar_counts is a 31x7 matrix which first column contains samples ID's and columns 2 to 7
# comprise counts for six regions gated for bacterial communities, named B1, B2, B3, B4, B5 and B6.

# flowCyBar::normalize
cybar=normalize(Cybar_counts[,-1],digits=2)

# Correcting problems for Vegan (class of cybar is "AsIs")
cybar$B1=as.numeric(cybar$B1)
cybar$B2=as.numeric(cybar$B2)
cybar$B3=as.numeric(cybar$B3)
cybar$B4=as.numeric(cybar$B4)
cybar$B5=as.numeric(cybar$B5)
cybar$B6=as.numeric(cybar$B6)

# Prepare to rename
rownames(cybar)<-Cybar_counts$X1

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)

# ID's CyBar

id.cybar_raw=rownames(cybar)
id.cybar_raw_regexpr=regexpr("[0-9]{2} dil", rownames(cybar), ignore.case = T)
id.cybar_clean=regmatches(id.cybar_raw,id.cybar_raw_regexpr, invert = F)

# Checking order
sum(match(id.cybar_clean, id.data_clean)==c(1:31))

# New data
DATA3=cbind("ID"=id.cybar_clean, DATA)

mat=vegdist(cybar)
#Matrix mat is converted to dissimilarity matrix ('mat.dist') to be used with metaMDS
mat.cybar<-as.matrix(mat)
# Renamed matrix
colnames(mat.cybar)<-rownames(mat.cybar)<-DATA3$Lago


# Heatmap

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



# Plot
myheatmap(mat.cybar[match(dgge$X1, rownames(mat.cybar)),match(dgge$X1, rownames(mat.cybar))])

