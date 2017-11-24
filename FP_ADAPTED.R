# Link to original code:
#http://www.ufz.de/index.php?en=38443

# Packages needed
library("vegan")
library("flowWorkspace")
library("ggplot2")
library("ggcyto")
library("plyr")
library("car")
library("corrplot")
library("pangaear")
library("ggbiplot")
library("psych")
library("cowplot")
library("RColorBrewer")
library("dunn.test")
library("seriation")
# Sourcing code (gates are based on FlowJo)
source("flowDiv_internals.R")
source("flowDiv_main.R")
# FlowFP
library("flowFP")
# The directory where FCS files are
# Parsing data
## MyWorkspace.wsp is a FlowJo workspace with 31 samples
## "Bact tot" os the gate of insterest
wksp=opc("MyWorkspace.wsp")
samples=nn(wksp,nod = "Bact tot", use.beads = F, nod2="NULL")
mydata=samples$nodesample[[1]]


#Defining a model on sample 1 using the parameters FL1-H and SSC-H and 6 recursions.
model <- flowFPModel(mydata[[1]], parameters=c("FL1-H", "SSC-H"), nRecursions=6)

#The model is shown
plot(model)


# Defining a flowSet from FCS (for FlowFP use)
pataFS=flowSet(mydata)
#Now, the model is applied to all samples.
fp <- flowFP(pataFS, model)

# From here we can use vegan
mat.fp0=counts(fp)
# Bray
mat.fp=vegdist(mat.fp0)
# Distance matrix
mat.fp=as.matrix(mat.fp)

# ID's flowDiv
id.flowdiv_raw=samples$names
id.flowdiv_raw_regexpr=regexpr("[0-9]{2} dil", id.flowdiv_raw, ignore.case = T)
id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)

# ID's DATA
id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)

# Checking order
sum(match(id.flowdiv_clean, id.data_clean)==c(1:31))

# New data
DATA3=cbind("ID"=id.flowdiv_clean, DATA)

# Renamed matrix
colnames(mat.fp)<-rownames(mat.fp)<-DATA3$Lago

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





myheatmap(mat.fp)

myheatmap(mat.fp[match(dgge$X1, rownames(mat.fp)),match(dgge$X1, rownames(mat.fp))])
