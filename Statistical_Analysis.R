## Required Packages
library("betapart")
library("reshape")
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
library("ggbiplot")
library("gvlma")
library("RVAideMemoire")
library("gplots")
## Sourcing flowDiv
source("flowDiv_internals.R")
source("flowDiv_main.R")
## Sourcing data
# DATA_IMPORT.R contains: dgge, env DATA and Cybar_counts
# dgge = DGGE band's information
# env = enviromental variables
# DATA = ID's and dilutions for samples
# Cybar_counts = gates counts used in FlowCybar pipeline
source("DATA_IMPORT.R")

# Removing Huapi and P.10 (no data avaiable)
env<-env[-4,]
dgge<-dgge[-c(4,10),]
# Correcting ID
env$`Water body`[24]="SL. Negra"
env$`Water body`[31]="L. Yehuín"
env$`Water body`[21]="SL. de los Cisnes"
#DATA
DATA$Lago[13]="P.15"
DATA$Lago[21]="L. Acigami"
DATA$Lago[28]="SL. Verde"
# dgge
dgge$X1[22]="SL. Victoria"
dgge$X1[27]="L. Yehuín"
dgge$X1[26]="SL. de los Cisnes"
# Correcting blank spaces
env$`Water body`=gsub("\\.","\\. ",gsub("\\. ", "\\.", env$`Water body`),env$`Water body`)
DATA$Lago=gsub("\\.","\\. ",gsub("\\. ", "\\.", DATA$Lago),DATA$Lago)
dgge$X1=gsub("\\.","\\. ",gsub("\\. ", "\\.", dgge$X1),dgge$X1)
# Visual checking 1
cbind(DATA$Lago, env$`Water body`)
# Matching sites
match(DATA$Lago, env$`Water body`)
# New ENV
env2=env[match(DATA$Lago, env$`Water body`),]
# Visual checking 2
cbind(DATA$Lago, env2$`Water body`)

# flowDiv
pata.fd=flowDiv2("../Romina_Bruno_1_10.wsp", "Bact tot", do.plot = T, static = F, dilutions=DATA$`x (flowDiv)`, transform = T, autotrans=T,psize = 1.08)


# Getting counts

wksp=opc("../Romina_Bruno_1_10.wsp")
samples=nn(wksp,nod = "Bact tot", use.beads = F, nod2="NULL")
counts=unlist(lapply(wksp[[1]], function(x)getTotal(x, "Bact tot")))
counts.corrected=counts*DATA$`x (flowDiv)`



# ID's flowDiv
id.flowdiv_raw=names(pata.fd$Alpha)
id.flowdiv_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.flowdiv_raw, ignore.case = T)
id.flowdiv_clean=regmatches(id.flowdiv_raw,id.flowdiv_raw_regexpr, invert = F)

# ID's DATA

id.data_raw=DATA$file
id.data_raw_regexpr=regexpr("[0-9]{2} dil [0-9]{1}_[0-9]{1}", id.data_raw, ignore.case = T)
id.data_clean=regmatches(id.data_raw,id.data_raw_regexpr, invert = F)
DATA2=cbind("ID"=id.data_clean, DATA)

# Visual checking 3
cbind(DATA$Lago, env2$`Water body`, DATA2$Lago)
# Changing names
rownames(pata.fd$Matrices)<-as.character(DATA2$Lago)

# Matriz renomeada
m1=pata.fd$Matrices

#Definindo os indices
rich=apply(m1, 1, function(x) specnumber(x, MARGIN = 1))
pielou=apply(m1, 1, function(x)diversity(x)/log(vegan::specnumber(x)))
shannon=apply(m1, 1, function(x) diversity(x,index = "shannon"))
simpson=apply(m1, 1, function(x) diversity(x,index = "simpson"))
invsimpson=apply(m1, 1, function(x) diversity(x,index = "invsimpson"))

# New data set

DATA3=cbind(DATA2, "FC_Richness"=rich, "FC_Pielou"=pielou, "FC_Shannon"=shannon, "FC_Simpson"=simpson, "FC_invSimpson"=invsimpson, "FC_counts"=counts.corrected)

# Mantel FD x DGGE
mantel(vegdist(m1[match(dgge$X1, rownames(m1)),]), vegdist(dgge[,-1]))

# FlowFP, CyBar, CHIC and Dalmation Plot results from files:
# CHIC_ADAPTED.R
# CYBAR_ADAPTED.R
# DALMATION_ADAPTED.R
# FP_ADAPTED.R
out<-c(21,22,23,25)
mantel(as.dist(mat.chic2[-out,-out]), vegdist(dgge[,-1]))
mantel(as.dist(mat.dp[-out,-out]), vegdist(dgge[,-1]))
mantel(as.dist(mat.cybar[-out,-out]), vegdist(dgge[,-1]))
mantel(as.dist(mat.fp[-out,-out]), vegdist(dgge[,-1]))

# FlowFP and flowDiv
mantel(mat.fp[-out,-out], vegdist(m1[match(dgge$X1, rownames(m1)),]))

# DP and flowDiv
mantel(mat.cybar[-out,-out], vegdist(m1[match(dgge$X1, rownames(m1)),]))


##
# Statistics (Dunn test)
dunn.test(DATA3$FC_Shannon, DATA3$`TROPHIC STATE`, method = "Bonferroni")
dunn.test(DATA3$FC_Pielou, DATA3$`TROPHIC STATE`, method = "Bonferroni")
dunn.test(DATA3$FC_Richness, DATA3$`TROPHIC STATE`, method = "Bonferroni")

# Statistics  (ANOVA)

# GVLMA
gvlma(lm(DATA3$FC_Shannon~DATA3$`TROPHIC STATE`)) #OK
gvlma(lm(DATA3$FC_Pielou~DATA3$`TROPHIC STATE`)) #OK
gvlma(lm(DATA3$FC_Richness~DATA3$`TROPHIC STATE`)) # NOT OK AT 0.05

# Log
gvlma(lm(log10(DATA3$FC_Richness)~DATA3$`TROPHIC STATE`)) # NOT OK

# ANOVA
aov.shannon=aov(lm(DATA3$FC_Shannon~DATA3$`TROPHIC STATE`)) # SIGNIFICANTE
aov.pielou=aov(lm(DATA3$FC_Pielou~DATA3$`TROPHIC STATE`)) # NOT SIGNIFICANTE

# Richness log for ANOVA
aov.rich=aov(lm(log10(DATA3$FC_Richness)~DATA3$`TROPHIC STATE`)) # SIGNIFICANTE
summary(aov.shannon)
summary(aov.pielou)
summary(aov.rich)
# Tukey
plot(TukeyHSD(aov.shannon))
plot(TukeyHSD(aov.rich))

# Adjusting labesl
label1=gsub("MESO", "MESOTROPHIC", DATA3$`TROPHIC STATE`)
label1=gsub("EUT", "EUTROPHIC", label1)
label1=gsub("OLIGO", "OLIGOTROPHIC", label1)
# Plot
mylabels=factor(label1)


# Boxplots

DATA3df=data.frame(DATA3)

DATA3df$Log_Richness=log10(DATA3$FC_Richness)

mybox<-function(x, lab){
  g0=ggplot(aes_string(y=x, x="TROPHIC.STATE"), data=DATA3df)+geom_boxplot(aes_string(fill=mylabels, alpha=.5))+ theme(legend.position="none")+ scale_x_discrete(name="")+geom_point()+
    theme(text = element_text(size=25),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
          axis.text.x=element_blank())+
    labs(y = lab)
  return(g0)
}


mybox("FC_Shannon", "Shannon Index")
mybox("FC_Pielou")
mybox("FC_Richness")
mybox("Log_Richness")

# Correlations

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, method = "spearman", exact = F)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value

    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

M4=cbind(FC_Shannon=DATA3$FC_Shannon,FC_Pielou=DATA3$FC_Pielou, FC_Richness=DATA3$FC_Richness,FC_Counts=DATA3$FC_counts, env2[,-c(1:3)])


corr.test(M4, method = "spearman")
res1 <- cor.mtest(M4, 0.95)
#corrplot(cor(M4, method = "spearman", use="complete.obs"), p.mat = res1[[1]], sig.level=.05, type = "lower", pch.col="red", method = "number", cl.cex=.8, number.cex = .7)

# Adjusting colnames
M4.1=M4
colnames(M4.1)<-c("Shannon","Pielou","Richness","Counts","Latitude","Longitude","Altitude","Area","Temperature","pH","Conductivity","DO","DIN","Kd","Chl a","Phosphate","DOC")
corrplot(cor(M4.1, method = "spearman", use="complete.obs"), p.mat = res1[[1]], sig.level=.05, type = "lower", pch.col="gray20", method = "square", cl.cex=.8, number.cex = .8,
         col=c("red3", "blue3"), tl.col = "black", pch.cex = 2.5)


# Scatterplots

M5=M4
colnames(M5)
colnames(M5)[10]="pH"
colnames(M5)[9]="Temperature"
colnames(M5)[15]="Chla"
M5$Chla = log10(M5$Chla)

mycor<-function(x,y){
    g0=ggplot(aes_string(x, "DATA3$FC_Shannon"), data=M5)+geom_text(aes(label=`Water body`))+geom_smooth(method='lm')+
      xlab(y)+ylab("Shannon Index")
  return(g0)
}

# Richness versus Shannon
ggplot(aes_string("FC_Richness", "FC_Shannon"), data=DATA3df)+geom_text(aes(label=Lago))+geom_smooth(method='lm')+
  xlab("FC_Richness")+ylab("Shannon Index")
# Other plots
mycor("pH", "pH")
mycor("Temperature", "Temp.(°C)")
mycor("Chla", "Log10 Chl a (μg L−1)")


# PCA

M5=M4
colnames(M5)[1:4]<-c("Shannon",  "Pielou",   "Richness", "Counts")
wine.pca <- prcomp(M5[,c(1:4)], scale. = TRUE)

ggplot(aes(PC1, PC2, colour=mylabels), data=pts)+
  geom_point(size=3)+
  stat_ellipse(level=0.95, geom = "polygon", alpha = .2, aes(fill = mylabels))+
  geom_segment(data=vects, mapping=aes(x=0, y=0, xend=PC1*3, yend=PC2*3), arrow=arrow(length = unit(0.2, "cm")), size=.5, color="black", linetype=1)+
  geom_text(data=vects, aes(x=PC1*3.5, y=PC2*3.5, label=rownames(vects)), size=5, vjust=0,hjust=c(0,0,0,0), colour="black")+
  theme_bw()+
  theme(panel.border = element_blank(),
        text = element_text(size=15),
        legend.title=element_blank(),
        legend.background = element_rect(colour="white", size=.5, linetype=3),
        legend.position = "bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  labs(x = "PC1 (57.59%)", y = "PC2(33.30%)")


meta2=metaMDS(m1, autotransform = F)

plot1=function(x){

  plot(x, type="n")
  text(x$points, labels=DATA3$Lago, col=as.numeric(DATA3$`TROPHIC STATE`), pch=17)

}

#plot1(meta1)
plot1(meta2)

# Heatmaps and Mantel

# New datasets
dgge1=data.frame(dgge)
rownames(dgge1) <- c(dgge1$X1)
# Distance matrices
dgge_dist=vegdist(dgge1[,-1])
m1.dist=vegdist(m1[match(dgge$X1, rownames(m1)),])
# Heatplots
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

myheatmap(as.matrix(dgge_dist))
myheatmap(as.matrix(m1.dist))


fit=envfit(meta2, M5, na.rm = T)

# Linear variables
ordisurf(meta2~M5$pH) # 10
ordisurf(meta2~M5$`Chl a (μg L−1)`) # 15
ordisurf(meta2~M5$`DOC (ppm)`) # 17


# Non-linear variables (not used)
ordisurf(meta2~M5$`Lat. (°S)`) # 5
ordisurf(meta2~M5$`Long. (°W)`) # 6
ordisurf(meta2~ M5$`Altitude (m a.s.l.)`) # 7
ordisurf(meta2~M5$`Area (km2)`) # 8
ordisurf(meta2~ M5$`Temp.(°C)`) # 9
ordisurf(meta2~M5$`Conduct. (μS cm−1)`) # no
ordisurf(meta2~ M5$`DO (mg L−1)`) # 12
ordisurf(meta2~M5$`DIN (mg L−1)`) # no
ordisurf(meta2~M5$`Kd (m−1)`) # no

fit2=data.frame(fit$vectors$arrows)
rownames(fit2)<-c("Shannon","Pielou","Richness","Counts","Lat.","Lon.","Altitude","Area","Temp.","pH","Conduc.","DO","DIN","Kd","Chla","Phos.","DOC")
# Visual adjust
fit3=fit2[c(5:10,15,17),]*c(rep(1.1, 6), 1.2, 1.2)
fit3[8,1]=fit3[8,1]+.6

ggplot(aes(NMDS1, NMDS2, colour=mylabels), data=data.frame(scores(meta2)))+
  geom_point(size=2)+
  stat_ellipse(level=0.95, geom = "polygon", alpha = .2, aes(fill = mylabels))+
  geom_segment(data=fit2[c(5:10,15,17),], mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2), arrow=arrow(length = unit(0.2, "cm")), size=.5, color="black", linetype=1)+
  geom_text(data=fit3, aes(x=NMDS1, y=NMDS2, label=rownames(fit2[c(5:10,15,17),])), size=5, vjust=0,hjust=1, colour="black")+
  theme_bw()+
  theme(panel.border = element_blank(),
        text = element_text(size=15),
        legend.title=element_blank(),
        legend.background = element_rect(colour="white", size=.5, linetype=3),
        legend.position = "bottom",
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


# Adonis
adonis(vegdist(m1)~DATA3$`TROPHIC STATE`)

# Beta disper # Not significant!

anova(betadisper(vegdist(m1), DATA3$`TROPHIC STATE`))

# Pairwise
pairwise.perm.manova(vegdist(m1),DATA3$`TROPHIC STATE`, p.method = "bonferroni")

# New variables
m2=vegdist(m1)
m3=as.matrix(m2)

meso= DATA3$`TROPHIC STATE`=="MESO"
eut= DATA3$`TROPHIC STATE`=="EUT"
oligo= DATA3$`TROPHIC STATE`=="OLIGO"


m=m3[meso,meso]
o=m3[oligo,oligo]
e=m3[eut,eut]

o0=o[upper.tri(o)]
e0=e[upper.tri(e)]
m0=m[upper.tri(m)]

all=stack(list("o"=c(o), "m"=c(m),"e"=c(e)))
all0=stack(list("o"=c(o0), "m"=c(m0),"e"=c(e0)))


boxplot(values~ind, all, ylim=c(0,1))
boxplot(values~ind, all0, ylim=c(0.3,1))

# Bray-Curtis Distance
ggplot(aes(y=values, x=ind), data=all0)+geom_boxplot(aes(fill=ind, alpha=.5))+
  theme(legend.position="none",
        text = element_text(size=25),
        axis.text=element_text(size=25),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  labs(y ="Bray-Curtis Distance")

# Variance

a=betadisper(vegdist(m1), DATA3$`TROPHIC STATE`)
df=data.frame(a$distances, a$group)

ggplot(aes(y=a.distances, x=a.group), data=df)+geom_boxplot(aes(fill=a.group, alpha=.5))+
  theme(legend.position="none",
        text = element_text(size=25),
        axis.text=element_text(size=25),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())+
  labs(y ="Distance to Centroid")


# Turnover x Nestdness
oligo=apply(m1[oligo,], 2, sum)/nrow(m1[oligo,])
eut=apply(m1[eut,], 2, sum)/nrow(m1[eut,])
meso=apply(m1[meso,], 2, sum)/nrow(m1[meso,])
all=rbind("oligo"=oligo, "eut"=eut, "meso"=meso)
mall=vegdist(all)
braym=bray.part(all)
nest=braym$bray.bal/braym$bray
turn=1-nest

df1=data.frame(cbind("contrast"=c("ExO","MxO", "ExM"), "TURNOVER"=c(turn), "NESTDNESS"=c(nest)))
df2=melt(df1, id.var="contrast")
df2$value=as.numeric(as.character(df2$value))

ggplot(df2, aes(x = contrast, y = value, fill = variable))+
  geom_bar(stat = "identity", alpha=0.9)+
  scale_y_continuous(labels = scales::percent)+
  theme(text = element_text(size=23),
        axis.line = element_blank(),
        axis.text=element_text(size=20),
        axis.title.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = "top")+
  labs(y = "Bray-Curtis Distance (%)")





