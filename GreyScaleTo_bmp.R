# Creation of 236 greyscale images (.bmp) from gated data

# Code used
source("flowDiv_internals.R")
source("flowDiv_main.R")
# Packages used
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

# Parsing data
## MyWorkspace.wsp is a FlowJo workspace with 31 samples
## "Bact tot" os the gate of insterest
wksp=opc("MyWorkspace.wsp")
samples=nn(wksp,nod = "Bact tot", use.beads = F, nod2="NULL")

# Getting informations on channels FL1 and SSC
## FL1
fl1=unlist(lapply(samples$nodesample[[1]], function(x)exprs(x$"FL1-H")))
 max(fl1)
 min(fl1)
## SSC
ssc=unlist(lapply(samples$nodesample[[1]], function(x)exprs(x$"SSC-H")))
 max(ssc)
 min(ssc)

#Generate all colors needed
myColors <- gray.colors(236)
#Reverse it as you want dark for high and white for low
myColors <- rev(myColors)
#Wrap it into a ramp function
myColRamp <- colorRampPalette(myColors)

# plotbmp function
plotbmp<-function(x,y){

flowdf=data.frame((exprs(x)))

# Cytograms in greyscales (236)
gplot=ggplot(flowdf, aes(y=FL1.H, x=SSC.H))+
  geom_hex(bins = 128, na.rm = T)+
   scale_fill_gradientn("", colours = rev(gray.colors(236, end = 5/6)))+
  theme_bw()+
  scale_x_continuous(limits= c(0.34,1.25)) +
  scale_y_continuous(limits= c(0.34, 1.25)) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position="none")
    ggsave(y,gplot, device = "bmp")
}


# Saving as bmp files
mapply(function(z,w)plotbmp(z,w), samples$nodesample[[1]], samples$names)
