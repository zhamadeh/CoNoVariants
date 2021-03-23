################################################
# Packages #
################################################
library(tidyverse)
library(GenomicRanges)
library(scales)

################################################
#PLOTTING WITH TRANSPARENT BACKGROUND
################################################

transparentBackground <- function(p,filename){
	if (!file.exists("Plots") ) { dir.create("Plots")}
	p <- p +
		theme(
			panel.background = element_rect(fill = "transparent"), # bg of the panel
			plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
			panel.grid.major = element_blank(), # get rid of major grid
			panel.grid.minor = element_blank(), # get rid of minor grid
			legend.background = element_rect(fill = "transparent"), # get rid of legend bg
			legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
		)

	ggsave(p, filename = paste0("Plots/",filename),  bg = "transparent")
}


################################################
# Plotting #
################################################

#READ IN SAVED DATASETS IF NOT ALREADY LOADED

cnvPerCell=read.table("DataSummary/cnvPerCell.txt",header=T)

cnvPerCellSummary=read.table("DataSummary/cnvPerCellSummary.txt",header=T)



# CLEAN UP A BIT
cnvPerCell$ploidy <- as.factor(cnvPerCell$ploidy)
cnvPerCell$seqnames <- factor(cnvPerCell$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" , "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
cnvPerCell$file=as.factor(cnvPerCell$file)


eventDistributionPerPloidy =ggplot(cnvPerCell)+ geom_bar(mapping=aes(type,fill=type))+
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
	theme_classic()+
	theme(
		legend.position = c(0.9,0.8),
		legend.title = element_blank(),
		text=element_text(size=20,face="bold"),
		axis.text.x  = element_blank(),
		axis.ticks.x = element_blank(),
		axis.title =element_text(size=28))+
	labs(x="TYPE OF CNV",y="EVENTS > 20Mb")


cnvPerChromsome  =ggplot(cnvPerCell)+geom_bar(aes(seqnames,fill=type))+
	theme_classic()+
	theme(legend.position = c(0.84,0.65),
		  legend.title = element_blank(),
		  text=element_text(size=15,face="bold"),
		  axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
	labs(x="CHROMSOME",y="EVENTS > 20Mb")

sizeOfEvents = ggplot()+geom_density(data = cnvPerCell,aes(x = width,fill=type),alpha=0.5)+
	theme_classic()+labs(x="SIZE",y="DENSITY")+
	scale_x_continuous(trans="log10",breaks=trans_breaks('log10',function(x) 10^x), labels=trans_format('log10',math_format(10^.x)))+
	theme(legend.title = element_blank(),
		  axis.text.y = element_blank(),
		  axis.ticks.y  = element_blank(),
		  text=element_text(size=22,face="bold"))




r  <- as.data.frame(cnvPerCell %>%  group_by(name)  %>% dplyr::summarize(n=n()))

r$name<-as.factor(r$name)
r$gene="LIBRARIES"
numOfEventsPerCellByPloidy =ggplot(r)+geom_col(aes(gene,n,fill=name))+
	scale_fill_manual(values =  c("0"="black","1"="purple",  "2"="#00EE76" , "3"="#CD0000",
								  "4"="#EEC900" , "5" ="#000080","6"="#FFFACD")) + theme_classic()+
	labs(y="# OF CNVs (>30Mb) PER CELL")+
	theme(
		legend.position = c(0.75,0.85),
		legend.title = element_blank(),
		text=element_text(size=20,face="bold"),
		axis.title.x = element_blank())



################################################
# Exporting #
################################################

transparentBackground(eventDistributionPerPloidy,"eventDistributionPerPloidy.png")
transparentBackground(cnvPerChromsome,"cnvPerChromsome.png")
transparentBackground(sizeOfEvents,"sizeOfEvents.png")
transparentBackground(numOfEventsPerCellByPloidy,"numOfEventsPerCellByPloidy.png")
