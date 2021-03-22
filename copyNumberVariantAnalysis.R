######################################################################
				####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(plyr)

######################################################################
			####   Assemble datasets   ####
######################################################################

## WRAPPER FUNCTION FOR COLLECTING GAINED AND LOST CHROMOSOME SEGMENTS > 20Mb


assembleTandemDuplicationDF <- function(CNV="Input/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz", metricsDir="Metrics/"){

	## IMPORT GRANGES LIST OBJECT (ANEUFINDER OUPUT)
	CNV<- import(CNV)

	## FUNCTION FOR ADDING METRICS
	collectMetrics <-  function(metricsDir){
		metrics <- data.frame()

		for (file in list.files("Metrics/",full.names = T)){
			met <- read.table(file,header = T,fill=T,row.names = NULL)
			if ("Postfiltering_reads_aligned" %in% colnames(met)){
				met <- met %>% select(Library,Postfiltering_reads_aligned,Mode_GC)
				met$file = as.factor(basename(file))
				colnames(met)=c("Library"="Library","Postfiltering_reads_aligned"="Reads", "Mode_GC" = "Mode_GC" ,"file"="file" )
			} else  {
				met <- met %>% select(Library,Reads_aligned_postfiltering,Mode_GC)
				met$file = as.factor(basename(file))
				colnames(met)=c("Library"="Library","Reads_aligned_postfiltering"="Reads", "Mode_GC" = "Mode_GC" ,"file"="file" )
			}
			met$Library=paste0(met$Library,".trimmed.mdup.bam")
			metrics <- rbind(met,metrics)
		}

		metrics$date=plyr::revalue(metrics$file, c("Aug18.txt" ="18-08-2020", "Aug28.txt"  ="28-08-2020" ,"Feb25Metrics.txt" ="25-02-2019",
												   "metrics_summary-jan18.txt"="18-01-2021", "oct23Metrics.txt"="23-10-2021","March6Metrics.txt"="06-03-2019"))

		message("\nThere are ",length(metrics$Library)," libraries in the base metrics file with ",ncol(metrics)," features\n")
		metrics$date <- as.Date(metrics$date,format="%d-%m-%Y")
		return(metrics)
	}

	## COLLECT METRICS
	metrics = collectMetrics(metricsDir)

	## RUN CORE FUNCTION
	countGainsAndLosses(CNV)
}




# CORE FUNCTION FOR COLLECTING GAINED AND LOST CHROMOSOME SEGMENTS > 20Mb


countGainsAndLosses <- function(CNV){


	cnvPerCell<-data.frame()
	cnvPerCellSummary=data.frame()


	for (i in 1:length(CNV)){


		## ASSIGN GENOTYPE BASED OFF FILE NAME
		ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
		file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
		date=as.Date(metrics[metrics$Library==file,]$date)


		message("Reading file: ",file," ... ",round((i/length(CNV))*100,2),"%")


		## GENOTYPE BRUTE FORCE
		for (j in ID){
			if ("blm" %in% tolower(ID) ){
				if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
					id <- "BLM/RECQL5"
				} else {
					id <- "BLM"
				}
			} else if ("recq5" %in% tolower(ID) ||"recql5" %in% tolower(ID) ){
				if (! "blm" %in% tolower(ID) ){
					id <- "RECQL5"
				}
			} else {id <- "WT" }
		}


		## RETRIEVE FIRST LIBRARY AND REMOVE SMALL SEGMENTS < 500Kb and HIGH PLOIDY STATES TO ASSIGN  PLOIDY
		tmp <- as.data.frame(CNV[i])
		tmp$name <- as.factor(tmp$name)
		tmp=tmp[tmp$width>500000,]
		tmp=tmp[tmp$seqnames!="chrY",]
		#HIGH PLOIDY STATES
		tmp <- tmp[tmp$name!="20-somy" & tmp$name!="19-somy" &tmp$name!="18-somy" &tmp$name!="17-somy" & tmp$name!="16-somy" &tmp$name!="15-somy" &tmp$name!="14-somy"&tmp$name!="13-somy" &tmp$name!="12-somy"&tmp$name!="11-somy"&tmp$name!="10-somy"&tmp$name!="zero-inflation",]
		tmp$name <-droplevels(tmp$name)


		## BUILD PLOIDY TABLE SUMMARY TO CLASSIFY NATIVE STATE AND GAINS & LOSSES
		ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
		ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
		if (ploidy==0){
			ploidy=1
		} else if (ploidy >2){
			next
		}


		## REMOVE SMALL SEGMENTS  AFFTER HAVING  ASSIGNED PLOIDY  STATE
		tmp=tmp[tmp$width>20000000,]
		ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))


		## TURN PLOIDY INTO NUMERIC AND CLASSIFY  GAINS/LOSSES
		ploidyTable=  separate(ploidyTable,col = name,sep="-",into = c("name"))
		ploidyTable$name<- as.numeric(ploidyTable$name)
		ploidyTable =  ploidyTable[ploidyTable$name!=ploidy,]
		ploidyTable$type =  ifelse(ploidyTable$name < ploidy, "loss", ifelse((ploidyTable$name > ploidy), "gain", "unclear"))
		## REPEAT FOR WHOLE DF
		tmp=  separate(tmp,col = name,sep="-",into = c("name"))
		tmp$name<- as.numeric(tmp$name)
		tmp =  tmp[tmp$name!=ploidy,]
		tmp$type =  ifelse(tmp$name < ploidy, "loss", ifelse((tmp$name > ploidy), "gain", "unclear"))


		## QUANTIFY GAINS/LOSSES
		totalGain= sum(filter(ploidyTable,type=="gain")$sum)
		totalGainSeg = sum(filter(ploidyTable,type=="gain")$`n()`)
		totalLoss= sum(filter(ploidyTable,type=="loss")$sum)
		totalLossSeg = sum(filter(ploidyTable,type=="loss")$`n()`)


		#SUMMARY STATS
		row <- data.frame(ID=id,ploidy=ploidy,gainSeg=totalGainSeg,totalGain=totalGain,lossSeg=totalLossSeg,totalLoss=totalLoss,file=file,date=date)
		cnvPerCellSummary <- rbind(row,cnvPerCellSummary)

		#SUMMARIZE INDIVIDUAL EVENTS
		df = tmp %>% select(c(seqnames,start,end,width,name,type))
		if (nrow(df)>0){
			df$file=file
			df$ploidy=ploidy
			df$gene=id
			cnvPerCell<- rbind(df,cnvPerCell)
		}

	}
	write.table(cnvPerCell,"DataSummary/cnvPerCell.txt",row.names = F,col.names = T,quote = F,sep="\t")
	write.table(cnvPerCellSummary,"DataSummary/cnvPerCellSummary.txt",row.names = F,col.names = T,quote = F,sep="\t")
}


######################################################################
			####   EXPORATORY ANALYSIS ON CNVs   ####
######################################################################

#READ IN SAVED DATASETS IF NOT ALREADY LOADED

cnvPerCell=read.table("DataSummary/cnvPerCell.txt",header=T)

cnvPerCellSummary=read.table("DataSummary/cnvPerCellSummary.txt",header=T)




