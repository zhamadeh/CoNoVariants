######################################################################
####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(plyr)

######################################################################
####   Assemble datasets   ####
######################################################################


assembleTandemDuplicationDF <- function(CNV="SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz",
                                        CNVgene = "SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_hmm_with_gene.txt",
                                        metricsDir="Metrics/",name="HMM"){
  
  #ANNOTATE GENE INFORMATION
  CNV<- import(CNV)
  CNVgene <- read.table(CNVgene,fill=T,header=T,sep="\t",comment.char="")
  
  ### METRICS  #####
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
  
  metrics = collectMetrics(metricsDir)
  
  #cnv <- merge(CNVgene,metrics, by.x="file",by.y="Library")
  
  #Core functtiton for finding TDs per library
  countTDs <- function(CNV){
    cnv_tds<-data.frame()
    tandemRepeats=data.frame()
    
    for (i in 1:length(CNV)){
      
      ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
      file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
      date=as.Date(metrics[metrics$Library==file,]$date)
      
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
      
      tmp <- as.data.frame(CNV[i])
      tmp$name <- as.factor(tmp$name)
      tmp=tmp[tmp$width>500000,]
      
      tmp <- tmp[tmp$name!="20-somy" & tmp$name!="19-somy" &tmp$name!="18-somy" &tmp$name!="17-somy" & tmp$name!="16-somy" &tmp$name!="15-somy" &tmp$name!="14-somy"&tmp$name!="13-somy" &tmp$name!="12-somy"&tmp$name!="11-somy"&tmp$name!="10-somy"&tmp$name!="zero-inflation",] 
      tmp$name <-droplevels(tmp$name)
      
      
      ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
      ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
      if (ploidy==0){
        ploidy=1
      } else if (ploidy >2){
        next
      }
      remove=c()
      
      for (row in 1:nrow(ploidyTable)){
        #print(row)
        #print(ploidyTable[row,1]$name)
        if (ploidyTable[row,1]$name=="zero-inflation"){
          remove <- append(row,remove)
        }
        else {
          level = as.numeric(strsplit(as.character(ploidyTable[row,1]$name),split = "[-]")[[1]][1])
          if (level <= ploidy){
            remove <- append(row,remove)
          }
        }
      }
      if (length(remove)>0){
        td=ploidyTable[-c(remove),]
        td <- sum(td$`n()`)
      } else {td <- sum(ploidyTable$`n()`)}
      
      for (row in 1:nrow(tmp)){
        
        if (tmp[row,]$name=="zero-inflation"){
          remove <- append(row,remove)
        }
        else {
          level = as.numeric(strsplit(as.character(tmp[row,]$name),split = "[-]")[[1]][1])
          if (level <= ploidy){
            remove <- append(row,remove)
          }
        }
      }
      if (length(remove)>0){
        tmp <-tmp[-c(remove),]
        td_count <- sum(tmp$width)
      } else {
        tmp <-tmp[-c(remove),]
        td_count <- sum(tmp$width)
      }
      
      row <- data.frame(ID=id,TD=td,ploidy=ploidy,td_sum=td_count,file=file,date=date)
      tandemRepeats <- rbind(row,tandemRepeats)
    }
    return(tandemRepeats)
  }
  #Core functtion for finding all TDs
  countTDsInCNV <- function(CNV){
    cnv_tds<-data.frame()
    tandemRepeats=data.frame()
    
    
    for (i in 1:length(CNV)){
      #print(i)
      ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
      file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
      date=as.Date(metrics[metrics$Library==file,]$date)
      
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
      
      tmp <- as.data.frame(CNV[i])
      tmp$name <- as.factor(tmp$name)
      
      
      ploidyTable <- tmp %>% group_by(name) %>% dplyr::summarize(n(),sum=sum(width))
      ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$sum),]$name),split = "[-]")[[1]][1])
      
      tmp <- tmp[tmp$width>500000,]
      if (ploidy==0){
        ploidy=1
      } else if (ploidy >2){
        next
      }
      tmp$ploidy <- ploidy
      remove=c()
      
      for (row in 1:nrow(tmp)){
        
        if (tmp[row,]$name=="zero-inflation"){
          remove <- append(row,remove)
        }
        else {
          level = as.numeric(strsplit(as.character(tmp[row,]$name),split = "[-]")[[1]][1])
          if (level <= ploidy){
            remove <- append(row,remove)
          }
        }
      }
      if (length(remove)>0){
        td_count <- sum(tmp[-c(remove),]$width)
      } else {td_count <- sum(tmp[-c(remove),]$width)}
      
      tmp <- as.data.frame(tmp[-c(remove),])
      tmp<-tmp %>% select(-c(group,group_name))
      if (nrow(tmp)>0){
        tmp$gene <- id
        tmp$file = as.factor(file)
      }
      
      tmp <- merge(tmp,metrics,by.x="file",by.y="Library")
      cnv_tds <- rbind(tmp,cnv_tds)
    }
    return(cnv_tds)
  }
  
  
  #Run above functions 
  tandemRepeatsPerCell <- countTDs(CNV)
  
  tandemRepeatsPerCell <- merge(tandemRepeatsPerCell,metrics,by.x="file",by.y="Library")
  
  tandemRepeatsTotal<- countTDsInCNV(CNV)
  
  tandemRepeatsPerCell$norm_td <- tandemRepeatsPerCell$TD/as.numeric(tandemRepeatsPerCell$ploidy)
  
  tandemRepeatsPerCell$norm <- tandemRepeatsPerCell$td_sum/as.numeric(tandemRepeatsPerCell$ploidy)
  tandemRepeatsPerCell$norm_td <- tandemRepeatsPerCell$TD/as.numeric(tandemRepeatsPerCell$ploidy)
  td_dna_content <- tandemRepeatsPerCell %>% group_by(ID) %>% dplyr::summarize(mean(norm))
  td_dna_content$perc <- ((td_dna_content$`mean(norm)`/3e+09)*100)
  td_dna_content$name=name
  
  tandemRepeatsTotal_filtHighPloidy <- tandemRepeatsTotal[tandemRepeatsTotal$name!="20-somy" & tandemRepeatsTotal$name!="19-somy" &tandemRepeatsTotal$name!="18-somy" &tandemRepeatsTotal$name!="17-somy" & tandemRepeatsTotal$name!="16-somy" &tandemRepeatsTotal$name!="15-somy" &tandemRepeatsTotal$name!="14-somy"&tandemRepeatsTotal$name!="13-somy" &tandemRepeatsTotal$name!="12-somy"&tandemRepeatsTotal$name!="11-somy"&tandemRepeatsTotal$name!="10-somy"&tandemRepeatsTotal$name!="zero-inflation"&tandemRepeatsTotal$name!="21-somy"&tandemRepeatsTotal$name!="22-somy"&tandemRepeatsTotal$name!="23-somy"&tandemRepeatsTotal$name!="24-somy"&tandemRepeatsTotal$name!="25-somy"&tandemRepeatsTotal$name!="28-somy",] 
  tandemRepeatsTotal_filtHighPloidy$name <-droplevels(tandemRepeatsTotal_filtHighPloidy$name)
  tandemRepeatsTotal_filtHighPloidy$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy$itemRgb)
  tandemRepeatsTotal_filtHighPloidy$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy$itemRgb)
  
  
  #separate into haploid and diploid TDs
  haploidPerCell <- filter(tandemRepeatsPerCell,ploidy==1)
  diploidPerCell <- filter(tandemRepeatsPerCell,ploidy==2)
  haploidTotal <- filter(tandemRepeatsTotal_filtHighPloidy,ploidy==1)
  diploidTotal <- filter(tandemRepeatsTotal_filtHighPloidy,ploidy==2)
  
  datapath=paste0("DATA/",name)
  if (!file.exists(datapath) ) { dir.create(datapath)}
  tablepath=paste0(datapath,"/TABLES/")
  if (!file.exists(tablepath) ) { dir.create(tablepath)}
  plotspath=paste0(datapath,"/PLOTS/")
  if (!file.exists(plotspath) ) { dir.create(plotspath)}
  
  write.table(CNVgene,paste0(tablepath,"/CNV_with_gene.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(haploidPerCell,paste0(tablepath,"/haploidPerCell.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(haploidTotal,paste0(tablepath,"/haploidTotal.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(diploidPerCell,paste0(tablepath,"/diploidPerCell.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(diploidTotal,paste0(tablepath,"/diploidTotal.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(tandemRepeatsTotal_filtHighPloidy,paste0(tablepath,"/TD_total_filtHighPloidy.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  write.table(tandemRepeatsPerCell,paste0(tablepath,"/tandemRepeatsPerCell.txt"),quote=F,row.names = F,col.names = T,sep="\t")
  
  plotting(plotspath,tandemRepeatsPerCell,tandemRepeatsTotal_filtHighPloidy,haploidPerCell,diploidPerCell,haploidTotal,diploidTotal)

}

assembleTandemDuplicationDF()

plotting <- function(plotspath,tandemRepeatsPerCell,tandemRepeatsTotal_filtHighPloidy,haploidPerCell,diploidPerCell,haploidTotal,diploidTotal){
  
  plotsPath=plotspath
  tandemRepeatsPerCell$ploidy<-as.factor(tandemRepeatsPerCell$ploidy)
  tandemRepeatsTotal_filtHighPloidy$ploidy<- as.factor(tandemRepeatsTotal_filtHighPloidy$ploidy)
  
  ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(ploidy,fill=gene))+
    facet_wrap(~gene)+
    theme_classic()+
    ggsave(paste0(plotsPath,"1.png"))
  
  tandemRepeatsTotal_filtHighPloidy$seqnames <- factor(tandemRepeatsTotal_filtHighPloidy$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                                                                                             "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
  ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(seqnames,fill=ploidy))+
    theme_classic()+facet_wrap(~gene)+ggsave(paste0(plotsPath,"3.png"))
  
  ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
    theme_classic() +labs(title = "Aneuploidy count/library")+
    scale_x_log10()+ ggsave(paste0(plotsPath,"4.png"))
  
  ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,fill=ID),alpha=0.7)+
    theme_classic() +labs(title = "Aneuploidy count/library")+
    facet_wrap(~ID)+
    scale_x_log10()+ ggsave(paste0(plotsPath,"2.png"))
  
  tandemRepeatsPerCell$ID<-factor(tandemRepeatsPerCell$ID,levels=c("WT","RECQL5", "BLM","BLM/RECQL5"))
  ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,group=ploidy,fill=ploidy),alpha=0.7)+
    theme_classic() +labs(title = "Aneuploidy count/library")+
    facet_wrap(~ID)+
    scale_x_log10()+ ggsave(paste0(plotsPath,"5.png"))
  
  ggplot()+geom_density(data = tandemRepeatsTotal_filtHighPloidy,aes(x = width,fill=gene),alpha=0.5)+
    theme_classic()+labs(title = "CNV segment size")+
    #scale_x_log10(limits=c(1e+05,1e+08))+# = c(1e+01,1e+04,1e+5,1e+6,1e+7,1e+8)) + 
    ggsave(paste0(plotsPath,"6.png"))
  
  ggplot(tandemRepeatsPerCell,aes(ID,TD))+geom_violin(aes(fill=ID))+
    geom_boxplot(width=0.05)+scale_y_log10()+ theme_classic()+
    labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
    ggsave(paste0(plotsPath,"7.png"))
  
  ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(gene,fill=name))+
    scale_fill_manual(values =  c("2-somy"="#00EE76" , "3-somy"="#CD0000",
                                  "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                  "8-somy" ="#458B00","7-somy"="#1E90FF",
                                  "6-somy"="#FFFACD")) + theme_classic()+
    ggsave(paste0(plotsPath,"8.png"))
}

nrow(tandemRepeatsTotal_filtHighPloidy)
tandemRepeatsTotal_filtHighPloidy_40Mb= filter(tandemRepeatsTotal_filtHighPloidy,width>20000000)
tandemRepeatsTotal_filtHighPloidy_40Mb_noY= filter(tandemRepeatsTotal_filtHighPloidy,seqnames!="chrY")



ggplot(tandemRepeatsPerCell)+geom_point(aes(TD,Reads,color=as.factor(ploidy)))+facet_wrap(~ID)+
  ggsave("TDvsREADS.png")

ggplot(tandemRepeatsPerCell)+geom_point(aes(TD,Mode_GC,color=as.factor(ploidy)))+facet_wrap(~ID)+
  ggsave("TDvsMODE_GC.png")



