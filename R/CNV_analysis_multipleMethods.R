######################################################################
                ####   Packages   ####
######################################################################

library(rtracklayer)
library(tidyverse)
library(plyr)

######################################################################
                  ####   Assemble datasets   ####
######################################################################


CNV_hmm<- import("SERVER-OUTPUT/Output-WTvarWidthRef//BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
CNV_dnacopy<- import("SERVER-OUTPUT/Output-WTvarWidthRef//BROWSERFILES/method-dnacopy//binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
CNV_edivisive<- import("SERVER-OUTPUT/Output-WTvarWidthRef//BROWSERFILES/method-edivisive//binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")


cnv_hmm <- as.data.frame(CNV_hmm)
cnv_dnacopy <- as.data.frame(CNV_dnacopy)
cnv_edivisive <- as.data.frame(CNV_edivisive)

length(cnv_hmm)
length(cnv_dnacopy)
length(cnv_edivisive)

cnv_hmm$name<- as.factor(cnv_hmm$name)
cnv_dnacopy$name<- as.factor(cnv_dnacopy$name)
cnv_edivisive$name<- as.factor(cnv_edivisive$name)
cnv_hmm$gene = NA
cnv_dnacopy$gene = NA
cnv_edivisive$gene = NA

for (row in 1:nrow(cnv_hmm)){
  ID <- strsplit(strsplit(cnv_hmm$group_name[row],split = " ")[[1]][4],split = "[-_.]")[[1]]
  
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
  
  cnv_hmm[row,]$gene=id
}
for (row in 1:nrow(cnv_dnacopy)){
  ID <- strsplit(strsplit(cnv_dnacopy$group_name[row],split = " ")[[1]][4],split = "[-_.]")[[1]]
  
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
  
  cnv_dnacopy[row,]$gene=id
}
for (row in 1:nrow(cnv_edivisive)){
  ID <- strsplit(strsplit(cnv_edivisive$group_name[row],split = " ")[[1]][4],split = "[-_.]")[[1]]
  
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
  
  cnv_edivisive[row,]$gene=id
}

cnv_hmm <- separate(cnv_hmm,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
cnv_dnacopy <- separate(cnv_dnacopy,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
cnv_edivisive <- separate(cnv_edivisive,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))

write.table(cnv_hmm,"SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_hmm_with_gene.txt",quote=F,row.names = F,col.names = T,sep="\t")
write.table(cnv_dnacopy,"SERVER-OUTPUT/Output-WTvarWidthRef//BROWSERFILES/CNV_dnacopy_with_gene.txt",quote=F,row.names = F,col.names = T,sep="\t")
write.table(cnv_edivisive,"SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_edivisive_with_gene.txt",quote=F,row.names = F,col.names = T,sep="\t")

cnv_hmm <- read.table("SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_hmm_with_gene.txt",fill=T,header=T,sep="\t",comment.char="")
cnv_dnacopy <- read.table("SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_dnacopy_with_gene.txt",fill=T,header=T,sep="\t",comment.char="")
cnv_edivisive <- read.table("SERVER-OUTPUT/Output-WTvarWidthRef/BROWSERFILES/CNV_edivisive_with_gene.txt",fill=T,header=T,sep="\t",comment.char="")


### METRICS  #####
metrics <- data.frame()
file=list.files("Metrics/",full.names = T)[4]
for (file in list.files("Metrics/",full.names = T)){
  met <- read.table(file,header = T,fill=T,row.names = NULL)
  if ("Library" %in% colnames(met)){
    met <- met %>% select(Library)
    met$file = as.factor(basename(file))
  }
  met$Library=paste0(met$Library,".trimmed.mdup.bam")
  metrics <- rbind(met,metrics)
}
message("\nThere are ",length(metrics$Library)," libraries in the base metrics file with ",ncol(metrics)," features\n")
metrics$date <- as.Date(metrics$date,format="%d-%m-%Y")


cnv_hmm <- merge(cnv_hmm,metrics, by.x="file",by.y="Library")
cnv_dnacopy <- merge(cnv_dnacopy,metrics, by.x="file",by.y="Library")
cnv_edivisive <- merge(cnv_edivisive,metrics, by.x="file",by.y="Library")


### METRICS  #####


#Core functtiton for finding TDs per library
countTDs <- function(CNV){
  cnv_tds<-data.frame()
  tandemRepeats=data.frame()
  
  for (i in 1:length(CNV)){
    
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
    date=as.Date(metrics[metrics$Library==file,3])
    
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
    tmp=tmp[tmp$width>250000,]
    
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
    print(i)
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
    date=as.Date(metrics[metrics$Library==file,3])
    
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
    
    tmp <- tmp[tmp$width>250000,]
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
    cnv_tds <- rbind(tmp,cnv_tds)
  }
  return(cnv_tds)
}


#Run above functions 
tandemRepeatsPerCell_hmm <- countTDs(CNV_hmm)
tandemRepeatsPerCell_dnacopy <- countTDs(CNV_dnacopy)
tandemRepeatsPerCell_edivisive <- countTDs(CNV_edivisive)

tandemRepeatsTotal_hmm<- countTDsInCNV(CNV_hmm)
tandemRepeatsTotal_dnacopy<- countTDsInCNV(CNV_dnacopy)
tandemRepeatsTotal_edivisive<- countTDsInCNV(CNV_edivisive)

tandemRepeatsPerCell_hmm$norm_td <- tandemRepeatsPerCell_hmm$TD/as.numeric(tandemRepeatsPerCell_hmm$ploidy)
tandemRepeatsPerCell_dnacopy$norm_td <- tandemRepeatsPerCell_dnacopy$TD/as.numeric(tandemRepeatsPerCell_dnacopy$ploidy)
tandemRepeatsPerCell_edivisive$norm_td <- tandemRepeatsPerCell_edivisive$TD/as.numeric(tandemRepeatsPerCell_edivisive$ploidy)

#summarize
tandemRepeatsPerCell_dnacopy%>% group_by(ID)%>%dplyr::summarize(n())
tandemRepeatsPerCell_dnacopy%>% group_by(ID)%>%dplyr::summarize(sum(TD))
tandemRepeatsPerCell_dnacopy%>% group_by(ID)%>%dplyr::summarize(mean(TD))
tandemRepeatsPerCell_dnacopy%>% group_by(ID)%>%dplyr::summarize(mean(norm_td))

tandemRepeatsPerCell_edivisive%>% group_by(ID)%>%dplyr::summarize(n())
tandemRepeatsPerCell_edivisive%>% group_by(ID)%>%dplyr::summarize(sum(TD))
tandemRepeatsPerCell_edivisive%>% group_by(ID)%>%dplyr::summarize(mean(TD))
tandemRepeatsPerCell_edivisive%>% group_by(ID)%>%dplyr::summarize(mean(norm_td))

tandemRepeatsPerCell_hmm%>% group_by(ID)%>%dplyr::summarize(n())
tandemRepeatsPerCell_hmm%>% group_by(ID)%>%dplyr::summarize(sum(TD))
tandemRepeatsPerCell_hmm%>% group_by(ID)%>%dplyr::summarize(mean(TD))
tandemRepeatsPerCell_hmm%>% group_by(ID)%>%dplyr::summarize(mean(norm_td))




tandemRepeatsPerCell_dnacopy$norm <- tandemRepeatsPerCell_dnacopy$td_sum/as.numeric(tandemRepeatsPerCell_dnacopy$ploidy)
tandemRepeatsPerCell_dnacopy$norm_td <- tandemRepeatsPerCell_dnacopy$TD/as.numeric(tandemRepeatsPerCell_dnacopy$ploidy)
td_dna_content_dnacopy <- tandemRepeatsPerCell_dnacopy %>% group_by(ID) %>% dplyr::summarize(mean(norm))
td_dna_content_dnacopy$perc <- ((td_dna_content_dnacopy$`mean(norm)`/3e+09)*100)
td_dna_content_dnacopy$name="DNA-COPY"

tandemRepeatsPerCell_edivisive$norm <- tandemRepeatsPerCell_edivisive$td_sum/as.numeric(tandemRepeatsPerCell_edivisive$ploidy)
tandemRepeatsPerCell_edivisive$norm_td <- tandemRepeatsPerCell_edivisive$TD/as.numeric(tandemRepeatsPerCell_edivisive$ploidy)
td_dna_content_edivisive <- tandemRepeatsPerCell_edivisive %>% group_by(ID) %>% dplyr::summarize(mean(norm))
td_dna_content_edivisive$perc <- ((td_dna_content_edivisive$`mean(norm)`/3e+09)*100)
td_dna_content_edivisive$name="EDIVISIVE"

tandemRepeatsPerCell_hmm$norm <- tandemRepeatsPerCell_hmm$td_sum/as.numeric(tandemRepeatsPerCell_hmm$ploidy)
tandemRepeatsPerCell_hmm$norm_td <- tandemRepeatsPerCell_hmm$TD/as.numeric(tandemRepeatsPerCell_hmm$ploidy)
td_dna_content_hmm <- tandemRepeatsPerCell_hmm %>% group_by(ID) %>% dplyr::summarize(mean(norm))
td_dna_content_hmm$perc <- ((td_dna_content_hmm$`mean(norm)`/3e+09)*100)
td_dna_content_hmm$name="HMM"






tandemRepeatsTotal_filtHighPloidy_dnacopy <- tandemRepeatsTotal_dnacopy[tandemRepeatsTotal_dnacopy$name!="20-somy" & tandemRepeatsTotal_dnacopy$name!="19-somy" &tandemRepeatsTotal_dnacopy$name!="18-somy" &tandemRepeatsTotal_dnacopy$name!="17-somy" & tandemRepeatsTotal_dnacopy$name!="16-somy" &tandemRepeatsTotal_dnacopy$name!="15-somy" &tandemRepeatsTotal_dnacopy$name!="14-somy"&tandemRepeatsTotal_dnacopy$name!="13-somy" &tandemRepeatsTotal_dnacopy$name!="12-somy"&tandemRepeatsTotal_dnacopy$name!="11-somy"&tandemRepeatsTotal_dnacopy$name!="10-somy"&tandemRepeatsTotal_dnacopy$name!="zero-inflation"&tandemRepeatsTotal_dnacopy$name!="21-somy"&tandemRepeatsTotal_dnacopy$name!="22-somy"&tandemRepeatsTotal_dnacopy$name!="23-somy"&tandemRepeatsTotal_dnacopy$name!="24-somy",] 
tandemRepeatsTotal_filtHighPloidy_dnacopy$name <-droplevels(tandemRepeatsTotal_filtHighPloidy_dnacopy$name)
tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb)
tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb)

tandemRepeatsTotal_filtHighPloidy_edivisive <- tandemRepeatsTotal_edivisive[tandemRepeatsTotal_edivisive$name!="20-somy" & tandemRepeatsTotal_edivisive$name!="19-somy" &tandemRepeatsTotal_edivisive$name!="18-somy" &tandemRepeatsTotal_edivisive$name!="17-somy" & tandemRepeatsTotal_edivisive$name!="16-somy" &tandemRepeatsTotal_edivisive$name!="15-somy" &tandemRepeatsTotal_edivisive$name!="14-somy"&tandemRepeatsTotal_edivisive$name!="13-somy" &tandemRepeatsTotal_edivisive$name!="12-somy"&tandemRepeatsTotal_edivisive$name!="11-somy"&tandemRepeatsTotal_edivisive$name!="10-somy"&tandemRepeatsTotal_edivisive$name!="zero-inflation"&tandemRepeatsTotal_edivisive$name!="21-somy"&tandemRepeatsTotal_edivisive$name!="22-somy"&tandemRepeatsTotal_edivisive$name!="23-somy"&tandemRepeatsTotal_edivisive$name!="24-somy"&tandemRepeatsTotal_edivisive$name!="25-somy"&tandemRepeatsTotal_edivisive$name!="28-somy",] 
tandemRepeatsTotal_filtHighPloidy_edivisive$name <-droplevels(tandemRepeatsTotal_filtHighPloidy_edivisive$name)
tandemRepeatsTotal_filtHighPloidy_edivisive$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy_edivisive$itemRgb)
tandemRepeatsTotal_filtHighPloidy_edivisive$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy_edivisive$itemRgb)

tandemRepeatsTotal_filtHighPloidy_hmm <- tandemRepeatsTotal_hmm[tandemRepeatsTotal_hmm$name!="20-somy" & tandemRepeatsTotal_hmm$name!="19-somy" &tandemRepeatsTotal_hmm$name!="18-somy" &tandemRepeatsTotal_hmm$name!="17-somy" & tandemRepeatsTotal_hmm$name!="16-somy" &tandemRepeatsTotal_hmm$name!="15-somy" &tandemRepeatsTotal_hmm$name!="14-somy"&tandemRepeatsTotal_hmm$name!="13-somy" &tandemRepeatsTotal_hmm$name!="12-somy"&tandemRepeatsTotal_hmm$name!="11-somy"&tandemRepeatsTotal_hmm$name!="10-somy"&tandemRepeatsTotal_hmm$name!="zero-inflation"&tandemRepeatsTotal_hmm$name!="21-somy"&tandemRepeatsTotal_hmm$name!="22-somy"&tandemRepeatsTotal_hmm$name!="23-somy"&tandemRepeatsTotal_hmm$name!="24-somy"&tandemRepeatsTotal_hmm$name!="25-somy"&tandemRepeatsTotal_hmm$name!="28-somy",] 
tandemRepeatsTotal_filtHighPloidy_hmm$name <-droplevels(tandemRepeatsTotal_filtHighPloidy_hmm$name)
tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb)
tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb)




#separate into haploid and diploid TDs
haploidPerCell_dnacopy <- filter(tandemRepeatsPerCell_dnacopy,ploidy==1)
diploidPerCell_dnacopy <- filter(tandemRepeatsPerCell_dnacopy,ploidy==2)
haploidTotal_dnacopy <- filter(tandemRepeatsTotal_filtHighPloidy_dnacopy,ploidy==1)
diploidTotal_dnacopy <- filter(tandemRepeatsTotal_filtHighPloidy_dnacopy,ploidy==2)

haploidPerCell_edivisive <- filter(tandemRepeatsPerCell_edivisive,ploidy==1)
diploidPerCell_edivisive <- filter(tandemRepeatsPerCell_edivisive,ploidy==2)
haploidTotal_edivisive <- filter(tandemRepeatsTotal_filtHighPloidy_edivisive,ploidy==1)
diploidTotal_edivisive <- filter(tandemRepeatsTotal_filtHighPloidy_edivisive,ploidy==2)

haploidPerCell_hmm <- filter(tandemRepeatsPerCell_hmm,ploidy==1)
diploidPerCell_hmm <- filter(tandemRepeatsPerCell_hmm,ploidy==2)
haploidTotal_hmm <- filter(tandemRepeatsTotal_filtHighPloidy_hmm,ploidy==1)
diploidTotal_hmm <- filter(tandemRepeatsTotal_filtHighPloidy_hmm,ploidy==2)



######################################################################
              #### PLOTTING ####
######################################################################

tandemRepeatsPerCell_dnacopy$name="DNA-COPY"
tandemRepeatsPerCell_edivisive$name="EDIVISIVE"
tandemRepeatsPerCell_hmm$name="HMM"

genomicContentGains <- rbind(tandemRepeatsPerCell_edivisive,tandemRepeatsPerCell_dnacopy,tandemRepeatsPerCell_hmm)
genomicContentGains$perc <- ((genomicContentGains$norm/3e+09)*100)

ggplot(genomicContentGains)+geom_violin(aes(ID,perc,fill=name))+
  facet_wrap(~name)+theme_classic()+
  theme(legend.position = "none")+
  labs(y="% Gain in Genomic Content/Hap Genome")+
  ggsave("Plots/genomicContentGainss.png")


plotting <- function(tandemRepeatsPerCell,tandemRepeatsTotal_filtHighPloidy,haploidPerCell,diploidPerCell,haploidTotal,diploidTotal){
  
  method = tandemRepeatsPerCell$name[1]
  plotsPath=paste0("Plots/WT_REF/",method,"/")
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
    scale_x_log10(limits=c(1e+05,1e+08))+# = c(1e+01,1e+04,1e+5,1e+6,1e+7,1e+8)) + 
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


plotting(tandemRepeatsPerCell_dnacopy,tandemRepeatsTotal_filtHighPloidy_dnacopy,haploidPerCell_dnacopy,diploidPerCell_dnacopy,haploidTotal_dnacopy,diploidTotal_dnacopy)
plotting(tandemRepeatsPerCell_edivisive,tandemRepeatsTotal_filtHighPloidy_edivisive,haploidPerCell_edivisive,diploidPerCell_edivisive,haploidTotal_edivisive,diploidTotal_edivisive)
plotting(tandemRepeatsPerCell_hmm,tandemRepeatsTotal_filtHighPloidy_hmm,haploidPerCell_hmm,diploidPerCell_hmm,haploidTotal_hmm,diploidTotal_hmm)


######################################################################
                    #### SAVING TABLES ####
######################################################################

for (i in c("hmm","edivisive","dnacopy")){
  for (j in c("tandemRepeatsPerCell","tandemRepeatsTotal_filtHighPloidy","haploidPerCell","diploidPerCell","haploidTotal","diploidTotal")){
    name <- paste0(j,"_",i)
    obj <- get(name)
    write.table(obj,paste0("DataSummary/",toupper(i),"/",name,".txt"),quote = F,col.names = T,row.names = F,sep="\t")
  }
}


######################################################################
              ####  EXPLORATORY ANALYSIS  ####
######################################################################
metricsfileDir<- "Metrics/"

collectMetrics <- function(metricsfileDir){
  
  metrics <- data.frame()

  for (file in list.files(metricsfileDir,full.names = T)){
    met <- read.table(file,header = T,fill=T,row.names = NULL)
    if ("Postfiltering_reads_aligned" %in% colnames(met)){
      met <- met %>% select(Library,Postfiltering_reads_aligned,Mode_GC)
      met$file = basename(file)
    }else if ("Reads_aligned_postfiltering" %in% colnames(met)){
      met <- met %>% select(Library,Reads_aligned_postfiltering,Mode_GC)
      met$Postfiltering_reads_aligned <- met$Reads_aligned_postfiltering
      met <- select(met,-c(Reads_aligned_postfiltering))
      met$file = basename(file)
    }
    met$Library=paste0(met$Library,".trimmed.mdup.bam")
    metrics <- rbind(met,metrics)
  }
  
  return(metrics)
}

metrics <- collectMetrics("Metrics/")
metrics$date=revalue(metrics$file, c("Aug18.txt" ="18-08-2020", "Aug28.txt"  ="28-08-2020" ,"Feb25Metrics.txt" ="25-02-2019",        
                                     "metrics_summary-jan18.txt"="18-01-2021", "oct23Metrics.txt"="23-10-2021","March6Metrics.txt"="06-03-2019"))


merge <- merge(metrics,tandemRepeatsPerCell_hmm,by.x="Library",by.y="file")

ggplot(merge)+geom_point(aes(Postfiltering_reads_aligned,norm_td),alpha=0.5)+facet_wrap(~ID)


merge$Mode_GC<- as.numeric(gsub('%', '', merge$Mode_GC))
ggplot(merge)+geom_point(aes(Mode_GC,norm_td),alpha=0.5)+facet_wrap(~ID)+
  theme_classic()+
  ggsave("Plots/HMM/GC.png")
ggplot(merge)+geom_point(aes(Postfiltering_reads_aligned,norm_td),alpha=0.5)+
  facet_wrap(~ID)+
  theme_classic()+
  ggsave("Plots/HMM/reads.png")






tandemRepeatsPerCell=tandemRepeatsPerCell_hmm
tandemRepeatsTotal_filtHighPloidy=tandemRepeatsTotal_filtHighPloidy_hmm
haploidPerCell=haploidPerCell_hmm
diploidPerCell=diploidPerCell_hmm
haploidTotal=haploidTotal_hmm
diploidTotal=diploidTotal_hmm
tandemRepeatsPerCell_hmm$ploidy<-as.factor(tandemRepeatsPerCell_hmm$ploidy)
merge$date.x
test <- merge %>% filter(Postfiltering_reads_aligned > 250000)
test2 <- test %>% filter(ploidy==1)
merge$date.y <- as.Date(merge$date.y,format="Y%-m%-d%")
ggplot(test)+geom_point(aes(date.y,TD,color=ploidy, size=Postfiltering_reads_aligned))+
  geom_smooth(aes(date.y,TD,color=ploidy),method = "lm",se=F)+
  facet_wrap(~ID) + labs(title = "No read count cut-off")+
  ggsave("NoReadCountCutOff.png")

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
  scale_x_log10(limits=c(1e+05,1e+08))+# = c(1e+01,1e+04,1e+5,1e+6,1e+7,1e+8)) + 
  ggsave(paste0(plotsPath,"6.png"))

ggplot(tandemRepeatsPerCell,aes(ID,TD))+geom_violin(aes(fill=ID))+
  geom_boxplot(width=0.05)+scale_y_log10()+ theme_classic()+
  labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
  ggsave(paste0(plotsPath,"7.png"))

ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(gene,fill=name))+
  scale_fill_manual(values =  c("2-somy"="#00EE76" , "3-somy"="#CD0000",
                                "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                "8-somy" ="#458B00","7-somy"="#1E90FF",
                                "6-somy"="#FFFACD")) + theme_classic()

