######################################################################
              ####CNV analysis of good libraries (n=317)####
######################################################################
#packages
library(tidyverse)
library(rtracklayer)

#data import
CNV<- import("SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
CNV_dnacopy<- import("SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/method-dnacopy//binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
CNV_edivisive<- import("SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/method-edivisive//binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")

#formating data
cnv <- as.data.frame(CNV)
cnv$file <- as.factor(cnv$file)
cnv$name<- as.factor(cnv$name)
cnv$gene = NA
cnv <- separate(cnv,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
for (row in 1:nrow(cnv)){
  ID <- strsplit(strsplit(cnv$group_name[row],split = " ")[[1]][4],split = "[-_.]")[[1]]
  
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
  
  cnv[row,]$gene=id
}
write.table(cnv,"SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/CNV_with_gene.txt",quote=F,row.names = F,col.names = T,sep="\t")


#import intead of running above lines
cnv <-read.table("SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/method-HMM/CNV_with_gene.txt",header=T,fill=T,sep="\t")



######################################################################
        #### Collect copy number gains (TDs) ####
######################################################################

#Core functtiton for finding TDs per library
countTDs <- function(CNV){
  cnv_tds<-data.frame()
  tandemRepeats=data.frame()

  for (i in 1:length(CNV)){

    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
    
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
    
    
    ploidyTable <- tmp %>% group_by(name) %>% summarize(n(),sum=sum(width))
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
    
    row <- data.frame(ID=id,TD=td,ploidy=ploidy,td_sum=td_count,file=file)
    tandemRepeats <- rbind(row,tandemRepeats)
  }
  return(tandemRepeats)
}
#Core functtion for finding all TDs
countTDsInCNV <- function(CNV){
  cnv_tds<-data.frame()
  tandemRepeats=data.frame()
  
  for (i in 1:length(CNV)){
    
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4]
    
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
    ploidyTable <- tmp %>% group_by(name) %>% summarize(n(),sum=sum(width))
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
    tmp$gene <- id
    tmp$file = as.factor(file)
    cnv_tds <- rbind(tmp,cnv_tds)
  }
  return(cnv_tds)
}


#Run above functions 
tandemRepeatsPerCell <- countTDs(CNV)
tandemRepeatsTotal<- countTDsInCNV(CNV)

#separate into haploid and diploid TDs
haploidPerCell <- filter(tandemRepeatsPerCell,ploidy==1)
diploidPerCell <- filter(tandemRepeatsPerCell,ploidy==2)
haploidTotal <- filter(tandemRepeatsTotal_filtHighPloidy,ploidy==1)
diploidTotal <- filter(tandemRepeatsTotal_filtHighPloidy,ploidy==2)


#summarize
tandemRepeatsPerCell%>% group_by(ID)%>%summarize(n())
tandemRepeatsPerCell%>% group_by(ID)%>%summarize(sum(TD))
tandemRepeatsPerCell%>% group_by(ID)%>%summarize(mean(TD))
tandemRepeatsPerCell%>% group_by(ID)%>%summarize(mean(norm_td))


#tandemRepeatsPerCell<-tandemRepeatsPerCell[tandemRepeatsPerCell$td_sum<5e+08,]
#tandemRepeatsTotal <- tandemRepeatsTotal[tandemRepeatsTotal$width>250000,]
tandemRepeatsPerCell$norm <- tandemRepeatsPerCell$td_sum/as.numeric(tandemRepeatsPerCell$ploidy)
tandemRepeatsPerCell$norm_td <- tandemRepeatsPerCell$TD/as.numeric(tandemRepeatsPerCell$ploidy)
td_dna_content <- tandemRepeatsPerCell %>% group_by(ID) %>% summarize(mean(norm))
td_dna_content$perc <- ((td_dna_content$`mean(norm)`/3e+09)*100)
td_dna_content 
#write.table(td_dna_content,"percentGainOfDNA.txt",quote = F,col.names = T,row.names = F,sep = "\t")


tandemRepeatsTotal_filtHighPloidy <- tandemRepeatsTotal[tandemRepeatsTotal$name!="20-somy" & tandemRepeatsTotal$name!="19-somy" &tandemRepeatsTotal$name!="18-somy" &tandemRepeatsTotal$name!="17-somy" & tandemRepeatsTotal$name!="16-somy" &tandemRepeatsTotal$name!="15-somy" &tandemRepeatsTotal$name!="14-somy"&tandemRepeatsTotal$name!="13-somy" &tandemRepeatsTotal$name!="12-somy"&tandemRepeatsTotal$name!="11-somy"&tandemRepeatsTotal$name!="10-somy"&tandemRepeatsTotal$name!="zero-inflation",] 
tandemRepeatsTotal_filtHighPloidy$name <-droplevels(tandemRepeatsTotal_filtHighPloidy$name)
tandemRepeatsTotal_filtHighPloidy$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy$itemRgb)
tandemRepeatsTotal_filtHighPloidy$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy$itemRgb)



######################################################################
                    ####  PLOTTING  ####
######################################################################

ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(ploidy))


ggplot()+geom_density(data = filter(tandemRepeatsPerCell,ID=="BLM"),aes(x = TD,color=ID))+
  geom_density(data = filter(tandemRepeatsPerCell,ID=="BLM/RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeatsPerCell,ID=="RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeatsPerCell,ID=="WT"),aes(x = TD,color=ID))+ 
  theme_bw() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()

tandemRepeatsTotal_filtHighPloidy$seqnames <- factor(tandemRepeatsTotal_filtHighPloidy$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                               "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(seqnames),fill="#4a826b")+
  theme_classic()+facet_wrap(~gene)+ggsave("Plots/CNV_per_chr.png")

ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library.png")

ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,fill=ID),alpha=0.7)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  facet_wrap(~ID)+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library-facet.png")

ggplot()+geom_density(data = tandemRepeatsTotal_filtHighPloidy,aes(x = width,color=gene),size=1.5)+
  theme_bw() +labs(title = "CNV segment size")+
  scale_x_log10(breaks = c(1e+5,1e+6,1e+7,1e+8)) + 
  ggsave("Plots/CNV_segment_size.png")

ggplot(tandemRepeatsPerCell,aes(ID,TD))+geom_violin(aes(fill=ID))+
  geom_boxplot(width=0.05)+scale_y_log10()+ theme_bw()+
  labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
  ggsave("Plots/CNV_per_cell.png")

ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(gene,fill=name))+
  scale_fill_manual(values =  c("2-somy"="#00EE76" , "3-somy"="#CD0000",
                                "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                "8-somy" ="#458B00","7-somy"="#1E90FF",
                                "6-somy"="#FFFACD")) + theme_bw()+
  ggsave("Plots/spreadOfPloidyAcrossGenes.png")

######################################################################
                    ####  PLOTTING HAP vs DIP ####
######################################################################

ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(ploidy))


ggplot()+geom_density(data = filter(tandemRepeatsPerCell,ID=="BLM"),aes(x = TD,color=ID))+
  geom_density(data = filter(tandemRepeatsPerCell,ID=="BLM/RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeatsPerCell,ID=="RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeatsPerCell,ID=="WT"),aes(x = TD,color=ID))+ 
  theme_bw() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()

tandemRepeatsTotal_filtHighPloidy$seqnames <- factor(tandemRepeatsTotal_filtHighPloidy$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                                                                                           "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
tandemRepeatsTotal_filtHighPloidy$ploidy<-as.factor(as.character(tandemRepeatsTotal_filtHighPloidy$ploidy))
ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(seqnames,fill=ploidy))+
  scale_fill_manual(values =  c("1"="#00EE76" , "2"="#CD0000"))+
  theme_classic()+facet_wrap(~gene)+ggsave("Plots/CNV_per_chr.png")

ggplot(filter(tandemRepeatsPerCell,ploidy==1))+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library.png")

ggplot(filter(tandemRepeatsPerCell,ploidy==2))+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library.png")

tandemRepeatsPerCell$ploidy<-as.factor(as.character(tandemRepeatsPerCell$ploidy))
ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,fill=ploidy),alpha=0.7)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  facet_wrap(~ID)+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library-facet.png")

ggplot()+geom_density(data = tandemRepeatsTotal_filtHighPloidy,aes(x = width,color=gene),size=1.5)+
  theme_bw() +labs(title = "CNV segment size")+
  scale_x_log10(breaks = c(1e+5,1e+6,1e+7,1e+8)) + 
  ggsave("Plots/CNV_segment_size.png")

ggplot(tandemRepeatsPerCell,aes(ID,TD))+geom_violin(aes(fill=ploidy))+
  geom_boxplot(width=0.05)+scale_y_log10()+ theme_bw()+
  labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
  ggsave("Plots/CNV_per_cell.png")

ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(gene,fill=name))+
  scale_fill_manual(values =  c("2-somy"="#00EE76" , "3-somy"="#CD0000",
                                "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                "8-somy" ="#458B00","7-somy"="#1E90FF",
                                "6-somy"="#FFFACD")) + theme_bw()+
  ggsave("Plots/spreadOfPloidyAcrossGenes.png")



  
######################################################################
                    #### EXPORTING ####
######################################################################


blm <- filter(cnv,gene=="BLM") %>% select(c(seqnames,start,end))
write.table(blm,"blmCNV.bed",quote = F,row.names = F,col.names = F,sep = "\t")
blm <- filter(cnv,gene=="BLM/RECQL5") %>% select(c(seqnames,start,end))
write.table(blm,"blm-recql5CNV.bed",quote = F,row.names = F,col.names = F,sep = "\t")
blm <- filter(cnv,gene=="RECQL5") %>% select(c(seqnames,start,end))
write.table(blm,"recql5CNV.bed",quote = F,row.names = F,col.names = F,sep = "\t")
blm <- filter(cnv,gene=="WT") %>% select(c(seqnames,start,end))
write.table(blm,"wtCNV.bed",quote = F,row.names = F,col.names = F,sep = "\t")




