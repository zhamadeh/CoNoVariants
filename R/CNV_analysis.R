######################################################################
              ####CNV analysis of good libraries (n=317)####
######################################################################
#packages
library(tidyverse)
library(rtracklayer)

#data
CNV<- import("SERVER-OUTPUT/Output-VarWidth/BROWSERFILES_varWidthRef/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
hotspots<- import("SERVER-OUTPUT/BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_breakpoint-hotspots.bed")
cnv <- as.data.frame(CNV)
cnv$gene = NA
cnv$name<- as.factor(cnv$name)
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

write.table(as.data.frame(levels(cnv$file)),"goodLibs.txt",quote = F,row.names = F,col.names = F,sep = "\t")
library(stringr)

# Open Connection to file
pathToFile <- path.expand("SERVER-OUTPUT-2/Output-VarWidth/BROWSERFILES_varWidthRef/CNV_with_gene.txt")
cnv <-read.table("SERVER-OUTPUT-2/Output-VarWidth/BROWSERFILES_varWidthRef/CNV_with_gene.txt",header=T,fill=T,sep="\t",comment.char="")

cnv <- separate(cnv,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
cnv$file <- as.factor(cnv$file)

cnv <- cnv[cnv$width>150000,]

#core functions
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
    tmp=tmp[tmp$width>550000,]
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
      #tmp <- tmp[tmp$width > 550000,]
      td_count <- sum(tmp$width)
    } else {
      tmp <-tmp[-c(remove),]
      #tmp <- tmp[tmp$width > 550000,]
      td_count <- sum(tmp$width)
      }
    
    row <- data.frame(ID=id,TD=td,ploidy=ploidy,td_sum=td_count,file=file)
    tandemRepeats <- rbind(row,tandemRepeats)
  }
  return(tandemRepeats)
}
tandemRepeats <- countTDs(CNV)

ggplot(tandemRepeats)+geom_density(aes(td_sum))
tandemRepeats<-tandemRepeats[tandemRepeats$td_sum<5e+08,]


countTDsInCNV <- function(CNV){
  cnv_tds<-data.frame()
  tandemRepeats=data.frame()
  
  for (i in 1:length(CNV)){
    
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    
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
    if (ploidy==0){
      ploidy=1
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
    tmp$gene <- id
    cnv_tds <- rbind(tmp,cnv_tds)
  }
  return(cnv_tds)
}
tandemRepeatsCNV <- countTDsInCNV(CNV)

tandemRepeatsCNV<- tandemRepeatsCNV[tandemRepeatsCNV$ploidy <= 2,]

#not informative because most cell are haploid
countDELs <- function(CNV){
  
  losses=data.frame()
  
  for (i in 1:length(CNV)){
    
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    
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
    if (ploidy==0){
      ploidy=1
    }
    remove=c()
    
    for (row in 1:nrow(ploidyTable)){
      #print(row)
      #print(ploidyTable[row,1]$name)
      if (ploidyTable[row,1]$name=="zero-inflation" || ploidyTable[row,1]$name=="0-somy"){
        remove <- append(row,remove)
      }
      else {
        level = as.numeric(strsplit(as.character(ploidyTable[row,1]$name),split = "[-]")[[1]][1])
        if (level >= ploidy){
          remove <- append(row,remove)
        }
      }
    }
    if (length(remove)>0){
      td <- sum(ploidyTable[-c(remove),]$`n()`)
    } else {td <- sum(ploidyTable$`n()`)}
    
    for (row in 1:nrow(tmp)){
      
      if (tmp[row,]$name=="zero-inflation"|| tmp[row,]$name=="0-somy"){
        remove <- append(row,remove)
      }
      else {
        level = as.numeric(strsplit(as.character(tmp[row,]$name),split = "[-]")[[1]][1])
        if (level >= ploidy){
          remove <- append(row,remove)
        }
      }
    }
    if (length(remove)>0){
      td_count <- sum(tmp[-c(remove),]$width)
    } else {td_count <- sum(tmp[-c(remove),]$width)}
    
    row <- data.frame(ID=id,DEL=td,ploidy=ploidy,td_sum=td_count)
    losses <- rbind(row,losses)
  }
  return(losses)
}
losses <- countDELs(CNV)


#summarize
counts <- tandemRepeats%>% group_by(ID)%>%summarize(n())
tandemRepeats$norm=NA
for (row in 1:nrow(tandemRepeats)){
  if (tandemRepeats[row,]$ID=="WT"){
    tandemRepeats[row,]$norm <- tandemRepeats[row,]$TD/counts[counts$ID=="WT",2]
  } else if (tandemRepeats[row,]$ID=="RECQL5"){
    tandemRepeats[row,]$norm <- tandemRepeats[row,]$TD/counts[counts$ID=="RECQL5",2]
  } else if (tandemRepeats[row,]$ID=="BLM"){
    tandemRepeats[row,]$norm <- tandemRepeats[row,]$TD/counts[counts$ID=="BLM",2]
  } else{
    tandemRepeats[row,]$norm <- tandemRepeats[row,]$TD/counts[counts$ID=="BLM/RECQL5",2]
  }
}
tandemRepeats$ID<-as.factor(tandemRepeats$ID)

m <- tandemRepeats%>% group_by(ID)%>%summarize(sum(TD))
sum(m$`sum(TD)`)

tandemRepeats%>% group_by(ID)%>%summarize(mean(TD))


td_dna_content <- tandemRepeats %>% group_by(ID) %>% summarize(mean(td_sum))
td_dna_content$perc <- ((td_dna_content$`mean(td_sum)`/3e+09)*100)
write.table(td_dna_content,"percentGainOfDNA.txt",quote = F,col.names = T,row.names = F,sep = "\t")

removeHighPloidy <- cnv[cnv$name!="20-somy" & cnv$name!="19-somy" &cnv$name!="18-somy" &cnv$name!="17-somy" & cnv$name!="16-somy" &cnv$name!="15-somy" &cnv$name!="14-somy"&cnv$name!="13-somy" &cnv$name!="12-somy"&cnv$name!="11-somy"&cnv$name!="10-somy"&cnv$name!="zero-inflation",] 
removeHighPloidy$name <-droplevels(removeHighPloidy$name)
removeHighPloidy$itemRgb <-as.factor(removeHighPloidy$itemRgb)
removeHighPloidy$itemRgb <-droplevels(removeHighPloidy$itemRgb)

removeHighPloidy2 <- tandemRepeatsCNV[tandemRepeatsCNV$name!="20-somy" & tandemRepeatsCNV$name!="19-somy" &tandemRepeatsCNV$name!="18-somy" &tandemRepeatsCNV$name!="17-somy" & tandemRepeatsCNV$name!="16-somy" &tandemRepeatsCNV$name!="15-somy" &tandemRepeatsCNV$name!="14-somy"&tandemRepeatsCNV$name!="13-somy" &tandemRepeatsCNV$name!="12-somy"&tandemRepeatsCNV$name!="11-somy"&tandemRepeatsCNV$name!="10-somy"&tandemRepeatsCNV$name!="zero-inflation",] 
removeHighPloidy2$name <-droplevels(removeHighPloidy2$name)
removeHighPloidy2$itemRgb <-as.factor(removeHighPloidy2$itemRgb)
removeHighPloidy2$itemRgb <-droplevels(removeHighPloidy2$itemRgb)



######################################################################
                    ####PLOTTING####
######################################################################
ggplot(tandemRepeatsCNV)+geom_bar(aes(ploidy))


ggplot()+geom_density(data = filter(tandemRepeats,ID=="BLM"),aes(x = TD,color=ID))+
  geom_density(data = filter(tandemRepeats,ID=="BLM/RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeats,ID=="RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeats,ID=="WT"),aes(x = TD,color=ID))+ 
  theme_bw() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()

cnv$seqnames <- factor(cnv$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                               "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
ggplot(cnv)+geom_bar(aes(seqnames),fill="#4a826b")+
  theme_classic()+facet_wrap(~gene)+ggsave("Plots/CNV_per_chr.png")

tandemRepeatsCNV$seqnames <- factor(tandemRepeatsCNV$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                               "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
ggplot(tandemRepeatsCNV)+geom_bar(aes(seqnames),fill="#4a826b")+
  theme_classic()+facet_wrap(~gene)+ggsave("Plots/CNV_per_chr.png")


ggplot(tandemRepeats)+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library.png")

ggplot(losses)+geom_density(mapping = aes(x = DEL,color=ID,fill=ID),alpha=0.3)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()

ggplot(tandemRepeats)+geom_density(mapping = aes(x = TD,fill=ID),alpha=0.7)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  facet_wrap(~ID)+
  scale_x_log10()+ ggsave("Plots/Aneuploidy-count-library-facet.png")

ggplot(losses)+geom_density(mapping = aes(x = DEL,fill=ID),alpha=0.7)+
  theme_classic() +labs(title = "Aneuploidy count/library")+
  facet_wrap(~ID)+
  scale_x_log10()

ggplot()+geom_density(data = tandemRepeats,aes(x = TD))+
  theme_bw() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()

ggplot(tandemRepeats)+geom_density(aes(x = td_sum,color=ID,group=ID))+
  theme_bw() +
  labs(title = "Aneuploidy count/library")+ ylim(c(0,15))+
  scale_x_log10() + ggsave("Plots/aneuploidy_amount.png")  

ggplot()+geom_density(data = cnv,aes(x = width,color=gene),size=1.5)+
  theme_bw() +labs(title = "CNV segment size")+
  scale_x_log10(breaks = c(1e+5,1e+6,1e+7,1e+8),
                limits = c(50000, 1e+8)) + 
  ggsave("Plots/CNV_segment_size.png")

ggplot(tandemRepeats,aes(ID,TD))+geom_violin(aes(fill=ID))+
  geom_boxplot(width=0.05)+scale_y_log10()+ theme_bw()+
  labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
  ggsave("Plots/CNV_per_cell.png")

ggplot(removeHighPloidy)+geom_bar(aes(gene,fill=name))+
  scale_fill_manual(values =  c("0-somy" ="#E5E5E5","1-somy"="#9A32CD" , "2-somy"="#00EE76" , "3-somy"="#CD0000",
                                "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                "8-somy" ="#458B00","7-somy"="#1E90FF",
                                "6-somy"="#FFFACD")) + theme_bw()+
  ggsave("Plots/spreadOfPloidyAcrossGenes.png")

ggplot(removeHighPloidy2)+geom_bar(aes(gene,fill=name))+
  scale_fill_manual(values =  c("0-somy" ="#E5E5E5","1-somy"="#9A32CD" , "2-somy"="#00EE76" , "3-somy"="#CD0000",
                                "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                "8-somy" ="#458B00","7-somy"="#1E90FF",
                                "6-somy"="#FFFACD")) + theme_bw()+
  labs(x="DNA Repair Deficiency",y="Number of TDs")+
  ggsave("Plots/spreadOfPloidy_TDsOnly.png")

ggplot(losses)+geom_bar(aes(gene,fill=name))+
   theme_bw()
  
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




