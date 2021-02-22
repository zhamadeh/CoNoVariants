bs_hmm <- import("SERVER-OUTPUT/Output-Niek_BS/BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")
bs_dnacopy <- import("SERVER-OUTPUT/Output-Niek_BS/BROWSERFILES/method-dnacopy/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed.gz")

BS_hmm <- as.data.frame(bs_hmm)
BS_hmm <- separate(BS_hmm ,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
BS_hmm <- separate(BS_hmm ,file,c("file","d","a"),"[.]") %>% select(-c(a,d))

BS_dnacopy <- as.data.frame(bs_dnacopy)
BS_dnacopy<- separate(BS_dnacopy,group_name,c("a","b","c","file")," ") %>% select(-c(a,b,c))
BS_dnacopy<- separate(BS_dnacopy,file,c("file","d","a"),"[.]") %>% select(-c(a,d))


metadata <- read.table("../Alignment Pipeline/Anxilliary/Metadata/E-MTAB-5976.sdrf.txt",sep="\t",header=T)
metadata <-  select(metadata,c(Comment.ENA_RUN., Extract.Name,Factor.Value.protocol., Characteristics.genotype., Characteristics.organism.))

bs_meta_hmm <- merge(BS_hmm,metadata,by.x="file",by.y="Comment.ENA_RUN.")
bs_meta_hmm$Characteristics.genotype.<- as.character(bs_meta_hmm$Characteristics.genotype.)
bs_meta_hmm$Characteristics.genotype.[bs_meta_hmm$Characteristics.genotype.==as.factor("homozygous Blm knockout")] <- "BLM"
bs_meta_hmm$Characteristics.genotype.[bs_meta_hmm$Characteristics.genotype.==as.factor("wild type genotype")] <- "WT"
bs_meta_hmm <- bs_meta_hmm %>% select(-c(group,strand,score,itemRgb,thick.start,thick.end,thick.width, Factor.Value.protocol., Characteristics.organism.))
bs_meta_hmm <- bs_meta_hmm %>% select(-c(Extract.Name))
bs_meta_hmm$file <- as.factor(bs_meta_hmm$file)
bs_meta_hmm$name<- as.factor(bs_meta_hmm$name)

bs_meta_dnacopy <- merge(BS_dnacopy,metadata,by.x="file",by.y="Comment.ENA_RUN.")
bs_meta_dnacopy$Characteristics.genotype.<- as.character(bs_meta_dnacopy$Characteristics.genotype.)
bs_meta_dnacopy$Characteristics.genotype.[bs_meta_dnacopy$Characteristics.genotype.==as.factor("homozygous Blm knockout")] <- "BLM"
bs_meta_dnacopy$Characteristics.genotype.[bs_meta_dnacopy$Characteristics.genotype.==as.factor("wild type genotype")] <- "WT"
bs_meta_dnacopy <- bs_meta_dnacopy %>% select(-c(group,strand,score,itemRgb,thick.start,thick.end,thick.width, Factor.Value.protocol., Characteristics.organism.))
bs_meta_dnacopy <- bs_meta_dnacopy %>% select(-c(Extract.Name))
bs_meta_dnacopy$file <- as.factor(bs_meta_dnacopy$file)
bs_meta_dnacopy$name<- as.factor(bs_meta_dnacopy$name)


i=300



countTDs <- function(CNV){
  cnv_tds<-data.frame()
  tandemRepeats=data.frame()
  
  for (i in 1:length(CNV)){
    
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[.]")[[1]][1]
    id  <- as.character(metadata$Characteristics.genotype.[metadata$Comment.ENA_RUN.==file])
    
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
    print(i)
    ID <- strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[-_.]")[[1]]
    file=strsplit(strsplit(CNV[i]@listData[[1]]@metadata$trackLine@description,split = " ")[[1]][4],split = "[.]")[[1]][1]
    id  <- as.character(metadata$Characteristics.genotype.[metadata$Comment.ENA_RUN.==file])
    
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
      message("tmp row:",row)
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
    } else {td_count <- sum(tmp$width)}
    print("*")
    if (length(remove)>0){
      tmp <- as.data.frame(tmp[-c(remove),])
      tmp<-tmp %>% select(-c(group,group_name))
      if (nrow(tmp)>0){
        tmp$gene <- id
        tmp$file = as.factor(file)
        cnv_tds <- rbind(tmp,cnv_tds)
      }
    }
  }
  return(cnv_tds)
}


#Run above functions 
tandemRepeatsPerCell_dnacopy <- countTDs(bs_dnacopy)
tandemRepeatsTotal_dnacopy<- countTDsInCNV(bs_dnacopy)

tandemRepeatsPerCell_hmm <- countTDs(bs_hmm)
tandemRepeatsTotal_hmm<- countTDsInCNV(bs_hmm)



tandemRepeatsTotal_filtHighPloidy_dnacopy <- tandemRepeatsTotal_dnacopy[tandemRepeatsTotal_dnacopy$name!="20-somy" & tandemRepeatsTotal_dnacopy$name!="19-somy" &tandemRepeatsTotal_dnacopy$name!="18-somy" &tandemRepeatsTotal_dnacopy$name!="17-somy" & tandemRepeatsTotal_dnacopy$name!="16-somy" &tandemRepeatsTotal_dnacopy$name!="15-somy" &tandemRepeatsTotal_dnacopy$name!="14-somy"&tandemRepeatsTotal_dnacopy$name!="13-somy" &tandemRepeatsTotal_dnacopy$name!="12-somy"&tandemRepeatsTotal_dnacopy$name!="11-somy"&tandemRepeatsTotal_dnacopy$name!="10-somy"&tandemRepeatsTotal_dnacopy$name!="zero-inflation"&tandemRepeatsTotal_dnacopy$name!="21-somy"&tandemRepeatsTotal_dnacopy$name!="34-somy" &tandemRepeatsTotal_dnacopy$name!="22-somy"&tandemRepeatsTotal_dnacopy$name!="24-somy",] 
tandemRepeatsTotal_filtHighPloidy_dnacopy$name <-droplevels(tandemRepeatsTotal_filtHighPloidy_dnacopy$name)
tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb)
tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy_dnacopy$itemRgb)


tandemRepeatsTotal_filtHighPloidy_hmm <- tandemRepeatsTotal_hmm[tandemRepeatsTotal_hmm$name!="20-somy" & tandemRepeatsTotal_hmm$name!="19-somy" &tandemRepeatsTotal_hmm$name!="18-somy" &tandemRepeatsTotal_hmm$name!="17-somy" & tandemRepeatsTotal_hmm$name!="16-somy" &tandemRepeatsTotal_hmm$name!="15-somy" &tandemRepeatsTotal_hmm$name!="14-somy"&tandemRepeatsTotal_hmm$name!="13-somy" &tandemRepeatsTotal_hmm$name!="12-somy"&tandemRepeatsTotal_hmm$name!="11-somy"&tandemRepeatsTotal_hmm$name!="10-somy"&tandemRepeatsTotal_hmm$name!="zero-inflation",] 
tandemRepeatsTotal_filtHighPloidy_hmm$name <-droplevels(tandemRepeatsTotal_filtHighPloidy_hmm$name)
tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb <-as.factor(tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb)
tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb <-droplevels(tandemRepeatsTotal_filtHighPloidy_hmm$itemRgb)



tandemRepeatsPerCell_dnacopy$name="niek_dnacopy"
tandemRepeatsPerCell_hmm$name="niek_hmm"

plotting <- function(tandemRepeatsPerCell,tandemRepeatsTotal_filtHighPloidy){
  
  method = tandemRepeatsPerCell$name[1]
  plotsPath=paste0("Plots/",method,"/")

  
  
  ggplot()+geom_density(tandemRepeatsPerCell,mapping=aes(x = TD,color=ID))+
    theme_bw() +labs(title = "Aneuploidy count/library")+
    scale_x_log10()+
    ggsave(paste0(plotsPath,"2.png"))
  
  tandemRepeatsTotal_filtHighPloidy$seqnames <- factor(tandemRepeatsTotal_filtHighPloidy$seqnames,levels = c("chr1" ,"chr2", "chr3" , "chr4" , "chr5" , "chr6" , "chr7",  "chr8" , "chr9" , "chr10" ,
                                                                                                             "chr11" ,"chr12", "chr13" ,"chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20" ,"chr21" ,"chr22",  "chrX" , "chrY" ))
  ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(seqnames),fill="#4a826b")+
    theme_classic()+facet_wrap(~gene)+ggsave(paste0(plotsPath,"3.png"))
  
  ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,color=ID,fill=ID),alpha=0.4)+
    theme_classic() +labs(title = "Aneuploidy count/library")+
    scale_x_log10()+ ggsave(paste0(plotsPath,"4.png"))
  
  ggplot(tandemRepeatsPerCell)+geom_density(mapping = aes(x = TD,fill=ID),alpha=0.7)+
    theme_classic() +labs(title = "Aneuploidy count/library")+
    facet_wrap(~ID)+
    scale_x_log10()+ ggsave(paste0(plotsPath,"5.png"))
  
  ggplot()+geom_density(data = tandemRepeatsTotal_filtHighPloidy,aes(x = width,color=gene),size=1.5)+
    theme_bw() +labs(title = "CNV segment size")+
    scale_x_log10(breaks = c(1e+5,1e+6,1e+7,1e+8)) + 
    ggsave(paste0(plotsPath,"6.png"))
  
  ggplot(tandemRepeatsPerCell,aes(ID,TD))+geom_violin(aes(fill=ID))+
    geom_boxplot(width=0.05)+scale_y_log10()+ theme_bw()+
    labs(x="DNA Repair Deficiency",y="Number of CNV segments/cell")+
    ggsave(paste0(plotsPath,"7.png"))
  
  ggplot(tandemRepeatsTotal_filtHighPloidy)+geom_bar(aes(gene,fill=name))+
    scale_fill_manual(values =  c("2-somy"="#00EE76" , "3-somy"="#CD0000",
                                  "4-somy"="#EEC900" , "5-somy" ="#000080" ,"9-somy" ="#458B00",
                                  "8-somy" ="#458B00","7-somy"="#1E90FF",
                                  "6-somy"="#FFFACD")) + theme_bw()+
    ggsave(paste0(plotsPath,"8.png"))
}

plotting(tandemRepeatsPerCell_dnacopy,tandemRepeatsTotal_filtHighPloidy_dnacopy)


