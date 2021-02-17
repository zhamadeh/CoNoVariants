######################################################################
                    ####CNV analysis of good libraries (n=317)####
######################################################################
#packages
library(tidyverse)
library(rtracklayer)
#data
CNV<- import("SERVER-OUTPUT/BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_CNV.bed")
CNV<- import("SERVER-OUTPUT/BROWSERFILES/method-HMM/binsize_1e+05_stepsize_1e+05_StrandSeq_breakpoint-hotspots.bed")

cnv <- GRanges(as.data.frame(CNV))


#core function
countTDs <- function(CNV){
  
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
    ploidyTable <- tmp %>% group_by(name) %>% summarize(n())
    ploidy = as.numeric(strsplit(as.character(ploidyTable[which.max(ploidyTable$`n()`),]$name),split = "[-]")[[1]][1])
    if (ploidy==0){
      ploidy=1
    }
    remove=c()
    
    for (row in 1:nrow(ploidyTable)){
      print(row)
      print(ploidyTable[row,1]$name)
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
      td <- sum(ploidyTable[-c(remove),]$`n()`)
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
      td_count <- sum(tmp[-c(remove),]$width)
    } else {td_count <- sum(tmp[-c(remove),]$width)}

    row <- data.frame(ID=id,TD=td,ploidy=ploidy,td_sum=td_count)
    tandemRepeats <- rbind(row,tandemRepeats)
  }
}







######################################################################
                    ####PLOTTING####
######################################################################
ggplot()+geom_density(data = filter(tandemRepeats,ID=="BLM"),aes(x = TD,color=ID))+
  geom_density(data = filter(tandemRepeats,ID=="BLM/RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeats,ID=="RECQL5"),aes(x = TD,color=ID))+ 
  geom_density(data = filter(tandemRepeats,ID=="WT"),aes(x = TD,color=ID))+ 
  theme_bw() +labs(title = "Aneuploidy count/library")+
  scale_x_log10()


ggplot(tandemRepeats)+geom_density(aes(x = td_sum,color=ID,group=ID))+
  theme_bw() +
  labs(title = "Aneuploidy count/library")+
  scale_x_log10()+ylim(c(0,15))


tandemRepeats%>% group_by(ID)%>%summarize(n())

