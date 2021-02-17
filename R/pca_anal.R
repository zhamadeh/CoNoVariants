
test <- get(load("SERVER-OUTPUT/method-HMM/blm_recq5-1-single-c14-w19_.trimmed.mdup.bam_binsize_1e+05_stepsize_1e+05_reads.per.bin_3.94_.RData"))
test2<- get(load("/Library/Frameworks/R.framework/Versions/3.6/Resources/library/AneuFinderData/extdata/primary-lung/hmms/AvdB150303_I_096.bam"))

class(test)
str(test2)
library(AneuFinder)
plot_pca(list.files(path = "SERVER-OUTPUT/method-HMM/",full.names = T))


lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
lung.files <- list.files(lung.folder, full.names=TRUE)
plot_pca(lung.files)
files=list.files(path = "SERVER-OUTPUT/method-HMM/",full.names = T)
for (file in files){
  tmp <- get(load(file))
  class(tmp)="aneuHMM"
  save(tmp,file=file)
}
plotHeterogeneity(hmms.list  = list.files(path = "SERVER-OUTPUT/method-HMM/",full.names = T))

cl <- get(load("SERVER-OUTPUT/PLOTS-1e5-goodLibs/clusterByQuality.RData"))
plot(cl$Mclust,what="classification")
legend(legend=c("1", "2","3","4","5","6","7","8","9"))

levels(as.factor(cl$Mclust$classification))
