library(karyoploteR)
regions <- createRandomRegions(nregions=10000, length.mean = 1e6, mask=NA, non.overlapping = FALSE)
kp <- plotKaryotype()
kpPlotDensity(kp, data=regions)


hmm <- read.table("DATA/HMM/TABLES/TD_total_filtHighPloidy.txt",header=T,fill=T,comment.char="")
HMM <- GRanges(hmm)
kp <- plotKaryotype()
p <-kpPlotDensity(kp, data=HMM)
save(file = "dist.png")

kp <- plotKaryotype(plot.type=2, chromosomes = "chr1")
kpPlotDensity(kp, data=HMM)
kpPlotRegions(kp, data=HMM, data.panel=2)

hmm <- read.table("DataSummary/DNACOPY/tandemRepeatsTotal_filtHighPloidy_dnacopy.txt",header=T,fill=T,comment.char="")
HMM <- GRanges(hmm)
filtHMM =HMM[seqnames(HMM)!="chr8" & seqnames(HMM)!="chr15" & seqnames(HMM)!="chrX"]
#HMM[seqnames(HMM)=="chr8" | seqnames(HMM)=="chr15" | seqnames(HMM)=="chrX"]
kp <- plotKaryotype()
kpPlotDensity(kp, data=filtHMM)

kp <- plotKaryotype(plot.type=2, chromosomes = "chr1")
kpPlotDensity(kp, data=HMM)
kpPlotRegions(kp, data=HMM, data.panel=2)
