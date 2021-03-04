library(karyoploteR)
regions <- createRandomRegions(nregions=10000, length.mean = 1e6, mask=NA, non.overlapping = FALSE)
kp <- plotKaryotype()
kpPlotDensity(kp, data=regions)


hmm <- read.table("DataSummary/HMM/tandemRepeatsTotal_filtHighPloidy_hmm.txt",header=T,fill=T,comment.char="")
HMM <- GRanges(hmm)
kp <- plotKaryotype()
kpPlotDensity(kp, data=HMM)

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
