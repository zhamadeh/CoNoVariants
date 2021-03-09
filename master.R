
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicRanges))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(GenomicAlignments))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(DNAcopy))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(doParallel))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(tiydyverse))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ggplot2))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(cowplot))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(caTools))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(gplots))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(ReorderCluster))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(BiocManager))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(BiocManager::install("AneuFinder"))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(AneuFinder))))
suppressMessages(suppressPackageStartupMessages(suppressWarnings(require(BSgenome.Hsapiens.UCSC.hg19))))

Aneufinder(inputfolder='../Sequencing_data/WT_composite/', outputfolder='aneufinder_output_compWT_only/',
           numCPU=5, method=c('dnacopy','HMM',"edivisive"),
           configfile=NULL, reuse.existing.files=TRUE, binsizes=1e5, stepsizes=1e5, #variable.width.reference="../Sequencing_data/WT_forCompositeRef/comp.bam", 
           reads.per.bin=NULL, pairedEndReads=TRUE,remove.duplicate.reads=TRUE, min.mapq=10, 
           use.bamsignals=FALSE, reads.store=FALSE, 
           correction.method="GC", GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19, 
           strandseq=TRUE, 
           R=10, sig.lvl=0.1, eps=0.01, max.time=60, max.iter=5000, num.trials=15, states=c('zero-inflation',paste0(0:10,'-somy')), 
           confint=NULL, refine.breakpoints=FALSE, hotspot.bandwidth=NULL, hotspot.pval=5e-2, cluster.plots=TRUE)

