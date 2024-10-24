setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/")

library("karyoploteR")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

# ===== Import bigwig file =====================================================
ATAC.WT <- "./ATAC/data/BigWig/SRR8361537.bw"
ATAC.KO <- "./ATAC/data/BigWig/SRR8361539.bw"
RNA.WT <- "./RNA/data/bigwig/SRR8361594.bw"
RNA.KO <- "./RNA/data/bigwig/SRR8361596.bw"

# ===== Visualize bigwig file ==================================================
Trop2.region <- toGRanges("chr1:58,575,065-58,579,300")

pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 1
pp$data1inmargin <- 5
pp$topmargin <- 10
pp$bottommargin <- 10
pp$data1height <- 200
pp$rightmargin <- 0.1
pp$leftmargin <- 0.25
pp$data1outmargin <- 5

pdf("./ATAC/figure/BigWig.pdf", width = 15, height = 12)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                    karyoplot = kp,
                                    plot.transcripts = TRUE,
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)

kp <- plotKaryotype(zoom = Trop2.region, cex = 2, plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 100, add.units = TRUE,
                 cex = 1, digits = 12)
kpPlotGenes(kp, data = genes.data, r0=0.0210, r1=0.0725, gene.name.cex = 2, 
            gene.name.position = "right")
kpAddLabels(kp, labels = "RNA-seq", r0=0.10, r1=0.525, cex=2.5, srt=90, 
            pos = 1, label.margin = 0.11)
kpAddLabels(kp, labels = "ATAC-seq", r0=0.575, r1=1, cex=2.5, srt=90, 
            pos = 1, label.margin = 0.11)

### RNA ARID1A KO
kpPlotBigWig(kp, data = RNA.KO, ymax = 2700, r0=0.10, r1=0.2875, col = "pink")
kpAxis(kp, ymin = 0, ymax = 2700, numticks = 2, r0=0.10, r1=0.2875, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A KO", r0=0.2675, r1=0.2875, cex=2, label.margin = -0.63)
### RNA ARID1A WT
kpPlotBigWig(kp, data = RNA.WT, ymax = 2700, r0=0.3375, r1=0.525, col = "steelblue")
kpAxis(kp, ymin = 0, ymax = 2700, numticks = 2, r0=0.3375, r1=0.525, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A WT", r0=0.505, r1=0.525, cex=2, label.margin = -0.63)
### ATAC ARID1A KO
kpPlotBigWig(kp, data = ATAC.KO, ymax = 2.5, r0=0.575, r1=0.7625, col = "pink")
kpAxis(kp, ymin = 0, ymax = 2.5, numticks = 2, r0=0.575, r1=0.7625, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A KO", r0=0.7425, r1=0.7625, cex=2, label.margin = -0.63)
### ATAC ARID1A WT
kpPlotBigWig(kp, data = ATAC.WT, ymax = 2.5, r0=0.8125, r1=1, col = "steelblue")
kpAxis(kp, ymin = 0, ymax = 2.5, numticks = 2, r0=0.8125, r1=1, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A WT", r0=0.98, r1=1, cex=2, label.margin = -0.63)

dev.off()
