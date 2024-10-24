setwd("/Users/nakyung/Desktop/Bioinformatics/Gilead/ATAC")

library("karyoploteR")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

# ===== Import bigwig file =====================================================
ATAC.WT <- "./data/BigWig/SRR8361537.bw"
ATAC.KO <- "./data/BigWig/SRR8361539.bw"

# ===== Visualize motif and bigwig file ========================================
peak.region <- toGRanges("chr1:58,577,000-58,577,700")

## Set plotting parameters
pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 1
pp$data1inmargin <- 5
pp$topmargin <- 10
pp$bottommargin <- 10
pp$data1height <- 200
pp$rightmargin <- 0.1
pp$leftmargin <- 0.25
pp$data1outmargin <- 5

## Get gene information
gene.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                   karyoplot = kp, plot.transcripts = TRUE,
                                   plot.transcripts.structure = TRUE)
gene.data <- addGeneNames(gene.data)
gene.data <- mergeTranscripts(gene.data)

## Plot the bigwig and region
pdf("./figure/Motif_bigwig.pdf", width = 13.5, height = 7.5)
kp <- plotKaryotype(zoom = peak.region, cex = 2, plot.params = pp)
kpPlotGenes(kp, data = gene.data, r0=0.02, r1=0.12, gene.name.position = "right",
            gene.name.cex = 2)
kpAddLabels(kp, labels = "ATAC-seq", r0=0.475, r1=0.9, cex=2, srt=90, 
            pos = 1, label.margin = 0.18, font = 2)
kpText(kp, chr = "chr1", x = 58577555, y = 0.88, labels = "chr1:58,577,000-58,577,700",
       col = "black", cex = 1.5)

### ATAC ARID1A KO ###
kpPlotBigWig(kp, data = ATAC.KO, ymax = 2.5, r0=0.475, r1=0.6625, col = "pink")
kpAxis(kp, ymin = 0, ymax = 2.5, numticks = 2, r0=0.475, r1=0.6625, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A KO", r0=0.475, r1=0.6625, cex=1.5, label.margin = 0.035)
### ATAC ARID1A WT ###
kpPlotBigWig(kp, data = ATAC.WT, ymax = 2.5, r0=0.7125, r1=0.9, col = "steelblue")
kpAxis(kp, ymin = 0, ymax = 2.5, numticks = 2, r0=0.7125, r1=0.9, cex = 1.5)
kpAddLabels(kp, labels = "ARID1A WT", r0=0.7125, r1=0.9, cex=1.5, label.margin = 0.035)

### KLF5 motif site ###
kpPlotRegions(kp, data = "chr1:58,577,351-58,577,360", col = "#FFEECC", border = "#FFCCAA", 
              r0 = 0.43, r1 = 0.46)
kpPlotRegions(kp, data = "chr1:58,577,344-58,577353", col = "#FFEECC", border = "#FFCCAA", 
              r0 = 0.39, r1 = 0.42)
kpPlotRegions(kp, data = "chr1:58,577,159-58,577,168", col = "#FFEECC", border = "#FFCCAA", 
              r0 = 0.41, r1 = 0.44)
kpAddLabels(kp, labels = "Motif binding site", r0=0.35, r1=0.46, cex=1.5, label.margin = 0.035)
kpAddLabels(kp, labels = "Motif sequence", r0=0.15, r1=0.35, cex=1.5, label.margin = 0.035)

### chr1:58577159-58577168
kpText(kp, chr = "chr1", x = 58577140, y = 0.31, labels = "CC  CC   CCCC", col = "blue", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577140, y = 0.31, labels = "  G", col = "orange", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577079, y = 0.31, labels = "T", col = "#009900", cex = 2.5, font = 2)
kpSegments(kp, chr = "chr1", x0=58577045, x1=58577159, y0=0.35, y1=0.41, col = "darkgrey")
kpSegments(kp, chr = "chr1", x0=58577168, x1=58577248, y0=0.41, y1=0.35, col = "darkgrey")

### chr1:58577344-58577353
kpText(kp, chr = "chr1", x = 58577280, y = 0.2, labels = "A", col = "red", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577300, y = 0.2, labels = "T", col = "#009900", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577345, y = 0.2, labels = "CCC", col = "blue", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577395, y = 0.2, labels = "G", col = "orange", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577455, y = 0.2, labels = "CCCC", col = "blue", cex = 2.5, font = 2)
kpSegments(kp, chr = "chr1", x0=58577280, x1=58577344, y0=0.24, y1=0.39, col = "darkgrey")
kpSegments(kp, chr = "chr1", x0=58577353, x1=58577483, y0=0.39, y1=0.24, col = "darkgrey")

### chr1:58577351-58577360
kpText(kp, chr = "chr1", x = 58577540, y = 0.35, labels = "CCCCC  CCCC", col = "blue", cex = 2.5, font = 2)
kpText(kp, chr = "chr1", x = 58577553, y = 0.35, labels = "T", col = "#009900", cex = 2.5, font = 2)
kpSegments(kp, chr = "chr1", x0=58577351, x1=58577410, y0=0.43, y1=0.38, col = "darkgrey")
kpSegments(kp, chr = "chr1", x0=58577360, x1=58577640, y0=0.43, y1=0.39, col = "darkgrey")

dev.off()