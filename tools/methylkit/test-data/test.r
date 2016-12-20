library("methylKit")

file.list = list("input_test1.myCpG.txt", "input_test2.myCpG.txt", "input_control1.myCpG.txt", "input_control2.myCpG.txt")

myobj=methRead( file.list,
                sample.id=list("test 1","test 2","control 1","control 2"),assembly="hg18",pipeline="amp",treatment=c(1,1,0,0))

pdf('output_statistics.pdf')
for (obj in myobj){
  getMethylationStats(obj,plot=TRUE,both.strands=FALSE)
  getCoverageStats(obj,plot=TRUE,both.strands=FALSE)
}
devname = dev.off()

# unite function
methidh = unite(myobj)

pdf("output_correlation.pdf")
getCorrelation(object = methidh, plot=TRUE, method = "pearson")
devname = dev.off()

# differential methylation
myDiff = calculateDiffMeth(methidh, overdispersion="none",
                           adjust="SLIM", effect="wmean", test="Chisq",
                           slim=FALSE, weighted.mean=FALSE)

bedgraph(myDiff, file.name="output_myDiff.bedgraph", col.name="meth.diff",
         unmeth=FALSE, log.transform=FALSE, negative=FALSE, add.on="")

MethPerChr = diffMethPerChr(myDiff, plot=FALSE,
                            qvalue.cutoff=0.01,
                            meth.cutoff=25)
write.table(MethPerChr, sep="\t", row.names=FALSE, quote=FALSE, file="output_MethPerChr.tsv")

MethylDiff = getMethylDiff(myDiff, difference=25,
                               qvalue=0.01, type="all")
bedgraph(MethylDiff, file.name="output_MethylDiff.bedgraph", col.name="meth.diff",
         unmeth=FALSE,log.transform=FALSE,negative=FALSE,add.on="")

pdf( "output_clustering.pdf" )
methClust = clusterSamples(methidh, dist="correlation", method="ward")
devname = dev.off()

pdf( "output_PCA.pdf" )
PCASamples(methidh)
devname = dev.off()

## methSeg works for methylRaw or methylDiff with resolution region,
## so methylBase has to be tiled before
tileRaw = tileMethylCounts(myobj[[1]])
tileBase = tileMethylCounts(methidh)
tileDiff = calculateDiffMeth(tileBase)

## methseg generates Granges
segRaw = methSeg(tileRaw, diagnostic.plot = FALSE)
segDiff = methSeg(tileDiff, diagnostic.plot = FALSE)

## and can be exported as BED
methSeg2bed(segments = segRaw, filename = "output_seg_raw.bed")
methSeg2bed(segments = segDiff, filename = "output_seg_diff.bed")
