# get max intron size

library(GenomicFeatures)
library(rtracklayer)


# GET THE GFF FILE
gtf <- makeTxDbFromGFF("~/genomes/waterhemp/male/WH_contig_annotation/Male_WH_contig_functional_genes.gff3") #change me!


#get intron length per transcript

introns <- intronsByTranscript(gtf)

# get summary
summary(sum(width(introns)))
