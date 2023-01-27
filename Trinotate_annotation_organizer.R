#!/usr/bin/env Rscript

### Organizer for Trinity annotation annotation #####

# if package not installed
list.of.packages <- c("devtools", "readr","optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

devtools::install_github("cstubben/trinotateR", quiet = T)

# Packages
library(optparse)
library(trinotateR)
library(readr)

# set arguments
option_list = list(
  make_option(c("-r", "--report"), type="character", default=NULL, 
              help="Trinotate xls report", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# if no argument is specified
if (is.null(opt$report)){
  print_help(opt_parser)
  stop("At least one argument must be supplied trinote_report.xls", call.=FALSE)
}


# load annotation report
x <- read_trinotate(opt$report)

# get summary of annotation
cat("Trinotate report summary:")
sumtrino <- summary_trinotate(x)
sumtrino
write.table(sumtrino, "summary_trinotate.txt", quote = F, sep = "\t")


# pfam
cat("Organizing Pfam report:")
pfam_annotation <- split_pfam(x)
cat("Saving Pfam report:")
write_csv(pfam_annotation, "pfam_annotation.csv")
cat("Summarizing Pfam report:")
summ_pfam <- summary_pfam(pfam_annotation)
head(summ_pfam)
cat("Saving Pfam summary:")
write_csv(summ_pfam, "summary_Pfam.csv")

# BlastX
cat("Saving BlastX report:")
blastX_annotation <- split_blast(x, hit = "sprot_Top_BLASTX_hit")
cat("Saving BlastX annotation report:")
write_csv(blastX_annotation, "blastX_annotation.csv")
cat("Summarizing BlastX annotation report:")
summ_blastX <-summary_blast(blastX_annotation)
head(summ_blastX)
write_csv(summ_blastX, "summary_blastX.csv")

# BlastP
cat("Saving BlastP report:")
blastP_annotation <- split_blast(x, hit = "sprot_Top_BLASTP_hit")

cat("Saving BlastP annotation report:")
write_csv(blastP_annotation, "blastP_annotation.csv")
cat("Summarizing BlastP annotation report:")
summ_blastp <- summary_blast(blastP_annotation)
head(summ_blastp)
write_csv(summ_blastp, "summary_blastP.csv")

# GOterm
cat("Saving GOterms based on BlastX, BlastP and Pfam: ")
GO_blastX <- split_GO(x, hit = "gene_ontology_BLASTX")
GO_blastP <- split_GO(x, hit = "gene_ontology_BLASTP")
GO_pfam <- split_GO(x, hit = "gene_ontology_Pfam")

write_csv(GO_blastX, "GOterm_blastX.csv")
write_csv(GO_blastP, "GOterm_blastP.csv")
write_csv(GO_pfam, "GOterm_Pfam.csv")

summ_GO_blastX <- summary_GO(GO_blastX)
summ_GO_blastP <- summary_GO(GO_blastP)
summ_GO_Pfam <- summary_GO(GO_pfam)

write_csv(summ_GO_blastX, "summary_GOterm_blastX.csv")
write_csv(summ_GO_blastP, "summary_GOterm_blastP.csv")
write_csv(summ_GO_blastP, "summary_GOterm_Pfam.csv")

# EggNog

data(cogs)
download.file("http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz", "NOG.annotations.tsv.gz")
system("gunzip NOG.annotations.tsv.gz")
egg <- read.table("NOG.annotations.tsv", sep="\t", stringsAsFactors=FALSE, quote="")
names(egg) <- c("db", "nog", "proteins", "species", "class", "description")

pdf(file = "eggnog_plot.pdf", width = 12, height = 9)
plot_NOGs(x)
dev.off()
