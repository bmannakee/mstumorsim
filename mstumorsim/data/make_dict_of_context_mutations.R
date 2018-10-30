library(tidyverse)
library(GenomicRanges)
library(VariantAnnotation)
library(SomaticSignatures)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(deconstructSigs)
library(rtracklayer)
library(liftOver)



############################################################################################################################
# The TCGA files used here are described in TCGA_manifest.txt
# and were downloaded from TCGA in August 2018
#
# The Pan-cancer analysis of whole genomes (pcawg) were downloaded from
# https://www.synapse.org/#!Synapse:syn11801870 which is described as
# "This is supplementary data and results from the ICGC Pan Cancer Analysis Mutational Signatures Working Group (PCAWG7)"
# These simple vcf files were downloaded August 2018
#
# This script along with the above files generates "mutations_with_contexts_GRCh38.tsv" containing 29323166 variants along
# with their context
############################################################################################################################

# MAFS from TCGA are hg38 coordinates. We will be working in GRCh38, so no need to liftover.
mafs <- fs::dir_ls('mafs',glob = '*.gz', recursive=T) %>%
  map_df(read_tsv, .id = 'file', col_types = cols(.default = 'c'), comment = "#") %>%
  dplyr::select(-file) %>% dplyr::filter(Variant_Type == 'SNP') %>%
  dplyr::distinct(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele1,Tumor_Seq_Allele2) %>%
  dplyr::select("seqnames" = "Chromosome",
                "start" = "Start_Position",
                "end" = "End_Position",
                "ref" = "Tumor_Seq_Allele1",
                "alt" = "Tumor_Seq_Allele2") %>%
  dplyr::filter(ref != alt & start == end & ref != '-' & alt != '-') # some have the same ref and alt, which chokes later on
gr_mafs <- GenomicRanges::makeGRangesFromDataFrame(mafs)
gr_mafs$ref <- mafs$ref
gr_mafs$alt <- mafs$alt

# pcawg simple variant files are in hg19 coordinates, and need to be lifted to hg38 for use
pcawg <- fs::dir_ls('pcawg_variants') %>%
  map_df(read_tsv, .id = 'file', col_types = cols(.default = 'c'), col_names = F) %>%
  dplyr::select("seqnames"="X6",
                "start" = "X7",
                "end" = "X8",
                "ref" = "X9",
                "alt" = "X10") %>%
  dplyr::mutate(seqnames = case_when(seqnames == "23" ~ "X",
                                     seqnames == "24" ~ "Y",
                                     seqnames == "MT" ~ "M",
                                     TRUE ~ seqnames)) %>%
  dplyr::filter(ref != alt & start == end & ref != '-' & alt != '-') # some have the same ref and alt, which chokes later on

gr_pcawg <- GenomicRanges::makeGRangesFromDataFrame(pcawg)
gr_pcawg$ref <- pcawg$ref
gr_pcawg$alt <- pcawg$alt
seqlevelsStyle(gr_pcawg) = "UCSC"
ch <- import.chain("./hg19ToHg38.over.chain")
gr38_pcawg <- liftOver(gr_pcawg,ch)
gr38_pcawg <- unlist(gr38_pcawg)

# Combine the manipulated variants into a single genomic ranges object
gr <- unique(c(gr_mafs,gr38_pcawg))
# We used UCSC chains for liftover, now need to translate back
# to NCBI coordinates so we can use the GRCh38 reference in the rest of the process.
seqlevelsStyle(gr) <- "NCBI"
vr <- VariantAnnotation::makeVRangesFromGRanges(gr)
vr <- vr[isSNV(vr)]
mfr <- SomaticSignatures::mutationContext(vr,ref = BSgenome.Hsapiens.NCBI.GRCh38) # The alignment for the simulated genomes is NCBI GRCh38
tmpfr <- as.data.frame(mfr) # loses contexts
mc <- as.data.frame(mcols(mfr)) # contexts
mc <- mc %>% as_tibble %>% tidyr::separate(context,into=c("before","after"),sep="\\.") %>%
  mutate(cref = stringr::str_sub(alteration,1L,1L),
         calt = stringr::str_sub(alteration,2L))
nrow(tmpfr)==nrow(mc) # They have the same length
context_fr <- dplyr::bind_cols(c(tmpfr,mc)) %>%
   dplyr::mutate(context = paste0(before,"[",cref,">",calt,"]",after),
                 start = as.character(start),
                 end = as.character(end)) %>% # Converted to integer but need character to re-merge
   dplyr::select("chr"="seqnames",start,end,ref,alt,context)

# make dictionary from context to integer
contexts_in_order <- colnames(deconstructSigs::signatures.cosmic)
context_ints <- seq(1,96,1)
names(context_ints) <- contexts_in_order
contexts_enum <- as_tibble(t(context_ints)) %>% gather(context,context_id,everything()) %>%
  mutate(context_id = as.integer(context_id))

final_fr <- context_fr %>% left_join(contexts_enum, by='context')
final_fr %>% write_tsv('mutations_with_contexts_GRCh38.tsv')

############## Label the variants that appear in AgilentV5exome.bed
library(fuzzyjoin) # Gives us the wonderful function genome_inner_join
bed_file <- 'Agilent-ExomeV5.bed'
col_spec <- cols_only(chr = "c", start = "n", end = "n")
bed_fr <- readr::read_tsv(bed_file,col_names = c('chr','start','end','gene'), col_types = col_spec)
# Write out the bed files (interval_list) for each chromosome to use for mutect, there is some purrr magic here.
bed_fr %>%
  split(.$chr) %>% # turns the data frame into a named list, one for each chromosome
  list(names(.)) %>% # pwalk takes two vectors as arguments, this one holds the names after the split
  pwalk(~ write_tsv(.x,paste0('./bed_by_chr/chr',.y,'.interval_list')))

# Do a conditional inner join to make a data frame where only the variants within exomes are included
final_fr <- final_fr %>% dplyr::mutate(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end))
exome_fr <- final_fr %>%
  genome_inner_join(bed_fr, by = c('chr','start','end')) %>% # 2,422,321 variants
  dplyr::select(-c(chr.y,start.y,end.y))
exome_fr %>% readr::write_tsv('AgilentV5-Exome_mutations_with_contexts_GRCh38.tsv')

