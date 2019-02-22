library(ChIPpeakAnno)
library(GenomicFeatures)
library(ChIPseeker)
library(GenomicRanges)
library(stringr)
library(ReactomePA)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)


#get files
list <- list.files(path='.', pattern='.narrowPeak')

#fix names
names(list) <- substr((str_replace_all(list, 'Group', '')),1,12)
   
#convert the list to GRanges
list_gr <- lapply(list, function(x) toGRanges(x, format=c('narrowPeak'), header=F))

#import gff
gff3_file <- ('Mus_musculus.GRCm38.95.gff3')

#make txdb
txdb <- makeTxDbFromGFF(gff3_file, format='gff', dataSource='ENSEMBL', organism='Mus musculus')


#make peak annotation list
peakannolist <- lapply(list, annotatePeak)

#######################################################################################################

#compare all of HO1
all_ol_HO1 <- findOverlapsOfPeaks(list_gr[['A_HO1_A_LH20']], list_gr[['A_HO1_B_H20_']], list_gr[['B_HO1_A_LH20']], list_gr[['B_HO1_B_H20_']])

#compare all of IgG
all_ol_IgG <- findOverlapsOfPeaks(list_gr[['A_IgG_A_LH20']], list_gr[['A_IgG_B_H20_']], list_gr[['B_IgG_A_LH20']], list_gr[['B_IgG_B_H20_']])


#calculate a more permissive estimate of overlaps using only the LH20 library:
permissive_ol_HO1 <- findOverlapsOfPeaks(list_gr[['A_HO1_A_LH20']], list_gr[['B_HO1_A_LH20']])
permissive_ol_IgG <- findOverlapsOfPeaks(list_gr[['A_IgG_A_LH20']], list_gr[['B_IgG_A_LH20']])



#find overlaps in 2 granges objects that represent the peaks present in all HO1 and all IgG chips
all_ol_HO1_peaklist <- all_ol_HO1$peaklist[['A_HO1_A_LH20///A_HO1_B_H20_///B_HO1_A_LH20///B_HO1_B_H20_']]
all_ol_IgG_peaklist <- all_ol_IgG$peaklist[['A_IgG_A_LH20///A_IgG_B_H20_///B_IgG_A_LH20///B_IgG_B_H20_']]
all_ol_HO1_vs_IgG <- findOverlapsOfPeaks(all_ol_HO1_peaklist, all_ol_IgG_peaklist)
all_ol_HO1_vs_IgG_venn_count<-all_ol_HO1_vs_IgG$venn_cnt

#fint the unique peaks...
stringent_HO1_peaks <- all_ol_HO1_vs_IgG$peaklist[['all_ol_HO1_peaklist']]

#calculate a more permissive estimate of overlaps using only the LH20 library:
LH20_HO1_ol <- findOverlapsOfPeaks(A_HO1_A_LH20, B_HO1_A_LH20)
LH20_IgG_ol <- findOverlapsOfPeaks(A_IgG_A_LH20, B_IgG_A_LH20)

#get the intersection peaklists
LH20_HO1_ol_peaklist <- LH20_HO1_ol$peaklist[['A_HO1_A_LH20///B_HO1_A_LH20']]
LH20_IgG_ol_peaklist <- LH20_IgG_ol$peaklist[['A_IgG_A_LH20///B_IgG_A_LH20']]

#find the more permissive overlaps...
permissive_HO1_peaks <- findOverlapsOfPeaks(LH20_HO1_ol_peaklist, LH20_IgG_ol_peaklist)

#######################################################################################################

#pull out just the peaks...
permissive_HO1_peaks_HO1_specific <- permissive_HO1_peaks$peaklist$LH20_HO1_ol_peaklist

#make peak annotation list
permissive_HO1_peaks_anno <- annotatePeak(permissive_HO1_peaks_HO1_specific, TxDb=txdb)

#covert permissive_HO1_peaks_anno to a GRanges object
permissive_HO1_peaks_anno_GRanges <- as.GRanges(permissive_HO1_peaks_anno)

###############

#make some figs
plotAnnoPie(permissive_HO1_peaks_anno)
ggsave('permissive_HO1_peaks_plotAnnoPie.png')

plotDistToTSS(permissive_HO1_peaks_anno,title="Distribution of transcription factor-binding loci relative to TSS")
ggsave('permissive_HO1_peaks_plotDistToTSS.png')

#make a plot of HO1-specific ChIP peaks over chromosomes
covplot(permissive_HO1_peaks_anno_GRanges)
ggsave('permissive_HO1_peaks_covplot.png')

############

#get a list of gene IDs
gene_list=(as.data.frame(permissive_HO1_peaks_anno)$geneId)

#convert gene IDs to entrez format
gene_list_formatted <- bitr(gene_list, fromType='ALIAS', toType='ENTREZID', OrgDb="org.Mm.eg.db")

#enrich...
#pathway1 <- enrichPathway(gene=gene_list_formatted$ENTREZID, organism='mouse')
pathway1 <- enrichPathway(gene=gene_list_formatted$ENTREZID, organism='mouse', pAdjustMethod='bonferroni')
dotplot(pathway1)

ggsave('permissive_HO1_peaks_enrichPathway.png')

#######################################################################################################
#print supplementary outputs
#print an output of the HO1-specific peak loci
permissive_HO1_peaks_HO1_specific_df <- as(permissive_HO1_peaks_anno_GRanges, 'data.frame')
write.table(permissive_HO1_peaks_HO1_specific_df, file='permissive_HO1_peaks.txt', row.names=F)

#######################################################################################################
#Make a plot of the significance value of each peak...

#find OLs betwn A_HO1_A_LH20 / A_HO1_A_LH20 and permissive_HO1_peaks_HO1_specific.. this should give me the 211 pvalues?
HO1_A_specific_Overlaps_log10_pval <- findOverlapsOfPeaks(permissive_HO1_peaks_HO1_specific, A_HO1_A_LH20)
HO1_B_specific_Overlaps_log10_pval <- findOverlapsOfPeaks(permissive_HO1_peaks_HO1_specific, B_HO1_A_LH20)

#pull out just the peaks...
HO1_A_specific_Overlaps_log10_pval_GRanges <- HO1_A_specific_Overlaps_log10_pval$overlappingPeaks$`permissive_HO1_peaks_HO1_specific///A_HO1_A_LH20`
HO1_B_specific_Overlaps_log10_pval_GRanges <- HO1_B_specific_Overlaps_log10_pval$overlappingPeaks$`permissive_HO1_peaks_HO1_specific///B_HO1_A_LH20`

#convert to dataframes...
HO1_A_specific_Overlaps_log10_pval_GRanges_df<- as(HO1_A_specific_Overlaps_log10_pval_GRanges, 'data.frame')
HO1_B_specific_Overlaps_log10_pval_GRanges_df<- as(HO1_B_specific_Overlaps_log10_pval_GRanges, 'data.frame')

#write them...
write.table(HO1_A_specific_Overlaps_log10_pval_GRanges_df, file='HO1_A_specific_Overlaps_log10_pval_GRanges_df.txt', row.names=F, sep='\t')
write.table(HO1_B_specific_Overlaps_log10_pval_GRanges_df, file='HO1_B_specific_Overlaps_log10_pval_GRanges_df.txt', row.names=F, sep='\t')

#######################################################################################################
HO1_B_volplot_dat <- data.frame('pValue'=HO1_B_specific_Overlaps_log10_pval_GRanges_df$pValue, 'signalValue'=HO1_B_specific_Overlaps_log10_pval_GRanges_df$signalValue)
B <- ggplot(data=HO1_B_volplot_dat, aes(HO1_B_volplot_dat$pValue)) + geom_histogram(binwidth=5) + xlab('-log10pValue') + geom_vline(xintercept=3, color='red', size=.25)
B
ggsave('HO1_B_-log10pValue_histogram.png', width=2, height=2)

HO1_A_volplot_dat <- data.frame('pValue'=HO1_A_specific_Overlaps_log10_pval_GRanges_df$pValue, 'signalValue'=HO1_A_specific_Overlaps_log10_pval_GRanges_df$signalValue)
A <- ggplot(data=HO1_A_volplot_dat, aes(HO1_A_volplot_dat$pValue)) + geom_histogram(binwidth=5) + xlab('-log10pValue') + geom_vline(xintercept=3, color='red', size=.25)
A
ggsave('HO1_A_-log10pValue_histogram.png', width=2, height=2)


