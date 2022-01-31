library(ChIPpeakAnno)
library(GenomicFeatures)
library(ChIPseeker)
library(GenomicRanges)
library(stringr)
library(ReactomePA)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)

#get files
list <- list.files(path='.', pattern='.narrowPeak')

#fix names
names(list) <- substr((str_replace_all(list, 'Group', '')),1,12)
   
#convert the list to GRanges
list_gr <- lapply(list, function(x) toGRanges(x, format=c('narrowPeak'), header=F))

#import gff
gff3_file <- ('Mus_musculus.GRCm38.95.gff3')

#make txdb
txdb <- makeTxDbFromGFF(gff3_file)

#make peak annotation list of all peaks
peakannolist <- lapply(list, annotatePeak)

#######################################################################################################

#calculate more strict overlaps
strict_ol_HO1_peaklist <- findOverlapsOfPeaks(list_gr[['A_HO1_A_LH20']], list_gr[['A_HO1_B_H20_']], list_gr[['B_HO1_A_LH20']], list_gr[['B_HO1_B_H20_']])
strict_ol_IgG_peaklist<- findOverlapsOfPeaks(list_gr[['A_IgG_A_LH20']], list_gr[['A_IgG_B_H20_']], list_gr[['B_IgG_A_LH20']], list_gr[['B_IgG_B_H20_']])

#calculate a more permissive estimate of overlaps using only the LH20 library:
permissive_ol_HO1 <- findOverlapsOfPeaks(list_gr[['A_HO1_A_LH20']], list_gr[['B_HO1_A_LH20']])
permissive_ol_IgG <- findOverlapsOfPeaks(list_gr[['A_IgG_A_LH20']], list_gr[['B_IgG_A_LH20']])

#get the strict overlapping peaks out...
strict_ol_HO1_peaklist <- strict_ol_HO1_peaklist$peaklist$`list_gr...A_HO1_A_LH20...///list_gr...A_HO1_B_H20_...///list_gr...B_HO1_A_LH20...///list_gr...B_HO1_B_H20_...`
strict_ol_IgG_peaklist <- strict_ol_IgG_peaklist$peaklist$`list_gr...A_IgG_A_LH20...///list_gr...A_IgG_B_H20_...///list_gr...B_IgG_A_LH20...///list_gr...B_IgG_B_H20_...`

#get the permissive overlapping peaks out...
permissive_ol_HO1_peaklist <- permissive_ol_HO1$peaklist$`list_gr...A_HO1_A_LH20...///list_gr...B_HO1_A_LH20...`
permissive_ol_IgG_peaklist <- permissive_ol_IgG$peaklist$`list_gr...A_IgG_A_LH20...///list_gr...B_IgG_A_LH20...`

#find the most strict overlaps...
strict_HO1_IgG_peaks_comparison <- findOverlapsOfPeaks(strict_ol_HO1_peaklist, strict_ol_IgG_peaklist)

#find the more permissive overlaps...
permissive_HO1_IgG_peaks_comparison <- findOverlapsOfPeaks(permissive_ol_HO1_peaklist, permissive_ol_IgG_peaklist)

#######################################################################################################
#just work w the permissive set now...
#pull out just the peaks that are in HO1A and B but not in both IgG libraries....
permissive_HO1_peaks_HO1_specific <- permissive_HO1_IgG_peaks_comparison$peaklist$permissive_ol_HO1_peaklist

#make a dataframe version too...
permissive_HO1_peaks_HO1_specific_df <- as(permissive_HO1_peaks_HO1_specific, 'data.frame')

#make peak annotation list
permissive_HO1_peaks_anno <- annotatePeak(permissive_HO1_peaks_HO1_specific, TxDb=txdb)

#covert permissive_HO1_peaks_anno to a GRanges object
permissive_HO1_peaks_anno_GRanges <- as.GRanges(permissive_HO1_peaks_anno)



###############

#make some figs



#Make a plot of the significance value of each peak...

#find OLs betwn A_HO1_A_LH20 / A_HO1_A_LH20 and permissive_HO1_peaks_HO1_specific.. really this is just to get the associated q-value data back out for each input peak
#note that the numbers of peaks might change slightly if the numbers of peaks in 2 given overlapping regions are different (e.g., groupA HO1 peak is spread out over 3 small peaks and groupB is 1 larger peak but the regions containing the peaks still overlap)
HO1_A_specific_Overlaps_metadata <- findOverlapsOfPeaks(permissive_HO1_peaks_HO1_specific, list_gr[['A_HO1_A_LH20']])
HO1_B_specific_Overlaps_metadata <- findOverlapsOfPeaks(permissive_HO1_peaks_HO1_specific, list_gr[['B_HO1_A_LH20']])

#pull out just the peaks...
HO1_A_specific_Overlaps_metadata_GRanges <- HO1_A_specific_Overlaps_metadata$overlappingPeaks$`permissive_HO1_peaks_HO1_specific///list_gr...A_HO1_A_LH20...`
HO1_B_specific_Overlaps_metadata_GRanges <- HO1_B_specific_Overlaps_metadata$overlappingPeaks$`permissive_HO1_peaks_HO1_specific///list_gr...B_HO1_A_LH20...`

#convert to dataframes...
HO1_A_specific_Overlaps_metadata_GRanges_df<- as(HO1_A_specific_Overlaps_metadata_GRanges, 'data.frame')
HO1_B_specific_Overlaps_metadata_GRanges_df<- as(HO1_B_specific_Overlaps_metadata_GRanges, 'data.frame')

#pull out just the qValue info and make the histo...
HO1_A_specific_Overlaps_metadata_GRanges_df_qvalue <- HO1_A_specific_Overlaps_metadata_GRanges_df[,'qValue', drop=FALSE]
HO1_A_specific_Overlaps_metadata_GRanges_df_qvalue_histo <- ggplot(data=HO1_A_specific_Overlaps_metadata_GRanges_df_qvalue, aes(HO1_A_specific_Overlaps_metadata_GRanges_df_qvalue$qValue)) + geom_histogram(binwidth=5) + xlab('-log10qValue') + geom_vline(xintercept=2, color='red', size=.25)

HO1_B_specific_Overlaps_metadata_GRanges_df_qvalue <- HO1_B_specific_Overlaps_metadata_GRanges_df[,'qValue', drop=FALSE]
HO1_B_specific_Overlaps_metadata_GRanges_df_qvalue_histo <- ggplot(data=HO1_B_specific_Overlaps_metadata_GRanges_df_qvalue, aes(HO1_B_specific_Overlaps_metadata_GRanges_df_qvalue$qValue)) + geom_histogram(binwidth=5) + xlab('-log10qValue') + geom_vline(xintercept=2, color='red', size=.25)

#pull out the qValue info for the permissive IgG peaks also...

permissive_IgG_metadata <- permissive_ol_IgG$all.peaks$list_gr...A_IgG_A_LH20...
write.table(permissive_IgG_metadata, 'permissive_IgG_metadata.txt', sep='\t', row.names=F)

permissive_IgG_metadata_uniquePeaks <- permissive_ol_IgG$uniquePeaks
write.table(permissive_IgG_metadata_uniquePeaks, 'permissive_IgG_metadata_uniquePeaks.txt', sep='\t', row.names=F)

#write the dfs
write.table(HO1_A_specific_Overlaps_metadata_GRanges_df, 'permissive_HO1_peaks_HO1_specific_A_peaklist_df.txt', sep='\t', row.names=F)
write.table(HO1_B_specific_Overlaps_metadata_GRanges_df, 'permissive_HO1_peaks_HO1_specific_B_peaklist_df.txt', sep='\t', row.names=F)

#Pie chart and bar chart of peak location in context of genomic features

plotAnnoPie(permissive_HO1_peaks_anno)
ggsave('permissive_HO1_peaks_plotAnnoPie.png')

plotDistToTSS(permissive_HO1_peaks_anno)
ggsave('permissive_HO1_peaks_plotDistToTSS.png')

#make a plot of HO1-specific ChIP peaks over chromosomes
covplot(permissive_HO1_peaks_anno_GRanges)
ggsave('permissive_HO1_peaks_covplot.png')


#Peaks per chromosome amd [eak width per chromosome

boxplot <-ggplot(permissive_HO1_peaks_HO1_specific_df, aes(x=seqnames, y=width)) + geom_boxplot() + xlab('chromosome')
boxplot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/jwalla12/Desktop/Scaffa_report/docs/assets/peak_width_per_chromosome_HO1_specific.png', height=3, width=4)

p <- ggplot(data=permissive_HO1_peaks_HO1_specific_df, aes(permissive_HO1_peaks_HO1_specific_df$seqnames)) + geom_histogram(binwidth=5, stat='count') + xlab('chromosome')
p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('/Users/jwalla12/Desktop/Scaffa_report/docs/assets/peaks_per_chromosome_HO1_specific.png', height=3, width=4)

#####

#get a list of gene IDs
gene_list=(as.data.frame(permissive_HO1_peaks_anno)$geneId)

#convert gene IDs to entrez format
gene_list_formatted <- bitr(gene_list, fromType='ALIAS', toType='ENTREZID', OrgDb="org.Mm.eg.db")

#enrich...
#pathway1 <- enrichPathway(gene=gene_list_formatted$ENTREZID, organism='mouse')
pathway1 <- enrichPathway(gene=gene_list_formatted$ENTREZID, organism='mouse', pAdjustMethod='bonferroni')
dotplot(pathway1)

ggsave('permissive_HO1_peaks_enrichPathway.png')

