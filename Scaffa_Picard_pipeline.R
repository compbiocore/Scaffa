setwd("/Users/jwalla12/Scaffa_project/Scaffa_picardtools_metrics")
library(tidyverse)

#import the _wgs_metrics file
wgs_metrics_files <- list.files(pattern = '_wgs_metrics.txt')
#import the metrics file, skip the first 7 rows and import only 1 row. Supply the column names and delimiter and specify no header.
wgs_metrics <- lapply(wgs_metrics_files, function(wgs_metrics_file) {
  wgs_metric_rows<-read.table(wgs_metrics_file, 
                              skip=7, 
                              nrows=1, 
                              col.names=c('GENOME_TERRITORY','MEAN_COVERAGE','SD_COVERAGE','MEDIAN_COVERAGE','MAD_COVERAGE','PCT_EXC_MAPQ','PCT_EXC_DUPE','PCT_EXC_UNPAIRED','PCT_EXC_BASEQ','PCT_EXC_OVERLAP','PCT_EXC_CAPPED','PCT_EXC_TOTAL','PCT_1X','PCT_5X','PCT_10X','PCT_15X','PCT_20X','PCT_25X','PCT_30X','PCT_40X','PCT_50X','PCT_60X','PCT_70X','PCT_80X','PCT_90X','PCT_100X','HET_SNP_SENSITIVITY','HET_SNP_Q'), header=FALSE, sep='\t')
  #add a column called library and populate it with the name of the file that was used to populate the rest of the row.
  wgs_metric_rows$library <- wgs_metrics_file
  wgs_metric_rows
}
)
#turn the information into a dataframe and write the outputs, specify row.names=false so that there isn't an index column returned and then column headers line up properly.
wgs_metrics_df<-map_df(wgs_metrics, ~as.data.frame(.x))





marked_dups_metrics_files <- list.files(pattern = 'marked_duplicates_metrics.txt')
marked_dups_metrics <- lapply(marked_dups_metrics_files, function(marked_dups_metrics_file) {
  marked_dups_metric_rows <- read.table(marked_dups_metrics_file, skip=7, nrows=1,col.names=c('LIBRARY','UNPAIRED_READS_EXAMINED','READ_PAIRS_EXAMINED','SECONDARY_OR_SUPPLEMENTARY_RDS','UNMAPPED_READS','UNPAIRED_READ_DUPLICATES','READ_PAIR_DUPLICATES','READ_PAIR_OPTICAL_DUPLICATES','PERCENT_DUPLICATION','ESTIMATED_LIBRARY_SIZE'), header=FALSE, sep='\t')
  marked_dups_metric_rows$library <- marked_dups_metrics_file
  marked_dups_metric_rows
}
)
marked_dups_metrics_df<-map_df(marked_dups_metrics, ~as.data.frame(.x))




#import the _alignment_summary_metrics file
alignment_summary_metrics_files <- list.files(pattern = '_alignment_summary_metrics.txt')
#import the metrics file, skip the first 9 rows and import only 1 row. Supply the column names and delimiter and specify no header.
alignment_summary_metrics <- lapply(alignment_summary_metrics_files, function(alignment_summary_metrics_file) {
  alignment_summary_metric_rows<-read.table(alignment_summary_metrics_file, 
                                            skip=9, 
                                            nrows=1, 
                                            col.names=c('CATEGORY', 'TOTAL_READS', 'PF_READS', 'PCT_PF_READS', 'PF_NOISE_READS', 'PF_READS_ALIGNED', 'PCT_PF_READS_ALIGNED', 'PF_ALIGNED_BASES', 'PF_HQ_ALIGNED_READS', 'PF_HQ_ALIGNED_BASES', 'PF_HQ_ALIGNED_Q20_BASES', 'PF_HQ_MEDIAN_MISMATCHES', 'PF_MISMATCH_RATE', 'PF_HQ_ERROR_RATE', 'PF_INDEL_RATE', 'MEAN_READ_LENGTH', 'READS_ALIGNED_IN_PAIRS', 'PCT_READS_ALIGNED_IN_PAIRS', 'PF_READS_IMPROPER_PAIRS', 'PCT_PF_READS_IMPROPER_PAIRS', 'BAD_CYCLES', 'STRAND_BALANCE', 'PCT_CHIMERAS', 'PCT_ADAPTER', 'SAMPLE', 'LIBRARY', 'READ_GROUP'), header=FALSE, sep='\t')
  #add a column called library and populate it with the name of the file that was used to populate the rest of the row.
  alignment_summary_metric_rows$library <- alignment_summary_metrics_file
  alignment_summary_metric_rows
}
)
#turn the information into a dataframe and write the outputs, specify row.names=false so that there isn't an index column returned and then column headers line up properly.
library(tidyverse)
alignment_summary_metrics_df<-map_df(alignment_summary_metrics, ~as.data.frame(.x))


marked_dups_metrics_df$library <- substr(marked_dups_metrics_df$library, 1, 10)
wgs_metrics_df$library <- substr(wgs_metrics_df$library, 1, 10)
alignment_metrics_df$library <- substr(wgs_metrics_df$library, 1, 10)

marked_dups_metrics_df$library <- factor(marked_dups_metrics_df$library, levels=c('GroupA_LH2','GroupA_IgG','GroupA_HO1','GroupB_H20','GroupB_IgG','GroupB_HO1'))
wgs_metrics_df$library <- factor(wgs_metrics_df$library, levels=c('GroupA_LH2','GroupA_IgG','GroupA_HO1','GroupB_H20','GroupB_IgG','GroupB_HO1'))
alignment_metrics_df$library <- factor(alignment_metrics_df$library, levels=c('GroupA_LH2','GroupA_IgG','GroupA_HO1','GroupB_H20','GroupB_IgG','GroupB_HO1'))


write.table(wgs_metrics_df, 'wgs_metrics_df.txt', sep='\t', row.names=FALSE)
write.table(marked_dups_metrics_df, 'marked_dups_metrics_df.txt', sep='\t', row.names=FALSE)
write.table(alignment_summary_metrics_df, 'alignment_summary_metrics_df.txt', sep='\t', row.names=FALSE)

