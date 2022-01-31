# Scaffa

* These experiments were performed on MLE-12 cells, a *Mus musculus* lung epithelial cell line.
* Chromatin was isolated from two biological replicates (GroupA and GroupB cultures) and immunoprecipitated using HO1, IgG, or a no-antibody control (LH20/H20), leading to 6 samples in total.
* Genewiz performed 150-BP PE sequencing with approximately 58 million reads per sample.
* All analyses were ran on Oscar using Bioflows and the CBC conda env, see the `Scaffa_report/docs/Report.html` file for more details on the analysis pipeline used.

In brief:

* QC analyses of raw and trimmed reads was performed in FastQC (V0.11.5)

* To remove adapters, low-quality sequences, and short sequences, trimming was performed in Trimmomatic (V0.36) using the following parameters:
```
- trimmomatic:
    subcommand: PE
    options:
      "ILLUMINACLIP:/gpfs/data/cbc/cbc_conda_v1/envs/cbc_conda/opt/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:5:6:true":
      "SLIDINGWINDOW:10:25 MINLEN:75":
```

* Trimmed reads were aligned against the Ensembl assembly GRCm38.p6 (`/gpfs/data/cbc/references/Ensembl_M_musculus.GRCm38.p6`).
* Reference indexing (`bwa index`)and read mapping (`bwa mem -M`) were performed in bwa (V0.7.15).
* Subsequent SAM conversion to BAM (`samtools view -Sbh`), as well as BAM sorting (`samtools sort`) and indexing (`samtools index`) were performed in Samtools (V1.9).
* Duplicates were marked (`MarkDuplicates`) and alignment statistics (`CollectAlignmentSummaryMetrics`, `CollectWGSMetrics`) were extracted using Picard tools (V2.9.2).

* MACS2 (V 2.1.2) for peak calling, used the input libraries (H20 and LH20) as controls, while the IgG and HO1 libraries were considered treatments.
* Each treatment was compared to both control files.
* The following additional flags were used: 
- `-f BAMPE` to specify that the libraries were paired-end
- `-g mm` to specify that the *Mus musculus* pre-built specs should be used to estimate genome size
- `-p 1e-3` a slightly relaxed p-value was used so that we would recover more peaks (with the understanding that we could filter out low p-values later if necessary)

* Using ChipPeakAnno (V3.16.1), we calculated a permissive and conservative estimate of the peaks that overlap between HO1 and IgG (`findOverlapsOfPeaks`).
* The conservative estimate is based on comparing all of the ChIP libraries (HO1 and IgG) to both input libraries (GroupA LH20 and GroupB H20) to make peak calls.
* It requires that a list of peaks be present in all HO1 peak calls but not both of the IgG peak calls.
* This convervative estimate finds 38 peaks specific to HO1.
* A more permissive estimate is based on all the GroupA and GroupB ChIP libraries but input control from GroupA only, as there were some issues with the GroupB input library having more unmapped reads.
* This more permissive estimate finds 211 peaks specific to HO1.
* All analyses from this point forward are based on this more permissive set.

* Using the more permissive set of peaks (Part 13), we annotated the peaks using ChipSeeker (V1.18.0) (`annotatePeak`).
* This required creating a TxDb which was done using GenomicFeatures (V1.34.3) (`makeTxDbFromGFF`) (based the TxDb on this [GFF3](ftp://ftp.ensembl.org/pub/release-95/gff3/mus_musculus).)

