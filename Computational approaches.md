# Details of our data analysis procedures in listed below;

**All the bioinformatics analyses were performed at IBG Bioinformatics Platform.**

**Alignment and Processing of Next-Generation Sequencing Data**
The ChIP-seq data were aligned to the human genome (hg38) using the nf-core/chipseq pipeline (https://nf-co.re/chipseq/1.2.1, accessed on 5 February 2021). Using this pipeline, 
KDM6A ChIP peaks were called, peak annotation was performed, and KDM6A signals on peaks were quantified in CPM.

**Peak Calling and Peak Overlap**
KDM6A peak calling was performed using MACS2 (embedded in the nf-core chip-seq pipeline) using the parameters: p < 10−5 and FDR < 0.1 and a filter with the narrow option.
For the set of peaks called for the T24, SV-HUC-1, and BdEC cell lines, the overlap among the peaks was checked using the “VennDiagram package” in R.

**Clustering of KDM6A Peaks**
Only the peaks which were located within a ±2 kb distance of the transcriptional start site (TSS) of genes and that had an fc value > 3 (value coming from MACS) were used
for the clustering analysis. Before clustering, CPM counts were normalized for the library sizes, and the log2 of the normalized values were used for plotting. Hierarchical 
clustering was performed using the “pheatmap package” (https://cran.r-project.org/web/packages/ pheatmap/index.html, accessed on 15 March 2021) in R. Clusters were determined 
and visualized using the cutree = 4 option in the package.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/47274420-97ad-46c2-95e2-807c2b2af2b4)


**Analysis of T24 KDM6A Peaks**
The TSS ±2 kb filtered peaks in the Venn Diagram were used for the analysis. Peaks were ranked according to their CPM normalized log2 values for each cluster. The pheatmap
R package was used for visualization.

**Transcription Factor Motif Finding**
Transcription factor motif finding for different KDM6A clusters was performed using the findMotifsGenome.pl command of HOMER (http://homer.ucsd.edu/homer/motif/, accessed on 
17 March 2021), with a size parameter of 300, using genomic positions.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/9803ebae-1b04-4486-8985-f2f2f325ecb4)


**Gene Ontology Analysis**
Gene Ontology (GO) term analysis was performed using ConsensusPathDB (http://cpdb.molgen.mpg.de/, accessed on 22 March 2021) with a GO Level of 3 and biological process (BP) 
options. The top 10 terms for ‘IH-MULTI’, ‘normal’, and ‘normal-immortal’ clusters were identified, and for those terms, GO Level 3 was visualized using the dotplot function 
of the ggplot2 package (https://cran.r-project.org/web/packages/ggplot2/index.html, accessed on 5 August 2021). Only the terms with a qval < 0.01 filter were visualized.

**Data Visualization**
ChIP-seq signal data visualization was performed using the Bioconductor Gviz package (https://bioconductor.org/packages/release/bioc/html/Gviz.html, accessed on 1 June 2021). 
“BW” files were obtained from CPM-normalized BAM and BAI files using deeptools “BamCoverage” functions’ default parameters. Visualization was performed with “BW” files using 
the ‘horizon’ type of GViz package. GENCODE hg38 Comprehensive Gene Annotation Version 30 data were used as the GeneRegionTrack for the annotation.

**String Protein Interaction Data Visualization**
For protein interaction data visualization, Cytoscape (version 3.8.2)’s default style was used with String data.

**HHEX and HES1 Gene Expression-Related Associations**
To obtain correlation heatmaps between the expressions of HHEX and HES1, primary bladder cancer Hiseq data from TCGA BLCA 2017 (n = 426) and bladder cancer cell line RNAseq 
data (n = 25) from CCLE (https://depmap.org/portal/download/all/(CCLE_RNAseq_genes_rpkm_20180929.gct, accessed on 17 December 2021) were used. Pearson correlation coefficients 
for the correlation of the expression of HHEX and HES1 and the genes associated with the term ‘regulation of the developmental process’ within the normalimmortal cluster (n = 34) 
were calculated using RPKM-normalized log2 expression values. The ggplot2 (version 3.3.5) and pheatmap (version 1.10.12) packages in R were used for visualization.

**Survival Analysis**
To generate Kaplan–Meier graphs showing patients survival times according to different expression levels of HHEX and HES1, gene expression RNAseq IlluminaHiSeq* (n = 426) data
and phenotype-curated survival data (n = 436) belonging to a TCGA Bladder Cancer (BLCA) cohort from Xena Browser Datasets [15] was used. Patients with survival data of OS.time > 120
(n = 382) were included in the Kaplan–Meier graphs. Low/high expression groups were formed by calculating the median values for HES1 and HHEX genes separately. Kaplan–Meier graphs 
were created using a log-rank test according to OS.time values. The coin (https://cran.r-project.org/web/packages/coin/index.html, accessed on 4 January 2022) (version 1.4-2) and
survival (https://cran.r-project.org/web/packages/survival/index.html, accessed on 4 January 2022) (version 3.2-7) R packages were used for the analysis.

![image](https://github.com/guldenozden/KDM6A/assets/147516133/735bc026-60d4-494b-86c9-9223fbfdd567)


**Comparison of T24 KDM6A Peaks with Published Data**
To compare the KDM6A peak profile in the T24 cell line to a wild-type KDM6A peak profile called in a bladder cancer cell line, published WT KDM6A-expressing UMUC-1 bladder cancer 
cell line’s narrowpeak data (Barrows et al., GSE157091) were used. Peak overlaps were performed using the ‘findOverlaps’ function and the rtracklayer (https://bioconductor.org/
packages/release/bioc/html/rtracklayer.html, accessed on 15 January 2023) and Genomic Ranges R packages (https://bioconductor.org/packages/release/bioc/ html/GenomicRanges.html, 
accessed on 15 January 2023), and subsequently visualized using a Venn diagram. A snapshot showing KDM6A occupancy at HES1 and TLE3 loci in wild-type KDM6A-expressing UMUC-1 cells
was created using the IGV genome browser(https://software.broadinstitute.org/software/igv/, accessed on 15 January 2023).
