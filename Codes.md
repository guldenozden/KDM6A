
**Code for generating the heatmap in Figure 1D**

[Usetwd("C:/Users/gulden.ozden/Desktop/Loreal_Chip-seq_data_22.01.21/narrowPeak/consensus/KDM6A")
peaks<-read.delim("KDM6A.consensus_peaks.boolean.annotatePeaks.txt")
featurecounts<-read.table("KDM6A.consensus_peaks.featureCounts.txt",header = T)

#b<- peaks[which(peaks$Distance.to.TSS < 2000 & peaks$Distance.to.TSS > -2000 & peaks$BdEC_KDM6A_IP_R1.qval > 1),]
#s<- peaks[which(peaks$Distance.to.TSS < 2000 & peaks$Distance.to.TSS > -2000 & peaks$SV.HUC.1_KDM6A_IP_R1.qval > 1),]
#t<-peaks[which(peaks$Distance.to.TSS < 2000 & peaks$Distance.to.TSS > -2000 & peaks$T24_KDM6A_IP_R1.qval > 1),]
#m<-rbind(b, s, t)
#m[10:12][is.na(m[10:12])] <- 0
#fc<-m[rowSums(m[10:12] > 3) >=1 ,]
#mcall<-featurecounts[featurecounts$Geneid%in%fc$interval_id,]

filt<- peaks[which(peaks[10:12]>1)>=1,]
filt1<-filt[which(filt$Distance.to.TSS<2000 & filt$Distance.to.TSS>-2000),]
filt1[10:12][is.na(filt1[10:12])] <- 0
fctest<-filt1[rowSums(filt1[10:12] > 3) >=1 ,]
mcalltest<-featurecounts[featurecounts$Geneid%in%fctest$interval_id,]


library(pheatmap)
df<-data.frame(row.names = mcalltest$Geneid, T24 = mcalltest$T24_KDM6A_IP_R1.mLb.clN.sorted.bam, BdEC = mcalltest$BdEC_KDM6A_IP_R1.mLb.clN.sorted.bam, SVHUC1 = mcalltest$SV.HUC.1_KDM6A_IP_R1.mLb.clN.sorted.bam)
#cpm normalization
cpm<- data.frame(row.names = rownames(df), T24 = df$T24/64002493*1000000, BdEC = df$BdEC/53300895*1000000, SVHUC1 = df$SVHUC1/52457026*1000000)
#log2
logx<- data.frame(row.names = rownames(cpm), T24 = log2(cpm$T24 +1), BdEC = log2(cpm$BdEC +1), SVHUC1 = log2(cpm$SVHUC1 +1))
set.seed(12345)
obj <- pheatmap(logx,show_rownames = F, cluster_cols = F, scale = "row", cutree_rows = 4)
cl <- cutree(obj$tree_row,4)
ann <- as.data.frame(cl)
rownames(ann) <- rownames(logx)
out<-pheatmap(logx,show_rownames = F, cluster_rows = T, cluster_cols = F, scale = "row", cutree_rows = 4 , annotation_row = ann, fontsize = 7)

clusDesignation <- cutree(as.hclust(out$tree_row), 4) #number of clusters
clusDesignation[clusDesignation==4] #cluster number you want
test<-names(clusDesignation[clusDesignation==4])
int<-peaks[peaks$interval_id%in%test ,]
clu <- data.frame(int$chr , int$start , int$end , int$interval_id , int$Gene.Name)
#write.table(clu, "tss2000cluster4anno.txt" , quote = F)
gene<-clu$int.Gene.Name
ugene<-unique(gene)
write.table(ugene, "tss2kb_fc3_np_cut4_cl4_t24genes.txt" , sep = "/t", quote = F , row.names = F)


###----for motif finding analysis----
library(rtracklayer)
clusDesignation <- cutree(as.hclust(out$tree_row), 4) #number of clusters
clusDesignation[clusDesignation==2] #cluster number you want
test<-names(clusDesignation[clusDesignation==2])
int<-peaks[peaks$interval_id%in%test ,]
region<-data.frame(int$chr , int$start , int$end)
write.table(region,"cluster1narrowpeakregion.txt" , quote = F , row.names = F, sep = "\t")
export(region, "cluster3.bed")


width<- mean(int$end - int$start)

###----duplicated interval_ids-----## cause of rbind ##

#test2 <- m$T24_KDM6A_IP_R1.fc
#test3 <- m$SV.HUC.1_KDM6A_IP_R1.fc
#int1 <- m$interval_id[grep(";",test2)]
#int2 <- m$interval_id[grep(";",test3)]
#allint <- c(int1,int2)
#summary(allint %in% fc$interval_id)
#m[which(m$interval_id=="Interval_7413"),]
#which(m$interval_id=="Interval_7413")
#allids <- m$interval_id
#length(unique(allids))
#dup <- duplicated(allids)
#length(unique(mcall$Geneid))


ploading narrowpeaks_heatmap.R…]()


