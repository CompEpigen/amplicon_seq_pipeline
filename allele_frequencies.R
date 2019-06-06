library(deepSNV)
library(ggplot2)

DKFZ_RUN<-FALSE


#DATA_DIR<-"/C010-datasets/Internal/GCTB/MiSeq/raw_data/GCTB_Miseq_01_2019/03_bwa"

if(DKFZ_RUN){
    SAMPLE_SHEETS<-c(
            "/C010-datasets/Internal/GCTB/MiSeq/raw_data/GCTB_Miseq_01_2019/190128_M00269_0388_000000000-D5N72/190128_M00269_0388_000000000-D5N72_meta.tsv",
            "/C010-datasets/Internal/GCTB/MiSeq/raw_data/DKFZ_run/180806_M00269_0346_000000000-BRDD3_meta.tsv"
    )
    
    DATA_DIRS<-c(
            "/C010-datasets/Internal/GCTB/MiSeq/raw_data/DKFZ_run/03_bwa",
            "/C010-datasets/Internal/GCTB/MiSeq/raw_data/GCTB_Miseq_01_2019/03_bwa"
    )
}else{
    SAMPLE_SHEET<-"/C010-datasets/Internal/GCTB/MiSeq/sample_info/full_sample_annotation.txt"
    DATA_DIRS<-"/C010-datasets/Internal/GCTB/MiSeq/raw_data/H3F3A_MiSeq_reads_all"
}


#### PREPARE SAMPLE INFORMATION DATA

sample_sheet_full<-read.table("/ngs_share/scratch/pavlo/gctb/data/sample_info/170525_GCTB_sample_annotation_AB.tsv", sep="\t", header=TRUE)

sample_sheet_miseq<-sample_sheet_full[sample_sheet_full[["MiSeq"]]=="+", ]
sample_sheet_miseq$SampleID<-gsub("\\s", "_",sample_sheet_miseq$SampleID)
rownames(sample_sheet_miseq)<-sample_sheet_miseq$SampleID


#### READ THE BAM DATA IN

#bams<-list.files(DATA_DIR, recursive=TRUE, pattern=".bam$")

bams<-lapply(DATA_DIRS, list.files, recursive=TRUE, pattern=".bam$", full.names=TRUE)
bams<-unlist(bams)

bam_names<-lapply(DATA_DIRS, list.files, recursive=TRUE, pattern=".bam$")
bam_names<-unlist(bam_names)

rbams<-list()

if(DKFZ_RUN){
    chr_name<-"chr1"
}else{
    chr_name<-"1"
}
for(bami in seq(bams)){
    rbams[[bam_names[bami]]]<-bam2R(bams[bami], chr=chr_name, start=226252052, stop=226252320)
}

if(DKFZ_RUN){
    rbams<-rbams[sprintf("%s_sorted.bam",amplicon)]
    names(rbams)<-sample_sheet$SAMPLE_ID
}else{
    rbams<-rbams[sprintf("%s/%s.bam", sh_primers,amplicon)]
    names(rbams)<-sample_sheet$Sample
}


rbams<-rbams[!sapply(rbams, is.null)]

fmats<-lapply(rbams, function(fm){
    print(dim(fm))
    fm<-fm[,1:5]
    #fm<-fm[rowSums(fm)>0,]
    
    fm<-fm/rowSums(fm)
})


library(pheatmap)
for(fmn in seq(fmats)){
pheatmap(t(fmats[[fmn]]), cluster_col=FALSE, cluster_row=FALSE, col=grey.colors(200, 1,0), 
		file=sprintf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/test_heatmap_%s_sb_run.pdf",names(fmats)[fmn]))
}
#### first base

mutated_position<-104
mut_base<-2

mut_freq<-sapply(fmats, "[", i=mutated_position, j=mut_base)
mut_freq<-mut_freq[order(names(mut_freq))]


pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/barplot_of_T_G34_pos1_sb_run.pdf")
par(oma=c(12,4,4,2))
barplot(mut_freq, las=2, ylab="freq. of T at pos 226252052 (G34)")
dev.off()

lapply(rbams, "[", i=103:109, j=1:5)

#### second base

mutated_position<-105
mut_base<-2

mut_freq2<-sapply(fmats, "[", i=mutated_position, j=mut_base)
mut_freq2<-mut_freq2[order(names(mut_freq2))]

pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/barplot_of_T_G34_pos2_sb_run.pdf")
par(oma=c(12,4,4,2))
barplot(mut_freq2, las=2, ylab="freq. of T at pos 226252053 (G34)")
dev.off()


#### third base

mutated_position<-106
mut_base<-2

mut_freq3<-sapply(fmats, "[", i=mutated_position, j=mut_base)
mut_freq3<-mut_freq3[order(names(mut_freq3))]

pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/barplot_of_T_G34_pos3_sb_run.pdf")
par(oma=c(12,4,4,2))
barplot(mut_freq3, las=2, ylab="freq. of T at pos 226252054 (G34)")
dev.off()

#### -1 base

#mutated_position<-103
#mut_base<-2
#
#mut_freq2<-sapply(fmats, "[", i=mutated_position, j=mut_base)
#mut_freq2<-mut_freq2[order(names(mut_freq2))]
#
#pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/barplot_of_T_G34_pos-1_last_run.pdf")
#par(oma=c(12,4,4,2))
#barplot(mut_freq2, las=2, ylab="freq. of T at pos 226252054 (G34)")
#dev.off()



##### Combined plot
library(ggplot2)
sample_sheet$MutFreq1<-mut_freq[match(sample_sheet$Sample, gsub("\\.T", "", names(mut_freq)))]
sample_sheet$MutFreq2<-mut_freq2[match(sample_sheet$Sample, gsub("\\.T", "", names(mut_freq2)))]


##### save the data frames for combining later
if(DKFZ_RUN){
    saveRDS(sample_sheet, file="/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/dkfz_run_df.RDS")
}else{
    saveRDS(sample_sheet, file="/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/sb_run_df.RDS")
}

sample_sheets<-lapply(list.files("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/", pattern="*_run_df.RDS", full.names=TRUE), readRDS)
cnames<-lapply(sample_sheets, colnames)
com_cnames<-Reduce("intersect", cnames)
sample_sheets_std<-lapply(sample_sheets, "[", com_cnames)
sample_sheet_combined<-do.call("rbind", sample_sheets_std)

sample_sheet_combined$Source<-as.character(sample_sheet_combined$Source)
sample_sheet_combined$Source[sample_sheet_combined$Source=="CL1"]<-"cells"
sample_sheet_combined$Source[sample_sheet_combined$Source=="CL2"]<-"new_cells"

sample_sheet_combined<-sample_sheet_combined[sample_sheet_combined$Primer=="Primer 1",]


sample_subset<-c("06", "AS", "AO", "AZ")
df2plot<-sample_sheet_combined[sample_sheet_combined$SampleName %in% sample_subset & sample_sheet_combined$Primer=="Primer 1" & sample_sheet_combined$MaterialType != "old",]

#df2plot<-sample_sheet

#df2plot<-sample_sheet[sample_sheet$SampleName %in% c("AO", "AZ") & sample_sheet$Primer=="Primer 1",]
#df2plot<-sample_sheet[sample_sheet$Primer=="Primer 1",]

#df2plot$Source[5]<-"CL1"


#df2plot<-df2plot[df2plot$Source%in%c("tissue", "CL2"),]

####df2plot$Source<-factor(gsub("CL2", "cell line", df2plot$Source), levels=c("tissue","cell line"))

df2plot$Source<-factor(df2plot$Source, levels=c("tissue","cells", "new_cells"))

df2plot$Source<-gsub("new_cells", "cells", df2plot$Source)
df2plot$Source<-factor(df2plot$Source, levels=c("tissue","cells"))
df2plot$SampleName<-factor(df2plot$SampleName, levels=sample_subset)



df2plot<-df2plot[!is.na(df2plot$MutFreq1),]

#### first pos
p<-ggplot(df2plot, aes(x=SampleName, y=MutFreq1))

p<-p+geom_bar(aes(fill=Source), stat="identity", #fun.y=log_mean, 
                position=position_dodge(preserve = "single"), width=0.75, color="black")

p<-p+ggpubr:::theme_pubr() + theme(axis.text.x=element_text(angle=45,vjust=c(1), hjust=c(1)))
#p<-p+facet_wrap("Primer")
p<-p + ylab("VAF, c.103G > T")
p<-p + xlab("Patient")

pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/MiSeq_barplot_pos1_new_all_runs_final.pdf", #height=2.5, width=3.5)
height=5,width=7)
print(p)
dev.off()

#### 2nd pos

p<-ggplot(df2plot, aes(x=SampleName, y=MutFreq2))

p<-p+geom_bar(aes(fill=Source), stat="identity", #fun.y=log_mean, 
        position=position_dodge(preserve = "single"), width=0.75, color="black")

p<-p+ggpubr:::theme_pubr()+ theme(axis.text.x=element_text(angle=45,vjust=c(1), hjust=c(1)))
#p<-p+facet_wrap("Primer")
p<-p + ylab("VAF, c.104G > T")
p<-p + xlab("Patient")

pdf("/ngs_share/scratch/pavlo/gctb/analysis/MiSeq/MiSeq_barplot_pos2_new_all_runs_final.pdf",# height=2.5, width=3.5)
        height=5,width=7)
print(p)
dev.off()


