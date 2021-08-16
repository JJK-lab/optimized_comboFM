# optimized_comboFM
#bin width=50
all_output=read.table("filter.snv_rm_exon.bed",sep="\t",stringsAsFactors=F)[,-c(2,4,7:8)]
colnames(all_output)[c(1:2,5)]=c("CHROM","POS","SAMPLE")
all_output$SAMPLE=sub("_WGS","",all_output$SAMPLE)
all_output$bin1=(ceiling(all_output$POS/50)-1)*50+1
all_output$bin2=(ceiling(all_output$POS/50))*50
all_output$chr_bin1_bin2=paste(all_output$CHROM,all_output$bin1,all_output$bin2,sep="_")

chr_bin=names(table(all_output$chr_bin1_bin2)[table(all_output$chr_bin1_bin2)>2])
dt_cluster=list()
for(i in 1:length(chr_bin)){
    dt_bin=all_output[all_output$chr_bin1_bin2==chr_bin[i],]
    dt_cluster[[i]]=data.frame(CHROM_bin=chr_bin[i],start=min(dt_bin$POS),end=max(dt_bin$POS),n_Sample=length(unique(dt_bin$SAMPLE)),Sample=paste(unique(dt_bin$SAMPLE),collapse=","),recurrent_percentage=length(unique(dt_bin$SAMPLE))/206,n_mutation_site=nrow(dt_bin),n_mutation_site_per_Sample=paste(paste(names(table(dt_bin$SAMPLE)),table(dt_bin$SAMPLE),sep=":"),collapse=";"),mutation_site_per_Sample=paste(paste(dt_bin$SAMPLE,dt_bin$POS,sep=":"),collapse="|"))
    }

library("dplyr")
output_cluster=bind_rows(dt_cluster)
write.table(output_cluster,"Output_cluster.xls",sep="\t",row.names=F,quote=F)

output_cluster=read.table("Output_cluster.xls",sep="\t",header=T,stringsAsFactors=F)
output_cluster$len=output_cluster$end-output_cluster$start+1

#calculate background mutation frequency
library(foreach)
arm=read.table("cytoband.data",sep="\t",stringsAsFactors = F)
chr=paste("chr",c(1:22,"X","Y"),sep="")
chr_length=foreach(i=1:length(chr),.combine=c) %do% max(as.numeric(arm[arm$V1==chr[i],3]))

library(maftools)
maf=read.maf("Somatic_mutation.filter.208.maf", removeDuplicatedVariants =F)
smt=subsetMaf(maf, includeSyn = T, tsb = NULL, genes = NULL,fields = NULL, query = NULL, mafObj = FALSE, isTCGA = FALSE)

smt$Tumor_Sample_Barcode=as.character(smt$Tumor_Sample_Barcode)
smt$Tumor_Sample_Barcode=sub("_WGS","",smt$Tumor_Sample_Barcode)
smt_tab=table(smt$Tumor_Sample_Barcode)
backgr_mt_rate=smt_tab/sum(chr_length)

output_cluster$backgr_mt_rate=foreach(i=1:nrow(output_cluster),.combine=c) %do% mean(backgr_mt_rate[strsplit(output_cluster$Sample[i],split=",")[[1]]])
output_cluster$p_value=foreach(i=1:nrow(output_cluster),.combine=c) %do% pnbinom(output_cluster$n_mutation_site[i], output_cluster$len[i], prob=output_cluster$backgr_mt_rate[i])
output_cluster$p.adj=p.adjust(output_cluster$p_value,method = "BH")
write.table(output_cluster,"hotspot_Output_cluster.xls",sep="\t",row.names=F,quote=F)

#annotation
id_file=read.table("gencode.v27.metadata.HGNC/data")
colnames(id_file)[1]="TRANS"

output_cluster=read.table("hotspot_Output_cluster.xls",sep="\t",header=T,stringsAsFactors=F)
output_cluster$CHROM=foreach(i=1:nrow(output_cluster),.combine=c) %do% strsplit(output_cluster$CHROM_bin[i],split="_")[[1]][1]

library(readxl)
library(foreach)
utr3=read.table("hg38_annotation/gencode.v27.3utr.bed/data",sep="\t",stringsAsFactors=F)
utr5=read.table("hg38_annotation/gencode.v27.5utr.bed/data",sep="\t",stringsAsFactors=F)
intron=read.table("hg38_annotation/gencode.v27.intron.bed/data",sep="\t",stringsAsFactors=F)
enhancer=read.table("hg38_annotation/genehancer.v4.7.bed/data",sep="\t",stringsAsFactors=F)
promoter=read.table("hg38_annotation/gencode.v27.promoterCore.bed/data",sep="\t",stringsAsFactors=F)
IGR=read.table("hg38_annotation/IGR2.1/data",sep="\t",stringsAsFactors=F)
colnames(IGR)=paste("V",c(1:3,7:9,4:6),sep="")

chr=paste("chr",c(1:22,"X","Y"),sep="")
dt=list()
fun=c("Promoter","3UTR","5UTR","Intron","Enhancer","IGR")
dat=list()
for(i in 1:length(chr)){
  dt[[i]]=output_cluster[output_cluster$CHROM==chr[i],]
  gn=list()
  gn[[1]]=promoter[promoter$V1==chr[i],]
  gn[[2]]=utr3[utr3$V1==chr[i],]
  gn[[3]]=utr5[utr5$V1==chr[i],]
  gn[[4]]=intron[intron$V1==chr[i],]
  gn[[5]]=enhancer[enhancer$V1==chr[i],]
  gn[[6]]=IGR[IGR$V1==chr[i],]
  for(j in 1:nrow(dt[[i]])){
    gl=foreach(k=1:length(gn)) %do% which(dt[[i]]$end[j]>=gn[[k]]$V2&dt[[i]]$start[j]<=gn[[k]]$V3)
    gn_rg=foreach(k=1:length(gn)) %do% ifelse(length(gl[[k]])==0,NA,paste(paste(gn[[k]]$V4[gl[[k]]],collapse="|"),fun[k],sep=":"))
    gn_rg[[5]]=ifelse(length(gl[[5]])==0,NA,paste(paste(gn[[5]]$V6[gl[[5]]],collapse="|"),fun[5],sep=":"))
    gn_lst=na.omit(unique(unlist(gn_rg)))
    dat[[j]]=cbind(dt[[i]][rep(j,length(gn_lst)),],Gene_region=gn_lst)
    }
   dt[[i]]=do.call(rbind,dat)
}
new_cluster=do.call(rbind,dt)
new_cluster=as.data.frame(new_cluster)
new_cluster$Region=unlist(foreach(j=1:nrow(new_cluster)) %do% strsplit(as.character(new_cluster$Gene_region[j]),split=":")[[1]][2])

cluster=list()
for(j in 1:nrow(new_cluster)){
gn=strsplit(strsplit(as.character(new_cluster$Gene_region[j]),split=":")[[1]][1],split="\\|")[[1]]
cluster[[j]]=cbind(new_cluster[rep(j,length(gn)),],Gene=gn)
}
cluster_gene=do.call(rbind,cluster)

gl=which(cluster_gene$Region=="Enhancer")
enhancer=list()
for(j in 1:length(gl)){
gn=strsplit(as.character(cluster_gene$Gene[gl[j]]),split=";")[[1]]
gn_sy=sub("connected_gene=","",gn[2])
gn_sc=sub("score=","",gn[3])
enhancer[[j]]=cbind(cluster_gene[rep(gl[j],length(gn_sy)),-ncol(cluster_gene)],data.frame(Gene=gn_sy,Enhancer_Score=gn_sc))
}
enhancer_gn=do.call(rbind,enhancer)
cluster_gene$Enhancer_Score=NA
cluster_gn=rbind(cluster_gene[-gl,],enhancer_gn)

cluster_gn$TRANS=cluster_gn$Gene
cluster_gnsy=merge(cluster_gn[,colnames(cluster_gn)!="Gene"],id_file,by="TRANS")
colnames(cluster_gnsy)[ncol(cluster_gnsy)]="Gene"
cluster_gnsy=rbind(cluster_gnsy,cluster_gn[!is.element(cluster_gn$TRANS,cluster_gnsy$TRANS),])

cluster_all=rbind(cluster_gn[cluster_gn$Region=="Enhancer",],cluster_gnsy[,c(2:19,1)])

write.table(unique(cluster_all),"hotspot_promoterCore.xls",sep="\t",row.names=F,quote=F)

#annotation 
cgc_smg=read.table("CGC_SMG_gene.xls",sep="\t",stringsAsFactors=F,header=T)
rownames(cgc_smg)=cgc_smg$Gene

file=c("hotspot_promoterCore.xls")

dat=foreach(i=1:length(file)) %do% read.table(file[i],sep="\t",stringsAsFactors=F,header=T)
for(i in 1:length(dat)){
dat[[i]]$smg_source=""
dat[[i]]$smg_source[is.element(dat[[i]]$Gene,rownames(cgc_smg))]=cgc_smg[dat[[i]]$Gene[is.element(dat[[i]]$Gene,rownames(cgc_smg))],"smg_source"]

dat[[i]]$cgc_smg=""
dat[[i]]$cgc_smg[is.element(dat[[i]]$Gene,rownames(cgc_smg))]=cgc_smg[dat[[i]]$Gene[is.element(dat[[i]]$Gene,rownames(cgc_smg))],"cgc_smg"]

write.table(dat[[i]],sub(".xls","_cgc_msg.xls",file[i]),sep="\t",row.names=F,quote=F)
}

#RNA expression test
  mRNA=read.table("DESeq2_normalizedCount_134.xls",header=T,stringsAsFactors = F)
  mRNA=mRNA[,-which(is.element(colnames(mRNA),c("T502","N502")))]
  t_sp=colnames(mRNA)[grep("T",colnames(mRNA))]
  mRNA$T_Mean=rowMeans(mRNA[,grep("T",colnames(mRNA))])
  mRNA$N_Mean=rowMeans(mRNA[,grep("N",colnames(mRNA))])
   
  file=c("hotspot_promoterCore_cgc_msg.xls") 
  ht=foreach(i=1:length(file)) %do% read.table(file[i],sep="\t",header=T,stringsAsFactors = F)
  for(i in 1:length(ht)){
  dt=ht[[i]]
  sp=foreach(j=1:nrow(dt)) %do% intersect(strsplit(gsub("_WGS","",dt$Sample[j]),split=",")[[1]],t_sp)
  n_sample=foreach(j=1:nrow(dt),.combine=c) %do% length(sp[[j]])
  mut_rna=foreach(j=1:nrow(dt)) %do% as.numeric(mRNA[dt$Gene[j],sp[[j]]])
  wt_rna=foreach(j=1:nrow(dt)) %do% as.numeric(mRNA[dt$Gene[j],setdiff(t_sp,sp[[j]])])
  dt$MUT_mRNA=unlist(foreach(j=1:nrow(dt)) %do% paste(mut_rna[[j]],collapse=";"))
  dt$WT_mRNA=unlist(foreach(j=1:nrow(dt)) %do% paste(wt_rna[[j]],collapse=";"))
  dt$MUT_mean=unlist(foreach(j=1:nrow(dt)) %do% mean(mut_rna[[j]]))
  dt$WT_mean=unlist(foreach(j=1:nrow(dt)) %do% mean(wt_rna[[j]]))
  dt$T_Mean=mRNA[dt$Gene,"T_Mean"]
  dt$N_Mean=mRNA[dt$Gene,"N_Mean"]
  dt$p.value=NA
  gl=which(n_sample>1&(dt$T_Mean>=1.5|dt$N_Mean>=1.5))
  dt$p.value[gl]=foreach(j=gl,.combine=c) %do% t.test(mut_rna[[j]],wt_rna[[j]])$p.value
  dt$q.value=p.adjust(dt$p.value,method="BH")
  #write.table(dt,"hotspot_promoter2.2k_assocaite_gene.xls",sep="\t",row.names=F,quote=F)
  write.table(dt,sub("cgc_smg","cgc_smg_associate",file[i]),sep="\t",row.names=F,quote=F)
  }

#associate gene expression
  mRNA=read.table("Annotated_lncRNA_FPKM_134.xls",sep="\t",header = T,stringsAsFactors = F)
  colnames(mRNA)=sub("_WTS","",colnames(mRNA))
  mRNA=mRNA[,-which(is.element(colnames(mRNA),c("T502","N502","T13","N13")))]
  t_sp=colnames(mRNA)[-c(1:4)][grep("T",colnames(mRNA)[-c(1:4)])]
  mRNA$T_Mean=rowMeans(mRNA[,t_sp])
  mRNA$N_Mean=rowMeans(mRNA[,sub("T","N",t_sp)])
  
  muti=names(table(mRNA$Official_Symbol))[table(mRNA$Official_Symbol)!=1]
  for(i in 1:length(muti)){
    gl=which(mRNA$Official_Symbol==muti[i])
    sgl=which.max(rowMeans(mRNA[gl,c("T_Mean","N_Mean")]))
    mRNA=mRNA[-gl[-sgl],]
  }
  rownames(mRNA)=mRNA$Official_Symbol
  mRNA=mRNA[,-c(1:4)]
  
 ht=foreach(i=1:length(file)) %do% read.table(sub("cgc_smg","cgc_smg_associate",file[i]),sep="\t",header=T,stringsAsFactors = F)
   for(i in 1:length(ht)){
  dt=ht[[i]]
  lnc=intersect(dt$Gene,rownames(mRNA))
  gl=which(is.element(dt$Gene,lnc))
  sp=foreach(j=gl) %do% intersect(strsplit(gsub("_WGS","",dt$Sample[j]),split=",")[[1]],t_sp)
  n_sample=foreach(j=1:length(sp),.combine=c) %do% length(sp[[j]])
  mut_rna=foreach(j=1:length(sp)) %do% as.numeric(mRNA[dt$Gene[gl[j]],sp[[j]]])
  wt_rna=foreach(j=1:length(sp)) %do% as.numeric(mRNA[dt$Gene[gl[j]],setdiff(t_sp,sp[[j]])])
  dt$MUT_mRNA[gl]=unlist(foreach(j=1:length(mut_rna)) %do% paste(mut_rna[[j]],collapse=";"))
  dt$WT_mRNA[gl]=unlist(foreach(j=1:length(mut_rna)) %do% paste(wt_rna[[j]],collapse=";"))
  dt$MUT_mean[gl]=unlist(foreach(j=1:length(mut_rna)) %do% mean(mut_rna[[j]]))
  dt$WT_mean[gl]=unlist(foreach(j=1:length(mut_rna)) %do% mean(wt_rna[[j]]))
  dt$T_Mean[gl]=mRNA[dt$Gene[gl],"T_Mean"]
  dt$N_Mean[gl]=mRNA[dt$Gene[gl],"N_Mean"]
  dt$p.value[gl]=NA
  gcl=which(n_sample>1&(dt$T_Mean[gl]>=1.5|dt$N_Mean[gl]>=1.5))
  dt$p.value[gl[gcl]]=foreach(j=gcl,.combine=c) %do% t.test(mut_rna[[j]],wt_rna[[j]])$p.value
  dt$q.value[gl[gcl]]=p.adjust(dt$p.value[gl[gcl]],method="BH")
   write.table(dt,sub("cgc_msg","cgc_smg_associate",file[i]),sep="\t",row.names=F,quote=F)
}
#
dir.create("figure1")
file.copy("hotspot_promoterCore_cgc_smg_associate.xls","figure1/hotspot_promoterCore_cgc_smg_associate.xls")

library(foreach)
setwd("figure1")
file=dir()
file=file[grep("associate.xls",file)]

mRNA=read.table("../DESeq2_normalizedCount_134.xls",header=T,stringsAsFactors = F)
  mRNA=mRNA[,-which(is.element(colnames(mRNA),c("T502","N502")))]
  t_sp=colnames(mRNA)[grep("T",colnames(mRNA))]
  mRNA$T_Mean=rowMeans(mRNA[,grep("T",colnames(mRNA))])
  mRNA$N_Mean=rowMeans(mRNA[,grep("N",colnames(mRNA))])

ht=foreach(i=1:length(file)) %do% read.table(file[i],sep="\t",header=T,stringsAsFactors = F)
for(i in 1:length(ht)){
  dt=ht[[i]]
  sp=foreach(j=1:nrow(dt)) %do% intersect(strsplit(gsub("_WGS","",dt$Sample[j]),split=",")[[1]],t_sp)
  n_sample=foreach(j=1:nrow(dt),.combine=c) %do% length(sp[[j]])
  mut_rna=foreach(j=1:nrow(dt)) %do% as.numeric(mRNA[dt$Gene[j],sp[[j]]])
  wt_rna=foreach(j=1:nrow(dt)) %do% as.numeric(mRNA[dt$Gene[j],setdiff(t_sp,sp[[j]])])
  dt$MUT_max=unlist(foreach(j=1:nrow(dt)) %do% max(mut_rna[[j]]))
  dt$WT_max=unlist(foreach(j=1:nrow(dt)) %do% max(wt_rna[[j]]))
  dt$MUT_min=unlist(foreach(j=1:nrow(dt)) %do% min(mut_rna[[j]]))
  dt$WT_min=unlist(foreach(j=1:nrow(dt)) %do% min(wt_rna[[j]]))
  dt$MUT_WT=""
  gl=which(dt$T_Mean>=1.5|dt$N_Mean>=1.5)
  dt$MUT_WT[gl]=foreach(j=gl,.combine=c) %do% ifelse(dt$MUT_min[gl]>dt$WT_max[gl],"+",ifelse(dt$MUT_max[gl]<dt$WT_min[gl],"-",""))
  #write.table(dt,"hotspot_promoter2.2k_assocaite_gene.xls",sep="\t",row.names=F,quote=F)
  write.table(dt,sub("cgc_smg_associate","cgc_smg_associate_min_max",file[i]),sep="\t",row.names=F,quote=F)
}

file=dir()
file=file[grep("associate_min_max.xls",file)]

mRNA=read.table("../Annotated_lncRNA_FPKM_134.xls",sep="\t",header = T,stringsAsFactors = F)
 colnames(mRNA)=sub("_WTS","",colnames(mRNA))
mRNA=mRNA[,-which(is.element(colnames(mRNA),c("T502","N502","T13","N13")))]
t_sp=colnames(mRNA)[-c(1:4)][grep("T",colnames(mRNA)[-c(1:4)])]
mRNA$T_Mean=rowMeans(mRNA[,t_sp])
mRNA$N_Mean=rowMeans(mRNA[,sub("T","N",t_sp)])

muti=names(table(mRNA$Official_Symbol))[table(mRNA$Official_Symbol)!=1]
for(i in 1:length(muti)){
  gl=which(mRNA$Official_Symbol==muti[i])
  sgl=which.max(rowMeans(mRNA[gl,c("T_Mean","N_Mean")]))
  mRNA=mRNA[-gl[-sgl],]
}
rownames(mRNA)=mRNA$Official_Symbol
mRNA=mRNA[,-c(1:4)]

ht=foreach(i=1:length(file)) %do% read.table(file[i],sep="\t",header=T,stringsAsFactors = F)
for(i in 1:length(ht)){
  dt=ht[[i]]
  lnc=intersect(dt$Gene,rownames(mRNA))
  gl=which(is.element(dt$Gene,lnc))
  sp=foreach(j=gl) %do% intersect(strsplit(gsub("_WGS","",dt$Sample[j]),split=",")[[1]],t_sp)
  n_sample=foreach(j=1:length(sp),.combine=c) %do% length(sp[[j]])
  mut_rna=foreach(j=1:length(sp)) %do% as.numeric(mRNA[dt$Gene[gl[j]],sp[[j]]])
  mut_rna[[which(n_sample==0)]]=NA
  wt_rna=foreach(j=1:length(sp)) %do% as.numeric(mRNA[dt$Gene[gl[j]],setdiff(t_sp,sp[[j]])])
  dt$MUT_max[gl]=unlist(foreach(j=1:length(mut_rna)) %do% max(mut_rna[[j]]))
  dt$WT_max[gl]=unlist(foreach(j=1:length(mut_rna)) %do% max(wt_rna[[j]]))
  dt$MUT_min[gl]=unlist(foreach(j=1:length(mut_rna)) %do% min(mut_rna[[j]]))
  dt$WT_min[gl]=unlist(foreach(j=1:length(mut_rna)) %do% min(wt_rna[[j]]))
  gcl=which((dt$T_Mean[gl]>=1.5|dt$N_Mean[gl]>=1.5)&n_sample!=0)
  dt$MUT_WT[gl[gcl]]=ifelse(dt$MUT_min[gl[gcl]]>dt$WT_max[gl[gcl]],"+",ifelse(dt$MUT_max[gl[gcl]]<dt$WT_min[gl[gcl]],"-",""))
  dt$MUT_mRNA[grep("NA",dt$MUT_mRNA)]=""
  dt$WT_mRNA[grep("NA",dt$WT_mRNA)]=""
  dt$MUT_min[dt$MUT_min==NA|dt$MUT_min==Inf]=""
  dt$MUT_max[dt$MUT_max==NA|dt$MUT_max==Inf]=""
  dt[is.na(dt)]=""
  write.table(dt,file[i],sep="\t",row.names=F,quote=F)
}

#function annotation
file=dir()
file=file[grep("hotspot_promoterCore_cgc_smg_associate_min_max",file)]
i=1
dt=read.table(file[i],sep="\t",header = T,stringsAsFactors = F)
dt=unique(dt)

func=read.table("../CancerGenesList.txt",sep="\t",header=T,stringsAsFactors=F,row.names=1)
func$OncoKB.OG[func$OncoKB.OG==""]=NA
func$OncoKB.TSG[func$OncoKB.TSG==""]=NA
func$OG_TSG=foreach(i=1:nrow(func),.combine = c) %do% paste(na.omit(c(func$OncoKB.OG[i],func$OncoKB.TSG[i])),collapse=",")
cgc_smg=read.table("../CGC_SMG_gene.xls",sep="\t",header=T,row.names=1,stringsAsFactors=F)
ccg <- as.data.frame(readr::read_csv("../Census_allTue Jan  2 12_17_09 2018.csv"))
cgc_gene=ccg$`Gene Symbol`
rownames(ccg)=ccg$`Gene Symbol`

pathw=as.data.frame(readxl::read_excel("../Pathway+gene+oncogene+TSG+curated2.xlsx")[,1:3])
colnames(pathw)[1]="Gene"
rownames(pathw)=pathw$Gene

dt$OncoKB.Annotated=""
dt[is.element(dt$Gene,rownames(func)),"OncoKB.Annotated"]=func[dt$Gene[is.element(dt$Gene,rownames(func))],2]
dt$OncoKB_OG_TSG=""
dt[is.element(dt$Gene,rownames(func)),"OncoKB_OG_TSG"]=func[dt$Gene[is.element(dt$Gene,rownames(func))],"OG_TSG"]

#CGC　anotation
dt$CGC_OG_TSG=""
dt[is.element(dt$Gene,rownames(ccg)),"CGC_OG_TSG"]=ccg[dt$Gene[is.element(dt$Gene,rownames(ccg))],"Role in Cancer"]

dt$pathwaylist=""
dt$OG_TSG=""
pthl=which(is.element(dt$Gene,unique(pathw$Gene)))
dt[pthl,c("OG_TSG","pathwaylist")]=pathw[dt$Gene[pthl],2:3]
write.table(dt,"hotspot_promoterCore_cgc_smg_associate_min_max_oncokb.xls",sep="\t",row.names=F,quote=F)
