suppressPackageStartupMessages({library(edgeR);library(limma);library(readxl);library(jsonlite)})
base <- normalizePath('audit_work_2026-10-06'); dl <- file.path(base,'downloads'); out <- file.path(base,'outputs'); summaries <- list(); matrices <- list(); annotations <- list(); models <- list()
qc <- function(label,counts,meta,is_array=FALSE){
 stopifnot(identical(colnames(counts),meta$sample), all(is.finite(counts)))
 raw <- counts; dup <- sum(duplicated(rownames(counts)))
 if(dup>0) stop(paste('Duplicate primary IDs',label))
 if(is_array){logmat<-counts;keep<-rep(TRUE,nrow(counts));y<-NULL} else {
  y<-DGEList(counts=counts);keep<-filterByExpr(y,group=meta$group);y<-y[keep,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);logmat<-cpm(y,log=TRUE,prior.count=2)
 }
 vars<-apply(logmat,1,var);pc<-prcomp(t(logmat[order(vars,decreasing=TRUE)[seq_len(min(2000,nrow(logmat)))],,drop=FALSE]),center=TRUE,scale.=FALSE)
 corrs<-cor(logmat); hc<-hclust(dist(t(logmat))); colors<-as.integer(factor(meta$group));pal<-c('#3B82F6','#F97316','#16A34A','#A855F7','#DC2626','#14B8A6')
 png(file.path(out,paste0(label,'_QC.png')),width=1800,height=1400,res=160);par(mfrow=c(2,2),mar=c(6,4,3,1))
 plot(pc$x[,1:2],col=pal[colors],pch=19,xlab=sprintf('PC1 (%.1f%%)',100*pc$sdev[1]^2/sum(pc$sdev^2)),ylab=sprintf('PC2 (%.1f%%)',100*pc$sdev[2]^2/sum(pc$sdev^2)),main=label);text(pc$x[,1:2],labels=sub('MHU-8_RTN_','',meta$sample),cex=.55,pos=3);legend('topright',legend=unique(meta$group),col=pal[match(unique(meta$group),levels(factor(meta$group)))],pch=19,cex=.65)
 boxplot(logmat,names=sub('MHU-8_RTN_','',colnames(logmat)),las=2,cex.axis=.55,col=pal[colors],main=if(is_array)'log2 RMA intensity' else 'TMM log2 CPM')
 image(seq_len(ncol(corrs)),seq_len(ncol(corrs)),corrs,col=hcl.colors(50,'YlOrRd',rev=TRUE),xlab='Sample index',ylab='Sample index',main=sprintf('Pearson correlation %.3f–%.3f',min(corrs[lower.tri(corrs)]),max(corrs[lower.tri(corrs)])))
 plot(hc,main='Unsupervised sample clustering',cex=.5,xlab='',sub='');dev.off()
 write.csv(data.frame(sample=meta$sample,group=meta$group,pc$x[,1:4,drop=FALSE]),file.path(out,paste0(label,'_PCA.csv')),row.names=FALSE)
 write.csv(corrs,file.path(out,paste0(label,'_correlation.csv')))
 nn<-sapply(seq_len(ncol(corrs)),function(i){v<-corrs[,i];v[i]<-NA;meta$sample[which.max(v)]}); ng<-meta$group[match(nn,meta$sample)]
 sm<-list(dim=dim(raw),missing=sum(is.na(raw)),duplicate_ids=dup,noninteger_fraction=mean(abs(raw-round(raw))>1e-7),raw_range=range(raw),library_size=if(!is_array)setNames(colSums(raw),colnames(raw)) else NULL,expressed=sum(keep),correlation_range=range(corrs[lower.tri(corrs)]),pc_variance=100*pc$sdev[1:4]^2/sum(pc$sdev^2),group=as.list(table(meta$group)),nearest_same_group=mean(ng==meta$group),sample_medians=setNames(apply(logmat,2,median),meta$sample))
 summaries[[label]]<<-sm;matrices[[label]]<<-rownames(raw)[keep]
 list(y=y,logmat=logmat,meta=meta,keep=keep)
}
fitcounts<-function(label,obj,design,contrasts){
 y<-estimateDisp(obj$y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 for(n in names(contrasts)){
  co<-makeContrasts(contrasts=contrasts[[n]],levels=design);res<-topTags(glmQLFTest(fit,contrast=co),n=Inf,sort.by='none')$table
  res$id<-rownames(res);write.csv(res,file.path(out,paste0(label,'_',n,'_DE.csv')),row.names=FALSE);models[[paste(label,n,sep='_')]]<<-res
  summaries[[label]][[n]]<<-list(tested=nrow(res),FDR05=sum(res$FDR<.05),FDR05_LFC1=sum(res$FDR<.05 & abs(res$logFC)>1),residual_df=nrow(design)-qr(design)$rank,top=res[order(res$FDR)[1:8],c('id','logFC','FDR')])
 }
}
# Original published spaceflight quantification; independent biological right retinas.
sc<-as.matrix(read.csv(file.path(dl,'GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv'),row.names=1,check.names=FALSE));sm<-read.csv(file.path(dl,'GLDS-664_rna-seq_SampleTable.csv'),row.names=1)
sm<-data.frame(sample=rownames(sm),group=make.names(sm$condition));sm<-sm[match(colnames(sc),sm$sample),];so<-qc('OSD758_original',sc,sm);group<-factor(sm$group);de<-model.matrix(~0+group);colnames(de)<-levels(group)
fitcounts('OSD758_original',so,de,list(MG_vs_HGC='FLT_0G-HGC',AG1_vs_HGC='FLT_1G-HGC',MG_vs_AG1='FLT_0G-FLT_1G',AG033_vs_HGC='FLT_0.33G-HGC',AG067_vs_HGC='FLT_0.67G-HGC'))
# New processing: keep original biological labels, current STAR counts only as sensitivity.
sc2<-as.matrix(read.csv(file.path(dl,'GLDS-664_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv'),row.names=1,check.names=FALSE));so2<-qc('OSD758_STAR',sc2,sm);fitcounts('OSD758_STAR',so2,de,list(MG_vs_HGC='FLT_0G-HGC'))
# P23H P90 published read counts, not FPKM; sex adjustment.
pf<-as.data.frame(read_excel(file.path(dl,'GSE156533_readcounts_organized_by_largest_to_smallest_count.xlsx')));cols<-c('P23H_6','WT_7','P23H_8','WT_9','P23H_11','WT_12','P23H_13');pc<-as.matrix(pf[,cols]);rownames(pc)<-pf$gene_id
pm<-data.frame(sample=cols,group=ifelse(grepl('P23H',cols),'P23H','WT'),sex=c('F','F','F','F','M','M','M'));po<-qc('GSE156533',pc,pm);group<-factor(pm$group,levels=c('WT','P23H'));sex<-factor(pm$sex);de<-model.matrix(~sex+group);colnames(de)<-make.names(colnames(de));fitcounts('GSE156533',po,de,list(P23H_vs_WT='groupP23H'))
annotations[['P23H']]<-data.frame(id=pf$gene_id,symbol=pf$gene_name)
# Glaucoma and optic nerve crush: retina only, never treating the other tissues as replicates.
gf<-read.csv(file.path(dl,'GSE241782_gene_count_matrix.csv.gz'),check.names=FALSE);gids<-sub('\\|.*','',gf$gene_id);gsym<-sub('^[^|]*\\|','',gf$gene_id)
gm<-read.csv(file.path(out,'GSE241782_samples.csv'),check.names=FALSE);gm<-gm[gm$tissue=='Retina',];gm$sample<-sub('.*_(R-[0-9]+)_genes.*','\\1',gm$supplementary);gm$group<-ifelse(gm$treatment=='Naïve','Naive',paste0(gm$treatment,'_',gsub(' ','',gm$time)));gm$sex<-gm$sex
# Some StringTie novel IDs share ENSMUSG identifiers? Primary ID check before collapse.
gc<-as.matrix(gf[,gm$sample]);rownames(gc)<-gf$gene_id;go<-qc('GSE241782_retina',gc,gm);group<-factor(gm$group);sex<-factor(gm$sex);de<-model.matrix(~0+group+sex);colnames(de)<-sub('group','',colnames(de));fitcounts('GSE241782_retina',go,de,list(Glaucoma3D='Glaucoma_3days-Naive',Glaucoma2W='Glaucoma_2weeks-Naive',Glaucoma6W='Glaucoma_6weeks-Naive',Crush3D='Crush_3days-Naive',Crush2W='Crush_2weeks-Naive'))
annotations[['Glaucoma']]<-data.frame(id=gf$gene_id,ensembl=gids,symbol=gsym)
# EAU: exclude sleep-deprived groups from the disease contrast.
ef<-read.delim(file.path(dl,'GSE284505_R.HTSeq.counts.txt.gz'),check.names=FALSE);ec<-as.matrix(ef[,-1]);rownames(ec)<-ef[,1];em<-data.frame(sample=colnames(ec),group=sub('[123]$','',colnames(ec)));eo<-qc('GSE284505_all',ec,em);sel<-em$group %in% c('CON','EAU');eobj<-qc('GSE284505_main',ec[,sel],em[sel,]);group<-factor(em$group[sel],levels=c('CON','EAU'));de<-model.matrix(~group);fitcounts('GSE284505_main',eobj,de,list(EAU_vs_CON='groupEAU'))
# Light injury: QC only until independence, age, sex, dose and time are verified.
lf<-read.delim(file.path(dl,'GSE314043_genecount.txt.gz'),check.names=FALSE);lc<-as.matrix(lf[,-1]);rownames(lc)<-lf[,1];lm<-data.frame(sample=colnames(lc),group=ifelse(grepl('^LD',colnames(lc)),'LD','CTRL'));lo<-qc('GSE314043',lc,lm)
# Arrays: native log2 RMA matrix; no count model, no cross-study expression merging.
readmatrix<-function(file){lines<-readLines(gzfile(file));a<-which(lines=='!series_matrix_table_begin');b<-which(lines=='!series_matrix_table_end');d<-read.delim(text=paste(lines[(a+1):(b-1)],collapse='\n'),check.names=FALSE);m<-as.matrix(d[,-1]);rownames(m)<-d[,1];m}
dc<-readmatrix(file.path(dl,'GSE87433_series_matrix.txt.gz'));dm<-read.csv(file.path(out,'GSE87433_samples.csv'),check.names=FALSE);dm$sample<-dm$accession;dm$group<-ifelse(dm$disease=='normoglycemic','NC',ifelse(dm[['beta cells transplantation']]=='yes','Tx','DC'));dm<-dm[match(colnames(dc),dm$sample),];do<-qc('GSE87433_all',dc,dm,is_array=TRUE);sel<-dm$group!='Tx';do2<-qc('GSE87433_main',dc[,sel],dm[sel,],is_array=TRUE);group<-factor(dm$group[sel],levels=c('NC','DC'));de<-model.matrix(~group);fit<-eBayes(lmFit(do2$logmat,de),trend=TRUE,robust=TRUE);dr<-topTable(fit,coef='groupDC',number=Inf,sort.by='none');dr$id<-rownames(dr);write.csv(dr,file.path(out,'GSE87433_DC_vs_NC_DE.csv'),row.names=FALSE);models[['GSE87433_DC_vs_NC']]<-dr;summaries[['GSE87433_main']]$DC_vs_NC<-list(tested=nrow(dr),FDR05=sum(dr$adj.P.Val<.05),FDR05_LFC1=sum(dr$adj.P.Val<.05 & abs(dr$logFC)>1),top=dr[order(dr$adj.P.Val)[1:8],c('id','logFC','adj.P.Val')])
# Annotation audit and conservative Ensembl overlap through supplied annotations.
dp<-read.delim(file.path(out,'GSE87433_platform.tsv'),check.names=FALSE);ens<-regmatches(dp$dbXrefs,regexpr('ENSMUSG[0-9]+',dp$dbXrefs));dp$ensembl<-ifelse(grepl('ENSMUSG',dp$dbXrefs),sapply(strsplit(dp$dbXrefs,'\\|'),function(v){z<-v[grepl('^Ensembl:ENSMUSG',v)];if(length(z))sub('Ensembl:','',z[1]) else NA}),NA)
spaceids<-matrices[['OSD758_original']]
getens<-function(ids){sub('\\|.*','',ids)}
annot<-list(P23H=list(named=sum(!is.na(pf$gene_name)&nzchar(pf$gene_name)),total=nrow(pf)),Glaucoma=list(ensembl=sum(grepl('^ENSMUSG',gids)),novel=sum(!grepl('^ENSMUSG',gids)),total=nrow(gf)),EAU=list(unique_symbols=nrow(ef),mapped_to_P23H=sum(ef[,1] %in% pf$gene_name)),Diabetes=list(platform_mapped=sum(rownames(dc)%in%as.character(dp$ID)),total=nrow(dc),ensembl_mapped=sum(!is.na(dp$ensembl))),Light=list(ensembl=sum(grepl('^ENSMUSG',rownames(lc))),total=nrow(lc)))
# EAU map only uniquely assigned symbols; do not assume ambiguous genes are independent.
tab<-table(pf$gene_name);valid<-!is.na(pf$gene_name) & tab[pf$gene_name]==1;symbolmap<-setNames(pf$gene_id[valid],pf$gene_name[valid]);eids<-unname(symbolmap[matrices[['GSE284505_main']]])
ov<-list(P23H=length(intersect(spaceids,matrices[['GSE156533']])),Glaucoma=length(intersect(spaceids,getens(matrices[['GSE241782_retina']]))),EAU=length(intersect(spaceids,eids)),Diabetes=length(intersect(spaceids,dp$ensembl[match(rownames(dc),as.character(dp$ID))])),Light=length(intersect(spaceids,matrices[['GSE314043']])))
write_json(list(qc=summaries,annotation=annot,expressed_overlap_with_OSD758=ov),file.path(out,'preliminary_results.json'),pretty=TRUE,auto_unbox=TRUE,na='null');saveRDS(list(summaries=summaries,models=models,matrices=matrices,annotation=annotations),file.path(out,'analysis_objects.rds'));cat(toJSON(list(qc=summaries,annotation=annot,overlap=ov),pretty=TRUE,auto_unbox=TRUE),sep='\n')
