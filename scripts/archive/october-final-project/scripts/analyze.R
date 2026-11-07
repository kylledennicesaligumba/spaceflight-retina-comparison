options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(edgeR);library(limma);library(readxl);library(jsonlite)})
dir.create('outputs',showWarnings=FALSE)
D <- 'downloads'; O <- 'outputs'; results <- list(); objects <- list()
readmeta <- function(a) fromJSON(file.path(O,paste0(a,'_metadata.json')),simplifyVector=FALSE)$samples
chars <- function(s) {
 x<-unlist(s$Sample_characteristics_ch1); k<-sub(':.*','',x); v<-sub('^[^:]+: *','',x)
 setNames(v,tolower(k))
}
writecsv <- function(x,n) write.csv(x,file.path(O,n),row.names=FALSE)
readcounts <- function(file,sep=',') {
 x<-read.table(file.path(D,file),sep=sep,header=TRUE,check.names=FALSE,quote='"',comment.char='')
 mat<-as.matrix(x[,-1,drop=FALSE]);storage.mode(mat)<-'numeric';rownames(mat)<-as.character(x[[1]]);mat
}
readseries <- function(f) {
 x<-readLines(gzfile(file.path(D,f)),warn=FALSE)
 b<-which(x=='!series_matrix_table_begin');e<-which(x=='!series_matrix_table_end')
 z<-read.delim(text=paste(x[(b+1):(e-1)],collapse='\n'),check.names=FALSE,quote='"')
 mat<-as.matrix(z[,-1,drop=FALSE]);storage.mode(mat)<-'numeric';rownames(mat)<-as.character(z[[1]]);mat
}
qc <- function(mat,meta,name,counts=TRUE,units=NULL) {
 stopifnot(identical(colnames(mat),meta$sample),all(is.finite(mat)),!anyDuplicated(rownames(mat)),nrow(meta)==ncol(mat))
 z<-list(rows=nrow(mat),columns=ncol(mat),missing=sum(is.na(mat)),duplicate_ids=sum(duplicated(rownames(mat))),
  fractional_entry_fraction=if(counts)mean(abs(mat-round(mat))>1e-6) else NA,
  groups=as.list(table(meta$group)),raw_value_range=range(mat))
 if(counts) {
  y<-DGEList(mat); k<-filterByExpr(y,group=meta$group); y<-calcNormFactors(y[k,,keep.lib.sizes=FALSE]);expr<-cpm(y,log=TRUE,prior.count=2)
  z$raw_library_totals<-as.list(setNames(colSums(mat),colnames(mat)))
  z$TMM_factors<-as.list(setNames(y$samples$norm.factors,colnames(mat)))
 } else {expr<-mat;y<-NULL;k<-rep(TRUE,nrow(mat))}
 co<-cor(expr);v<-apply(expr,1,var);top<-order(v,decreasing=TRUE)[seq_len(min(2000,sum(v>0)))]
 p<-prcomp(t(expr[top,,drop=FALSE]),center=TRUE,scale.=FALSE);pct<-p$sdev^2/sum(p$sdev^2)*100
 z$tested_features<-nrow(expr);z$correlation_range<-range(co[upper.tri(co)]);z$PC_variance_percent<-head(pct,5)
 g<-factor(meta$group); colors<-setNames(hcl.colors(nlevels(g),'Dark 3'),levels(g));cl<-colors[as.character(g)]
 labels<-if('label'%in%names(meta))meta$label else meta$sample
 png(file.path(O,paste0(name,'_QC.png')),width=1900,height=1500,res=150)
 par(mfrow=c(2,2),mar=c(7,4,3,1))
 boxplot(expr,names=labels,las=2,col=cl,outline=FALSE,cex.axis=.65,main=paste(name,if(counts)'TMM log2 CPM' else if(is.null(units))'processed log2 intensities' else units),ylab='log2 expression')
 par(mar=c(5,4,3,1))
 xr<-range(p$x[,1]);yr<-range(p$x[,2]);xr<-xr+c(-1,1)*diff(xr)*.25;yr<-yr+c(-1,1)*diff(yr)*.25
 plot(p$x[,1:2],col=cl,pch=if('sex'%in%names(meta))ifelse(meta$sex=='Male',17,16) else 16,
      xlim=xr,ylim=yr,xlab=sprintf('PC1 (%.1f%%)',pct[1]),ylab=sprintf('PC2 (%.1f%%)',pct[2]),main='PCA: 2,000 most variable features')
 text(p$x[,1],p$x[,2],labels=labels,pos=3,cex=.55)
 legend('topright',legend=levels(g),col=colors,pch=16,cex=.60,bty='n')
 par(mar=c(7,7,3,1));ord<-hclust(as.dist(1-co),method='average')$order
 image(seq_len(ncol(expr)),seq_len(ncol(expr)),co[ord,ord],col=hcl.colors(50,'YlOrRd',rev=TRUE),axes=FALSE,xlab='',ylab='',main=sprintf('Pearson correlation (%.3f–%.3f)',min(co[upper.tri(co)]),max(co[upper.tri(co)])))
 axis(1,at=seq_along(ord),labels=labels[ord],las=2,cex.axis=.55);axis(2,at=seq_along(ord),labels=labels[ord],las=2,cex.axis=.55)
 par(mar=c(7,4,3,1));h<-hclust(as.dist(1-co),method='average');h$labels<-labels
 plot(h,main='Clustering: 1 – Pearson correlation',xlab='',sub='',cex=.6)
 dev.off()
 writecsv(meta,paste0(name,'_analysis_samples.csv'));write.csv(co,file.path(O,paste0(name,'_correlations.csv')))
 writecsv(cbind(meta,p$x[,1:3]),paste0(name,'_PCA.csv'))
 results[[name]]<<-z;objects[[name]]<<-list(raw=mat,meta=meta,y=y,expr=expr,keep=k)
 cat(name,nrow(mat),'features',ncol(mat),'samples',nrow(expr),'tested\n')
 invisible(objects[[name]])
}
fitcounts <- function(name,design,contrasts) {
 z<-objects[[name]];y<-estimateDisp(z$y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
 results[[name]]$design_columns<<-colnames(design);results[[name]]$residual_df<<-nrow(design)-qr(design)$rank
 out<-list()
 for(n in colnames(contrasts)) {
  tt<-topTags(glmQLFTest(fit,contrast=contrasts[,n]),n=Inf,sort.by='none')$table
  tt$gene_id<-rownames(tt);tt<-tt[,c('gene_id',setdiff(names(tt),'gene_id'))]
  writecsv(tt,paste0(name,'_',n,'_DE.csv'))
  out[[n]]<-list(BH05=sum(tt$FDR<.05),BH05_abslog2FC1=sum(tt$FDR<.05 & abs(tt$logFC)>1),
   top=tt[order(tt$FDR),][1:5,])
  objects[[name]][[n]]<<-tt
 }
 results[[name]]$DE<<-out;objects[[name]]$design<<-design
}

# Spaceflight: retain every available retina sample and all artificial-gravity groups.
space <- readcounts('GLDS-664_rna-seq_RSEM_Unnormalized_Counts.csv')
sm <- read.csv(file.path(D,'GLDS-664_rna-seq_MHU-8_RTN_metadata.csv'))
sm<-sm[match(colnames(space),sm$sample),];sm$group<-c('AG033','AG067','MG','AG1','HGC')[match(sm$condition,c('FLT_0.33G','FLT_0.67G','FLT_0G','FLT_1G','HGC'))]
sm$label<-sub('.*_(FLT[0-9]+|HGC[0-9]+)$','\\1',sm$sample)
writecsv(data.frame(gene_id=rownames(space)[!grepl('^ENSMUSG',rownames(space))]),'OSD758_nonmouse_features.csv')
space<-space[grepl('^ENSMUSG',rownames(space)),,drop=FALSE]
qc(space,sm,'OSD758_RSEM')
design<-model.matrix(~0+group,data=sm);colnames(design)<-sub('group','',colnames(design))
ct<-makeContrasts(MG_HGC=MG-HGC,AG1_HGC=AG1-HGC,MG_AG1=MG-AG1,AG033_HGC=AG033-HGC,AG067_HGC=AG067-HGC,levels=design)
fitcounts('OSD758_RSEM',design,ct)
star<-readcounts('GLDS-664_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv');star<-star[grepl('^ENSMUSG',rownames(star)),sm$sample,drop=FALSE]
qc(star,sm,'OSD758_STAR');fitcounts('OSD758_STAR',design,ct[,'MG_HGC',drop=FALSE])

# Adult P23H, individual samples with sex as a covariate.
rp<-as.data.frame(read_excel(file.path(D,'GSE156533_counts.xlsx')))
rpcols<-grep('^(P23H|WT)_',names(rp),value=TRUE)
x<-as.matrix(rp[,rpcols]);rownames(x)<-rp$gene_id
meta<-data.frame(sample=rpcols,group=ifelse(grepl('^P23H',rpcols),'P23H','WT'),
 sex=ifelse(rpcols%in%c('P23H_11','P23H_13','WT_12'),'Male','Female'),label=rpcols)
meta$genotype<-factor(meta$group,levels=c('WT','P23H'))
qc(x,meta,'GSE156533');des<-model.matrix(~sex+genotype,data=meta)
ct<-matrix(0,ncol(des),1,dimnames=list(colnames(des),'P23H_WT'));ct['genotypeP23H',1]<-1
fitcounts('GSE156533',des,ct)
writecsv(rp[,c('gene_id','gene_name','gene_biotype')],'P23H_gene_annotation.csv')

# Glaucoma and crush: only retina, one male pool and one female pool per condition.
gc<-readcounts('GSE241782_gene_count_matrix.csv.gz');gc<-gc[,grep('^R-',colnames(gc)),drop=FALSE]
symbol<-sub('^[^|]*\\|','',rownames(gc));rownames(gc)<-sub('\\|.*','',rownames(gc))
stopifnot(!anyDuplicated(rownames(gc)))
gm<-readmeta('GSE241782');rows<-lapply(names(gm),function(gsm){
 s<-gm[[gsm]];ch<-chars(s);if(!identical(unname(ch['tissue']),'Retina'))return(NULL)
 file<-unlist(s$Sample_supplementary_file_1);id<-regmatches(file,regexpr('R-[0-9]+',file))
 gr<-if(ch['treatment']=='Naïve')'Naive' else paste0(ch['treatment'],ifelse(ch['time']=='3 days','3d',ifelse(ch['time']=='2 weeks','2w','6w')))
 data.frame(sample=id,accession=gsm,group=unname(gr),sex=unname(ch['sex']),label=paste0(id,ifelse(ch['sex']=='Male','M','F')))
});gm<-do.call(rbind,rows);gm<-gm[match(colnames(gc),gm$sample),]
qc(gc,gm,'GSE241782_retina');des<-model.matrix(~0+group+sex,data=gm);colnames(des)<-sub('group','',colnames(des))
ct<-makeContrasts(Glaucoma3d_Naive=Glaucoma3d-Naive,Glaucoma2w_Naive=Glaucoma2w-Naive,Glaucoma6w_Naive=Glaucoma6w-Naive,Crush3d_Naive=Crush3d-Naive,Crush2w_Naive=Crush2w-Naive,levels=des)
fitcounts('GSE241782_retina',des,ct)
writecsv(data.frame(gene_id=rownames(gc),symbol=symbol),'GSE241782_gene_annotation.csv')

# Uveitis: exploratory inference explicitly conditional on independent libraries.
ea<-readcounts('GSE284505_counts.txt.gz',sep='\t')
em<-data.frame(sample=colnames(ea),group=sub('[0-9]+$','',colnames(ea)),label=colnames(ea))
qc(ea,em,'GSE284505_all')
sel<-em$group%in%c('CON','EAU');em<-em[sel,];em$condition<-factor(em$group,levels=c('CON','EAU'))
qc(ea[,sel],em,'GSE284505_main');des<-model.matrix(~condition,data=em)
ct<-matrix(c(0,1),2,1,dimnames=list(colnames(des),'EAU_CON'));fitcounts('GSE284505_main',des,ct)
results$GSE284505_main$inference_condition<-'Assumes each deposited library is an independent animal or independent pool; pooling is not established.'

# Light injury: processed-data QC only until its design is verified.
li<-readcounts('GSE314043_counts.txt.gz',sep='\t')
lm<-data.frame(sample=colnames(li),group=ifelse(grepl('^LD',colnames(li)),'LD','CTRL'),label=colnames(li))
qc(li,lm,'GSE314043');results$GSE314043$DE_status<-'Not run: age, sex, timing and independence unresolved.'

# Diabetes: already RMA-normalized log2 array values; no second log or RNA-seq model.
di<-readseries('GSE87433_series_matrix.txt.gz');dm<-readmeta('GSE87433')
mt<-do.call(rbind,lapply(colnames(di),function(gsm){ch<-chars(dm[[gsm]]);data.frame(sample=gsm,
 group=if(ch['beta cells transplantation']=='yes')'Tx' else if(ch['disease']=='normoglycemic')'NC' else 'DC',label=sub('GSM233','',gsm))}))
qc(di,mt,'GSE87433_all',FALSE);sel<-mt$group!='Tx';mt<-mt[sel,];mt$condition<-factor(mt$group,levels=c('NC','DC'))
qc(di[,sel],mt,'GSE87433_main',FALSE);des<-model.matrix(~condition,data=mt)
fit<-eBayes(lmFit(di[,sel],des),trend=TRUE,robust=TRUE);tt<-topTable(fit,coef='conditionDC',number=Inf,sort.by='none')
tt$gene_id<-rownames(tt);writecsv(tt,'GSE87433_main_DC_NC_DE.csv');objects$GSE87433_main$DC_NC<-tt
results$GSE87433_main$DE<-list(DC_NC=list(BH05=sum(tt$adj.P.Val<.05),BH05_abslog2FC1=sum(tt$adj.P.Val<.05 & abs(tt$logFC)>1),top=tt[order(tt$adj.P.Val),][1:5,]))
results$GSE87433_main$residual_df<-10

# Identifier coverage and expressed overlap are kept distinct from array detectability.
spgenes<-rownames(objects$OSD758_RSEM$expr)
for(n in c('GSE156533','GSE241782_retina','GSE314043'))results[[n]]$expressed_Ensembl_overlap<-length(intersect(spgenes,rownames(objects[[n]]$expr)))
sym<-rp$gene_name;valid<-!is.na(sym)&nzchar(sym)&sym!='-'&!duplicated(sym)&!duplicated(sym,fromLast=TRUE)
mapping<-setNames(rp$gene_id[valid],sym[valid]);eg<-rownames(objects$GSE284505_main$expr)
results$GSE284505_main$unique_symbol_mapping_coverage<-sum(!is.na(mapping[eg]))
results$GSE284505_main$expressed_Ensembl_overlap<-length(intersect(spgenes,unname(mapping[eg])))
results$GSE156533$nonmissing_symbol_annotations<-sum(!is.na(rp$gene_name)&rp$gene_name!='-'&nzchar(rp$gene_name))
gp<-read.delim(file.path(O,'GSE87433_platform.tsv'),check.names=FALSE,quote='',comment.char='',fill=TRUE)
writecsv(data.frame(platform_columns=names(gp)),'GSE87433_platform_columns.csv')
ic<-grep('^ID$',names(gp),value=TRUE);ec<-grep('Ensembl|ENSEMBL',names(gp),value=TRUE)[1]
if(is.na(ec))ec<-'dbXrefs'
ids<-as.character(gp[[ic]]);vals<-as.character(gp[[ec]]);allids<-rownames(di)
results$GSE87433_main$platform_ID_coverage<-sum(allids%in%ids)
mapped<-vals[match(allids,ids)];ens<-unlist(regmatches(mapped,gregexpr('ENSMUSG[0-9]+',mapped)))
results$GSE87433_main$array_assayed_overlap_not_expression_calls<-length(intersect(spgenes,ens))
results$GSE87433_main$features_with_Ensembl_annotation<-sum(grepl('ENSMUSG',mapped))
write_json(results,file.path(O,'preliminary_results.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
saveRDS(objects,file.path(O,'analysis_objects.rds'))
writeLines(capture.output(sessionInfo()),file.path(O,'R_session_info.txt'))
cat('Preliminary analyses completed\n')
