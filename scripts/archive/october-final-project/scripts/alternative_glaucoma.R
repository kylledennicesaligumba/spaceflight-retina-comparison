# Reuse the QC and I/O function definitions, without running the main analyses.
code<-readLines('scripts/analyze.R');eval(parse(text=code[seq_len(grep('^# Spaceflight:',code)-1)]))
arr<-readseries('GSE26299_series_matrix.txt.gz')
pm<-fromJSON('outputs/GSE26299_retina_pairing_inferred.json')
pm$group<-ifelse(grepl('D2-Gpnmb',pm$title),'Control',ifelse(grepl('Moderate',pm$title),'Moderate',ifelse(grepl('Severe',pm$title),'Severe',ifelse(grepl('early 1',pm$title),'NoEarly1','NoEarly2'))))
pm$animal_inferred<-as.character(pm$animal_inferred);pm$array_date_inferred<-factor(pm$array_date_inferred)
pm$label<-paste0(pm$eye,pm$animal_inferred)
pm<-pm[match(intersect(colnames(arr),pm$sample),pm$sample),];arr<-arr[,pm$sample,drop=FALSE]
qc(arr,pm,'GSE26299_retina',FALSE)
pm$group<-factor(pm$group,levels=c('Control','NoEarly1','NoEarly2','Moderate','Severe'))
des<-model.matrix(~0+group,data=pm);colnames(des)<-sub('group','',colnames(des))
stopifnot(qr(des)$rank==ncol(des))
du<-duplicateCorrelation(arr,des,block=pm$animal_inferred)
fit<-lmFit(arr,des,block=pm$animal_inferred,correlation=du$consensus.correlation)
fit<-eBayes(contrasts.fit(fit,makeContrasts(Moderate_Control=Moderate-Control,Severe_Control=Severe-Control,levels=des)),trend=TRUE,robust=TRUE)
dd<-list()
for(n in colnames(fit$coefficients)) {
 t<-topTable(fit,coef=n,number=Inf,sort.by='none');t$probe_id<-rownames(t)
 writecsv(t,paste0('GSE26299_retina_',n,'_DE.csv'))
 dd[[n]]<-list(BH05_probes=sum(t$adj.P.Val<.05),BH05_abslog2FC1_probes=sum(t$adj.P.Val<.05&abs(t$logFC)>1))
}
results$GSE26299_retina$DE<-dd
results$GSE26299_retina$inference_condition<-'Mouse blocks inferred from primary CEL filenames; confirm before final inference. Dates in filenames have an unresolved role and are not adjusted in the primary model. Results are probe-level.'
results$GSE26299_retina$block_correlation<-du$consensus.correlation
results$GSE26299_retina$unique_animals_inferred_per_group<-as.list(sapply(split(pm$animal_inferred,pm$group),function(x)length(unique(x))))
results$GSE26299_retina$filename_date_by_group_unresolved_role<-as.data.frame(table(pm$group,pm$array_date_inferred))
gp<-read.delim('outputs/GSE26299_platform.tsv',quote='',comment.char='',check.names=FALSE,fill=TRUE)
results$GSE26299_retina$platform_ID_coverage<-sum(rownames(arr)%in%gp$ID)
ent<-as.character(gp$ENTREZ_GENE_ID);one<-grepl('^[0-9]+$',ent)
results$GSE26299_retina$probes_with_unambiguous_Entrez<-sum(one)
results$GSE26299_retina$unique_Entrez<-length(unique(ent[one]))
rgp<-read.delim('outputs/GSE87433_platform.tsv',quote='',comment.char='',check.names=FALSE,fill=TRUE)
enscol<-grep('Ensembl|ENSEMBL',names(rgp),value=TRUE)[1]
if(is.na(enscol))enscol<-'dbXrefs'
vals<-rgp[[enscol]][match(unique(ent[one]),as.character(rgp$ID))]
ens<-unlist(regmatches(vals,gregexpr('ENSMUSG[0-9]+',vals)))
main<-readRDS('outputs/analysis_objects.rds')
results$GSE26299_retina$array_assayed_overlap_not_expression_calls<-length(intersect(rownames(main$OSD758_RSEM$expr),ens))
writecsv(pm,'GSE26299_retina_analysis_samples.csv')
write_json(results,'outputs/alternative_glaucoma_results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Alternative glaucoma analysis complete\n')
