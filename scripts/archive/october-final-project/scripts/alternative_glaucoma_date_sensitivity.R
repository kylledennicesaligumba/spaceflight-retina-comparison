suppressPackageStartupMessages({library(limma);library(jsonlite)})
code<-readLines('scripts/analyze.R');eval(parse(text=code[seq_len(grep('^# Spaceflight:',code)-1)]))
x<-readseries('GSE26299_series_matrix.txt.gz');m<-read.csv('outputs/GSE26299_retina_analysis_samples.csv')
x<-x[,m$sample,drop=FALSE];m$group<-factor(m$group);m$array_date_inferred<-factor(m$array_date_inferred)
d<-model.matrix(~0+group+array_date_inferred,data=m);colnames(d)<-sub('group','',colnames(d))
du<-duplicateCorrelation(x,d,block=m$animal_inferred)
fit<-eBayes(contrasts.fit(lmFit(x,d,block=m$animal_inferred,correlation=du$consensus.correlation),makeContrasts(Moderate_Control=Moderate-Control,Severe_Control=Severe-Control,levels=d)),trend=TRUE,robust=TRUE)
z<-list(block_correlation=du$consensus.correlation,interpretation='Dates are filename codes whose role is unverified; this is sensitivity analysis, not a verified batch correction.')
for(n in colnames(fit$coefficients)) {
 t<-topTable(fit,coef=n,number=Inf,sort.by='none');t$probe_id<-rownames(t);old<-read.csv(paste0('outputs/GSE26299_retina_',n,'_DE.csv'))
 write.csv(t,paste0('outputs/GSE26299_retina_',n,'_filename_date_sensitivity_DE.csv'),row.names=FALSE)
 i<-match(t$probe_id,old$probe_id);z[[n]]<-list(BH05_probes=sum(t$adj.P.Val<.05),BH05_abslog2FC1_probes=sum(t$adj.P.Val<.05&abs(t$logFC)>1),
  Spearman_logFC_with_unadjusted=cor(t$logFC,old$logFC[i],method='spearman'))
}
v<-apply(x,1,var);p<-prcomp(t(x[order(v,decreasing=TRUE)[1:2000],]),center=TRUE,scale.=FALSE)
a<-summary(lm(p$x[,1]~m$group));b<-summary(lm(p$x[,1]~m$array_date_inferred));c<-summary(lm(p$x[,1]~m$group+m$array_date_inferred))
z$PC1_R_squared<-list(severity_group=a$r.squared,filename_date=b$r.squared,both=c$r.squared)
write_json(z,'outputs/alternative_glaucoma_date_sensitivity.json',pretty=TRUE,auto_unbox=TRUE)
cat('Filename-date sensitivity complete\n')
