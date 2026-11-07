# Fit glaucoma and optic nerve crush separately so dispersion and the sex effect
# need not be shared between biologically different injury models.
code<-readLines('scripts/analyze.R');eval(parse(text=code[seq_len(grep('^# Spaceflight:',code)-1)]))
results<-fromJSON('outputs/preliminary_results.json',simplifyVector=FALSE)
objects<-readRDS('outputs/analysis_objects.rds');gc<-objects$GSE241782_retina
for(injury in c('Glaucoma','Crush')) {
 sel<-gc$meta$group=='Naive'|grepl(paste0('^',injury),gc$meta$group)
 m<-gc$meta[sel,];n<-paste0('GSE241782_',tolower(injury))
 qc(gc$raw[,sel,drop=FALSE],m,n)
 d<-model.matrix(~0+group+sex,data=m);colnames(d)<-sub('group','',colnames(d))
 if(injury=='Glaucoma')ct<-makeContrasts(Glaucoma3d_Naive=Glaucoma3d-Naive,Glaucoma2w_Naive=Glaucoma2w-Naive,Glaucoma6w_Naive=Glaucoma6w-Naive,levels=d)
 else ct<-makeContrasts(Crush3d_Naive=Crush3d-Naive,Crush2w_Naive=Crush2w-Naive,levels=d)
 fitcounts(n,d,ct)
 results[[n]]$expressed_Ensembl_overlap<-length(intersect(rownames(objects$OSD758_RSEM$expr),rownames(objects[[n]]$expr)))
}
results$GSE241782_retina$DE_model_status<-'Combined model is a methodological sensitivity check. Main inference uses separate glaucoma and crush models.'
write_json(results,'outputs/preliminary_results.json',pretty=TRUE,auto_unbox=TRUE,na='null')
saveRDS(objects,'outputs/analysis_objects.rds')
cat('Separate models completed\n')
