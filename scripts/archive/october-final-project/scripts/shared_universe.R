suppressPackageStartupMessages(library(jsonlite))
o<-readRDS('outputs/analysis_objects.rds');ann<-read.csv('outputs/P23H_gene_annotation.csv')
good<-!duplicated(ann$gene_name)&!duplicated(ann$gene_name,fromLast=TRUE)&!is.na(ann$gene_name)&ann$gene_name!='-'
mp<-setNames(ann$gene_id[good],ann$gene_name[good])
ea<-data.frame(symbol=rownames(o$GSE284505_main$expr));ea$Ensembl<-unname(mp[ea$symbol])
write.csv(ea,'outputs/EAU_unique_symbol_mapping.csv',row.names=FALSE)
g<-read.delim('outputs/GSE87433_platform.tsv',quote='',comment.char='',fill=TRUE,check.names=FALSE)
g<-g[match(rownames(o$GSE87433_main$expr),as.character(g$ID)),]
es<-regmatches(g$dbXrefs,gregexpr('ENSMUSG[0-9]+',g$dbXrefs))
v<-vapply(es,function(x){x<-unique(x);if(length(x)==1)x else NA_character_},character(1))
unique_target<-!is.na(v)&!duplicated(v)&!duplicated(v,fromLast=TRUE)
dm<-data.frame(Entrez=as.character(g$ID),Ensembl=v,unambiguous_one_to_one=unique_target)
write.csv(dm,'outputs/diabetes_Entrez_Ensembl_mapping.csv',row.names=FALSE)
sets<-list(space=rownames(o$OSD758_RSEM$expr),RP=rownames(o$GSE156533$expr),
 glaucoma=rownames(o$GSE241782_glaucoma$expr),EAU=na.omit(ea$Ensembl),diabetic_assayed=dm$Ensembl[dm$unambiguous_one_to_one])
shared<-Reduce(intersect,sets)
z<-list(interpretation='RNA-seq features passed within-study expression filtering; array targets are assayed and mapped, without expression-detection calls. Mapping uses deposited annotations, not a newly curated gene release.',
 per_study_features=sapply(sets,length),diabetes_one_to_one_assayed_overlap=length(intersect(sets$space,sets$diabetic_assayed)),all_four_diseases_and_reference_shared=length(shared))
write.csv(data.frame(Ensembl=shared),'outputs/core_panel_shared_gene_universe.csv',row.names=FALSE)
write_json(z,'outputs/shared_gene_universe_summary.json',pretty=TRUE,auto_unbox=TRUE)
print(z)
