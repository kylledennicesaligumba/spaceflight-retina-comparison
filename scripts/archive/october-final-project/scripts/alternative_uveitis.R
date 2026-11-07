code<-readLines('scripts/analyze.R');eval(parse(text=code[seq_len(grep('^# Spaceflight:',code)-1)]))
cp<-as.data.frame(read_excel('downloads/GSE144168_normalized_CPM.xlsx'))
cols<-grep('^(NR|DR)[0-9]+$',names(cp),value=TRUE)
x<-as.matrix(cp[,cols,drop=FALSE]);rownames(x)<-cp[[1]]
duplicate_symbols<-sum(duplicated(rownames(x)))
ambiguous<-duplicated(rownames(x))|duplicated(rownames(x),fromLast=TRUE)
writecsv(data.frame(label=rownames(x)[ambiguous]),'GSE144168_ambiguous_identifier_rows.csv')
# The repository says log2 CPM; keep the deposited scale and exclude only all-zero rows for exploratory QC.
allrows<-nrow(x);x<-x[!ambiguous & rowSums(abs(x))>0,,drop=FALSE]
m<-data.frame(sample=cols,group=ifelse(grepl('^NR',cols),'Naive','EAU'),label=cols)
qc(x,m,'GSE144168_total_live_retinal_cells',FALSE,units='deposited log2 CPM')
results$GSE144168_total_live_retinal_cells$deposited_rows<-allrows
results$GSE144168_total_live_retinal_cells$duplicate_label_entries_in_deposit<-duplicate_symbols
results$GSE144168_total_live_retinal_cells$ambiguous_rows_excluded_from_QC<-sum(ambiguous)
results$GSE144168_total_live_retinal_cells$DE_status<-'Not run: technical sort replicates and independent mouse pool mapping unresolved; normalized log2 CPM cannot substitute for raw counts.'
main<-readRDS('outputs/analysis_objects.rds');ann<-read.csv('outputs/P23H_gene_annotation.csv')
good<-!duplicated(ann$gene_name)&!duplicated(ann$gene_name,fromLast=TRUE)
map<-setNames(ann$gene_id[good],ann$gene_name[good])
results$GSE144168_total_live_retinal_cells$symbol_mapped_assayed_overlap<-length(intersect(rownames(main$OSD758_RSEM$expr),unname(map[rownames(x)])))
write_json(results,'outputs/alternative_uveitis_results.json',pretty=TRUE,auto_unbox=TRUE)
