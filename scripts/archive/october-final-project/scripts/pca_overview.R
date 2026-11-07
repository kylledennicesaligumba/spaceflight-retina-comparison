o<-readRDS('outputs/analysis_objects.rds')
ns<-c('OSD758_RSEM','GSE156533','GSE241782_glaucoma','GSE284505_main','GSE87433_main','GSE314043')
titles<-c('Spaceflight: original RSEM','Adult P23H: P90','Experimental glaucoma','Autoimmune uveitis: untreated','Diabetic retina: untreated','Light injury: vehicle groups')
png('outputs/PCA_overview.png',width=2100,height=1450,res=170)
par(mfrow=c(2,3),mar=c(4.2,4.2,3.0,1),oma=c(4,0,2,0))
for(i in seq_along(ns)){
 z<-o[[ns[i]]];v<-apply(z$expr,1,var);p<-prcomp(t(z$expr[order(v,decreasing=TRUE)[1:min(2000,sum(v>0))],]),center=TRUE,scale.=FALSE)
 pc<-p$sdev^2/sum(p$sdev^2)*100;g<-factor(z$meta$group);co<-setNames(hcl.colors(nlevels(g),'Dark 3'),levels(g))
 if(ns[i]=='OSD758_RSEM')co<-c(AG033='#999999',AG067='#999999',AG1='#999999',HGC='#24749C',MG='#CC5B28')
 xr<-range(p$x[,1]);yr<-range(p$x[,2]);xr<-xr+c(-1,1)*diff(xr)*.17;yr<-yr+c(-1,1)*diff(yr)*.17
 plot(p$x[,1:2],col=co[as.character(g)],pch=16,cex=1.1,xlim=xr,ylim=yr,
  xlab=sprintf('PC1 (%.1f%%)',pc[1]),ylab=sprintf('PC2 (%.1f%%)',pc[2]),main=titles[i],cex.main=.95)
 if(ns[i]=='OSD758_RSEM')legend('topright',c('MG (n=5)','HGC (n=12)','AG groups (n=18)'),col=c(co['MG'],co['HGC'],'#999999'),pch=16,bty='n',cex=.7)
 else legend('topright',paste0(levels(g),' (n=',as.integer(table(g)),')'),col=co,pch=16,bty='n',cex=.62)
}
mtext('Each PCA is computed within its own study from the 2,000 most variable features.',outer=TRUE,side=3,line=.5,cex=.9)
mtext('n denotes libraries/arrays; all available main-contrast samples are retained.',outer=TRUE,side=1,line=1.5,cex=.75)
mtext('Uveitis and light-injury independence remains unverified; separation is not a disease-similarity score.',outer=TRUE,side=1,line=2.5,cex=.75)
dev.off()
