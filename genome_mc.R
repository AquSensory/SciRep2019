setwd("~/Desktop")

genome = read.table("AquPCgenes_CAP_ave.txt",header=T,row.names=1,sep="\t")
Ngenome = length(rownames(genome))
GOI=list("Aqu2.1.37536_001","Aqu2.1.39742_001","Aqu2.1.18863_001","Aqu2.1.22337_001",
         "Aqu2.1.13963_001","Aqu2.1.31831_001","Aqu2.1.41903_001","Aqu2.1.27832_001")
iGOI = which(row.names(genome) %in% GOI)
GOC=sample(rownames(genome),length(GOI),replace=F)
iGOC = which(row.names(genome) %in% GOC)
nGOIclades = c()
nGOCclades = c()
#Part 1 generate runs with GOI
for (r in c(1:200))
{ iG1 = sample(setdiff(1:Ngenome,iGOI),119-length(iGOI),replace=F)
  genesubset1 = genome[union(iGOI,iG1),]
  hm1<-pheatmap(as.matrix(genesubset1), 
                 cluster_cols=F, cluster_rows=T, scale = "row", 
                 clustering_method = "complete",cutree_rows = 10,silent=T)
  out1<-sort(cutree(hm1$tree_row, k=10))
  GOI_clades = out1[which(names(out1) %in% GOI)]
  GOI_clades = GOI_clades[order(names(GOI_clades))]
  nGOIclades[r] = length(unique(GOI_clades))
  print(r)  
}
#Part 2 generate runs with GOC
for (r in c(1:200))
{ iG2 = sample(setdiff(1:Ngenome,iGOC),119-length(iGOC),replace=F)
  genesubset2 = genome[union(iGOC,iG2),]
  hm2<-pheatmap(as.matrix(genesubset2), 
                cluster_cols=F, cluster_rows=T, scale = "row", 
                clustering_method = "complete",cutree_rows = 10,silent=T)
  out2<-sort(cutree(hm2$tree_row, k=10))
  GOC_clades = out2[which(names(out2) %in% GOC)]
  GOC_clades = GOC_clades[order(names(GOC_clades))]
  nGOCclades[r] = length(unique(GOC_clades))
  print(r)
}
as.data.frame(table(nGOIclades))
as.data.frame(table(nGOCclades))
par(mfrow=c(1,2)) 
hist(nGOIclades)
hist(nGOCclades)
t.test(nGOIclades,nGOCclades) 
#for unequal variances, so df not necess a integer - conventional to round down 
#to the nearest integer before consulting standard t tables