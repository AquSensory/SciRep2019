###Monte Carlo sampling 10000 dendrograms to obtain nclades GOIs fall in, contrast against randomly selected (same no. of) GOCs###

genome = read.table("AquPCgenes_BLIND_ave.txt",header=T,row.names=1,sep="\t")
Ngenome = length(rownames(genome))

#Part 1 generate runs with GOI
GOI=list("Aqu2.1.18457_001","Aqu2.1.34242_001","Aqu2.1.36578_001","Aqu2.1.41533_001","Aqu2.1.41520_001","Aqu2.1.13963_001","Aqu2.1.31831_001")
iGOI = which(row.names(genome) %in% GOI)
nGOIclades = c()
for (r in c(1:10000))
{ iG1 = sample(setdiff(1:Ngenome,iGOI),123-length(iGOI),replace=F)
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
as.data.frame(table(nGOIclades))

#Part 2 generate runs with GOC
GOC=sample(rownames(genome),length(GOI),replace=F)
iGOC = which(row.names(genome) %in% GOC)
nGOCclades = c()
for (r in c(1:10000))
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
as.data.frame(table(nGOCclades))

par(mfrow=c(1,2)) 
hist(nGOIclades,xlim=c(0,10),xlab="# clades",main="genes of interest")
hist(nGOCclades,xlim=c(0,10),xlab="# clades",main="randomly selected genes")

#Part 3 Kruskal-Wallis & Dunn's test (after running part1 & part2 x 5)
clades=c(nGOIclades,nGOCclades1,nGOCclades2,nGOCclades3,nGOCclades4,nGOCclades5)
genetype=c(rep("GOI",10000),rep("Ct1",10000),rep("Ct2",10000),rep("Ct3",10000),
           rep("Ct4",10000),rep("Ct5",10000))
type=c(rep("GOI",10000),rep("control",50000))
boxplot(clades~genetype,ylim=c(0,10),ylab="# clades genes fall in",
        main="CDH, CLTB, CLTC, CTNNA, CTTN, DNM1, DNM2 | developmental")
#Kruskal-Wallis
genetype=as.factor(genetype)
type=as.factor(type)
kruskal.test(clades~genetype)
kruskal.test(clades~type)
#posthoc / pairwise differences
library(FSA)
dunnTest(clades~genetype,method="bh") 
