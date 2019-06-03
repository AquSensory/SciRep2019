library("igraph")
p = 0.05
G=suppressWarnings(read.graph("~/Desktop/synaptome.graphml",format="graphml"))
G=as.undirected(G,mode="collapse")
G_genes=V(G)

G_metric=function(G) gsize(G) #DEFINE CONNECTIVITY METRIC FOR COMPARISON eg. NO. OF EDGES
GOI=c("Aqu2.1.08624_001","Aqu2.1.39742_001","Aqu2.1.18863_001","Aqu2.1.26149_001","Aqu2.1.41520_001","Aqu2.1.13963_001","Aqu2.1.31831_001","Aqu2.1.41903_001","Aqu2.1.43783_001","Aqu2.1.24262_001","Aqu2.1.27832_001","Aqu2.1.39924_001")
G_GOI_subgraph=induced.subgraph(G,GOI)

# MC SAMPLING OF RANDOM SUBGRAPHS
MC_rounds = 10000
MC_samples = c()
message( sprintf("Perform %i rounds of mc sampling...", MC_rounds ) )
MC_progressbar=txtProgressBar(min=0, max=MC_rounds, style=3)

for (r in 1:MC_rounds)
{ G_subgraph= induced.subgraph(G,sample(G_genes,length(GOI),replace=F))
  MC_samples[r]= G_metric(G_subgraph)
  setTxtProgressBar(MC_progressbar,r)
}
close(MC_progressbar)
message("")

MC_quantile= quantile(MC_samples,1-p)
m_G_GOI_subgraph= G_metric(G_GOI_subgraph)
if(m_G_GOI_subgraph>MC_quantile)
  {message(sprintf("Induced subgraph is more connected (%i>%i)",m_G_GOI_subgraph,MC_quantile))
} else {
  message(sprintf("Induced subgraph is more connected (%i<%i)",m_G_GOI_subgraph,MC_quantile))
}
ecdf(MC_samples)(m_G_GOI_subgraph)
1-(ecdf(MC_samples)(m_G_GOI_subgraph))

pValue= 1
d_pValue= 0.0001
while(pValue>0)
{ pValue_=max(pValue-d_pValue,0)
  if(m_G_GOI_subgraph>quantile(MC_samples,1-pValue_))
  {pValue= pValue_;
  } else{
    break;
  }
}
message(sprintf("P-value:%f",pValue))
