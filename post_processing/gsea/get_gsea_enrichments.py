import pandas as pd
import gseapy as gp
import sys
import scipy.stats as st
import numpy as np
data_file=sys.argv[1]
data=pd.read_csv(data_file,header=0,sep='\t')
data['z']=np.sign(data['logFC'])*st.norm.ppf(1-data['adj.P.Val']*.5)
gsea_inputs=data[['Gene','z']].sort_values('z',ascending=False)
result = gp.prerank(gsea_inputs,
                    gene_sets="/data/GxE_gsea_and_go_gene_sets/msigdb_v7.4/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt",
                    outdir=None,
                    min_size=5,
                    max_size=300,
                    permutation_num=100,
                    weighted_score_type=0,
                    ascending=False,
                    no_plot=True,
                    verbose=True,
                    processes=10,
                    seed=0).res2d
result.to_csv(data_file+'.GSEA.tsv',sep='\t',index=True,header=True)

