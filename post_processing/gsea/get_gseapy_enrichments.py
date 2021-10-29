import pandas as pd
import gseapy as gp
import sys 
import scipy.stats as st
import numpy as np
from statsmodels.stats.multitest import multipletests
min_set_size=10
num_permutations=10
#note: we run single-sample GSEA on wilcoxon auc values
data_file=sys.argv[1]
data=pd.read_csv(data_file,header=0,sep='\t')
data['z']=np.sign(data['logFC'])*st.norm.ppf(1-data['adj.P.Val']*.5)
gsea_inputs=data[['Gene','z']].sort_values('z',ascending=False)
result_list=[]
gsea_weight=1
result = gp.prerank(gsea_inputs,
                    gene_sets="/mnt/data/annotations/msigdb_v7.4/msigdb_v7.4_GMTs/msigdb.v7.4.symbols.gmt",
                    outdir=None,
                    min_size=min_set_size,
                    max_size=1e6,
                    weighted_score_type=gsea_weight,
                    permutation_num=num_permutations,
                    ascending=False,
                    no_plot=True,
                    verbose=True,
                    processes=40,
                    seed=0).res2d
result.index += f"_w{gsea_weight}"
result["gsea_weight"] = gsea_weight
result = result.rename(columns={"pval": "pv", "fdr": "fdr_gseapy"})
result["pv"] = (result["pv"] * num_permutations + 1) / (num_permutations + 1)
result["fdr_bh"] = multipletests(result["pv"], method="fdr_bh")[1]
result["fwer_bf"] = np.clip(result["pv"] * result.shape[0], 0, 1)
result.to_csv(data_file+'.GSEA.tsv',sep='\t',index=True,header=True)
