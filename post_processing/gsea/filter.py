import sys
import pandas as pd
import pdb 
data=pd.read_csv(sys.argv[1],header=0,sep='\t')
data_sig=data[(data.fdr_gseapy<0.001)]
data_sig=data_sig[abs(data_sig.nes)>1]
with pd.option_context('mode.use_inf_as_na', True):
        data_sig = data_sig.dropna(subset=['nes'], how='all')
go_terms=data_sig[data_sig.Term.str.startswith('GO')]
gobp_terms=data_sig[data_sig.Term.str.startswith('GOBP')]

reactome_terms=data_sig[data_sig.Term.str.startswith('REACTOME')]
go_terms.to_csv(sys.argv[1]+'.sig.GO.tsv',header=True,index=False,sep='\t')
gobp_terms.to_csv(sys.argv[1]+'.sig.GOBP.tsv',header=True,index=False,sep='\t')
reactome_terms.to_csv(sys.argv[1]+'.sig.REACTOME.tsv',header=True,index=False,sep='\t')


