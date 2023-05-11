import pandas as pd

de_files=["oxid.ConditionNPCASDDM_vs_ConditionNPCASDN.de.sig.tsv",
          "oxid.ConditionNPCASDDM_vs_ConditionNPCTDN.de.sig.tsv",
          "oxid.ConditionNPCASDN_vs_ConditionNPCTDN.de.sig.tsv"]
comparisons=['ASDDM_vs_ASDN','ASDDM_vs_TDN','ASDN_vs_TDN']
de_genes={}
for i in range(len(de_files)):
    print(i) 
    data=open(de_files[i],'r').read().strip().split('\n')
    cur_comparison=comparisons[i] 
    for line in data[1::]:
        tokens=line.split('\t')
        gene=tokens[0]
        lfc=float(tokens[1])
        if gene not in de_genes:
            de_genes[gene]={}
        de_genes[gene][cur_comparison]=lfc
de_df=pd.DataFrame.from_dict(de_genes).transpose()
de_df.to_csv("oxid.de.summary.txt",header=True,index=True,sep='\t')


