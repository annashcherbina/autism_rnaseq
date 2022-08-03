import pandas as pd
import pdb

comparisons=["ConditionASDDM_vs_ConditionTDN",
             "ConditionASDN_vs_ConditionTDN",
             "ConditionASDDM_vs_ConditionASDN"]
for comparison in comparisons:
    print(comparison)
    clin_sig_first=pd.read_csv(f"spearman.{comparison}.first.tsv",header=0,sep='\t').to_dict(orient='index')
    clin_sig_last=pd.read_csv(f"spearman.{comparison}.first.tsv",header=0,sep='\t').to_dict(orient='index')
    #merge taking the max
    merged_clin=clin_sig_first
    for gene in clin_sig_last:
        if gene not in merged_clin:
            merged_clin[gene]={}
        for feat in clin_sig_last[gene]:
            if feat not in merged_clin[gene]:
                merged_clin[gene][feat]=clin_sig_last[gene][feat]
            else:
                merged_clin[gene][feat]=max([clin_sig_last[gene][feat],merged_clin[gene][feat]])
                
    merged_clin_sig_df=pd.DataFrame.from_dict(merged_clin,orient='index')
    
    #get the de metrics
    de=pd.read_csv(f"{comparison}.de.sig.varfilter.tsv",header=0,sep='\t',index_col=['Gene'])
    merged=merged_clin_sig_df.join(de,how='left')

    #add WGCNA module membership
    modules=pd.read_csv("modules.txt",header=0,sep='\t',index_col=['ProbeID'])
    merged=merged.join(modules,how='left')
    
    
    #add the GSEA metrics
    gsea=pd.read_csv(f"{comparison}.de.tsv.GSEA.tsv",header=0,sep='\t')
    gsea=gsea[['Term','ledge_genes']]
    gsea_dict={}
    for index,row in gsea.iterrows():
        term=row['Term']
        genes=row['ledge_genes'].split(';')
        for gene in genes:
            if gene not in gsea_dict:
                gsea_dict[gene]=[term]
            else:
                gsea_dict[gene].append(term)
    for gene in gsea_dict:
        gsea_dict[gene]=';'.join(gsea_dict[gene])
    gsea_df=pd.DataFrame.from_dict(gsea_dict,orient='index',columns=['GSEA'])
    
    merged=merged.join(gsea_df,how='left')

    #add gProfiler
    gprofiler=pd.read_csv(f"{comparison}.gProfiler.tsv",header=0,sep='\t')
    gprofiler=gprofiler[['term.name','intersection']]
    gprofiler_dict={}
    for index,row in gprofiler.iterrows():
        term=row['term.name']
        genes=row['intersection'].split(',')
        for gene in genes:
            if gene not in gprofiler_dict:
                gprofiler_dict[gene]=[term]
            else:
                gprofiler_dict[gene].append(term)
    for gene in gprofiler_dict:
        gprofiler_dict[gene]=';'.join(gprofiler_dict[gene])
    gprofiler_df=pd.DataFrame.from_dict(gprofiler_dict,orient='index',columns=['gprofiler'])
    merged=merged.join(gprofiler_df,how='left')

    merged.to_csv(f"merged.{comparison}.tsv",header=True,index=True,sep='\t')
    
        
