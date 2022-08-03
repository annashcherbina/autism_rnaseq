from scipy.stats import hypergeom
import pandas as pd

modules=open('modules.txt','r').read().strip().split('\n')
module_dict={}
for line in modules[1::]:
    tokens =line.split('\t')
    gene=tokens[0]
    module=tokens[1]
    if module not in module_dict:
        module_dict[module]=[gene]
    else:
        module_dict[module].append(gene)
comparisons=['ASDDM_ASDN','ASDDM_TDN','ASDN_TDN']
expressed_genes={}
diff_genes={}
wgcna_ledges={}
module_to_pval={}
module_to_wgcna={}
outf_genes=open('module.gene.enrichments.tsv','w')
outf_gsea=open('module.gsea.enrichments.tsv','w')
outf_genes.write('Comparison\tModule\tP.val\tCommonGenes\n')
outf_gsea.write('Comparison\tModule\tGOTerm\tP.val\tCommonGenes\n')

for comparison in comparisons:
    
    expressed_genes[comparison]=open(f'{comparison}.expressed.tsv','r').read().strip().split('\n') 
    diff_comparison=pd.read_csv(f'{comparison}.diff.tsv',header=0,sep='\t')
    diff_genes[comparison]= set(diff_comparison[(diff_comparison['adj.P.Val']<0.01) & (abs(diff_comparison['logFC'])>2)]['Gene'].tolist())
    wgcna=pd.read_csv(f'{comparison}.diff.tsv.GSEA.tsv.sig.GOBP.tsv',header=0,sep='\t')
    
    module_to_pval[comparison]={}
    module_to_wgcna[comparison]={} 


    
    for module in module_dict:
        module_genes=set(module_dict[module])
        common_genes=module_genes.intersection(diff_genes[comparison])
        n_module_genes=len(module_genes)
        n_diff_genes=len(diff_genes[comparison])
        #get module enrichments
        cur_pval=hypergeom.sf(len(common_genes)-1,
                              len(expressed_genes[comparison]),
                              n_diff_genes,
                              n_module_genes)
        outf_genes.write(comparison+'\t'+module+'\t'+str(cur_pval)+'\t'+str(','.join(common_genes))+'\n')
        
        if cur_pval < 0.1:
            module_to_wgcna[comparison][module]=[]
            #check for hypergeometric enrichment of WGCNA
            for index,row in wgcna.iterrows():
                term=row['Term']
                ledge_genes=set(row['ledge_genes'].split(';'))
                common_module_gsea=ledge_genes.intersection(module_genes)
                if len(common_module_gsea)==0:
                    continue
                cur_pval_gsea=hypergeom.sf(len(common_module_gsea)-1,
                                      len(expressed_genes[comparison]),
                                      len(ledge_genes),
                                      n_module_genes)
                outf_gsea.write(comparison+'\t'+module+'\t'+term+'\t'+str(cur_pval_gsea)+'\t'+str(','.join(set(diff_genes[comparison]).intersection(common_module_gsea)))+'\n')
outf_genes.close()
outf_gsea.close()
    
