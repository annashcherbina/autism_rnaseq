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
measurements=['first','last']
de_genes={}

outf_genes=open('module.enrichments.clinical.assoc.genes.tsv','w')
outf_genes.write('Comparison\tClinFeat\tModule\tP.val\tCommonGenes\n')

for comparison in comparisons:
    print(comparison) 
    de_genes[comparison]=open(f'all.sig.{comparison}.tsv','r').read().strip().split('\n')
    n_de_genes=len(de_genes[comparison])
    print(f"n_de_genes:{n_de_genes}")
    clin_cor_first=pd.read_csv(f'spearman.{comparison}.first.tsv',header=0,sep='\t')
    clin_cor_last=pd.read_csv(f'spearman.{comparison}.last.tsv',header=0,sep='\t')
    clin_cor=pd.concat([clin_cor_first,clin_cor_last])
    for feature in clin_cor.columns:
        print(feature) 
        cor_genes=set(clin_cor.index[clin_cor[feature]!=0].tolist())
        for module in module_dict:
            module_genes=set(module_dict[module]).intersection(set(de_genes[comparison]))
            common_genes=module_genes.intersection(cor_genes)
            n_common_genes=len(common_genes)
            if n_common_genes < 2:
                continue 
            n_module_genes=len(module_genes)
            n_cor_genes=len(cor_genes)
            #get module enrichments
            cur_pval=hypergeom.sf(n_common_genes-1,
                                  n_de_genes,
                                  n_cor_genes,
                                  n_module_genes)
            outf_genes.write(comparison+'\t'+feature+'\t'+module+'\t'+str(cur_pval)+'\t'+str(','.join(common_genes))+'\n')
outf_genes.close()

