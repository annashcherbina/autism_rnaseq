#venn diagram inputs

asdn_tdn_genes=set(open('pval.lt.0.001.lfc.gt.2.NPC_ASDN_TDN.diff.tsv','r').read().strip().split('\n'))
asddm_tdn_genes=set(open('pval.lt.0.001.lfc.gt.2.NPC_ASDDM_TDN.diff.tsv','r').read().strip().split('\n'))
asddm_asdn_genes=set(open('pval.lt.0.001.lfc.gt.2.NPC_ASDDM_ASDN.diff.tsv','r').read().strip().split('\n'))
asddm_asdn_unique=list(asddm_asdn_genes - asddm_tdn_genes - asdn_tdn_genes)
print(len(asddm_asdn_unique))
#we want to restrict to the pink, tan, magenta modules
outf=open('NPC_ASDDM_ASDN_unique.txt','w')
outf.write('\n'.join(asddm_asdn_unique))
