ASDDM_TDN=read.table("merged.ConditionASDDM_vs_ConditionTDN.tsv",header=T,sep='\t',row.names=1)
ASDN_TDN=read.table("merged.ConditionASDN_vs_ConditionTDN.tsv",header=T,sep='\t',row.names=1)

set1=row.names(ASDDM_TDN)
set2=row.names(ASDN_TDN)
asddm_unique=setdiff(set1,set2)

ASDDM_TDN_filtered=ASDDM_TDN[asddm_unique,]
write.table(ASDDM_TDN_filtered,file="merged.ConditionASDDM_vs_ConditionTDN.UNIQUE.tsv",row.names = T,col.names=T,sep='\t')