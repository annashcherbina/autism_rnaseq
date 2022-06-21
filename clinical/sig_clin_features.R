rm(list=ls())

#conditions
conditions=read.table("subject_to_condition.txt",header=T,sep='\t')
conditions$MajorGroup='ASD'
conditions$MajorGroup[conditions$Condition=='TD-N']='TD'

#first
first=as.data.frame(scale(read.table('first_vals.txt',header=T,sep='\t',row.names=1)))
features=colnames(first)

first$Subject=rownames(first)
first=merge(first,conditions,by='Subject')

#last 
last=as.data.frame(scale(read.table('last_vals.txt',header=T,sep='\t',row.names=1)))
last$Subject=rownames(last)
last=merge(last,conditions,by="Subject")

#merge them 
#merged=rbind(first,last)
merged=last
anova_effects=NULL
wilcox_effects=NULL
for(feature in features)
{
  print(feature)
  #anova 
  df_subset=merged[,c(feature,'Condition','MajorGroup')]
  colnames(df_subset)=c('feature','Condition','MajorGroup')
  anova.out=as.data.frame(TukeyHSD(aov(feature ~ Condition,data=df_subset))$Condition)
  anova.out$feature=feature 
  if(is.null(anova_effects)){
    anova_effects=anova.out 
  }else{
    anova_effects=rbind(anova_effects,anova.out)
  }
  #wilcox-test between ASD & TD
  asd_group=merged[merged$MajorGroup=='ASD',c(feature)]
  asd_group<-asd_group[!is.na(asd_group)]
  td_group=merged[merged$MajorGroup=='TD',c(feature)]
  td_group<-td_group[!is.na(td_group)]
  
  if((length(asd_group)>0) & (length(td_group>0)))
  {
    wilcox_out=wilcox.test(asd_group,td_group,paired=FALSE,exact=FALSE)
    wilcox_out_df=t(data.frame(c(wilcox_out$statistic[['W']], wilcox_out$p.value, feature)))
    rownames(wilcox_out_df)=NULL
    colnames(wilcox_out_df)=c("W",'p', 'feature')
    if(is.null(wilcox_effects)){
      wilcox_effects=wilcox_out_df
    }else{
      wilcox_effects=rbind(wilcox_effects,wilcox_out_df)
    }
    
  }
}
# write to output 
write.table(wilcox_effects,'wilcox_effects.last.tsv',row.names = T,col.names = T,sep='\t')
write.table(anova_effects,'anova_effects.last.tsv',row.names = T,col.names = T,sep='\t')
