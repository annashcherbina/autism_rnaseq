import pandas as pd
import sys
data=open(sys.argv[1],'r').read().strip().split('\n')
header=data[0].split('\t')
header=['ENSGID','Gene']+header
outf=open(sys.argv[1]+'.genes.txt','w')
outf.write('\t'.join(header)+'\n')

for line in data[1::]:
    tokens=line.split('\t')
    gene=tokens[0].split('-')
    ensgid=gene[0]
    genename=''.join(gene[1::])
    outf.write(ensgid+'\t'+genename+'\t'+'\t'.join(tokens[1::])+'\n')
outf.close()

