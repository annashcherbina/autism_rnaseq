import argparse
import pandas as pd
from os import listdir
from os.path import  join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--gobp_id_to_term",default="/data/GxE_gsea_and_go_gene_sets/msigdb_v7.4/map_go_to_description_msigdb/GO_ID_TO_DESCRIPTION/GOBP_ID_TO_DESCRIPTION_MSIGDB.txt")
    parser.add_argument("--gsea_sig_term_inputs", nargs="+")
    parser.add_argument("--revigo_annotations", nargs="+")
    return parser.parse_args()

def main():
    args=parse_args()
    # get map of go term name to id
    gobp_id_to_term=open(args.gobp_id_to_term,'r').read().strip().split('\n')
    gobp_id_to_term_dict={}
    for line in gobp_id_to_term:
        tokens=line.split('\t')
        gobp_id_to_term_dict[tokens[0]]=tokens[1]        
    print("made term-id dict")
    header = ['Term', 'es', 'nes', 'pv', 'fdr_gseapy', 'geneset_size', 'matched_size','genes','ledge_genes','go','cluster','parent','parentSimScore','score','size','term','parentTerm','ClusterRepresentative']
    for i in range(len(args.revigo_annotations)):
        rows=[]
        cur_revigo_annot=args.revigo_annotations[i]
        go_to_info=dict()
        for line in open(cur_revigo_annot,'r').read().strip().split('\n')[1::]: 
            tokens=line.split('\t')
            term=tokens[0].strip()
            #print(term)
            go_to_info[term]=[token.strip() for token in tokens]
        # outf=open(cur_revigo_annot+".joined.txt",'w')
        # outf.write(header+'\n')
        cur_gsea=args.gsea_sig_term_inputs[i]
        for line in open(cur_gsea,'r').read().strip().split('\n')[1::]: 
            tokens=line.split('\t')
            go_term=gobp_id_to_term_dict[tokens[0]]
            if go_term in go_to_info:
                revigo_vals=go_to_info[go_term]
                cluster_rep=revigo_vals[-1]==revigo_vals[-2]
                rows.append(tokens+revigo_vals+[cluster_rep])
        out_df=pd.DataFrame(rows,columns=header)
        # sort by ClusterRepresentative and then by cluster
        out_df.sort_values(by=['ClusterRepresentative','cluster'],
                           ascending=[False, True],
                           inplace=True)
        # save to outfile
        out_df.to_csv(cur_revigo_annot+".joined.txt",sep='\t',index=False,header=True)
    
    
if __name__=='__main__':
    main()

