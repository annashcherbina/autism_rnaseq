import argparse
from os import listdir
from os.path import  join

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--gobp_id_to_term",default="/data/GxE_gsea_and_go_gene_sets/msigdb_v7.4/map_go_to_description_msigdb/GO_ID_TO_DESCRIPTION/GOBP_ID_TO_DESCRIPTION_MSIGDB.txt")
    parser.add_argument("--gsea_sig_term_dir",default="filtered")
    parser.add_argument("--out_dir",default="revigo")
    return parser.parse_args()

def main():
    args=parse_args()
    #get map of go term name to id
    gobp_id_to_term=open(args.gobp_id_to_term,'r').read().strip().split('\n')
    gobp_id_to_term_dict={}
    for line in gobp_id_to_term:
        tokens=line.split('\t')
        gobp_id_to_term_dict[tokens[0]]=tokens[1]

    gsea_files = [f for f in listdir(args.gsea_sig_term_dir) if f.endswith('GSEA.tsv')]
    for f in gsea_files:
        outf=open(join(args.out_dir,f),'w')
        outf.write("GOterm\tenrichment_P-value\n")
        data=open(join(args.gsea_sig_term_dir,f),'r').read().strip().split('\n')
        for line in data[1::]:
            tokens=line.split('\t')
            term_name=tokens[0]
            term_id=gobp_id_to_term_dict[term_name]
            fdr=float(tokens[4])+1e-300
            outf.write(term_id+'\t'+str(fdr)+'\n')
        outf.close() 
            
if __name__=='__main__':
    main()

