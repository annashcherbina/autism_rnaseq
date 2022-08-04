#aggregate RSEM gene & transcript files across all age replicates
#convert ENS id's to known gene names
import argparse
import pdb 
def parse_args():
    parser=argparse.ArgumentParser(description="aggregate RSEM gene and transcript files across all age replicates ")
    parser.add_argument("--rsem_file_list",default="rsem.gene.cd47.files.txt")
    parser.add_argument("--ensembl_id_to_gene_name",default="hg38.gencode.v1.gid_to_name.txt")
    parser.add_argument("--outf_prefix",default="rna.counts.txt")
    return parser.parse_args()

def main():
    args=parse_args()
    #generate dictionary mapping ensembl id's to gene names
    id_to_name=dict()
    map_file=open(args.ensembl_id_to_gene_name,'r').read().strip().split('\n')
    map_dict=dict() 
    for line in map_file:
        tokens=line.split('\t')
        ensembl_id=tokens[0]
        gene_name=tokens[1]
        map_dict[ensembl_id]=gene_name
    print("generated id to gene name dictionary")

    expected_count_dict=dict()
    tpm_dict=dict()
    fpkm_dict=dict()
    samples=[] 
    
    rsem_file_list=open(args.rsem_file_list,'r').read().strip().split('\n')
    for line in rsem_file_list:
        tokens=line.split('\t')
        sample_id=tokens[0]
        samples.append(sample_id) 
        expected_count_dict[sample_id]=dict()
        
        rsem_file=open(tokens[1],'r').read().strip().split('\n')
        header=rsem_file[0].split('\t')
        gene_id_index=header.index('gene_id')
        expected_count_index=header.index('expected_count')
        tpm_index=header.index('TPM')
        fpkm_index=header.index('FPKM')
        for entry in rsem_file[1::]:
            entry=entry.split('\t')
            cur_gene_id=entry[gene_id_index]
            cur_expected_count=entry[expected_count_index]
            cur_tpm=entry[tpm_index]
            cur_fpkm=entry[fpkm_index]
            expected_count_dict[sample_id][cur_gene_id]=[cur_expected_count,cur_tpm,cur_fpkm]
            #pdb.set_trace()
    outf_expected_count=open(args.outf_prefix+'.expected_count','w')
    outf_expected_count.write('GeneID\tGeneName\t'+'\t'.join(samples)+'\n')
    for gene_id in map_dict:
        outf_expected_count.write(gene_id+'\t'+map_dict[gene_id])
        for sample in samples:
            try:
                outf_expected_count.write('\t'+expected_count_dict[sample][gene_id][0])
            except:
                outf_expected_count.write('\t0') 
        outf_expected_count.write('\n')
        
    outf_tpm=open(args.outf_prefix+'.tpm','w')
    outf_tpm.write('GeneID\tGeneName\t'+'\t'.join(samples)+'\n')
    for gene_id in map_dict:
        outf_tpm.write(gene_id+'\t'+map_dict[gene_id])
        for sample in samples:
            try:
                outf_tpm.write('\t'+expected_count_dict[sample][gene_id][1])
            except:
                outf_tpm.write('\t0') 
        outf_tpm.write('\n')
        
    outf_fpkm=open(args.outf_prefix+'.fpkm','w')
    outf_fpkm.write('GeneID\tGeneName\t'+'\t'.join(samples)+'\n')
    for gene_id in map_dict:
        outf_fpkm.write(gene_id+'\t'+map_dict[gene_id])
        for sample in samples:
            try:
                outf_fpkm.write('\t'+expected_count_dict[sample][gene_id][2])
            except:
                outf_fpkm.write('\t0')
        outf_fpkm.write('\n')
        
        
    
    
    
    
if __name__=="__main__":
    main()
    
    
    
