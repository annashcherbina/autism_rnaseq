# aggregate RSEM gene & transcript files across all age replicates
# convert ENS id's to known gene names
import argparse
from functools import reduce

import pandas as pd

pd.options.mode.chained_assignment = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="aggregate RSEM gene and transcript files across all age replicates "
    )
    parser.add_argument(
        "--rsem_file_list",
        help="2 columns, first is sample name, second is sample path; file should be tab delimited \
        without a header",
        default="isoform.samples.txt"
    )
    parser.add_argument(
        "--ensembl_id_to_gene_name",
        default=(
            "hg38.gencode.v1.gid_to_name.txt"
        ),
        help=" 2 column file with ENSGID in column 1 and Gene Symbol\
        in column 2; no header; tab-delimited",
    )
    parser.add_argument(
        "--out_dir", help="directory (either local or on s3 to save the aggregated rsem outputs",default="/srv/scratch/annashch/autism_rnaseq/post_processing/merged_rsem"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    # generate dictionary mapping ensembl id's to gene names
    try:
        map_file = pd.read_csv(args.ensembl_id_to_gene_name, header=None, sep="\t", index_col=0)
    except Exception as e:
        print(
            'make sure you have pandas installed with boto3 support:\
            pip install boto3 pandas "s3fs<=0.4"'
        )
        print(f"Unexpected error:{e}")
        raise
    map_file.index=[i.split('.')[0] for i in map_file.index]
    map_dict = map_file.to_dict()[1]
    print("generated ENSGID -> Gene Symbol dictionary")
    print(map_dict) 
    tpm_vals = []
    fpkm_vals = []
    expected_count_vals = []
    iso_pct_vals=[]
    samples = []

    rsem_file_list = pd.read_csv(
        args.rsem_file_list, header=None, sep="\t", index_col=0
    ).to_dict()[1]
    for sample_id in rsem_file_list:
        print(f"processing:{sample_id}")
        samples.append(sample_id)
        rsem_data = pd.read_csv(rsem_file_list[sample_id], header=0, sep="\t")
        rsem_tpm_sub = rsem_data[["transcript_id","gene_id","TPM"]]
        rsem_tpm_sub.rename(columns={"TPM": sample_id}, inplace=True)
        tpm_vals.append(rsem_tpm_sub)

        rsem_fpkm_sub = rsem_data[["transcript_id", "gene_id", "FPKM"]]
        rsem_fpkm_sub.rename(columns={"FPKM": sample_id}, inplace=True)
        fpkm_vals.append(rsem_fpkm_sub)

        rsem_expected_count_sub = rsem_data[["transcript_id", "gene_id", "expected_count"]]
        rsem_expected_count_sub.rename(columns={"expected_count": sample_id}, inplace=True)
        expected_count_vals.append(rsem_expected_count_sub)

        rsem_iso_pct_sub = rsem_data[["transcript_id", "gene_id", "IsoPct"]]
        rsem_iso_pct_sub.rename(columns={"IsoPct": sample_id}, inplace=True)
        iso_pct_vals.append(rsem_iso_pct_sub)

        

    # merge the data frame lists
    tpm_df = reduce(lambda x, y: pd.merge(x, y, on=["transcript_id", "gene_id"]), tpm_vals)
    fpkm_df = reduce(lambda x, y: pd.merge(x, y, on=["transcript_id","gene_id"]), fpkm_vals)
    expected_count_df = reduce(
        lambda x, y: pd.merge(x, y, on=["transcript_id", "gene_id"]), expected_count_vals
    )
    
    iso_pct_df = reduce(
        lambda x, y: pd.merge(x, y, on=["transcript_id", "gene_id"]), iso_pct_vals
    )
    
    print("merged sample data frames")

    tpm_df["Symbol"] = [
        map_dict[gene_id.split('.')[0]] if gene_id.split('.')[0] in map_dict else gene_id
        for gene_id in tpm_df["gene_id"]
    ]
    fpkm_df["Symbol"] = [
        map_dict[gene_id.split('.')[0]] if gene_id.split('.')[0] in map_dict else gene_id
        for gene_id in fpkm_df["gene_id"]
    ]
    expected_count_df["Symbol"] = [
        map_dict[gene_id.split('.')[0]] if gene_id.split('.')[0] in map_dict else gene_id
        for gene_id in expected_count_df["gene_id"]
    ]

    iso_pct_df["Symbol"] = [
        map_dict[gene_id.split('.')[0]] if gene_id.split('.')[0] in map_dict else gene_id
        for gene_id in iso_pct_df["gene_id"]
    ]

    tpm_df.to_csv("/".join([args.out_dir, "tpm.isoform.txt"]), sep="\t", header=True, index=False)
    fpkm_df.to_csv("/".join([args.out_dir, "fpkm.isoform.txt"]), sep="\t", header=True, index=False)
    expected_count_df.to_csv(
        "/".join([args.out_dir, "expected_count.isoform.txt"]), sep="\t", header=True, index=False
    )
    iso_pct_df.to_csv(
        "/".join([args.out_dir, "iso_pct.isoform.txt"]), sep="\t", header=True, index=False
    )

if __name__ == "__main__":
    main()

