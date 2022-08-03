"""
This script parses GSEA results from multiple analyses, filters by nes and fdr_gseapy_thresholds.
It can optionally aggregate across multiple GSEA runs to generate a GO term x sample matrix \
of GSEA output fields (i.e. like a heatmap)
"""

import argparse
import os
from os import listdir
from os.path import join

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="aggregate results from GSEA analysis")
    parser.add_argument("--term_prefix", default="")
    parser.add_argument("--fdr_gseapy_thresh", type=float, default=0.01)
    parser.add_argument("--nes_thresh", type=float, default=0.2)
    parser.add_argument("--filter_RP_driven_terms",action="store_true")
    parser.add_argument("--RP_driven_thresh",type=float,default=0.25,help="fraction of ledge genes in a GO term that are RPL or RPS")
    parser.add_argument("--min_set_size", type=int, default=50)
    parser.add_argument("--max_set_size", type=int, default=500)
    parser.add_argument(
        "--source_dir", help="directory with source GSEA results to aggregate", default=None
    )
    parser.add_argument(
        "--source_files",
        nargs="+",
        help="list of GSEA files to aggregate, alternative input to --source_dir",
        default=None,
    )
    parser.add_argument("--out_dir", help="directory to write filtered GSEA results to")
    parser.add_argument(
        "--keep_fields",
        nargs="+",
        default="all",
        help="fields to parse out and keep from the gseapy resultt",
    )
    parser.add_argument("--aggregate", action="store_true", default=False)
    parser.add_argument(
        "--aggregation_fields",
        nargs="+",
        default=["fdr_thresholded_nes"],
        help="aggregated terms x samples into a matrix with values specified by aggregation_fields. \
        fdr_thresholded_nes will store nes output if fdr < thresh; else None",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.aggregate is True:
        terms = set([])
        aggregated_dfs = {}

    # identify files to aggregate
    if os.path.exists(args.out_dir) is False:
        os.makedirs(args.out_dir)
    if args.source_dir:
        source_files = [join(args.source_dir, f) for f in listdir(args.source_dir)]
    else:
        assert args.source_files is not None
        source_files = args.source_files

    # filter all GSEA out files by nes & term prefix (i.e. GOBP) & nes score.
    for f in source_files:
        data = pd.read_csv(f, header=0,sep='\t')
        data_filtered = data[
            (data.Term.str.startswith(args.term_prefix))
            & (data.fdr <= args.fdr_gseapy_thresh)
            & (abs(data.nes) >= args.nes_thresh)
            & (data.geneset_size >= args.min_set_size)
            & (data.geneset_size <= args.max_set_size)
        ]
        if args.keep_fields != "all":
            data_filtered = data_filtered[args.keep_fields]

        if args.filter_RP_driven_terms == True:
            # remove GO terms driven by RPL/RPS genes
            terms_to_keep=[]
            for index,row in data_filtered.iterrows():
                ledge=row['ledge_genes'].split(';')
                rp_tally=[i for i in ledge if (i.startswith('RPL') or i.startswith('RPS')) ]
                rp_fract=len(rp_tally)/len(ledge)
                if rp_fract < args.RP_driven_thresh:
                    terms_to_keep.append(row['Term'])
            data_filtered=data_filtered[data_filtered.Term.isin(terms_to_keep)]
            
        if args.aggregate is True:
            terms = terms.union(set(data_filtered.Term))
            aggregated_dfs[f] = data
        else:
            data_filtered.to_csv(join(args.out_dir, f.split('/')[-1]), sep="\t", header=True, index=False)
        print("filtered:" + str(f))

    # if aggregate flag set, generate Term x Sample matrix of aggregation fields
    # (i.e. nes, fdr_gseapy, etc)
    if args.aggregate is True:
        terms = list(terms)
        for aggregation_field in args.aggregation_fields:
            # aggregate GO terms x samples, with value = aggregation field
            outf = open(
                join(
                    args.out_dir,
                    f"merged.filtered.{args.fdr_gseapy_thresh}"
                    f".{args.nes_thresh}.{aggregation_field}.tsv",
                ),
                "w",
            )
            outf.write("Term\t" + "\t".join(source_files) + "\n")
            for term in terms:
                outf.write(term)
                for source_file in source_files:
                    df = aggregated_dfs[source_file]
                    row = df[df.Term == term]
                    if aggregation_field == "fdr_thresholded_nes":
                        if float(row["fdr"]) < args.fdr_gseapy_thresh:
                            outf.write("\t" + str(round(float(row["nes"]), 2)))
                        else:
                            outf.write("\tNA")
                    else:
                        outf.write("\t" + str(str(row[aggregation_field])))
                outf.write("\n")
            outf.close()


if __name__ == "__main__":
    main()
