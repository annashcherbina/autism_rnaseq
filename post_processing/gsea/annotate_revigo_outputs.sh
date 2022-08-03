for comparison in polarizationM0_vs_polarizationM1 polarizationM0_vs_polarizationM2 polarizationM1_vs_polarizationM2
do
    revigo_file=$comparison.clustered.go.tsv
    filtered_gsea_file=gsea/filtered/$comparison.gsea.tsv
    python annotate_revigo_outputs.py --gsea_sig_term_inputs $filtered_gsea_file --revigo_annotations $revigo_file
done



