import pandas as pd

def gene_info(
    gencode_row,
):
    """Extract gene_id, gene_name, gene_type for a single entry in gencode table;
    if any of these are not present in the entry, store a None value
    gencode data frame follows format:
    chr1    HAVANA    gene    11869    14409    .    +    .    gene_id "ENSG00000223972.5";
    gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id
    "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
    """
    gencode_row_entries = [
        entry.strip('"').split('"') for entry in gencode_row.strip(";").replace(" ", "").split(";")
    ]
    # parse the row to find instances of gene_id, gene_name, hgnc_id, gene_type.
    # verify that there is 0 or 1 entries of each type

    # each attribute defaults to None (i.e. absent from the row)
    g_id = None
    t_id = None
    g_name = None
    g_type = None
    hgnc_id = None

    gene_ids = []
    transcript_ids = []
    gene_names = []
    gene_types = []
    hgnc_ids = []
    for annotation in gencode_row_entries:
        if annotation[0] == "gene_id":
            gene_ids.append(annotation[1])
        if annotation[0] == "transcript_id":
            transcript_ids.append(annotation[1])
        if annotation[0] == "gene_name":
            gene_names.append(annotation[1])
        if annotation[0] == "gene_type":
            gene_types.append(annotation[1])
        if annotation[0] == "hgnc_id":
            hgnc_ids.append(int(annotation[1].split(":")[1]))
    try:
        assert (
            (len(gene_ids) < 2)
            & (len(transcript_ids) < 2)
            & (len(gene_names) < 2)
            & (len(gene_types) < 2)
            & (len(hgnc_ids) < 2)
        )
    except AssertionError as err:
        logging.exception(
            f"Row has > 1 entries for versioned_gene_ids:{versioned_gene_ids},"
            f" gene_names:{gene_names}, gene_types:{gene_types}, or hgnc_ids:{hgnc_ids}"
        )
        raise err
    if len(gene_ids) == 1:
        g_id = gene_ids[0]
    if len(transcript_ids)==1:
        t_id=transcript_ids[0]
    if len(gene_names) == 1:
        g_name = gene_names[0]
    if len(gene_types) == 1:
        g_type = gene_types[0]
    if len(hgnc_ids) == 1:
        hgnc_id = hgnc_ids[0]

    return g_id, t_id, g_name, g_type, hgnc_id

data=pd.read_csv("gencode.v41.annotation.gtf",header=0,sep='\t',skiprows=5)
outf=open("filtered.gencode.v41.annotation.gtf",'w')
outf.write('\t'.join(['gene_id','transcript_id','gene_name','gene_type','hgnc_id'])+'\n')
print("loaded gencode") 
for index,row in data.iterrows(): 
    g_id,t_id,g_name,g_type,hgnc_id=gene_info(row[8]) 
    outf.write('\t'.join([str(i) for i in [g_id,t_id,g_name,g_type,hgnc_id]])+'\n')
outf.close() 


