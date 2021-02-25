import pandas as pd

# related to processing splitseq

def get_bc1_matches():
	# from spclass.py - barcodes and their well/primer type identity
	bc_file = '/Users/fairliereese/mortazavi_lab/bin/pacbio-splitpipe/barcodes/bc_8nt_v2.csv'
	bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
	bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
	bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]

	# pivot on well to get df that matches bcs with one another from the same well
	bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
	bc_df = bc_df.rename_axis(None, axis=1).reset_index()
	bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)

	return bc_df

def get_illumina_metadata(fname):
    
    # we want info about primer type as well
    bc_df = get_bc1_matches()

    # read in illumina bcs
    df = pd.read_csv(fname, sep='\t')
    cols = ['sample', 'umi_count', 'gene_count', 'bc']
    df = df[cols]
    df.rename({'umi_count':'ill_umi_count', 'gene_count': 'ill_gene_count'},
              axis=1, inplace=True)
    df['bc3'] = df.bc.str.slice(start=0, stop=8)
    df['bc2'] = df.bc.str.slice(start=8, stop=16)
    df['bc1'] = df.bc.str.slice(start=16, stop=24)

    # merge bc1 df with illumina data 
    df = df.merge(bc_df, how='left', left_on='bc1', right_on='bc1_dt')
    df.rename({'bc1_randhex': 'Random hexamer', 'bc1_dt':'Poly dT'}, axis=1, inplace=True)

    # melt df to duplicate entries b/w dt and randhex primers for each cell
    id_vars = ['sample', 'ill_umi_count', 'ill_gene_count',\
               'bc', 'bc3', 'bc2', 'bc1', 'well']
    value_vars = ['Poly dT', 'Random hexamer']

    temp = pd.melt(df, id_vars=id_vars, value_vars=value_vars)
    temp.rename({'variable': 'primer_type', 'value': 'raw_bc1', 'bc': \
                 'merged_bc'}, axis=1, inplace=True)
    temp['raw_bc'] = temp.bc3+temp.bc2+temp.raw_bc1
    
    return temp 

def find_randhex_polydt_datasets(df, bc_df):
    dt_bc1s = bc_df.bc1_dt.tolist()
    randhex_bc1s = bc_df.bc1_randhex.tolist()
    
    dt_datasets = df.loc[df.bc1.isin(dt_bc1s), 'dataset'].tolist()
    randhex_datasets = df.loc[df.bc1.isin(randhex_bc1s), 'dataset'].tolist()
    
    return dt_datasets, randhex_datasets

def add_read_annot_metadata(df, bulk=False):
    fname = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/c2c12_paper_2020/sc_pacbio/illumina_cell_metadata_full_bcs.tsv'
    i_df = get_illumina_metadata(fname)
#     i_df = i_df[['raw_bc', 'primer_type']]
    df = df.merge(i_df, how='left', left_on='dataset', right_on='raw_bc')
    
    if bulk:
        sample_df = df['dataset'].to_frame()
        sample_df.drop_duplicates(inplace=True)
        
        mini_idf = i_df[['raw_bc', 'sample']]
        
        sample_df['experiment'] = sample_df.apply(lambda x: 'bulk' if 'PB' in x.dataset else 'sc', axis=1)
        sample_df = sample_df.merge(mini_idf, how='left', left_on='dataset', right_on='raw_bc')
        sample_df.loc[sample_df.dataset.isin(['PB154', 'PB155']), 'sample'] = 'MB'
        sample_df.loc[sample_df.dataset.isin(['PB213', 'PB214']), 'sample'] = 'MT'
        sample_df['tech'] = 'temp'
        sample_df.loc[sample_df['sample'].str.contains('nuclei'), 'tech'] = 'Single-nucleus'
        sample_df.loc[sample_df['sample'].str.contains('cells'), 'tech'] = 'Single-cell'
        sample_df.loc[sample_df.experiment == 'bulk', 'tech'] = 'Bulk'
        sample_df.drop(['sample', 'raw_bc'], axis=1, inplace=True)
        
        df = df.merge(sample_df, how='left', on='dataset')
    
    return df

def add_bcs_df(df):
    df['bc3'] = df.dataset.str.slice(start=0, stop=8)
    df['bc2'] = df.dataset.str.slice(start=8, stop=16)
    df['bc1'] = df.dataset.str.slice(start=16, stop=24)
    return df

# related to talon processing

def read_whitelist(fname):
    df = pd.read_csv(fname, header=None, names=['gid', 'tid'])
    whitelist = df.tid.tolist()
    return whitelist
