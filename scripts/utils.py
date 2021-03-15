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

def get_illumina_metadata():
    
    fname = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/c2c12_paper_2020/sc_pacbio/illumina_cell_metadata_full_bcs.tsv'
    
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
    df.rename({'bc1_randhex': 'Random hexamer', 'bc1_dt':'Oligo dT'}, axis=1, inplace=True)

    # melt df to duplicate entries b/w dt and randhex primers for each cell
    id_vars = ['sample', 'ill_umi_count', 'ill_gene_count',\
               'bc', 'bc3', 'bc2', 'bc1', 'well']
    value_vars = ['Oligo dT', 'Random hexamer']

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

def get_sample_df(datasets):
    
    
    sample_df = pd.DataFrame(data=datasets, columns=['dataset'])
    sample_df.drop_duplicates(inplace=True)
        
    i_df = get_illumina_metadata()
    mini_idf = i_df[['raw_bc', 'sample']]    
    sample_df['experiment'] = sample_df.apply(lambda x: 'bulk' if 'PB' in x.dataset else 'sc', axis=1)
    
    
    sample_df = sample_df.merge(mini_idf, how='left', left_on='dataset', right_on='raw_bc')
    
    sample_df.loc[sample_df.dataset.isin(['PB154', 'PB155']), 'sample'] = 'MB'
    sample_df.loc[sample_df.dataset.isin(['PB213', 'PB214']), 'sample'] = 'MT'
    
    sample_df['tech'] = 'temp'
    sample_df.loc[sample_df['sample'].str.contains('nuclei'), 'tech'] = 'Single-nucleus'
    sample_df.loc[sample_df['sample'].str.contains('cells'), 'tech'] = 'Single-cell'
    sample_df.loc[sample_df.experiment == 'bulk', 'tech'] = 'Bulk'
    
    return sample_df
    
def add_read_annot_metadata(df, bulk=False):
    i_df = get_illumina_metadata()
    df = df.merge(i_df, how='left', left_on='dataset', right_on='raw_bc')
    
    if bulk:
        datasets = df.dataset.unique().tolist()
        sample_df = get_sample_df(datasets)
        sample_df.drop(['sample', 'raw_bc'], axis=1, inplace=True)
        
        df = df.merge(sample_df, how='left', on='dataset')
        df.loc[df.dataset.isin(['PB154', 'PB155']), 'sample'] = 'MB'
        df.loc[df.dataset.isin(['PB213', 'PB214']), 'sample'] = 'MT'
    
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

# from a talon abundance file, get a list of columns that correspond to the datasets
def get_dataset_names(df):
    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                       'annot_transcript_id', 'annot_gene_name',
                       'annot_transcript_name', 'n_exons', 'length',
                       'gene_novelty', 'transcript_novelty', 'ISM_subtype', 'experiment']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]
    return dataset_cols

def get_gtf_info(fname, kind='gene'):
    
    df = pd.read_csv(fname, sep='\t', comment='#', usecols=[0,2,3,4,8], header=None)
    df.columns = ['chr', 'entry_type', 'start', 'stop', 'fields']
    if kind == 'gene':
        df = df.loc[df.entry_type == 'gene']
        
        name_pat = 'gene_name "'
        id_pat = 'gene_id "'
        
        df['gene_type'] = df.fields.str.split(pat='gene_type "', n=1, expand=True)[1]
        df['gene_type'] = df.gene_type.str.split(pat='"', n=1, expand=True)[0]
        
    elif kind == 'transcript':
        df = df.loc[df.entry_type == 'transcript']
        
        name_pat = 'transcript_name "'
        id_pat = 'transcript_id "'
    
    else:
        raise ValueError('Only genes or transcript u dumb')
        
    df['name'] = df.fields.str.split(pat=name_pat, n=1, expand=True)[1]
    df['name'] = df.name.str.split(pat='"', n=1, expand=True)[0]

    df['id'] = df.fields.str.split(pat=id_pat, n=1, expand=True)[1]
    df['id'] = df.id.str.split(pat='"', n=1, expand=True)[0]

    df['len'] = df.start - df.stop
    df['len'] = df.len.abs()
    
    df.drop('fields', axis=1, inplace=True)
        
    return df    

def make_counts_table(bulk, sc, sample_df, \
                      gtf, kind='gene', novelty='Known'):

    # generate merged T/F table
    df = group_table(bulk, sc, sample_df, kind=kind)
    
    # subset on novelty
    if kind == 'gene':
        nov_col = 'gene_novelty'
    elif kind == 'transcript':
        nov_col = 'transcript_novelty'
    df = df.loc[df[nov_col] == novelty]
    df_copy = df.copy(deep=True)
    
    # count T/F occurrences on the different dataset categories
    if kind == 'gene':
        id_col = 'annot_gene_id'
    elif kind == 'transcript':
        id_col = 'annot_transcript_id'
        
    # groupby and count the bois
    df = df[['Bulk MB', 'Bulk MT', 'sc MB', 'sn MB', 'sn MT', id_col]].groupby(['Bulk MB', 'Bulk MT', 'sc MB', 'sn MB', 'sn MT']).count()
    
    # add annot info
    if novelty == 'Known' and gtf:
        if kind == 'gene':
            gtf_annot = get_gtf_info(gtf, kind='gene')
        elif kind == 'transcript':
            gtf_annot = get_gtf_info(gtf, kind='transcript')
        df_copy = df_copy.merge(gtf_annot, how='left', left_on=id_col, right_on='id')
        
    df_copy.set_index(['Bulk MB', 'Bulk MT', 'sc MB', 'sn MB', 'sn MT'], inplace=True)
        
    return df, df_copy

def group_table(bulk, sc, sample_df, kind='gene'):
        
    # what columns do we need?
    bulk_datasets = sample_df.loc[sample_df.experiment=='bulk', 'dataset'].tolist()
    sc_datasets = sample_df.loc[sample_df.experiment=='sc', 'dataset'].tolist()
    cols = []
    cols = ['annot_gene_id', 'annot_transcript_id', 'gene_novelty', 'transcript_novelty']
    bulk_cols = cols+bulk_datasets
    sc_cols = cols+sc_datasets
        
    # subset on cols we need
    bulk = bulk[bulk_cols]
    sc = sc[sc_cols]

    # merge columns to get datasets all in one df
    df = bulk.merge(sc, how='outer', on=cols)
    df.fillna(value=0, inplace=True)

    # if gene, groupby on gene id and sum counts
    if kind == 'gene':
        cols = ['annot_gene_id', 'gene_novelty']
        agg_dict = {}
        for d in bulk_datasets+sc_datasets:
            agg_dict[d] = 'sum'
        df = df.groupby(cols).agg(agg_dict).reset_index()
    
    # get a dict of dataset ids to category
    dataset_dict = {}
    dataset_dict['Bulk MB'] = sample_df.loc[(sample_df.experiment=='bulk')&(sample_df['sample']=='MB'), 'dataset'].tolist()
    dataset_dict['Bulk MT'] = sample_df.loc[(sample_df.experiment=='bulk')&(sample_df['sample']=='MT'), 'dataset'].tolist()
    dataset_dict['sc MB'] = sample_df.loc[(sample_df.experiment=='sc')&(sample_df['sample']=='MB_cells'), 'dataset'].tolist()
    dataset_dict['sn MB'] = sample_df.loc[(sample_df.experiment=='sc')&(sample_df['sample']=='MB_nuclei'), 'dataset'].tolist()
    dataset_dict['sn MT'] = sample_df.loc[(sample_df.experiment=='sc')&(sample_df['sample']=='MT_nuclei'), 'dataset'].tolist()

    # add T/F columns to denote presence/absence of gene/transcript for each dataset group
    for key in dataset_dict:
        datasets = dataset_dict[key]
        df[key] = df[datasets].any(axis=1)

    # reduce to just the columns we care about
    if kind == 'gene':
        cols = ['annot_gene_id', 'gene_novelty']
    elif kind == 'transcript':
        cols = ['annot_transcript_id', 'transcript_novelty']
    cols = cols+[k for k in dataset_dict]
    df = df[cols]
    
    return df

def make_cond_map(groups, group_names):
    cond_map = dict()
    for group, group_name in zip(groups, group_names):
        for c in group:
            cond_map[c] = group_name
    return cond_map

# calculate the normalized average or sum of TSS expression 
# per cell from the TSS anndata object
def calc_exp(adata, groups, group_names, how='tss', cpm=False):
    
    try:
        adata.var.reset_index(inplace=True)
    except:
        pass
    
    if how == 'tss':
        id_col = 'tss_id'
    elif how == 'iso':
        id_col = 'transcript_id'
        
    # conditions map
    cond_map = make_cond_map(groups, group_names)
    col = 'condition'
    adata.obs[col] = adata.obs.leiden.map(cond_map)
    
    # make df that we can groupby
    colnames = adata.var[id_col].tolist()
    rownames = adata.obs.merged_bc.tolist()    
    raw = adata.X
    gene_names = adata.var.gene_name.tolist()
    df = pd.DataFrame(data=raw, index=rownames, columns=colnames)
    df.reset_index(inplace=True)
    df.rename({'index':'merged_bc'}, axis=1, inplace=True)
    samp = adata.obs[['merged_bc', col]]
    df = df.merge(samp, how='left', on='merged_bc')
    
    # limit to only the cells that we want in this condition
    df[col] = df[col].astype('str')
    df = df.loc[df[col].isin(group_names)]
        
    # groupby sample type and sum over gen
    df.drop('merged_bc', axis=1, inplace=True)
    df = df.groupby(col).sum().reset_index()
    
    if cpm:
        # since these values haven't been normalized yet, do that
        # CPM : (counts/total_counts)* 1**6
        # Note : ATAC values were pre-normalized
        df.set_index(col, inplace=True)
        df = df.transpose()
        for c in group_names:
            total_counts = df[c].sum()
            df[c] = (df[c]/total_counts)*(1^6)
        df = df.transpose()
        df.reset_index(inplace=True)
    
    # melty boi
    tss_cols = df.columns.tolist()[1:]
    df = df.melt(id_vars=col, value_vars=tss_cols)
    
    # rename some cols
    df.rename({'variable':id_col,'value':'counts'}, axis=1, inplace=True)
            
    # add gene name
    if how == 'tss':
        temp = adata.var[[id_col, 'gene_name']]
    df = df.merge(temp, how='left', on=id_col)
    
    return df  

def calc_iso_exp(df, groups, group_names, cpm=False):
    
    counts_cols = [i for j in groups for i in j]
    df = df[counts_cols]
    colnames = df.columns
    rownames = df.index.tolist()
    raw = df.to_numpy()
    
    df = pd.DataFrame(data=raw, index=rownames, columns=colnames)
    df = df.transpose()   
    df.reset_index(inplace=True)
    df.rename({'index':'cluster_id'}, axis=1, inplace=True)
    
    cond_map = make_cond_map(groups, group_names)
    df['condition'] = df.cluster_id.map(cond_map)   
    df.set_index('cluster_id', inplace=True)
    df = df.groupby('condition').sum()
    
    df = df.transpose()    
    if cpm:
        # since these values haven't been normalized yet, do that
        # CPM : (counts/total_counts)* 1**6
        # Note : ATAC values were pre-normalized
        for c in group_names:
            total_counts = df[c].sum()
            df[c] = (df[c]/total_counts)*(1^6)
        
    df = df.rename_axis(None, axis=1)
    df.reset_index(inplace=True)
    df.rename({'index': 'tid'}, axis=1, inplace=True)
    
    # get tss ID
    df['tss_id'] = df['tid'].str.split('_', n=1, expand=True)[1]
    df['tid'] = df['tid'].str.split('_', n=1, expand=True)[0]
    
    df.set_index(['tss_id', 'tid'], inplace=True)
    
    # sort by expression
    df['total'] = df.sum(axis=1)
    df = df.sort_values(by=['tss_id', 'total'], ascending=[True, False])
    df.drop('total', axis=1, inplace=True)
    df.reset_index(inplace=True)
    
    return df   

