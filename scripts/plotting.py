import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import upsetplot
from scipy import stats
from .utils import *

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def get_talon_nov_colors():
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Genomic': '#F0E442'}
    order = ['Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Genomic']
    
    return c_dict, order

def get_priming_colors():
    blue = '#0072B2'
    red_orange = '#D55E00'
    c_dict = {'Random hexamer': red_orange, 'Oligo dT': blue}
    order = ['Oligo dT', 'Random hexamer']
    
    return c_dict, order

def get_lr_cluster_colors():
    purple = '#8DA0CB'
    yellow = '#FFD92F'
    green = '#A6D854'
    pink = '#E78AC3'
    beige = '#E5C494'
    salmon = '#FC8D62'
    teal = '#66C2A5'
    
    c_dict = {'1': yellow, '2': beige, '3': teal,
              '4': salmon, '5': purple, '6': pink, '7': green}
    order = ['1', '2', '3', '4', '5', '6', '7']
    
    return c_dict, order

def get_sample_colors(samples=None):
    sample_green = '#019f73'
    sample_blue = '#57b4e9'
    sample_pink = '#cb79a7'
    c_dict = {'MB_cells': sample_pink, 'MB_nuclei': sample_blue, 'MT_nuclei': sample_green}
    order = ['MB_cells', 'MB_nuclei', 'MT_nuclei']
    
    if samples:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in samples:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in samples]            
    
    return c_dict, order

def get_condition_colors():
    sample_pink = '#cb79a7'
    sample_green = '#019f73'
    sample_dark_green = '#066b4b'
    c_dict = {'MB': sample_pink, 'MNC': sample_green, 'MT': sample_dark_green}
    order = ['MB', 'MNC', 'MT']
    
    return c_dict, order

def get_tech_colors():
    sample_green = '#019f73'
    sample_blue = '#57b4e9'
    sample_pink = '#cb79a7'
    c_dict = {'Bulk': '#EBC046', \
                       'Single-cell': sample_pink,
                       'Single-nucleus': sample_blue}
    order = ['Bulk', 'Single-cell', 'Single-nucleus']
    
    return c_dict, order

def get_2_tech_colors():

    yellow = '#EBC046'
    magenta = 'mediumvioletred'

    c_dict = {'Bulk': yellow, \
              'Single-cell': magenta}
    order = ['Bulk', 'Single-cell']
    
    return c_dict, order

def get_cat_cmap():
    
    # get condition and cluster colors
    clust_cdict, _ = get_lr_cluster_colors()
    cond_cdict, _ = get_condition_colors()
    
    # concatenate dictionaries
    c_dict = {**clust_cdict, **cond_cdict}
    
    # raw values 
    clusters = ['1', '2', '3', '4', '5', '6', '7']
    conditions = ['MB', 'MB', 'MB', 'MNC', 'MNC', 'MT', 'MT']
    data = np.transpose([clusters,conditions])
    df = pd.DataFrame(data=data, 
                      columns=['cluster', 'condition'])
 
    # assign arbitrary numbers to each category (cluster, condition)
    cats = clusters+['MB', 'MNC', 'MT']
    cat_dict = dict([(cat, i) for i, cat in enumerate(cats)])
    df['cluster'] = df.cluster.map(cat_dict)
    df['condition'] = df.condition.map(cat_dict)
    df = df.transpose()
    
    colors = [c_dict[cat] for cat in cats]
    cat_cmap = ListedColormap(colors)
    
    return df, cat_cmap
    
def plot_read_len_kde(df, hue, c_dict, order, opref, common_norm=True):
    sns.set_context('paper', font_scale=2)    
    
    ax = sns.displot(data=df, x='read_length', hue=hue,
                 palette=c_dict, kind='kde', hue_order=order, linewidth=3, 
                 common_norm=common_norm)
    ax.set(xlabel='Read length', ylabel='KDE of reads',
          title='Length distribution of Reads', xlim=(0,7500),
          xticks=[0, 2500, 5000, 7500])
    plt.savefig('{}_read_length_kde.pdf'.format(opref), dpi=300, bbox_inches='tight')
    
# plot proportion of reads per nov. cat per cell
# when grouped by another var
def plot_reads_per_cell_nov(df, hue, c_dict, order, opref):
    sns.set_context('paper', font_scale=2)    
    
    _, nov_order = get_talon_nov_colors()
    
    # groupby
    temp = df[['dataset', hue, 'transcript_novelty', 'read_name']].groupby(['dataset', hue, 'transcript_novelty']).count()
    temp.reset_index(inplace=True)
    temp.rename({'read_name':'counts'}, axis=1, inplace=True)
    
    # calculate proportion nov. cat per cell
    total_temp = df[['read_name', 'dataset']].groupby('dataset').count()
    total_temp.reset_index(inplace=True)
    total_temp.rename({'read_name':'counts'}, axis=1, inplace=True)

    temp['Proportion of reads'] = temp.apply(lambda x: x.counts/total_temp.loc[total_temp.dataset==x.dataset, 'counts'].tolist()[0], axis=1)
    
    g = sns.boxplot(data=temp, x='transcript_novelty', y='Proportion of reads', hue=hue,
                hue_order=order, order=nov_order, saturation=1, palette=c_dict)
    g.set_xticklabels(g.get_xticklabels(),rotation=90)
    g.spines['right'].set_visible(False)
    g.spines['top'].set_visible(False)
    g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('{}_transcript_nov_box.pdf'.format(opref), dpi=300, bbox_inches='tight')

def add_perc(ax, data, feature):
    total = data[feature].sum()
    ylim = ax.get_ylim()[1]
    for p in ax.patches:
        percentage = '{:.1f}%'.format(100 * p.get_height()/total)
        x = p.get_x() + p.get_width() / 2 - 0.45
        y = p.get_y() + p.get_height() + ylim*0.00625
        ax.annotate(percentage, (x, y), size = 12)
        
def plot_read_novelty(df, opref, c_dict, order,
                      ylim=None, title=None, 
                      datasets='all'):
    sns.set_context("paper", font_scale=1.6)
    
    temp = df.copy(deep=True)
    
    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]        
    
    # count number of reads per cat
    temp = temp[['transcript_novelty', 'read_name']].groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'read_name':'counts'}, axis=1, inplace=True)
    print(temp)
    
    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar', 
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Reads')
    g.set_xlabels('Transcript novelty')
    
    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')
    
    if ylim:
        g.set(ylim=(0,ylim))
    
    # add title
    if not title:
        g.fig.suptitle('Reads per novelty category')
    else:
        g.fig.suptitle('{} reads per novelty category'.format(title))
        
    # save figure
    fname = '{}_read_novelty'.format(opref)
    g.savefig(fname+'.pdf', dpi=300)
    
def plot_transcript_novelty(df, oprefix, c_dict, order, \
                            ylim=None, title=None,
                            whitelist=None, datasets='all', save_type='pdf'):
    sns.set_context('paper', font_scale=1.6)
    
    temp = df.copy(deep=True)
    
    # remove transcripts that are not on whitelist
    if whitelist:
        temp = temp.loc[temp.transcript_ID.isin(whitelist)]
    
    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]        
    
    # count number of isoforms per cat
    temp = temp[['transcript_ID', 'transcript_novelty', 'read_name']].groupby(['transcript_ID', 'transcript_novelty']).count()
    temp.reset_index(inplace=True)
    temp.drop('read_name', axis=1, inplace=True)
    temp = temp.groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'transcript_ID': 'counts'}, axis=1, inplace=True)
    print(temp)
    
    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar', 
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Isoforms')
    g.set_xlabels('Transcript novelty')
    
    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')
    
    if ylim:
        g.set(ylim=(0,ylim))
    
    # add title
    if not title:
        g.fig.suptitle('Transcript models per novelty category')
    else:
        g.fig.suptitle('{} transcript models per novelty category'.format(title))
        
    # save figure
    fname = '{}_isoform_novelty'.format(oprefix)
    if save_type == 'png':
        g.savefig(fname+'.png', dpi=300)        
    elif save_type == 'pdf':
        g.savefig(fname+'.pdf', dpi=300)        
    
    plt.show()
    plt.clf()
    
def plot_short_long_det(df, c_dict, order, opref, \
                    xlim, ylim, how='gene'):
    
    sns.set_context('paper', font_scale=2)
    
    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'ill_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'ill_umi_count'
    
    ax = sns.jointplot(data=df, x=c1, y=c2,
                     hue='sample', palette=c_dict,
                     xlim=(0,xlim), ylim=(0,ylim), 
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint
    
    # plot regression lines and equation of regression lines
    # https://stackoverflow.com/questions/48145924/different-colors-for-points-and-line-in-seaborn-regplot/68135585#68135585
    # https://stackoverflow.com/questions/45902739/seaborn-annotate-the-linear-regression-equation
    # https://stackoverflow.com/questions/62705904/add-entry-to-matplotlib-legend-without-plotting-an-object
    lines = []
    labels = []
    for s in df['sample'].unique().tolist():
        temp = df.loc[df['sample'] == s]
        color = c_dict[s]
        line_color = adjust_lightness(color, 0.5)
        
        # get coeffs of linear fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(temp[c1],temp[c2])
        lines += [mpl.lines.Line2D([0], [0], color=line_color)]
        labels += ['m={0:.1f}'.format(slope)]

        print('Slope of {} correlation: {}'.format(s, slope))
        
        sns.regplot(data=temp, x=c1, y=c2,
                    scatter=False, ax=ax, color=color)
        sns.regplot(data=temp, x=c1, y=c2,
            scatter=False, ax=ax, color=color, ci=0,
            line_kws={'color':line_color,
                      'linestyle':'-',
                      'label':"m={0:.1f}".format(slope)})
    
    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    
    
def plot_detection_venn(bulk, sc, opref, \
                        gene_nov=None, \
                        transcript_nov=None, \
                        sample='MB', cell_part='cell'):
    
    sns.set_context('paper', font_scale=1.8)
    
    # colors
    known_out_green = '#90D6C3'
    known_int_green = '#009E73'
    nnc_out_gold = '#F5DFAE'
    nnc_int_gold = '#E69F00'
    nic_out_orange = '#DEA67A'
    nic_int_orange = '#D55E00'
    mt_green = '#019f73'
    mb_pink = '#cb79a7'
    mb_blue = '#57b4e9'
    
    sc_datasets = get_dataset_names(sc)
    bulk_datasets = get_dataset_names(bulk)
    sample_df = get_sample_df(sc_datasets+bulk_datasets)
    
    if not gene_nov and not transcript_nov:
        print('must choose one!')
        return
    elif gene_nov:
        bulk = bulk.loc[bulk.gene_novelty == gene_nov]
        sc = sc.loc[sc.gene_novelty == gene_nov]
        id_col = 'annot_gene_id'
        out_color = known_out_green
        int_color = known_int_green
        
    elif transcript_nov:
        bulk = bulk.loc[bulk.transcript_novelty == transcript_nov]
        sc = sc.loc[sc.transcript_novelty == transcript_nov]
        id_col = 'annot_transcript_id'
        
        if transcript_nov == 'Known':
            out_color = known_out_green
            int_color = known_int_green
        elif transcript_nov == 'NNC':
            out_color = nnc_out_gold
            int_color = nnc_int_gold
        elif transcript_nov == 'NIC':
            out_color = nic_out_orange
            int_color = nic_int_orange
        
    # which genes/transcripts are in MB vs. MT?
    if sample == 'MB' and cell_part == 'cell':
        bulk_sample = 'MB'
        sc_sample = 'MB_cells'
        celltype_color = mb_pink
    elif sample == 'MB' and cell_part == 'nucleus':
        bulk_sample = 'MB'
        sc_sample = 'MB_nuclei'
        celltype_color = mb_blue
    elif sample == 'MT' and cell_part == 'nucleus':
        bulk_sample = 'MT'
        sc_sample = 'MT_nuclei'
        celltype_color = mt_green
    else:
        print('No cell myotubes in dataset')
        return
        
    bulk_datasets = sample_df.loc[(sample_df['sample'] == bulk_sample)&(sample_df.experiment == 'bulk'), 'dataset'].tolist()
    sc_datasets = sample_df.loc[(sample_df['sample'] == sc_sample)&(sample_df.experiment == 'sc'), 'dataset'].tolist()
    
    # filter abundance files based on which datasets we care about 
    bulk = bulk.loc[bulk[bulk_datasets].any(axis=1)]
    sc = sc.loc[sc[sc_datasets].any(axis=1)]
    
    # get list of unique gene or transcript ids from bulk/sc
    bulk_ids = bulk[id_col].unique().tolist()
    sc_ids = sc[id_col].unique().tolist()
    
    # intersection
    intersection = list(set(bulk_ids) & set(sc_ids))
    n_bulk_sc = len(intersection)
    n_bulk = len(bulk_ids) - n_bulk_sc
    n_sc = len(sc_ids) - n_bulk_sc
    counts = [n_bulk, n_sc, n_bulk_sc]
    log_counts = [np.log2(n) for n in counts]
    log_counts = tuple(counts)
    
    v = venn2(subsets=log_counts, set_labels=('',''))
    v.get_patch_by_id('10').set_color(out_color)
    v.get_patch_by_id('01').set_color(out_color)
    v.get_patch_by_id('11').set_color(int_color)
    v.get_patch_by_id('10').set_edgecolor(celltype_color)
    v.get_patch_by_id('01').set_edgecolor(celltype_color)
    v.get_patch_by_id('11').set_edgecolor(celltype_color)
    v.get_patch_by_id('10').set_linewidth(5)
    v.get_patch_by_id('01').set_linewidth(5)
    v.get_patch_by_id('11').set_linewidth(5)
    v.get_patch_by_id('10').set_alpha(1)
    v.get_patch_by_id('01').set_alpha(1)
    v.get_patch_by_id('11').set_alpha(1)
    v.get_label_by_id('10').set_text(counts[0])
    v.get_label_by_id('01').set_text(counts[1])
    v.get_label_by_id('11').set_text(counts[2])
    
    # save da ting
    fname = '{}'.format(opref)
    fname += '_{}_{}'.format(sample, cell_part)
    if gene_nov:
        fname += '_{}_genes'.format(gene_nov)
    if transcript_nov:
        fname += '_{}_transcripts'.format(transcript_nov)
    fname += '_venn.pdf'
    plt.savefig(fname, dpi=300, bbox_inches='tight') 

def plot_upset_plot(bulk, sc, opref, gtf=None, \
                    kind='gene', novelty='Known'):
    
    sns.set_context('paper', font_scale=1)    
    
    sc_datasets = get_dataset_names(sc)
    bulk_datasets = get_dataset_names(bulk)
    sample_df = get_sample_df(sc_datasets+bulk_datasets)

    # colors
    known_green = '#009E73'
    nnc_gold = '#E69F00'
    nic_orange = '#D55E00'
    
    if novelty == 'Known':
        color = known_green
    elif novelty == 'NNC':
        color = nnc_gold
    elif novelty == 'NIC':
        color = nic_orange

    # df is table with counts, df_copy is df right before groupby and counting
    if novelty != 'Known':
        gtf = None
    df, df_copy = make_counts_table(bulk, sc, sample_df, gtf, kind=kind, novelty=novelty)
    
    # what column?
    if kind == 'gene':
        id_col = 'annot_gene_id'
        len_col = 'Gene length'
        df_copy.rename({'len': len_col}, axis=1, inplace=True)  
    elif kind == 'transcript':
        id_col = 'annot_transcript_id'
        len_col = 'Transcript length'
        df_copy.rename({'len': len_col}, axis=1, inplace=True)
    if novelty == 'Known':
        nov = 'known'
    elif novelty == 'NNC':
        nov = 'NNC'
    elif novelty == 'NIC':
        nov = 'NIC'
        
    ylab = 'Number of {} {}s'.format(nov, kind)
          
    # plot de plot
    blot = upsetplot.UpSet(df_copy, subset_size='auto', show_counts='%d', sort_by='cardinality')
    
    if novelty == 'Known':
        blot.add_catplot(value=len_col, kind='box', color=color, fliersize=1, linewidth=1)
        
    ax_dict = blot.plot()
    ax_dict['intersections'].set_ylabel(ylab)
    
    if novelty == 'Known':
#         ax_dict['extra1'].set_ylim((-10000,200000))
        ax_dict['extra1'].set_yscale('log')

#         blot.plot()

    f = plt.gcf()
    fname = '{}_{}_{}_detection_upset.pdf'.format(opref, novelty, kind)
    f.savefig(fname, dpi=300, bbox_inches='tight')  
        
    return(df, df_copy)

def plot_iso_complexity_vs_reads(adata, gene_adata, obs_col, xlim, ylim, c_dict, order, opref, how='iso'):
    sns.set_context('paper', font_scale=2)
    
    # first format the data such that there's information about how 
    # many genes have more than one iso/TSS/TES per gene
    df = get_multiple_iso_genes(adata, gene_adata)

    c1 = 'gene_n_counts'
    c2 = 'n_genes_multiple_iso'

    ax = sns.jointplot(data=df, x=c1, y=c2,
                     hue=obs_col, palette=c_dict,
                     xlim=(-100,xlim), ylim=(-10,ylim), 
                     marginal_kws={'common_norm':False},
                     joint_kws={'data':df, 's':10, 'alpha':.5})
    ax = ax.ax_joint
    
    lines = []
    labels = []
    for s in df[obs_col].unique().tolist():
        temp = df.loc[df[obs_col] == s]
        color = c_dict[s]
        line_color = adjust_lightness(color, 0.5)

        c1 = 'gene_n_counts'
        c2 = 'n_genes_multiple_iso'

        # get coeffs of linear fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(temp[c1],temp[c2])
        lines += [mpl.lines.Line2D([0], [0], color=line_color)]
        labels += ['m={0:.1f}'.format(slope)]

        print('Slope of {} correlation: {}'.format(s, slope))

        sns.regplot(data=temp, x=c1, y=c2,
                    scatter=False, ax=ax, color=color)
        sns.regplot(data=temp, x=c1, y=c2,
            scatter=False, ax=ax, color=color, ci=0,
            line_kws={'color':line_color,
                      'linestyle':'-',
                      'label':"m={0:.1f}".format(slope)})

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if how == 'iso':
        _ = ax.set(xlabel='# reads/cell', ylabel='# Genes/cell with >1 isoform')
    elif how == 'tss':
        _ = ax.set(xlabel='# reads/cell', ylabel='# Genes/cell with >1 TSS')
    elif how == 'tes':
        _ = ax.set(xlabel='# reads/cell', ylabel='# Genes/cell with >1 TES')
       
    fname = '{}_{}_n_{}_genes_read_depth.pdf'.format(opref, obs_col, how)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
    return df

def plot_ends_iso_cell(df, tss_df, opref, kind='tss', xlim=None, ylim=None):
    sns.set_context("paper", font_scale=1.6)
    
    i_df = get_illumina_metadata()
    df = df.merge(i_df, how='left', left_on='dataset', right_on='raw_bc')

    # merge read annot df with ends
    df = df.merge(tss_df, how='left', left_on='read_name', right_on='read_name')

    # only known transcripts
    df = df.loc[df.transcript_novelty == 'Known']

    # get isoforms per gene
    df.rename({"merged_bc":"cell_ID"}, axis=1, inplace=True)
    ts_id_col = 'peak_id'
    iso_no_dups = df[["cell_ID","gene_ID","transcript_ID"]].drop_duplicates()
    n_isoforms = iso_no_dups.groupby(["cell_ID","gene_ID"]).count().reset_index()
    n_isoforms = n_isoforms.rename(columns = {"transcript_ID":"n_isoforms"})
    n_isoforms = pd.merge(n_isoforms, df[["gene_ID","annot_gene_name"]].drop_duplicates(),
                          how = "left")

    # get TSSs/TESs per gene
    TS_no_dups = df[["cell_ID","gene_ID",ts_id_col]].drop_duplicates()
    n_TS = TS_no_dups.groupby(["cell_ID","gene_ID"]).count().reset_index()
    n_TS = n_TS.rename(columns = {ts_id_col:"n_sites"})

    # merge
    all_data = pd.merge(n_isoforms, n_TS, on = ['gene_ID','cell_ID'], 
                        how = 'left')

    # count number of isoforms and sites 

    # first remove everything with 0 TSSs
    df = all_data.loc[all_data.n_sites > 0]
    parent_df = df.copy(deep=True)
    df = df[['gene_ID', 'n_isoforms', 'n_sites']].groupby(['n_isoforms', 'n_sites']).count()
    df.reset_index(inplace=True)
    df.rename({'gene_ID': 'counts'}, axis=1, inplace=True)

    # log 
    df['log2counts'] = np.log10(df.counts)
    
    if kind == 'tss':
        color = '#56B4E9'
        title = 'TSS'
    elif kind == 'tes':
        color = '#E69F00'
        title = 'TES'
       
    # jointgrid dotplot
    g = sns.JointGrid(xlim=(0.25, xlim+0.5), ylim=(0.5, ylim+0.5), marginal_ticks=True)
    palette = sns.color_palette("dark:"+color, as_cmap=True).reversed()
    sns.scatterplot(data=df, x='n_isoforms', y='n_sites',\
                    hue='log2counts', palette=palette,
                    size='log2counts', sizes=(50,300),\
                    ax=g.ax_joint)

    sns.histplot(data=parent_df, x='n_isoforms', ax=g.ax_marg_x,\
                 binwidth=1, discrete=True, color=color, edgecolor=None)
    sns.histplot(data=parent_df, y='n_sites', ax=g.ax_marg_y,\
                 binwidth=1, discrete=True, color=color, edgecolor=None)

    g.ax_joint.set_xticks(range(1, xlim+1))
    g.ax_joint.set_yticks(range(1, ylim+1))

    g.ax_marg_x.set_yscale('log')
    g.ax_marg_x.set_yticks([1, 100, 10000])
    g.ax_marg_x.set_ylabel('')
    g.ax_marg_y.set_xscale('log')
    g.ax_marg_y.set_xticks([1, 100, 10000])
    g.ax_marg_y.set_xlabel('')

    plt.setp(g.ax_marg_y.get_xticklabels(), fontsize=12)
    plt.setp(g.ax_marg_x.get_yticklabels(), fontsize=12)

    g.set_axis_labels('x', 'y', fontsize=16)
    g.set_axis_labels("Known splice isoforms detected per gene per cell", "{} per gene per cell".format(title))
    g.ax_joint.legend(bbox_to_anchor=(1.6, 1.035), title='log2(counts)')

    fname = '{}_{}_iso_cell.pdf'.format(opref, kind)
    plt.savefig(fname, dpi=300, bbox_inches='tight')


def plot_len_boxplot(sc, bulk, gtf, opref, novelty='Known', kind='gene'):
    
    sns.set_context('paper', font_scale=2)
    
    # colors
    known_green = '#009E73'
    
    sc_datasets = get_dataset_names(sc)
    bulk_datasets = get_dataset_names(bulk)
    sample_df = get_sample_df(sc_datasets+bulk_datasets)
    color = known_green
    
    df_copy, df = make_counts_table(bulk, sc, sample_df, gtf, kind=kind, novelty=novelty)
    
    try:
        df.reset_index(inplace=True)
    except:
        pass
    
    df['kind'] = np.nan
    order = ['Bulk only', 'sc only', 'sn only']
    df.loc[(df['sc MB']==False)&(df['sn MB']==False)&(df['sn MT']==False), 'kind'] = 'Bulk only'
    df.loc[(df['Bulk MB']==False)&(df['Bulk MT']==False)&(df['sn MB']==False)&(df['sn MT']==False), 'kind'] = 'sc only'
    df.loc[(df['Bulk MB']==False)&(df['Bulk MT']==False)&(df['sc MB']==False), 'kind'] = 'sn only'
      
    if kind == 'gene':
        y = 'Gene length'
    elif kind == 'transcript':
        y = 'Transcript length'
    df.rename({'len':y}, axis=1, inplace=True)
    
    ax = sns.boxplot(data=df, x='kind', y=y, color=color, saturation=1, order=order)
    ax.set_yscale('log')
    ax.set_xlabel('')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    fname = '{}_{}_{}_detection_bulk_sc_boxplot.pdf'.format(opref, novelty, kind)
    plt.savefig(fname, dpi=300, bbox_inches='tight')  
    
    plt.show()
    
def plot_exon_hist(df, opref, xlim=30):

    sns.set_context('paper', font_scale=2)
    c_dict, order = get_tech_colors()
    
    # loop through each technology
    for t in df.tech.unique():
        temp = df.loc[df.tech==t]
        ax = sns.displot(data=temp, x='n_exons', hue='tech',
                         palette=c_dict, kind='hist',
                         hue_order=order, discrete=True,
                         linewidth=0, alpha=1)
        plt.xlim(0,xlim)
        ax.set(xlabel='Number of exons', ylabel='Number of reads',
              title='')
        plt.savefig('{}_{}_n_exons_bulk_hist.pdf'.format(opref, t), dpi=300, bbox_inches='tight')
        
def plot_umis_v_barcodes(df, opref):	
    sns.set_context("paper", font_scale=1.8)
    bc_cols = ['bc1', 'bc2', 'bc3', 'state']

    # only want unique bc/umi combos
    temp = df[bc_cols+['umi']].drop_duplicates()

    # get the number of unique bc/umi combos
    temp = temp[bc_cols+['umi']].groupby(bc_cols).count()
    temp.reset_index(inplace=True)
    temp.rename({'umi':'counts'}, axis=1, inplace=True)
    temp.sort_values(by=['state', 'counts'], ascending=False, inplace=True)
    
    # add x coordinate by count and state
    temp['x_coord'] = np.nan
    temp.loc[temp.state == 'Pre-correction', 'x_coord'] = range(len(temp.loc[temp.state=='Pre-correction'].index))
    temp.loc[temp.state == 'Post-correction', 'x_coord'] = range(len(temp.loc[temp.state=='Post-correction'].index))

    # plot
    ax = sns.lineplot(data=temp, x='x_coord', y='counts',
                      hue='state', linewidth=3)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlabel('Ranked cells by # UMIs (logscale)')
    ax.set_ylabel('# UMIs (logscale)')

    plt.tight_layout()
    
    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    fname = '{}_umis_v_barcodes.pdf'.format(opref)
    plt.savefig(fname, dpi=300)
                
def plot_umis_per_cell(df, opref):	   
    bc_cols = ['bc1', 'bc2', 'bc3', 'state']

    # get the number of reads per barcode
    temp1 = df[bc_cols+['umi']].groupby(bc_cols).count()
    temp1.reset_index(inplace=True)
    temp1.rename({'umi':'reads'}, axis=1, inplace=True)
    temp1.sort_values(by='reads', ascending=False, inplace=True)

    # get the number of unique umis per barcode
    temp2 = df[bc_cols+['umi']].drop_duplicates()
    temp2 = temp2[bc_cols+['umi']].groupby(bc_cols).count()
    temp2.reset_index(inplace=True)
    temp2.rename({'umi':'umis'}, axis=1, inplace=True)
    temp2.sort_values(by='umis', ascending=False, inplace=True)

    # merge on barcode
    temp = temp1.merge(temp2, on=bc_cols)

    bins = [i for i in range(0, temp.reads.max(), 1000)]
    temp_reads = temp.reads.values.tolist()
    temp_bins = np.digitize(temp_reads, bins)
    temp['bin'] = temp_bins

    bins.append(temp.reads.max())
    bin_dict = dict([(bin, bin_num) for bin, bin_num in zip([i for i in range(1,len(bins))], bins)])
    temp['bin_total'] = temp['bin'].map(temp['bin'].value_counts())
    temp['bin_reads'] = temp.bin.map(bin_dict)

    # groupby the bin and get the median of number of umis
    temp = temp[['state', 'bin_total', 'bin_reads', 'umis']].groupby(['state', 'bin_total', 'bin_reads']).median()
    temp.reset_index(inplace=True)
    temp.rename({'umis': 'median_umis'}, axis=1, inplace=True)
    temp.sort_values(by='bin_reads', inplace=True)
    
    print(temp.head())

    # plot de plot
    ax = sns.lineplot(x='bin_reads', y='median_umis', alpha=0.5, linewidth=3,
                      hue='state', marker='o', data=temp)
    ax.set_ylabel('Median UMIs per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')
    
    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_title('')
    plt.draw()

    plt.tight_layout()

    fname = '{}_median_umis_v_reads.pdf'.format(opref)
    plt.savefig(fname, dpi=300)
    
def plot_ends_iso(df, tss_df, opref, kind='tss', xlim=None, ylim=None):
    
    # merge read annot df with TSSs
    df = df.merge(tss_df, how='left', left_on='read_name', right_on='read_name')

    # only known
    df = df.loc[df.transcript_novelty == 'Known']


    ts_id_col = 'peak_id'

    # get isoforms per gene
    iso_no_dups = df[["gene_ID","transcript_ID"]].drop_duplicates()
    n_isoforms = iso_no_dups.groupby(["gene_ID"]).count().reset_index()
    n_isoforms = n_isoforms.rename(columns = {"transcript_ID":"n_isoforms"})
    n_isoforms = pd.merge(n_isoforms, df[["gene_ID","annot_gene_name"]].drop_duplicates(),
                          how = "left")

    # get TSSs/TESs per gene
    TS_no_dups = df[["gene_ID",ts_id_col]].drop_duplicates()
    n_TS = TS_no_dups.groupby(["gene_ID"]).count().reset_index()
    n_TS = n_TS.rename(columns = {ts_id_col:"n_sites"})

    # merge
    all_data = pd.merge(n_isoforms, n_TS, on = ['gene_ID'], 
                        how = 'left')

    # count number of isoforms and sites 

    # first remove everything with 0 TSSs
    df = all_data.loc[all_data.n_sites > 0]
    parent_df = df.copy(deep=True)

    df = df[['gene_ID', 'n_isoforms', 'n_sites']].groupby(['n_isoforms', 'n_sites']).count()
    df.reset_index(inplace=True)
    df.rename({'gene_ID': 'counts'}, axis=1, inplace=True)

    # log 
    df['log2counts'] = np.log10(df.counts)
    
    # plot dot plot
    if kind == 'tss':
        color = '#56B4E9'
        title = 'TSS'
    elif kind == 'tes':
        color = '#E69F00'
        title = 'TES'
       
    # jointgrid dotplot
    g = sns.JointGrid(xlim=(0.25, xlim+0.5), ylim=(0.5, ylim+0.5), marginal_ticks=True)
    palette = sns.color_palette("dark:"+color, as_cmap=True).reversed()
    sns.scatterplot(data=df, x='n_isoforms', y='n_sites',\
                    hue='log2counts', palette=palette,
                    size='log2counts', sizes=(50,300),\
                    ax=g.ax_joint)

    sns.histplot(data=parent_df, x='n_isoforms', ax=g.ax_marg_x,\
                 binwidth=1, discrete=True, color=color, edgecolor=None)
    sns.histplot(data=parent_df, y='n_sites', ax=g.ax_marg_y,\
                 binwidth=1, discrete=True, color=color, edgecolor=None)

    g.ax_joint.set_xticks(range(1, xlim+1))
    g.ax_joint.set_yticks(range(1, ylim+1))

    g.ax_marg_x.set_yscale('log')
    g.ax_marg_x.set_yticks([1, 100, 10000])
    g.ax_marg_x.set_ylabel('')
    g.ax_marg_y.set_xscale('log')
    g.ax_marg_y.set_xticks([1, 100, 10000])
    g.ax_marg_y.set_xlabel('')

    plt.setp(g.ax_marg_y.get_xticklabels(), fontsize=12)
    plt.setp(g.ax_marg_x.get_yticklabels(), fontsize=12)



    g.set_axis_labels('x', 'y', fontsize=16)
    g.set_axis_labels("Known splice isoforms detected per gene", "{}s per gene".format(title))
    g.ax_joint.legend(bbox_to_anchor=(1.6, 1.035), title='log2(counts)')

    fname = '{}_{}_iso_cell.pdf'.format(opref, kind)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_tss_iso_heatmap(adata, sg_t_df, groups, group_names, gname, opref):
    
    # calculate TSS expression per condition
    tss_df = calc_exp(adata, groups, group_names, how='tss')
    
    # subset by gene and calculate DPI per gene
    tss_df = tss_df.loc[tss_df.gene_name == gname]
    tss_df.drop(['gene_name'], axis=1, inplace=True)
    tss_df = tss_df.pivot(index='tss_id', columns='condition', values='counts')
    tss_df = tss_df.div(tss_df.sum(axis=0), axis=1)
    
    # subset by gene and calculate iso expression per condition
    iso_df = sg_t_df.loc[sg_t_df.gname == gname]
    iso_df = calc_iso_exp(iso_df, groups, group_names, cpm=True)
    
    # get categorical colormap
    mini_obs, cat_cmap = get_cat_cmap()
    
    # plot the figure
    sns.set(rc={'figure.figsize':(12,7)})
    sns.set_context('paper', font_scale=1.5)
    fig = plt.figure()
    
    vmax = iso_df[group_names].max().max()
    n_tss = len(iso_df.tss_id.unique())

    # complicated subplot stuff
    tss_ax = plt.subplot2grid((n_tss+1,2), loc=(0,0), rowspan=n_tss)

    axes = []
    for i in range(1,n_tss+1):
        axes.append(fig.add_subplot(n_tss+1, 2, i*2))

    # fig, axes = plt.subplots(nrows=4)
    fig.subplots_adjust(hspace=0.00)
    fig.subplots_adjust(wspace=0.05)

    # plot tss only plot
    sns.heatmap(tss_df, cmap='magma', ax=tss_ax, cbar=False)
    tss_ax.set_yticklabels(tss_ax.get_yticklabels(), rotation=0)
    tss_ax.set_ylabel('')
    tss_ax.set_xlabel('')

    # plots for each iso per TSS
    for tss, ax in zip(iso_df.tss_id.unique().tolist(),axes):
        temp = iso_df.loc[iso_df.tss_id == tss]
        temp.drop('tss_id', axis=1, inplace=True)
        temp.set_index('tid', inplace=True)

        sns.heatmap(temp, yticklabels=True, cmap='viridis',
                vmin=0, vmax=5,
                cbar=False, ax=ax)

        plt.yticks(rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
        ax.yaxis.set_label_position("right")
        ax.tick_params(right=True, labelright=True, rotation=0)
        ax.tick_params(left=False, labelleft=False, rotation=0)

        ax.set_ylabel('')
        ax.set_xlabel('')

    # plot sample labels
    tss_colorbar_ax = fig.add_subplot((n_tss+1)*2,2,((n_tss+1)*4)-3)
    sns.heatmap(mini_obs, cmap=cat_cmap,
                ax=tss_colorbar_ax, cbar=False)
    tss_colorbar_ax.set_ylabel('')
    tss_colorbar_ax.set_xlabel('')
    tss_colorbar_ax.tick_params(left=False, labelleft=False, rotation=0)
    tss_colorbar_ax.tick_params(right=False, labelright=False, rotation=0)
    tss_colorbar_ax.set_xticklabels('')

    iso_colorbar_ax = fig.add_subplot((n_tss+1)*2,2,((n_tss+1)*4)-2)
    sns.heatmap(mini_obs, cmap=cat_cmap,
                ax=iso_colorbar_ax, cbar=False)
    iso_colorbar_ax.set_ylabel('')
    iso_colorbar_ax.set_xlabel('')
    iso_colorbar_ax.tick_params(left=False, labelleft=False, rotation=0)
    iso_colorbar_ax.set_xticklabels('')

    # plot colorbars
    tss_colorbar_ax = fig.add_subplot((n_tss+1)*5,2,((n_tss+1)*10)-1)
    cmap = plt.get_cmap('magma')
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cb = mpl.colorbar.ColorbarBase(tss_colorbar_ax, cmap=cmap,
                                  norm=norm, orientation='horizontal')
    cb.set_label('Proportion TSS usage')

    iso_colorbar_ax = fig.add_subplot((n_tss+1)*5,2,((n_tss+1)*10))
    cmap = plt.get_cmap('viridis')
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(iso_colorbar_ax, cmap=cmap,
                                  norm=norm, orientation='horizontal')
    cb.set_label('TPM')
    
    fname = '{}_{}_heatmap.pdf'.format(opref, gname)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    