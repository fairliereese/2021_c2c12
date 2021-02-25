import matplotlib.pyplot as plt
import seaborn as sns

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
    c_dict = {'Random hexamer': red_orange, 'Poly dT': blue}
    order = ['Poly dT', 'Random hexamer']
    
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

def get_tech_colors():
    sample_green = '#019f73'
    sample_blue = '#57b4e9'
    sample_pink = '#cb79a7'
    c_dict = {'Bulk': '#EBC046', \
                       'Single-cell': sample_pink,
                       'Single-nucleus': sample_blue}
    order = ['Bulk', 'Single-cell', 'Single-nucleus']
    
    return c_dict, order

def plot_read_len_kde(df, hue, c_dict, order, opref, common_norm=True):
    sns.set_context('paper', font_scale=2)    
    
    ax = sns.displot(data=df, x='read_length', hue=hue,
                 palette=c_dict, kind='kde', hue_order=order, linewidth=3, 
                 common_norm=common_norm)
    ax.set(xlabel='Read length', ylabel='KDE of reads',
          title='Length distribution of Reads')
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
                            whitelist=None, datasets='all', save_type='png'):
    
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