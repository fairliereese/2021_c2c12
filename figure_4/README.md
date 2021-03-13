# Figure 4

#### For the TSS calling for both single-cell and bulk long-read
* BAM files output from minimap2 (from the processing steps outlined in [figure 1](https://github.com/fairliereese/2021_c2c12/tree/master/figure_1)) were filtered for reads that were assigned the known, NIC, NNC, or prefix-ISM novelty category by TALON that were identified via the read_annot TALON output, also from the steps in [figure 1](https://github.com/fairliereese/2021_c2c12/tree/master/figure_1)
* TSS peak regions were called on the filtered BAM with the [ENCODE PacBio TSS caller](https://github.com/ENCODE-AWG/tss-annotation/blob/master/long_read/pacbio_to_tss.py)
* The output TSS peak regions were filtered for those that were supported by more than 1 read
* Reads in the corresponding TALON-output read_annot files were associated with a TSS using bedtools intersect. Reads that did not intersect a TSS peak were not assigned a TSS
* Using the gene IDs present in the TALON read_annot files, each TSS was filtered on the gene level where each TSS was required to have a number of reads >10% (>5% for bulk) of the number of reads that supported the most highly-expressed TSS within the same genes

#### For the TES calling for both single-cell and bulk long-read
* BAM files output from minimap2 (from the processing steps outlined in [figure 1](https://github.com/fairliereese/2021_c2c12/tree/master/figure_1)) were filtered for reads that were assigned the known, NIC, NNC, or suffix-ISM novelty category by TALON that were identified via the read_annot TALON output, also from the steps in [figure 1](https://github.com/fairliereese/2021_c2c12/tree/master/figure_1)
* TES peak regions were called on the filtered BAM with the [ENCODE PacBio TSS caller](https://github.com/ENCODE-AWG/tss-annotation/blob/master/long_read/pacbio_to_tss.py)
* The output TES peak regions were filtered for those that were supported by more than 1 read
* Reads in the corresponding TALON-output read_annot files were associated with a TES using bedtools intersect. Reads that did not intersect a TES peak were not assigned a TES
* Using the gene IDs present in the TALON read_annot files, each TES was filtered on the gene level where each TES was required to have a number of reads >50% of the number of reads that supported the most highly-expressed TES within the same genes

#### Creating a Scanpy AnnData object for TSS expression
* After associating each read with a filtered TSS, reads per TSS per cell were summed up from the read_annot file
* Cell metadata was loaded in from the gene-level AnnData object and used along with the TSS expression to make a TSS-level AnnData object

```python
import pandas as pd
import sys
import os
import scanpy as sc
import swan_vis as swan

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from scripts.utils import *
from scripts.plotting import *
```

```python
# read in the data relevant for this figure

# output from TALON
def get_sc_data():
    fname = '../processing/talon/sc_talon_read_annot.tsv'

    df = pd.read_csv(fname, sep='\t')
    
    return df

# output from TSS calling script
def get_tss_bed():

    fname = '../processing/ends/sc_tss.bed'
    df = pd.read_csv(fname, sep='\t', header=None, usecols=[3,9])
    df.columns = ['peak_id', 'read_name']
    
    return df

# output from calculating TSS expression per cell 
# and adding it to scanpy
def get_tss_raw_adata():
    fname = '../processing/scanpy/sc_tss_raw.h5ad'
    adata = sc.read(fname)
    
    return adata

# output from Swan to get expression of TSS AND transcript
# note: this could have also been done from the adata file
# which is why Swan is not mentioned in the methods
# it just happens that I wrote the plotting and data formatting
# code for some of these figures using the Swan output
def get_swan_data():
    fname = '../processing/swan/sc_tss_iso_swan.p'
    sg = swan.SwanGraph(fname)
    t_df = sg.t_df
    cell_type_map = {'Cluster_1_counts': '1', 
                 'Cluster_2_counts': '2', 
                 'Cluster_3_counts': '3', 
                 'Cluster_4_counts': '4', 
                 'Cluster_5_counts': '5', 
                 'Cluster_6_counts': '6', 
                 'Cluster_7_counts': '7'}
    t_df.rename(cell_type_map, axis=1, inplace=True)
    return t_df
```

### Panel 4C


```python
df = get_sc_data()
tss = get_tss_bed()

ylim = 11
xlim = 5

opref = 'figures/sc_tss'

plot_ends_iso_cell(df, tss, opref, kind='tss', xlim=xlim, ylim=ylim)
```


    
![png](figures/output_3_0.png)
    


### Panel 4E


```python
groups = ['1', '2', '3', '4', '5', '6', '7']
group_names = groups
adata = get_tss_raw_adata()
t_df = get_swan_data()

plot_tss_iso_heatmap(adata, t_df, groups, group_names, 'Tnnt2', 'figures/tss_iso')
```

    Graph from ../processing/swan/sc_tss_iso_swan.p loaded


    /Users/fairliereese/miniconda3/lib/python3.7/site-packages/pandas/core/frame.py:4169: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      errors=errors,



    
![png](figures/output_5_2.png)
    

