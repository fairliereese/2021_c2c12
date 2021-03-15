# Figure 1

#### For the long-read single-cell data
* CCS reads were generated from raw PacBio reads usin the CCS softare from the SMRT analysis software suite (parameters: --skip-polish --min-length=10 --min-passes=3 --min-rq=0.9 --min-snr=2.5) 
* The Split-seq adapters were identified and removed using Lima (v2.0.0) (parameters: --ccs --min-score 0  --min-end-score 0 --min-signal-increase 0 --min-score-lead 0) 
* Full-length non-chimeric reads were generated with IsoSeq3's Refine (v3.4.0)
* Reads were demultiplexed for their Split-seq barcodes using a custom script (https://github.com/fairliereese/pacbio-splitpipe)
* Demultiplexed reads were filtered for those where the same combination of barcodes were also detected in the correspoding short-read experiment
* Reads were mapped to the mm10 reference genome using Minimap2 (v2.17-r94) (-ax splice:hq -uf --MD)
* Long-read sequencing artifacts were corrected with TranscriptClean (--canonOnly --primaryOnly)
* Reads were annotated to their transcripts of origin using TALON (https://github.com/mortazavilab/TALON/tree/cb_tag) (--cb) and the GENCODE vM21 annotation. Each cell is its own dataset in the output TALON database.
* NIC and NCC transcripts were filtered for those that were seen in 4 or more half-cells (as at this point the data is still represented by its barcode which could either be oligo dT primed or randomly primed)
* Random hexamer and oligo dT cell halves were merged 
* Scanpy was used to filter, create a UMAP, cluster, and normalize the single-cell data 

#### For the long-read bulk data
* Bulk PacBio data was processed following the [ENCODE Long-read RNA-seq analysis protocol](https://www.encodeproject.org/documents/a84b4146-9e2d-4121-8c0c-1b6957a13fbf) for Mouse Samples (v1.0) for CCS, Lima, Refine, and TranscriptClean steps. 
* TALON was run (--cov 0.9 --identity 0.8)
* Novel transcript models were filtered (--minCount 5 --minDatasets 2 --maxFracA 0.5)
* Transcript abundances were quantified using the talon_abundance module


## Figures made in Python

The following figures were made in Python and the relevant notebook used to generate them can be seen in [figure_1.ipynb](https://github.com/fairliereese/2021_c2c12/blob/master/figure_1/figure_1.ipynb)

```python
import sys
import os
import scanpy as sc
import pandas as pd

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from scripts.utils import *
from scripts.plotting import *
```


```python
# read in the data relevant for this figure

# output from the TALON run
def get_sc_data():
    fname = '../processing/talon/sc_talon_read_annot.tsv'

    df = pd.read_csv(fname, sep='\t')
    df = add_read_annot_metadata(df)
    
    return df

# output from transforming the unfiltered TALON abundanca matrix into 
# a scanpy AnnData object
def get_sc_adata():
    fname = '../processing/scanpy/sc_gene.h5ad'
    adata = sc.read(fname)
    
    return adata
    
# output from TALON filtering
def get_sc_whitelist():
    fname = '../processing/talon/sc_whitelist.csv'
    whitelist = read_whitelist(fname)
    
    return whitelist

def get_sc_bulk_data():
    fname = '../processing/talon/bulk_sc_talon_read_annot.tsv'

    df = pd.read_csv(fname, sep='\t')
    df = add_read_annot_metadata(df, bulk=True)
    
    return df
```

### Panel 1B


```python
df = get_sc_data()
opref = 'figures/dt_v_randhex'
c_dict, order = get_priming_colors()
plot_read_len_kde(df, 'primer_type', c_dict, order, opref)
```


    
![png](figures/output_3_0.png)
    


### Panel 1C


```python
df = get_sc_data()
opref = 'figures/dt_v_randhex'
c_dict, order = get_priming_colors()
plot_reads_per_cell_nov(df, 'primer_type', c_dict, order, opref)
```


    
![png](figures/output_5_0.png)
    


### Panel 1D


```python
opref = 'figures/sample'
adata = get_sc_adata()
c_dict, order = get_sample_colors()
plot_short_long_det(adata.obs, c_dict, order, opref, 
                    xlim=7000, ylim=11000, how='gene')
```


    
![png](figures/output_7_0.png)
    


### Panel 1E


```python
opref = 'figures/sample'
adata = get_sc_adata()
c_dict, order = get_sample_colors()
plot_short_long_det(adata.obs, c_dict, order, opref, 
                    xlim=27000, ylim=175000, how='read')
```


    
![png](figures/output_9_0.png)
    


### Panel 1F


```python
opref = 'figures/cells_v_nuclei'
df = get_sc_data()

# remove genomic reads and limit to only MB cells and nuclei
df = df.loc[df['sample'] != 'MT_nuclei']
df = df[df.transcript_novelty!='Genomic']

c_dict, order = get_sample_colors(samples=['MB_cells', 'MB_nuclei'])
plot_read_len_kde(df, 'sample', c_dict, order, opref)
```


    
![png](figures/output_11_0.png)
    


### Panel 1G


```python
df = get_sc_data()

# limit to only MB cells and nuclei
df = df.loc[df['sample'] != 'MT_nuclei']
c_dict, order = get_sample_colors(samples=['MB_cells', 'MB_nuclei'])

plot_reads_per_cell_nov(df, 'sample', c_dict, order, opref)
```


    
![png](figures/output_13_0.png)
    


### Panel 1H


```python
opref = 'figures/bulk_vs_sc'
df = get_sc_bulk_data()

# myoblasts only, no genomic reads, no single-nucleus
df = df.loc[df['sample'].isin(['MB', 'MB_cells'])]
df = df.loc[(df.primer_type == 'Oligo dT')|(df.tech=='Bulk')]
df = df.loc[df.transcript_novelty != 'Genomic']

c_dict, order = get_2_tech_colors()
plot_read_len_kde(df, 'tech', c_dict, order, opref, common_norm=False)
```
    
![png](figures/output_15_1.png)
    


### Panel 1I


```python
opref = 'figures/bulk'
df = get_sc_bulk_data()
c_dict, order = get_talon_nov_colors()
d = df.loc[df.experiment == 'bulk', 'dataset'].unique().tolist()
plot_read_novelty(df, opref, c_dict, order, title='Bulk', \
                  datasets=d, ylim=7200000)
```

      transcript_novelty   counts
    0          Antisense   259540
    1            Genomic   448898
    2                ISM  1102134
    3         Intergenic    40109
    4              Known  6654958
    5                NIC   531944
    6                NNC   284313



    
![png](figures/output_17_1.png)
    


### Panel 1J


```python
opref = 'figures/sc'
df = get_sc_bulk_data()
c_dict, order = get_talon_nov_colors()
d = df.loc[df.experiment == 'sc', 'dataset'].unique().tolist()
plot_read_novelty(df, opref, c_dict, order, title='Single cell/nucleus', \
                  datasets=d)
```

      transcript_novelty  counts
    0          Antisense   70402
    1            Genomic  709033
    2                ISM  386720
    3         Intergenic   25641
    4              Known  317621
    5                NIC   20503
    6                NNC   15492



    
![png](figures/output_19_1.png)
    


### Panel 1K


```python
opref = 'figures/sc_filtered'
df = get_sc_data()
whitelist = get_sc_whitelist()
c_dict, order = get_talon_nov_colors()
plot_transcript_novelty(df, opref, c_dict, order, title='Filtered',
                        whitelist=whitelist, ylim=27000)
```

      transcript_novelty  counts
    0          Antisense     927
    1                ISM   16207
    2         Intergenic    1042
    3              Known   22279
    4                NIC     393
    5                NNC     134



    
![png](figures/output_21_1.png)
    

