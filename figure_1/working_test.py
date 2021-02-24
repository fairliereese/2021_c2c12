import sys
import os
import seaborn as sns

print(sys.path)
# p = os.path.dirname(os.getcwd())+'/scripts/'
p = os.path.dirname(os.getcwd())
print(p)
sys.path.append(p)

print(sys.path)

from scripts.utils import *
from scripts.plotting import *

sns.set_context('paper', font_scale=2)

# read in talon read_annot file from sc data
fname = '../processing/talon/sc_talon_read_annot.tsv'
fig_o = 'figures/dt_v_randhex/'
df = pd.read_csv(fname, sep='\t')

fname = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/c2c12_paper_2020/sc_pacbio/illumina_cell_metadata_full_bcs.tsv'
i_df = get_illumina_metadata(fname)
i_df = i_df[['raw_bc', 'primer_type']]
print(i_df.head())