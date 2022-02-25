import os
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

marker_genes_dict = ["Tmem119", "P2ry12", "Lcp1", "Cx3cr1", "Irf8",
                        "Sall1", "Mpeg1", "Csf1r", "Mrc1"]

# hammond lpc graph
adata1 = sc.read("hammond_lpc.h5ad")

vp1= sc.pl.stacked_violin(adata1, marker_genes_dict, groupby='cell_ident',swap_axes=True, dendrogram=False, return_fig = True, figsize = (6,6))
vp1.style(row_palette = 'husl', yticklabels=True, ylim = (0,5), x_padding=0, y_padding=0, linewidth= 1.5).legend(show=False).show()
vp1.savefig('hammond_lpc_vln_plot')