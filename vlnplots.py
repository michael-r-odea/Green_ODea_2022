import os
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import rcParams

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# hammond graph (figure 3C)
adata1 = sc.read("hammond_microglia.h5ad")

marker_genes_dict = ["Tmem119", "P2ry12", "Lcp1", "Cx3cr1", "Irf8",
                        "Sall1", "Mpeg1", "Csf1r", "Mrc1"]

vp = sc.pl.stacked_violin(adata1, marker_genes_dict, groupby='cell_ident',swap_axes=True, dendrogram=False, return_fig = True, figsize = (6,6))
vp.style(row_palette = 'husl', yticklabels=True, ylim = (0,5), x_padding=0, y_padding=0, linewidth= 1.5).legend(show=False).show()
vp.savefig('hammond_vln_plot')

# zhong graph (Figure 3G)
adata2 = sc.read("zhong_microglia.h5ad")

human_marker_genes_dict = ["TMEM119", "P2RY12", "LCP1", "CX3CR1", "IRF8",
                        "SALL1", "MPEG1", "CSF1R", "MRC1"]

vp2 = sc.pl.stacked_violin(adata2, human_marker_genes_dict, groupby='cell_ident',swap_axes=True, dendrogram=False, return_fig = True, figsize = (6,6))
vp2.style(row_palette = 'husl', yticklabels=True, ylim = (0,6), x_padding=0, y_padding=0, linewidth= 1.5).legend(show=False).show()
vp2.savefig('zhong_vln_plot')

# hammond lpc graph (figure s6E)
adata3 = sc.read("hammond_lpc.h5ad")

vp3= sc.pl.stacked_violin(adata3, marker_genes_dict, groupby='cell_ident',swap_axes=True, dendrogram=False, return_fig = True, figsize = (6,6))
vp3.style(row_palette = 'husl', yticklabels=True, ylim = (0,5), x_padding=0, y_padding=0, linewidth= 1.5).legend(show=False).show()
vp3.savefig('hammond_lpc_vln_plot')