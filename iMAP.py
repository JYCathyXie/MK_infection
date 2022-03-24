import scanpy as sc
import imap
import pandas as pd
import numpy as np
import os

os.chdir('/media/user/sdd')
adata = sc.read_loom('/media/user/sdd/Merge_all_round.loom',sparse=False)
adata.obs['nCount_RNA'] = adata.obs['nCount_RNA'].astype("int64")
adata = imap.stage1.data_preprocess(adata, 'batch', n_top_genes=2000, min_genes=0, min_cells=0)
adata.obs['n_counts'] = adata.obs['n_counts'].astype("int64")
adata = imap.stage1.data_preprocess(adata, 'batch', n_top_genes=2000, min_genes=0, min_cells=0)
EC, ec_data = imap.stage1.iMAP_fast(adata, key="batch", n_epochs=100)
output_results = imap.stage2.integrate_data(adata, ec_data, n_epochs=100)

import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import umap
def data2umap(data, n_pca=0):
    if n_pca > 0:
        pca = PCA(n_components=n_pca)
        embedding = pca.fit_transform(data)
    else:
        embedding = data
    embedding_ = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric='cosine',
        n_components = 2,
        learning_rate = 1.0,
        spread = 1.0,
        set_op_mix_ratio = 1.0,
        local_connectivity = 1,
        repulsion_strength = 1,
        negative_sample_rate = 5,
        angular_rp_forest = False,
        verbose = False
    ).fit_transform(embedding)
    return embedding_
  
embedding_ = data2umap(output_results, n_pca=30)
test = pd.DataFrame(embedding_, columns=['UMAP_1', 'UMAP_2'])
batch_info = np.array(adata.obs['batch'])
celltype_info = np.array([adata.obs['cell_type'][item] for item in adata.obs_names])
test['Label1'] = batch_info
test['Label2'] = celltype_info
test['Label3'] = adata.obs_names
test.to_csv('top5000_nep100_npca30_umap.csv',index= True ,header=True)


