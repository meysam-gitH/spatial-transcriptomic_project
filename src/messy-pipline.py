import scanpy as sc
import squidpy as sq

# Load data
adata = sc.datasets.visium_sge()

# Fix duplicate gene names
adata.var_names_make_unique()


# 2. QC
sc.pp.calculate_qc_metrics(adata, inplace=True)

# 3. NORMALIZATION
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 4. PCA + UMAP + CLUSTERING
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, flavor="igraph", directed=False, n_iterations=2)

adata.obs['cell_type'] = adata.obs['leiden'].map({
    '0': 'immune',
    '1': 'epithelial',
    '2': 'stromal',
})
#mapping correction
adata.obs['cell_type'] = adata.obs['leiden'].replace({
    '0': 'immune_myeloid',
    '1': 'immune_lymphoid',
    '2': 'epithelial',
})
sc.pl.umap(adata, color='cell_type')
import squidpy as sq
sq.pl.spatial_scatter(adata, color='cell_type')



adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)
adata.obs['cell_type'] = adata.obs['cell_type'].replace({
    '3': 'unknown',
    '4': 'unknown',
    '5': 'unknown',
    '6': 'unknown',
    '7': 'unknown',
    '8': 'unknown',
    '9': 'unknown',
    '10': 'unknown',
    '11': 'unknown',
    '12': 'unknown',
    '13': 'unknown',
    '14': 'unknown',
    '15': 'unknown',
    '16': 'unknown',
    '17': 'unknown',
})




# 5. NOW THIS WORKS
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10)
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=5,
    standard_scale='var'
)


#Identify clusters using marker genes
sc.pl.umap(
    adata,
    color=['SLC17A7', 'OLIG1', 'GFAP'],
    cmap='viridis'
)

#Step — Extract top genes per cluster
sc.get.rank_genes_groups_df(adata, group='0').head(10)
for i in ['0', '1', '2']:
    print(f"\nCluster {i}")
    print(sc.get.rank_genes_groups_df(adata, group=i).head(5))
