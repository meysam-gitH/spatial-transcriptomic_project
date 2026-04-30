import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt


def load_data():
    adata = sc.datasets.visium_sge()
    adata.var_names_make_unique()
    return adata


def preprocess(adata):
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def run_analysis(adata):
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)


def annotate(adata):
    adata.obs['cell_type'] = adata.obs['leiden'].replace({
        '0': 'immune_myeloid',
        '1': 'immune_lymphoid',
        '2': 'epithelial',
    })


def plot_results(adata):
    sc.pl.umap(adata, color='cell_type', save="_umap.png")
    sq.pl.spatial_scatter(adata, color='cell_type', save="_spatial.png")


def run_pipeline():
    adata = load_data()
    preprocess(adata)
    run_analysis(adata)
    annotate(adata)
    plot_results(adata)


if __name__ == "__main__":
    run_pipeline()
