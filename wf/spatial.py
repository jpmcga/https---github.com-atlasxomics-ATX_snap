import anndata
import squidpy as sq


def add_spatial(
    adata: anndata.AnnData, x_key: str = "xcor", y_key: str = "ycor"
) -> anndata.AnnData:
    """Add move x and y coordinates from .obs to .obsm["spatial"] for squidpy.
    """
    adata.obsm["spatial"] = adata.obs[[y_key, x_key]].values

    return adata


def plot_spatial(
    adata: anndata.AnnData, run_id: str, metric: str, dot_size: int = 50
):
    """Extracts cells for run_id, creates squidpy spatial_scatter plot with
    color as metric from .obs.
    """

    adata = adata[adata.obs["run_id"] == run_id]
    plot = sq.pl.spatial_scatter(
        adata, color=metric, size=dot_size, shape=None, library_id=run_id
    )

    return plot
