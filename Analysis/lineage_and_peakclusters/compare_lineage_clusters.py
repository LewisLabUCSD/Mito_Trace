import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join
import click


def load_lineage(indir):
    """Load the cell indices for each cluster and each sample"""
    return


def load_clusters(indirs):
    """Load Flt3l and control cluster assignments"""
    for ind in indirs:
        ### Load clusters
        clusters_df = pd.read_csv(ind)
    return


def cell_labels_extract_clusters():
    # Combine the samples, create new cluster IDs (will merge another time)
    return


def cell_labels_extract_lineage():
    # Extract the cell names in each cluster, seeing which sample it comes
    # from
    return


def remap(clusters, lineages, cluster_labels,  lineage_labels):
    return


def count_overlap_lineage_cluster(df, out_f):
    """

    df: Columns: ["Cell", "Cluster", "Donor", "Lineage", "Sample]
    """
    df["ind"] = df["Lineage"].astype(str) + "_" + df["Donor"].astype(str)
    df["cluster"] = df["Cluster"].astype(str) + "_" + df["Sample"].astype(str)
    sns.countplot(y="ind", hue="cluster", data=df)
    plt.savefig(out_f, dpi=300)
    return


def plot_overlap_lineage_cluster():
    return


@click.command()
@click.argument("aggregate_in", type=click.Path(exists=True))
@click.argument("lineage_in", type=click.Path(exists=True))
@click.argument("outdir", type=click.Path(exists=False))
def main(aggregate_in, lineage_in, outdir):
    #aggregate_out
    return


if __name__ == "__main__":
    main()

