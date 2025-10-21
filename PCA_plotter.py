import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")  #use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

def plot_pca(full_matrix, sample_groups, output_file="PCA_plot.jpeg", title="PCA of Mutational Signatures"):
    """
    full_matrix: pd.DataFrame, index=mutation types, columns=samples
    sample_groups: list of group labels corresponding to columns in full_matrix
    output_file: filename to save the PCA plot
    """
    X = full_matrix.T  #transposing(rows=samples, columns=features)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"])
    pca_df["Group"] = sample_groups
    pca_df["Sample"] = X.index

    plt.figure(figsize=(8,6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Group", s=100, palette="Set1")

    for i in range(pca_df.shape[0]):
        plt.text(pca_df.PC1[i]+0.2, pca_df.PC2[i]+0.2, pca_df.Sample[i], fontsize=8)

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file, dpi=300)
    print(f"PCA plot saved as {output_file}")
    plt.close()
