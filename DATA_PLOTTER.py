import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def Bar_plotter(exposures_df, output_file="BAR_PLOT.jpeg", title="BAR PLOT of Mutational Signatures"):
    exposures_clean = exposures_df.replace([np.inf, -np.inf], np.nan).fillna(0) #Cleaning the data
    exposures_clean = exposures_clean.loc[(exposures_clean != 0).any(axis=1)] #dropping zero value signatures

    colors = sns.color_palette("Set2", n_colors=len(exposures_clean))
    fig, ax = plt.subplots(figsize=(9, 6))
    exposures_clean.T.plot(kind="bar", stacked=True, color=colors, edgecolor="black", linewidth=0.7,ax=ax)

    ax.set_ylabel("Signature Contribution", fontsize=12)
    ax.set_xlabel("Samples", fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_title(title, fontsize=14, weight="bold")

    legend = ax.legend(title="Signature", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, fontsize=9)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Bar plot saved as {output_file}")
    plt.close(fig)

def heatmap_plotter(exposures_df, output_file="HEAT_MAP.jpeg", title="HEAT MAP of Mutational Signatures"):
    exposures_clean = exposures_df.replace([np.inf, -np.inf], np.nan).fillna(0) #Cleaning the data
    exposures_clean = exposures_clean.loc[(exposures_clean != 0).any(axis=1)] #dropping zero value signatures

    plt.figure(figsize=(7,6))
    sns.heatmap(exposures_clean, cmap="mako", annot=False, cbar_kws={"label": "Exposure"}, linewidths=0.5, linecolor="white")

    plt.title(title, fontsize=14, weight="bold")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Heatmap saved as {output_file}")
    plt.close()

    if exposures_clean.shape[0] > 1 and exposures_clean.shape[1] > 1:
        cluster_grid = sns.clustermap(exposures_clean, cmap="mako", metric="cosine", method="average", standard_scale=0, figsize=(7,6), cbar_kws={"label": "Exposure"})
        cluster_grid.fig.suptitle("Clustered Heatmap of Signature Exposures", fontsize=14, weight="bold", y=1.02)
        cluster_grid.savefig("Clustered_" + output_file, dpi=300)
        print(f"Clustered heatmap saved as Clustered_{output_file}")
        plt.close()
    else:
        print(" Not enough non zero rows/columns for clustering.")
