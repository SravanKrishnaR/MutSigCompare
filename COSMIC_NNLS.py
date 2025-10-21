import pandas as pd
import numpy as np
from scipy.optimize import nnls
from numpy.linalg import norm
from DATA_PLOTTER import Bar_plotter, heatmap_plotter

def cosmic_nnlss(cosmic_data,matrix_rel):
    cosmic_data=pd.read_csv(cosmic_data, sep="\t",index_col=0)

    X = cosmic_data.values
    exposures = {}

    for sample in matrix_rel.columns:
        y = matrix_rel[sample].values
        coeffs, _ = nnls(X, y)
        exposures[sample] = coeffs

    exposures_df = pd.DataFrame(exposures, index=cosmic_data.columns)
    Bar_plotter(exposures_df)
    heatmap_plotter(exposures_df)

if __name__ == "__main__":
    cosmic_data = pd.read_csv(
        "DATA/COSMIC_v3.4_SBS_GRCh37.txt", sep="\t", index_col=0
    )

    dummy_matrix = pd.DataFrame(
        np.random.randint(50, 200, size=(cosmic_data.shape[0], 2)),
        index=cosmic_data.index,
        columns=["Sample1", "Sample2"]
    )

    cosmic_nnlss("DATA/COSMIC_v3.4_SBS_GRCh37.txt", dummy_matrix)
