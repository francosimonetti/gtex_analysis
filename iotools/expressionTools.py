import numpy as np
from sklearn.decomposition import PCA
import pandas as pd

def normalize_expr(Y):
    if isinstance(Y, pd.DataFrame):
        Y_cent = (Y.values - np.mean(Y.values, axis = 1).reshape(-1, 1)) / np.std(Y.values, axis = 1).reshape(-1, 1)
        Y_cent = pd.DataFrame(Y_cent, index=Y.index, columns=Y.columns)
        Y_cent.index.name = Y.index.name
    else:
        Y_cent = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return Y_cent

def quant_normalize_expr(Y):
    from sklearn.preprocessing import normalize
    from sklearn.preprocessing import QuantileTransformer

    Y_quant = QuantileTransformer(output_distribution='normal').fit_transform(Y.T).T
    if isinstance(Y, pd.DataFrame):
        Y_quant_df = pd.DataFrame(Y_quant, index=Y.index, columns=Y.columns)
        Y_quant_df.index.name = Y.index.name
        return Y_quant_df
    else:
        return Y_quant


def knn_correction(expr, dosage, K, f=1):
    pca = PCA(n_components=min(expr.shape[0], expr.shape[1]) )
    # pca = PCA(n_components=30 )
    pca.fit(expr) # requires N x G
    expr_pca = pca.transform(expr)

    def gene_distance(a, b):
        return np.linalg.norm(a - b)

    nsample = expr.shape[0]
    distance_matrix = np.zeros((nsample, nsample))
    for i in range(nsample):
        for j in range(i+1, nsample):
            dist = gene_distance(expr_pca[i,:], expr_pca[j,:])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist

    kneighbor = K
    gx_knn = np.zeros_like(expr)
    if dosage is not None:
        gt_knn = np.zeros_like(dosage)
    neighbor_list = list()

    for i in range(nsample):
        neighbors = np.argsort(distance_matrix[i, :])[:kneighbor + 1][1:]
        gx_knn[i, :] = expr[i, :] - np.mean(expr[neighbors, :], axis = 0)
        # noisy_neighbors = np.random.choice(neighbors, size = int(2 * kneighbor / 3), replace = False)
        # noisy_neighbors = np.random.choice(neighbors, size = kneighbor, replace = True )
        noisy_neighbors = neighbors
        if dosage is not None:
            gt_knn[:, i] = dosage[:, i] - np.mean(dosage[:, noisy_neighbors], axis = 1)
        neighbor_list.append(neighbors)

    if dosage is not None:
        return gx_knn, gt_knn
    else:
        return gx_knn