import sys, os, warnings
warnings.filterwarnings("ignore")
try:
    import scanpy as sc
    import numpy as np
    from sklearn.neighbors import kneighbors_graph
    files = sys.argv[1:]
    out = []
    for f in files[:3]:
        try:
            ad = sc.read_h5ad(f)
            has_spatial = "spatial" in ad.obsm and ad.obsm["spatial"].shape[1] >= 2
            if has_spatial:
                coords = ad.obsm["spatial"][:, :2]
                coord_type = "spatial"
            elif "X_umap" in ad.obsm:
                coords = ad.obsm["X_umap"][:, :2]
                coord_type = "UMAP"
            elif "X_pca" in ad.obsm:
                coords = ad.obsm["X_pca"][:, :2]
                coord_type = "PCA"
            else:
                out.append(f"{os.path.basename(f)}|{ad.shape[0]}|{ad.shape[1]}|NA|no_coords")
                continue
            
            y = None
            if "NDUFB7" in ad.var_names:
                y = ad[:, "NDUFB7"].X.toarray().flatten() if hasattr(ad[:, "NDUFB7"].X, "toarray") else ad[:, "NDUFB7"].X.flatten()
            
            if y is not None and len(y) > 20:
                knn = kneighbors_graph(coords, n_neighbors = min(10, len(y)-1), mode = "connectivity")
                W = knn.astype(float)
                row_sums = W.sum(axis=1).A1
                W = W.tocsr()
                for i in range(W.shape[0]):
                    if row_sums[i] > 0:
                        W.data[W.indptr[i]:W.indptr[i+1]] /= row_sums[i]
                z = y - np.mean(y)
                if np.sum(z**2) > 0:
                    Iz = W.dot(z)
                    I = np.sum(z * Iz) / np.sum(z**2)
                    out.append(f"{os.path.basename(f)}|{ad.shape[0]}|{ad.shape[1]}|{I:.4f}|{coord_type}")
                else:
                    out.append(f"{os.path.basename(f)}|{ad.shape[0]}|{ad.shape[1]}|NA|zero_variance")
            else:
                out.append(f"{os.path.basename(f)}|{ad.shape[0]}|{ad.shape[1]}|NA|NDUFB7_missing")
        except Exception as e:
            out.append(f"{os.path.basename(f)}|NA|NA|ERROR|{str(e)}")
    print("\n".join(out))
except ImportError as e:
    print(f"SCANPY_MISSING|{str(e)}")

