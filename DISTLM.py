import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

# ================================
# 1) Cargar datos
# ================================
df = pd.read_csv("fatty.csv")
FA_cols = [
    'sumSAFA', 'sumMUFA', 'sumPUFA',
    'c18.2w6', 'c18.3w6', 'c18.3w3', 'c18.4w3', 'c18.5w3',
    'c20.4w6', 'c20.5w3', 'c22.6w3'
]
predictor_cols = ["Group", "Nutrient", "Light", "Temp", "Salinity"]

# ================================
# 2) Separar FA y predictores
# ================================
FA = df[FA_cols].copy()
pred = df[predictor_cols].copy()

# Detectar columnas en % y dividir entre 100
for col in FA_cols:
    if FA[col].max() > 1.2:
        FA[col] = FA[col] / 100.0

# Clip a [0,1] y arcsin–sqrt
FA_tr = np.arcsin(np.sqrt(FA.clip(0, 1)))

# ================================
# 3) Eliminar filas con NAs en predictores o FA
# ================================
combined = pd.concat([FA_tr, pred], axis=1)
combined.dropna(inplace=True)
FA_tr = combined[FA_cols].values
pred = combined[predictor_cols]

print(f"Filas finales después de limpieza: {FA_tr.shape[0]}")

# ================================
# 4) Distancia Bray–Curtis
# ================================
dist_matrix = squareform(pdist(FA_tr, metric="braycurtis"))


# ================================
# 5) Funciones DISTLM marginal
# ================================
def compute_ss_reg(G, X):
    X = np.array(X, dtype=float)
    Hx = X @ np.linalg.pinv(X.T @ X) @ X.T
    G_reg = Hx @ G @ Hx
    SS_reg = np.trace(G_reg)
    df_reg = np.linalg.matrix_rank(X)
    return SS_reg, df_reg


def distlm_marginal(G, X, permutations=1000):
    SS_reg, df_reg = compute_ss_reg(G, X)
    SST = np.trace(G)
    F_obs = (SS_reg / df_reg) / ((SST - SS_reg) / (G.shape[0] - df_reg - 1))

    perm_F = []
    for _ in tqdm(range(permutations), desc="Permutations"):
        perm_idx = np.random.permutation(G.shape[0])
        G_perm = G[perm_idx][:, perm_idx]
        SS_perm, _ = compute_ss_reg(G_perm, X)
        F_perm = (SS_perm / df_reg) / ((SST - SS_perm) / (G.shape[0] - df_reg - 1))
        perm_F.append(F_perm)
    pval = np.mean(np.array(perm_F) >= F_obs)
    prop_var = SS_reg / SST
    return F_obs, pval, prop_var


# ================================
# 6) Ejecutar DISTLM marginal para cada predictor
# ================================
G = dist_matrix
for col in predictor_cols:
    X = pd.get_dummies(pred[col], drop_first=True)
    F, p, prop_var = distlm_marginal(G, X.values, permutations=1000)
    print(f"{col:10} | F = {F:.4f} | p = {p:.4f} | prop.var = {prop_var:.4f}")
