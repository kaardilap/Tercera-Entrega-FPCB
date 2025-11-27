import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

df = pd.read_csv("fatty.csv")

# SOLO ácidos individuales (los sum no deben usarse en PERMANOVA real)
FA_cols = [
    'c18.2w6','c18.3w6','c18.3w3','c18.4w3','c18.5w3',
    'c20.4w6','c20.5w3','c22.6w3'
]

FA = df[FA_cols].copy()
group = df["Group"].astype("category").cat.codes

# Detectar % vs proporciones
percentage_cols = []
proportion_cols = []

for col in FA_cols:
    if FA[col].max() > 1.2:
        percentage_cols.append(col)
    else:
        proportion_cols.append(col)

print("\nColumnas (%) :", percentage_cols)
print("Columnas (proporción):", proportion_cols)

# Corregir negativos
FA[FA < 0] = 0

# Convertir % → proporciones
FA[percentage_cols] = FA[percentage_cols] / 100.0

# Transformación
FA_tr = np.arcsin(np.sqrt(FA.clip(0, 1) + 1e-12))

# Validación
print("\nVALIDACIÓN POST-TRANSFORMACIÓN")
print("Min =", FA_tr.min().min())
print("Max =", FA_tr.max().max())
print("NaN =", FA_tr.isna().sum().sum())
print("-------------------------------------")

if FA_tr.isna().sum().sum() > 0:
    FA_tr = FA_tr.fillna(FA_tr.mean())

# Distancia
dist_matrix = squareform(pdist(FA_tr, metric="braycurtis"))

# PERMANOVA
def permanova(distance_matrix, groups, permutations=999):
    groups = np.array(groups)
    n = len(groups)

    A = -0.5 * (distance_matrix ** 2)
    I = np.eye(n)
    H = I - np.ones((n, n)) / n
    G = H @ A @ H

    SST = np.trace(G)
    grand_mean = G.mean()

    SSB = 0
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        Gg = G[np.ix_(idx, idx)]
        SSB += len(idx) * (Gg.mean() - grand_mean) ** 2

    SSW = SST - SSB

    df_b = len(np.unique(groups)) - 1
    df_w = n - len(np.unique(groups))

    F = (SSB / df_b) / (SSW / df_w)

    perm_F = []
    for _ in tqdm(range(permutations)):
        perm = np.random.permutation(groups)
        SSBp = 0
        for g in np.unique(perm):
            idx = np.where(perm == g)[0]
            Gg = G[np.ix_(idx, idx)]
            SSBp += len(idx) * (Gg.mean() - grand_mean) ** 2
        SSWp = SST - SSBp
        perm_F.append((SSBp / df_b) / (SSWp / df_w))

    p_value = np.mean(np.array(perm_F) >= F)
    return F, p_value


F, p = permanova(dist_matrix, group, permutations=999)

print("\nPERMANOVA FINAL")
print("-----------------------------")
print(f"F = {F:.4f}")
print(f"p = {p:.4f}")
print("-----------------------------")
