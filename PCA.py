import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# ================================
# 1. Cargar dataset
# ================================
df = pd.read_csv("fatty.csv")

# ================================
# 2. Filtrar grupos incluidos (S2 Table)
# ================================
groups_included = [
    "Chlorophyta",
    "Cryptophyta",
    "Cyanobacteria",
    "Ochrophyta",     # diatoms
    "Dinophyta",
    "Haptophyta"
]

df = df[df["Group"].isin(groups_included)].reset_index(drop=True)

# ================================
# 3. Seleccionar columnas FA
# ================================
FA_cols = [
    'sumSAFA','sumMUFA','sumPUFA',
    'c18.2w6','c18.3w6','c18.3w3','c18.4w3','c18.5w3',
    'c20.4w6','c20.5w3','c22.6w3'
]

FA = df[FA_cols].copy()

# ================================
# 4. Convertir % → proporción
# ================================
for col in FA_cols:
    if FA[col].max() > 1.2:
        FA[col] = FA[col] / 100.0

# ================================
# 5. Transformación arcsin–sqrt
# ================================
FA_tr = np.arcsin(np.sqrt(FA.clip(0, 1)))

# ================================
# 6. Escalamiento (PCA requiere datos centrados)
# ================================
scaler = StandardScaler()
FA_scaled = scaler.fit_transform(FA_tr)

# ================================
# 7. PCA
# ================================
pca = PCA(n_components=3)
pc = pca.fit_transform(FA_scaled)

explained = pca.explained_variance_ratio_

print("\nVARIANZA EXPLICADA:")
print(f"PC1: {explained[0]*100:.2f}%")
print(f"PC2: {explained[1]*100:.2f}%")
print(f"PC3: {explained[2]*100:.2f}%")
print(f"PC1 + PC2 = {(explained[0]+explained[1])*100:.2f}%")

# ================================
# 8. Gráfica PCA
# ================================
plt.figure(figsize=(9,7))

groups = df["Group"].astype(str)
unique_groups = groups.unique()
colors = plt.cm.tab10(range(len(unique_groups)))

for g, c in zip(unique_groups, colors):
    idx = groups == g
    plt.scatter(pc[idx,0], pc[idx,1], label=g, alpha=0.7, color=c)

plt.xlabel(f"PC1 ({explained[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({explained[1]*100:.1f}%)")
plt.title("PCA de Ácidos Grasos (réplica del artículo)")
plt.legend()
plt.tight_layout()
plt.show()

# ================================
# 9. Cargas (contribución de cada FA)
# ================================
loadings = pd.DataFrame(
    pca.components_.T,
    columns=["PC1","PC2","PC3"],
    index=FA_cols
)

print("\nCARGAS (loadings) DE LOS ÁCIDOS GRASOS:")
print(loadings.round(3))
