import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Cargar datos
df = pd.read_csv("fatty.csv")  # reemplaza con tu archivo

# Convertir categorías a valores numéricos
df["N_numeric"] = df["Nutrient"].replace({"Replete": 10, "Limited": 1})
df["L_numeric"] = df["Light"]

# Seleccionar la variable dependiente
df["G_obs"] = df["sumPUFA"]

# Limpiar datos: eliminar NaN o inf
df_clean = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["N_numeric", "L_numeric", "G_obs"])

N = df_clean["N_numeric"].values
L = df_clean["L_numeric"].values
G_obs = df_clean["G_obs"].values

# Definir modelo Monod extendido con término basal
def modelo_fit(X, Gmax, KN, KL, Gmin):
    N, L = X
    return Gmin + (Gmax - Gmin) * (N / (KN + N)) * (L / (KL + L))

# Valores iniciales y límites (biológicamente razonables)
p0 = [30, 1, 1, 0.5]  # [Gmax, KN, KL, Gmin]
bounds = ([0, 0.1, 0.1, 0], [100, 5, 10, 10])

# Ajuste
popt, pcov = curve_fit(modelo_fit, (N, L), G_obs, p0=p0, bounds=bounds, maxfev=5000)

print("Parámetros ajustados:")
print(f"Gmax = {popt[0]:.2f}, KN = {popt[1]:.2f}, KL = {popt[2]:.2f}, Gmin = {popt[3]:.2f}")

# Graficar ajuste 3D
N_lin = np.linspace(min(N), max(N), 30)
L_lin = np.linspace(min(L), max(L), 30)
N_grid, L_grid = np.meshgrid(N_lin, L_lin)
G_fit = modelo_fit((N_grid, L_grid), *popt)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(N, L, G_obs, color='r', label='Datos')
ax.plot_surface(N_grid, L_grid, G_fit, alpha=0.5, cmap='viridis')
ax.set_xlabel("Nutrient")
ax.set_ylabel("Light")
ax.set_zlabel("sumPUFA")
plt.title("Ajuste del modelo biológico (Monod extendido)")
plt.legend()
plt.show()
