import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Cargar datos
df = pd.read_csv("fatty.csv")

# Convertir categorías a valores numéricos
df["N_numeric"] = df["Nutrient"].replace({"Replete": 10, "Limited": 1})
df["L_numeric"] = df["Light"]

# Seleccionar la variable dependiente
df["G_obs"] = df["sumPUFA"]

# Eliminar filas con NaN o inf
df_clean = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["N_numeric", "L_numeric", "G_obs"])

N = df_clean["N_numeric"].values
L = df_clean["L_numeric"].values
G_obs = df_clean["G_obs"].values

# Definir modelo Monod extendido
def modelo_fit(X, Gmax, KN, KL):
    N, L = X
    return Gmax * (N / (KN + N)) * (L / (KL + L))

# Valores iniciales y límites para evitar KN o KL = 0
p0 = [30, 1, 1]
bounds = ([0, 0.01, 0.01], [np.inf, np.inf, np.inf])

# Ajuste
popt, pcov = curve_fit(modelo_fit, (N, L), G_obs, p0=p0, bounds=bounds, maxfev=5000)

print("Parámetros ajustados:")
print(f"Gmax = {popt[0]:.2f}, KN = {popt[1]:.2f}, KL = {popt[2]:.2f}")

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
plt.title("Ajuste del modelo biológico")
plt.legend()
plt.show()
