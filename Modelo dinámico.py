import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parámetros del modelo
mu_max = 0.8       # tasa máxima de crecimiento (1/día)
K_N = 0.5          # constante de saturación para N
K_P = 0.3          # constante de saturación para P
K_C = 0.4          # constante de saturación para C
d = 0.05           # tasa de pérdida de biomasa
Y_N = 0.2          # consumo de N por unidad de biomasa
Y_P = 0.1          # consumo de P por unidad de biomasa
Y_C = 0.3          # consumo de C por unidad de biomasa

# Modelo de ODE
def modelo(y, t):
    B, N, P, C = y
    # tasa de crecimiento Monod extendida a tres nutrientes
    mu = mu_max * (N/(K_N + N)) * (P/(K_P + P)) * (C/(K_C + C))
    dBdt = mu * B - d * B
    dNdt = -Y_N * dBdt
    dPdt = -Y_P * dBdt
    dCdt = -Y_C * dBdt
    return [dBdt, dNdt, dPdt, dCdt]

# Condiciones iniciales
B0 = 0.1   # biomasa inicial
N0 = 10.0  # concentración inicial de N
P0 = 5.0   # concentración inicial de P
C0 = 8.0   # concentración inicial de C
y0 = [B0, N0, P0, C0]

# Vector de tiempo
t = np.linspace(0, 20, 200)  # 0 a 20 días

# Integrar el sistema de ODEs
sol = odeint(modelo, y0, t)

# Graficar resultados
plt.figure(figsize=(10,6))
plt.plot(t, sol[:,0], label='Biomasa B')
plt.plot(t, sol[:,1], label='Nitrógeno N')
plt.plot(t, sol[:,2], label='Fósforo P')
plt.plot(t, sol[:,3], label='Carbono C')
plt.xlabel('Tiempo (días)')
plt.ylabel('Concentración / Biomasa')
plt.title('Dinámica microbiana limitada por tres nutrientes')
plt.legend()
plt.grid()
plt.show()
