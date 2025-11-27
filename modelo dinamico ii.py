import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ipywidgets import interact, FloatSlider


# Definir el modelo dinámico
def sistema_biologico(y, t, rB, rN, rP, KN, KP, KC):
    B, N, P, C = y
    dBdt = rB * B * (N / (KN + N)) * (P / (KP + P)) * (C / (KC + C))
    dNdt = -rN * dBdt
    dPdt = -rP * dBdt
    dCdt = -dBdt
    return [dBdt, dNdt, dPdt, dCdt]


# Función para simular y graficar
def simular_dinamica(rB=0.5, rN=0.3, rP=0.2, KN=2.0, KP=1.0, KC=5.0):
    y0 = [1.0, 10.0, 5.0, 20.0]  # condiciones iniciales: B, N, P, C
    t = np.linspace(0, 50, 500)
    sol = odeint(sistema_biologico, y0, t, args=(rB, rN, rP, KN, KP, KC))

    plt.figure(figsize=(10, 6))
    plt.plot(t, sol[:, 0], label='Biomasa (B)', linewidth=2)
    plt.plot(t, sol[:, 1], label='Nitrógeno (N)', linewidth=2)
    plt.plot(t, sol[:, 2], label='Fósforo (P)', linewidth=2)
    plt.plot(t, sol[:, 3], label='Carbono (C)', linewidth=2)
    plt.xlabel('Tiempo')
    plt.ylabel('Concentración')
    plt.title('Dinámica del sistema biológico con tres nutrientes')
    plt.legend()
    plt.grid(True)
    plt.show()


# Crear sliders interactivos
interact(
    simular_dinamica,
    rB=FloatSlider(min=0.1, max=1.0, step=0.05, value=0.5, description='rB'),
    rN=FloatSlider(min=0.1, max=1.0, step=0.05, value=0.3, description='rN'),
    rP=FloatSlider(min=0.1, max=1.0, step=0.05, value=0.2, description='rP'),
    KN=FloatSlider(min=0.5, max=10.0, step=0.5, value=2.0, description='KN'),
    KP=FloatSlider(min=0.5, max=10.0, step=0.5, value=1.0, description='KP'),
    KC=FloatSlider(min=1.0, max=20.0, step=1.0, value=5.0, description='KC')
)
