# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 19:29:30 2026

@author: Leandro
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.integrate import solve_ivp
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


# =========================
# Modelo Lotka-Volterra
# =========================

def lotka_volterra(t, z, alpha, beta, delta, gamma):
    x, y = z

    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y

    return [dxdt, dydt]


def solve_model(alpha, beta, delta, gamma):
    sol = solve_ivp(
        lotka_volterra,
        t_span=(0, T),
        y0=[x0, y0],
        t_eval=t,
        args=(alpha, beta, delta, gamma),
        rtol=1e-8,
        atol=1e-10
    )

    return sol.y[0], sol.y[1]


# =========================
# Parámetros iniciales y solver
# =========================

LABEL_FONTSIZE = 14
TICK_FONTSIZE = 12

alpha0 = 1.0
beta0  = 0.1
delta0 = 0.075
gamma0 = 1.5

x0 = 10.0
y0 = 5.0

T = 50
t = np.linspace(0, T, 2000)

x, y = solve_model(alpha0, beta0, delta0, gamma0)


# =========================
# Figura y layout
# =========================

fig = plt.figure(figsize=(10, 6.5))

gs = GridSpec(
    2, 2,
    height_ratios=[3, 1.7],
    width_ratios=[1, 1],
    hspace=0.55,
    wspace=0.35
)

ax_time = fig.add_subplot(gs[0, :])
ax_phase = fig.add_subplot(gs[1, 1])

slider_gs = GridSpecFromSubplotSpec(
    5, 1,
    subplot_spec=gs[1, 0],
    hspace=0.9
)

ax_alpha = fig.add_subplot(slider_gs[0])
ax_beta  = fig.add_subplot(slider_gs[1])
ax_delta = fig.add_subplot(slider_gs[2])
ax_gamma = fig.add_subplot(slider_gs[3])
ax_reset = fig.add_subplot(slider_gs[4])


# =========================
# Gráficos
# =========================

line_x, = ax_time.plot(t, x, lw=2, label="Presa")
line_y, = ax_time.plot(t, y, lw=2, label="Depredador")

ax_time.set_xlabel("Tiempo", fontsize=LABEL_FONTSIZE)
ax_time.set_ylabel("Población", fontsize=LABEL_FONTSIZE)
ax_time.tick_params(axis="both", labelsize=TICK_FONTSIZE)

ax_time.legend(frameon=False, fontsize=12)
ax_time.grid(alpha=0.3)



line_phase, = ax_phase.plot(x, y, lw=2)

ax_phase.set_xlabel("Presa", fontsize=LABEL_FONTSIZE)
ax_phase.set_ylabel("Depredador", fontsize=LABEL_FONTSIZE)
ax_phase.tick_params(axis="both", labelsize=TICK_FONTSIZE)

ax_phase.grid(alpha=0.3)


# =========================
# Sliders
# =========================

def add_slider_title(ax, text):
    ax.text(
        0.0, 1.35,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        clip_on=False
    )


slider_alpha = Slider(
    ax_alpha,
    "",
    0.0, 3.0,
    valinit=alpha0,
    valstep=0.01
)

add_slider_title(ax_alpha, r"$\alpha$ : crecimiento de la presa")


slider_beta = Slider(
    ax_beta,
    "",
    0.0, 1.0,
    valinit=beta0,
    valstep=0.005
)

add_slider_title(ax_beta, r"$\beta$ : depredación")


slider_delta = Slider(
    ax_delta,
    "",
    0.0, 1.0,
    valinit=delta0,
    valstep=0.005
)

add_slider_title(ax_delta, r"$\delta$ : conversión presa $\rightarrow$ depredador")


slider_gamma = Slider(
    ax_gamma,
    "",
    0.0, 3.0,
    valinit=gamma0,
    valstep=0.01
)

add_slider_title(ax_gamma, r"$\gamma$ : muerte del depredador")


# =========================
# Actualización 
# =========================

def update(val):
    alpha = slider_alpha.val
    beta  = slider_beta.val
    delta = slider_delta.val
    gamma = slider_gamma.val

    x, y = solve_model(alpha, beta, delta, gamma)

    line_x.set_ydata(x)
    line_y.set_ydata(y)

    line_phase.set_xdata(x)
    line_phase.set_ydata(y)

    ymax_time = max(np.max(x), np.max(y))
    ax_time.set_ylim(0, 1.1 * ymax_time)

    ax_phase.set_xlim(0, 1.1 * np.max(x))
    ax_phase.set_ylim(0, 1.1 * np.max(y))

    fig.canvas.draw_idle()


slider_alpha.on_changed(update)
slider_beta.on_changed(update)
slider_delta.on_changed(update)
slider_gamma.on_changed(update)



button_reset = Button(ax_reset, "Reset")

def reset(event):
    slider_alpha.reset()
    slider_beta.reset()
    slider_delta.reset()
    slider_gamma.reset()

button_reset.on_clicked(reset)


# =========================
# Ajuste inicial
# =========================

update(None)

plt.show()
