import numpy as np
import matplotlib.pyplot as plt

filename = "simulation_L_0.09.txt"

# leer archivo y separar bloques por tiempo
times = []
positions = []

with open(filename, "r") as f:
    lines = f.readlines()

block = []
for line in lines:
    line = line.strip()
    if line.startswith("t"):  # nuevo bloque
        if block:
            positions.append(np.array(block))
            block = []
        times.append(float(line[1:]))
    elif line:
        values = [float(x) for x in line.split()]
        block.append(values[:2])  # solo x, y
if block:
    positions.append(np.array(block))

times = np.array(times)
N_particles = positions[0].shape[0]

# posiciones iniciales
x0 = positions[0][:,0]
y0 = positions[0][:,1]

# calcular MSD y desviación estándar
MSD = []
MSD_std = []
for pos in positions:
    dx = pos[:,0] - x0
    dy = pos[:,1] - y0
    msd_per_particle = dx**2 + dy**2
    MSD.append(np.mean(msd_per_particle))
    MSD_std.append(np.std(msd_per_particle))

MSD = np.array(MSD)
MSD_std = np.array(MSD_std)

# restar tiempo inicial estacionario t0 = 60 s
t0_stationary = 60
times_stationary = times - t0_stationary

# seleccionar solo tiempos entre 0 y 400 s desde estacionario (realmente t=60 a t=460)
mask_plot = (times_stationary >= 0) & (times_stationary <= 400)
times_plot = times_stationary[mask_plot]
MSD_plot = MSD[mask_plot]
MSD_std_plot = MSD_std[mask_plot]

# tomar un punto cada 40 s
sample_interval = 40
sample_mask = np.arange(0, len(times_plot), sample_interval)
times_sampled = times_plot[sample_mask]
MSD_sampled = MSD_plot[sample_mask]
MSD_std_sampled = MSD_std_plot[sample_mask]

# ajuste lineal solo hasta 200 s desde estacionario (realmente t=260)
mask_fit = times_sampled <= 200
times_fit = times_sampled[mask_fit]
MSD_fit = MSD_sampled[mask_fit]

coeffs = np.polyfit(times_fit, MSD_fit, 1)
slope = coeffs[0]
intercept = coeffs[1]
D = slope / 4  # 2D

# gráfico con puntos muestreados, línea uniendo puntos y ajuste lineal
plt.figure(figsize=(8,6))
plt.errorbar(times_sampled, MSD_sampled, yerr=MSD_std_sampled,
             fmt='o', markersize=6, color='tab:blue',
             ecolor='gray', elinewidth=1, capsize=3, label='Datos ± std')
plt.plot(times_sampled, MSD_sampled, 'b-', linewidth=1)
plt.plot(times_fit, np.poly1d(coeffs)(times_fit), 'r--', linewidth=2,
         label=f'Ajuste lineal, D={D:.2e} m²/s')

plt.xlabel("Tiempo desde estado estacionario (s)", fontsize=12)
plt.ylabel("<r²> (m²)", fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
