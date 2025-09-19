import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

m = 1.0
dt_bin = 1.0
t_estacionario = 60.0

def perimetro_compartment_0(L):
    return 0.09 * 4 - L

def perimetro_compartment_1(L):
    return 0.09 * 2 + L

def total_perimetro(L):
    return perimetro_compartment_0(L) + perimetro_compartment_1(L)

def process_simulation(filename, L):
    df = pd.read_csv(
        filename,
        names=["time", "compartment", "normal_velocity", "radius"],
        dtype={"time": float, "compartment": int, "normal_velocity": float, "radius": float},
        skiprows=1
    )

    # bin de tiempo
    df["bin"] = (df["time"] // dt_bin).astype(int)
    df["impulso"] = 2.0 * m * df["normal_velocity"]

    # impulsos por bin y compartimento
    presiones = df.groupby(["bin", "compartment"])["impulso"].sum().reset_index()
    presiones["tiempo"] = presiones["bin"] * dt_bin

    # perímetro según compartimento
    def calc_perimetro(comp):
        return perimetro_compartment_0(L) if comp == 0 else perimetro_compartment_1(L)

    presiones["perimetro"] = presiones["compartment"].apply(calc_perimetro)
    presiones["P"] = presiones["impulso"] / (presiones["perimetro"] * dt_bin)

    # presión del sistema completo por bin (suma de ambos compartimentos)
    system = presiones.groupby("bin")["impulso"].sum().reset_index()
    system["tiempo"] = system["bin"] * dt_bin
    system["P_system"] = system["impulso"] / (total_perimetro(L) * dt_bin)

    # estado estacionario: usar >60s si hay, si no, último 20%
    t_max = system["tiempo"].max()
    t_cut = t_estacionario if t_max > t_estacionario else t_max * 0.8
    estacionario = system[system["tiempo"] > t_cut]
    P_system_prom = estacionario["P_system"].mean() if len(estacionario) > 0 else np.nan

    return presiones, system, P_system_prom

def create_plots(presiones, L):
    """Presión vs tiempo por compartimento."""
    plt.figure(figsize=(10,6))
    for comp, color, label in [(0, "blue", "Área Fija"), (1, "red", "Área Variable")]:
        subset = presiones[presiones["compartment"] == comp]
        plt.plot(subset["tiempo"], subset["P"], color=color, label=label)
    plt.xlabel("Tiempo (s)", fontsize=16)
    plt.ylabel("Presión (kg/s²)", fontsize=16)
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.show()

def plot_pressure_vs_inv_area(resultados):
    """P promedio del sistema vs 1/A total."""
    Ls = np.array([r[0] for r in resultados])
    Ps = np.array([r[1] for r in resultados])

    Areas = 0.09 * 0.09 + 0.09 * Ls
    invA = 1.0 / Areas

    mask = ~np.isnan(Ps)
    Ls, Ps, invA = Ls[mask], Ps[mask], invA[mask]

    coeffs = np.polyfit(invA, Ps, 1)
    m_fit, b_fit = coeffs
    fit_fn = np.poly1d(coeffs)
    y_pred = fit_fn(invA)
    R2 = 1 - np.sum((Ps - y_pred)**2) / np.sum((Ps - Ps.mean())**2)

    cmap = plt.get_cmap("tab10")
    plt.figure(figsize=(8,6))
    for i, (x, y, L) in enumerate(zip(invA, Ps, Ls)):
        plt.scatter(x, y, label=f"L = {L:.2f}", s=100, edgecolor="k", color=cmap(i))

    x_lin = np.linspace(invA.min()*0.95, invA.max()*1.05, 200)
    plt.plot(x_lin, fit_fn(x_lin), "k--", linewidth=2, label="Ajuste del Modelo Lineal")

    plt.xlabel("1 / Área total (1/m²)", fontsize=16)
    plt.ylabel("Presión (kg/s²)", fontsize=16)
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.show()

def main():
    simulation_files = [
        ("wall_collisions_L_0.03.txt", 0.03),
        ("wall_collisions_L_0.05.txt", 0.05),
        ("wall_collisions_L_0.07.txt", 0.07),
        ("wall_collisions_L_0.09.txt", 0.09),
    ]

    resultados = []
    for filename, L in simulation_files:
        presiones, system, P_system_prom = process_simulation(filename, L)
        create_plots(presiones, L)
        print(f"{filename} -> L={L}, t_max={system['tiempo'].max():.1f}s, P_prom={P_system_prom}")
        resultados.append((L, P_system_prom))

    plot_pressure_vs_inv_area(resultados)

if __name__ == "__main__":
    main()
