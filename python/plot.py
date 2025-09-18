import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# Parámetros del sistema
# ----------------------------
m = 1.0       # masa de cada partícula
dt_bin = 1.0  # tamaño de ventana temporal para promediar

def process_simulation(filename, L):
    """
    Procesa un archivo de simulación y devuelve:
    - dataframe con presiones vs tiempo
    - presión promedio (estado estacionario)
    """
    df = pd.read_csv(
        filename,
        names=["time", "compartment", "normal_velocity", "radius"],
        dtype={"time": float, "compartment": int, "normal_velocity": float, "radius": float},
        skiprows=1
    )
    
    # Bin de tiempo
    df["bin"] = (df["time"] // dt_bin).astype(int)
    
    # impulso transferido
    df["impulso"] = 2 * m * df["normal_velocity"]
    
    # presión
    presiones = df.groupby(["bin", "compartment"])["impulso"].sum().reset_index()
    presiones["tiempo"] = presiones["bin"] * dt_bin
    presiones["P"] = presiones["impulso"] / (L * dt_bin)
    
    # presión promedio (estado estacionario → después de t > 20 s por ejemplo)
    estacionario = presiones[presiones["tiempo"] > 20]
    P_prom = estacionario.groupby("compartment")["P"].mean().to_dict()
    
    return presiones, P_prom


def create_plots(presiones, L):
    """Genera gráfico de presión vs tiempo para un L dado."""
    plt.figure(figsize=(10,6))
    for comp, color, label in [(0, "blue", "Cuadrado"), (1, "red", "Rectángulo")]:
        subset = presiones[presiones["compartment"] == comp]
        plt.plot(subset["tiempo"], subset["P"], color=color, label=label)
    
    plt.xlabel("Tiempo (s)")
    plt.ylabel("Presión (kg/s²)")
    plt.title(f"Gráfico de la Presión respecto al tiempo\nL = {L:.2f} m")
    plt.legend()
    plt.grid(True)
    plt.show()


def main():
    simulation_files = [
        ("wall_collisions_L_0.03.txt", 0.03),
        ("wall_collisions_L_0.05.txt", 0.05),
        ("wall_collisions_L_0.07.txt", 0.07),
        ("wall_collisions_L_0.09.txt", 0.09),
    ]
    
    resultados = []  # para guardar presiones promedio
    
    for filename, L in simulation_files:
        presiones, P_prom = process_simulation(filename, L)
        create_plots(presiones, L)
        
        resultados.append((L, P_prom))
    
    print("Presiones promedio en estado estacionario:")
    for L, P_prom in resultados:
        print(f"L = {L}: {P_prom}")


if __name__ == "__main__":
    main()
