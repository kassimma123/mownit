import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
import os

warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13

# ==========================================
# 1. PARAMETRY I FUNKCJA
# ==========================================
m_param = 5.0
k_param = 0.5
INTERVAL = (-5.0, 5.0)
X_DENSE = np.linspace(INTERVAL[0], INTERVAL[1], 1000)
SAVE_DIR = "wykresy_sprawozdanie"

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

def f(x):
    return x**2 - m_param * np.cos((np.pi * x) / k_param)

def get_nodes(n):
    x = np.linspace(INTERVAL[0], INTERVAL[1], n)
    return x, f(x)

def poly_approximation(x_nodes, y_nodes, degree):
    m = degree
    if len(x_nodes) <= m: return None 

    x_powers = [np.sum(x_nodes**k) for k in range(2 * m + 1)]
    gram_matrix = np.zeros((m + 1, m + 1))
    for j in range(m + 1):
        for k in range(m + 1):
            gram_matrix[j, k] = x_powers[j + k]

    moment_vector = np.zeros(m + 1)
    for j in range(m + 1):
        moment_vector[j] = np.sum(y_nodes * (x_nodes**j))

    try:
        coeffs = np.linalg.solve(gram_matrix, moment_vector)
    except np.linalg.LinAlgError:
        coeffs = np.zeros(m + 1) 

    def approx_func(x):
        return sum(coeffs[i] * x**i for i in range(m + 1))
    return approx_func

def calculate_mse(func1, func2, x_vals):
    y_approx = np.nan_to_num(func2(x_vals), nan=1e20, posinf=1e20, neginf=-1e20)
    return np.mean((func1(x_vals) - y_approx)**2)

def calculate_max_error(func1, func2, x_vals):
    y_approx = np.nan_to_num(func2(x_vals), nan=1e20, posinf=1e20, neginf=-1e20)
    return np.max(np.abs(func1(x_vals) - y_approx))

# ==========================================
# 2. WYZNACZANIE WYKRESÓW DLA STAŁEGO N
# ==========================================
def draw_single_isolated_panel(ax, n, m, color):
    x_nodes, y_nodes = get_nodes(n)
    ax.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginał', alpha=0.5, linewidth=1.5)
    
    m_size = max(5, 100 // n)
    m_alpha = 1.0 if n <= 20 else 0.4
    ax.scatter(x_nodes, y_nodes, color='black', s=m_size, alpha=m_alpha, label=f'Węzły', zorder=5)
    
    pf = poly_approximation(x_nodes, y_nodes, m)
    if pf:
        y_approx = pf(X_DENSE)
        y_approx = np.clip(y_approx, -20, 50) 
        ax.plot(X_DENSE, y_approx, color=color, label=f'Aproks. (m={m})', linewidth=2.5)
    else:
         ax.text(0, 15, "Zbyt mało węzłów!", ha='center', color='red', fontsize=12)
        
    ax.set_ylim(-10, 40)
    ax.set_title(f"Stopień m={m}")
    ax.legend(loc='lower center', fontsize=9, ncol=3)

def generate_isolated_grid(n, m_list, filename):
    cols = 3
    rows = (len(m_list) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
    if rows == 1: axes = axes.reshape(1, -1)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']
    for idx, m in enumerate(m_list):
        r, c = idx // cols, idx % cols
        draw_single_isolated_panel(axes[r, c], n, m, colors[idx % len(colors)])
        
    for idx in range(len(m_list), rows * cols):
        fig.delaxes(axes[idx // cols, idx % cols])
        
    plt.suptitle(f"Analiza odizolowana: Stała liczba węzłów n={n}", fontsize=18, y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, filename), bbox_inches='tight', dpi=150)
    plt.close()

print("Generowanie analizy dla stałego n...")
generate_isolated_grid(n=10, m_list=[2, 3, 4, 5, 6, 7], filename="fig_n10_stale.png")
generate_isolated_grid(n=20, m_list=[8, 10, 12, 14, 16, 18], filename="fig_n20_stale.png")
generate_isolated_grid(n=50, m_list=[10, 15, 20, 25, 30, 40], filename="fig_n50_stale.png")
generate_isolated_grid(n=100, m_list=[20, 30, 40, 50, 60, 80], filename="fig_n100_stale.png")


# ==========================================
# 3. WYZNACZANIE WYKRESÓW DLA STAŁEGO M
# ==========================================
def plot_panel_fixed_m(ax, m, n_list, title_prefix=""):
    ax.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginalna', alpha=0.6)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, n in enumerate(n_list):
        x_nodes, y_nodes = get_nodes(n)
        m_size = max(5, 100 // n)
        if idx == 0: ax.scatter(x_nodes, y_nodes, color='black', s=m_size, alpha=0.5, label='Węzły')
            
        pf = poly_approximation(x_nodes, y_nodes, m)
        if pf:
            y_approx = np.clip(pf(X_DENSE), -20, 50)
            ax.plot(X_DENSE, y_approx, color=colors[idx % len(colors)], label=f'n={n}', linewidth=1.5)
            
    ax.set_ylim(-10, 40)
    ax.set_title(f"Stały stopień m={m}")
    ax.legend(loc='lower center', fontsize=9, ncol=2)

print("Generowanie analizy dla stałego m...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
plot_panel_fixed_m(axes[0,0], 5, [10, 20, 50])
plot_panel_fixed_m(axes[0,1], 10, [15, 30, 60])
plot_panel_fixed_m(axes[1,0], 20, [25, 50, 100])
plot_panel_fixed_m(axes[1,1], 40, [45, 70, 100])
plt.suptitle("Aproksymacja dla stałego stopnia wielomianu (m)", fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "fig_stale_m.png"), bbox_inches='tight', dpi=150)
plt.close()

# ==========================================
# 4. WYKRES SPADKU BŁĘDÓW
# ==========================================
def run_error_vs_nodes_plot():
    n_list = list(range(10, 105, 5))
    m_tests = [4, 8, 12] 
    
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'red', 'green']
    
    for m, col in zip(m_tests, colors):
        errors = []
        ns = []
        for n_test in n_list:
            if m < n_test:
                xn, yn = get_nodes(n_test)
                pf = poly_approximation(xn, yn, m)
                if pf:
                    mse = calculate_mse(f, pf, X_DENSE)
                    errors.append(mse)
                    ns.append(n_test)
        plt.plot(ns, errors, marker='o', linestyle='-', color=col, label=f'Stopień m={m}')
        
    plt.title('Spadek błędu MSE w miarę zagęszczania węzłów (n)')
    plt.xlabel('Liczba węzłów (n)')
    plt.ylabel('Błąd średniokwadratowy MSE (Skala Log)')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(SAVE_DIR, "blad_vs_wezly.png"), dpi=150)
    plt.close()

print("Generowanie wykresu spadku błędów...")
run_error_vs_nodes_plot()
# ==========================================
# 5. NOWE MAPY CIEPŁA (DRASTYCZNIE POPRAWIONA SKALA)
# ==========================================
def generate_heatmaps_and_table():
    n_values = list(range(10, 105, 10))
    m_values = list(range(2, 41, 2))
    
    mse_results = np.full((len(n_values), len(m_values)), np.nan)
    max_results = np.full((len(n_values), len(m_values)), np.nan)
    table_data = []

    for i, n_val in enumerate(n_values):
        xn, yn = get_nodes(n_val)
        for j, m_val in enumerate(m_values):
            if m_val < n_val:
                pf = poly_approximation(xn, yn, m_val)
                if pf:
                    mse = calculate_mse(f, pf, X_DENSE)
                    max_err = calculate_max_error(f, pf, X_DENSE)
                    
                    mse_results[i, j] = mse
                    max_results[i, j] = max_err
                    
                    if n_val in [20, 50, 100] and m_val in [10, 18, 30]:
                        table_data.append({"n": n_val, "m": m_val, "MSE": mse, "E_max": max_err})

    print("Rysowanie mapy MSE...")
    plt.figure(figsize=(18, 10))
    
    # KLUCZOWA ZMIANA: vmin=1e-1 (0.1), vmax=1e2 (100). 
    # To rozciągnie kolory IDEALNIE na dobre wyniki. Reszta będzie żółta.
    sns.heatmap(mse_results, xticklabels=m_values, yticklabels=n_values, 
                annot=True, fmt=".1e", annot_kws={"size": 8}, cmap="viridis_r", 
                norm=plt.matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e2))
                
    plt.title("Błąd Średniokwadratowy (MSE) - Skala obcięta do błędu 100", fontsize=16)
    plt.xlabel("Stopień wielomianu (m)", fontsize=14)
    plt.ylabel("Liczba węzłów (n)", fontsize=14)
    plt.savefig(os.path.join(SAVE_DIR, "heatmap_mse.png"), bbox_inches='tight', dpi=150)
    plt.close()
    
    print("Rysowanie mapy Max Error...")
    plt.figure(figsize=(18, 10))
    
    # KLUCZOWA ZMIANA: vmin=1e0 (1.0), vmax=1e3 (1000).
    sns.heatmap(max_results, xticklabels=m_values, yticklabels=n_values, 
                annot=True, fmt=".1e", annot_kws={"size": 8}, cmap="plasma", 
                norm=plt.matplotlib.colors.LogNorm(vmin=1e0, vmax=1e3))
                
    plt.title("Błąd Maksymalny (E_max) - Skala obcięta do błędu 1000", fontsize=16)
    plt.xlabel("Stopień wielomianu (m)", fontsize=14)
    plt.ylabel("Liczba węzłów (n)", fontsize=14)
    plt.savefig(os.path.join(SAVE_DIR, "heatmap_max.png"), bbox_inches='tight', dpi=150)
    plt.close()

    df = pd.DataFrame(table_data)
    df['MSE'] = df['MSE'].map('{:.3e}'.format)
    df['E_max'] = df['E_max'].map('{:.3e}'.format)
    df.to_csv(os.path.join(SAVE_DIR, "tabela_bledow.csv"), index=False)

generate_heatmaps_and_table()

# ==========================================
# 6. WARIANT OPTYMALNY Z MATEMATYCZNYM OPISEM
# ==========================================
def best_fit_plot():
    n = 51
    m = 20
    x_nodes, y_nodes = get_nodes(n)
    pf = poly_approximation(x_nodes, y_nodes, m)
    
    plt.figure(figsize=(10, 6))
    plt.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginalna', alpha=0.6)
    plt.scatter(x_nodes, y_nodes, color='black', s=15, label=f'Węzły (n={n})')
    
    if pf:
        y_approx = pf(X_DENSE)
        mse = calculate_mse(f, pf, X_DENSE)
        plt.plot(X_DENSE, y_approx, 'b-', label=f'Aproksymacja optymalna')
        
        # Opis umieszczony bezpośrednio na wykresie
        plt.text(-4.8, 25, f"Zoptymalizowano do MSE: {mse:.3f}\n\nUzasadnienie (Kompromis):\nZwiększanie stopnia 'm' do 20 minimalizuje\nbłąd samej metody (lepsze dopasowanie).\nPrzekroczenie m=20 powoduje wybuch\nbłędów zaokrągleń macierzy Grama.", 
                 bbox=dict(facecolor='white', alpha=0.9, edgecolor='black', boxstyle='round,pad=0.5'))
                 
        plt.title(f"Aproksymacja optymalna (Zbalansowana): n={n}, m={m}")
    
    plt.legend()
    plt.ylim(-10, 40)
    plt.savefig(os.path.join(SAVE_DIR, "najlepsze_dopasowanie.png"), dpi=150)
    plt.close()

print("Generowanie wariantu optymalnego...")
best_fit_plot()
print("Gotowe! Pliki zapisane w folderze:", SAVE_DIR)