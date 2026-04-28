import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
import os

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11

# ==========================================
# 1. PARAMETRY I FUNKCJA
# ==========================================
m_param = 5.0
k_param = 0.5
INTERVAL = (-5.0, 5.0)
X_DENSE = np.linspace(INTERVAL[0], INTERVAL[1], 1000)
SAVE_DIR = "wykresy_trygonometryczna"

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

def f(x):
    return x**2 - m_param * np.cos((np.pi * x) / k_param)

def get_nodes(n, include_endpoint=True):
    # Generowanie węzłów z lub bez prawego krańca przedziału
    x = np.linspace(INTERVAL[0], INTERVAL[1], n, endpoint=include_endpoint)
    return x, f(x)

# ==========================================
# 2. SILNIK APROKSYMACJI TRYG. (MACIERZ GRAMA)
# ==========================================
def trig_approximation(x_nodes, y_nodes, m):
    n = len(x_nodes)
    if n <= 2 * m: # Warunek: n > 2m (lub n >= 2m+1)
        return None

    # Mapowanie [-5, 5] na [0, pi] (półokres - omija zjawisko Gibbsa na brzegach)
    v_nodes = np.pi * (x_nodes - INTERVAL[0]) / (INTERVAL[1] - INTERVAL[0])

    # Baza trygonometryczna: [1, cos(v), sin(v), cos(2v), sin(2v), ...]
    B = [np.ones_like(v_nodes)]
    for k in range(1, m + 1):
        B.append(np.cos(k * v_nodes))
        B.append(np.sin(k * v_nodes))
    
    B = np.column_stack(B)
    
    # Macierz Grama i wektor wyrazów wolnych
    G = B.T @ B
    b = B.T @ y_nodes

    try:
        # Rozwiązanie układu metodą rozkładu LU
        coeffs = np.linalg.solve(G, b)
    except np.linalg.LinAlgError:
        return None

    def approx_func(x):
        v = np.pi * (x - INTERVAL[0]) / (INTERVAL[1] - INTERVAL[0])
        res = coeffs[0]
        idx = 1
        for k in range(1, m + 1):
            res += coeffs[idx] * np.cos(k * v) + coeffs[idx+1] * np.sin(k * v)
            idx += 2
        return res

    return approx_func

def calculate_mse(func1, func2, x_vals):
    y_approx = np.nan_to_num(func2(x_vals), nan=1e10, posinf=1e10, neginf=-1e10)
    return np.mean((func1(x_vals) - y_approx)**2)

def calculate_max_error(func1, func2, x_vals):
    y_approx = np.nan_to_num(func2(x_vals), nan=1e10, posinf=1e10, neginf=-1e10)
    return np.max(np.abs(func1(x_vals) - y_approx))

# ==========================================
# 3. ANALIZA: WĘZEŁ KOŃCOWY
# ==========================================
def plot_endpoint_comparison():
    n, m = 45, 20
    x_inc, y_inc = get_nodes(n, include_endpoint=True)
    x_exc, y_exc = get_nodes(n, include_endpoint=False)
    
    approx_inc = trig_approximation(x_inc, y_inc, m)
    approx_exc = trig_approximation(x_exc, y_exc, m)
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    
    ax[0].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.3)
    ax[0].scatter(x_inc, y_inc, c='red', s=20, label='Węzły')
    if approx_inc: ax[0].plot(X_DENSE, approx_inc(X_DENSE), 'r', label="Aproksymacja")
    ax[0].set_title(f"Z węzłem końcowym x=5 (n={n})")
    ax[0].legend()
    
    ax[1].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.3)
    ax[1].scatter(x_exc, y_exc, c='blue', s=20, label='Węzły')
    if approx_exc: ax[1].plot(X_DENSE, approx_exc(X_DENSE), 'b', label="Aproksymacja")
    ax[1].set_title(f"Bez węzła końcowego x=5 (n={n})")
    ax[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "trig_wezel_koncowy.png"), dpi=150)
    plt.close()

# ==========================================
# 4. MOMENT DETEKCJI I ZAŁAMANIA
# ==========================================
def draw_panel(ax, n, m, color):
    x_n, y_n = get_nodes(n)
    ax.plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.4)
    # Dynamiczny rozmiar punktów, by się nie zlewały
    ax.scatter(x_n, y_n, color='black', s=max(3, 100//n), alpha=0.5, zorder=5)
    
    pf = trig_approximation(x_n, y_n, m)
    if pf:
        y_app = np.clip(pf(X_DENSE), -20, 50)
        ax.plot(X_DENSE, y_app, color=color, linewidth=2, label=f'm={m}')
    else:
        ax.text(0, 15, "Brak spełnienia\n n >= 2m+1", ha='center', color='red')
        
    ax.set_ylim(-10, 40)
    ax.set_title(f"Stopień m={m}")
    ax.legend(loc='lower center')

def generate_evolution(n, m_list, filename, title):
    cols = 3
    rows = (len(m_list) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4.5*rows))
    axes = axes.flatten()
    colors = sns.color_palette("husl", len(m_list))
    
    for i, m in enumerate(m_list):
        draw_panel(axes[i], n, m, colors[i])
        
    plt.suptitle(title, fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, filename), dpi=150)
    plt.close()

# ==========================================
# 5. MAPY CIEPŁA Z DANYMI (Ścisła Skala)
# ==========================================
def plot_heatmaps():
    n_vals = list(range(10, 105, 10))
    m_vals = list(range(2, 41, 2))
    
    mse_res = np.full((len(n_vals), len(m_vals)), np.nan)
    max_res = np.full((len(n_vals), len(m_vals)), np.nan)

    for i, n in enumerate(n_vals):
        x_n, y_n = get_nodes(n)
        for j, m in enumerate(m_vals):
            if n > 2*m:
                pf = trig_approximation(x_n, y_n, m)
                if pf:
                    mse_res[i, j] = calculate_mse(f, pf, X_DENSE)
                    max_res[i, j] = calculate_max_error(f, pf, X_DENSE)

    # Obcięcie skali vmin/vmax by uwydatnić drobne różnice i "wybuch" układu
    plt.figure(figsize=(16, 9))
    sns.heatmap(mse_res, xticklabels=m_vals, yticklabels=n_vals, annot=True, 
                fmt=".1e", annot_kws={"size": 7}, cmap="viridis_r", 
                norm=plt.matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e2))
    plt.title("Mapa Ciepła MSE (Skala uwydatniająca próg uwarunkowania)")
    plt.xlabel("Stopień szeregu (m)")
    plt.ylabel("Liczba węzłów (n)")
    plt.savefig(os.path.join(SAVE_DIR, "heatmap_trig_mse.png"), bbox_inches='tight', dpi=150)
    plt.close()

# ==========================================
# 6. NAJLEPSZE DOPASOWANIE
# ==========================================
def plot_best_fit():
    n, m = 100, 20
    x_n, y_n = get_nodes(n)
    pf = trig_approximation(x_n, y_n, m)
    
    plt.figure(figsize=(10, 6))
    plt.plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label="Oryginał")
    plt.scatter(x_n, y_n, c='black', s=5, alpha=0.5, label=f"Węzły (n={n})")
    
    if pf:
        mse = calculate_mse(f, pf, X_DENSE)
        plt.plot(X_DENSE, pf(X_DENSE), 'g-', linewidth=2, label=f"Aproks. (m={m})")
        plt.text(-4.8, 30, f"MSE: {mse:.4e}\nCzęstotliwość funkcji bazowej ulega\nidealnemu zrównaniu z 20-tą\nharmoniczną bazy trygonometrycznej.", 
                 bbox=dict(facecolor='white', alpha=0.9))
                 
    plt.title(f"Aproksymacja optymalna: n={n}, m={m}")
    plt.ylim(-10, 40)
    plt.legend()
    plt.savefig(os.path.join(SAVE_DIR, "trig_najlepsze.png"), dpi=150)
    plt.close()

# ==========================================
# WIZUALIZACJA DLA STAŁEGO STOPNIA m (Zmienne n)
# ==========================================
def plot_fixed_m_trig():
    m_values = [5, 10, 15, 20]
    # Zgodnie z warunkiem n >= 2m+1, n musi być odpowiednio duże
    n_lists = [
        [12, 15, 20, 50],    # dla m=5  (min n=11)
        [22, 30, 50, 70],    # dla m=10 (min n=21)
        [32, 40, 60, 90],    # dla m=15 (min n=31)
        [42, 50, 80, 100]    # dla m=20 (min n=41)
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for idx, m in enumerate(m_values):
        ax = axes[idx // 2, idx % 2]
        ax.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginalna', alpha=0.5)
        
        for n_idx, n in enumerate(n_lists[idx]):
            x_nodes, y_nodes = get_nodes(n, include_endpoint=True)
            if n_idx == 0:
                ax.scatter(x_nodes, y_nodes, color='black', s=10, alpha=0.5, label=f'Węzły (min)')
            
            pf = trig_approximation(x_nodes, y_nodes, m)
            if pf:
                ax.plot(X_DENSE, pf(X_DENSE), color=colors[n_idx], label=f'n={n}', linewidth=1.5)
                
        ax.set_ylim(-10, 40)
        ax.set_title(f"Stopień m={m}")
        ax.legend(loc='lower center', fontsize=9, ncol=2)

    plt.suptitle("Wpływ zagęszczania siatki (n) przy stałym stopniu (m)", fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "trig_stale_m.png"), dpi=150)
    plt.close()

if __name__ == "__main__":
    print("Start...")
    plot_endpoint_comparison()
    
    # Krokowe wyłapywanie momentu działania (Rygorystyczny warunek m <= (n-1)//2)
    
    # n=7  -> max m=3
    generate_evolution(n=7, m_list=[1, 2, 3], filename="trig_detekcja_n7.png", title="Zbyt mało węzłów (max m=3) dla n=7")
    
    # n=10 -> max m=4
    generate_evolution(n=10, m_list=[1, 2, 3, 4], filename="trig_detekcja_n10.png", title="Osiąganie granicy uwarunkowania dla n=10")
    
    # n=12 -> max m=5
    generate_evolution(n=12, m_list=[1, 2, 3, 4, 5], filename="trig_detekcja_n12.png", title="Próby detekcji fali dla n=12")
    
    # n=15 -> max m=7
    generate_evolution(n=15, m_list=[2, 3, 4, 5, 6, 7], filename="trig_detekcja_n15.png", title="Stopniowa detekcja do max m=7 dla n=15")
    
    # n=20 -> max m=9
    generate_evolution(n=20, m_list=[4, 5, 6, 7, 8, 9], filename="trig_detekcja_n20.png", title="Ewolucja dopasowania do max m=9 dla n=20")
    
    # n=30 -> max m=14
    generate_evolution(n=30, m_list=[8, 10, 11, 12, 13, 14], filename="trig_detekcja_n30.png", title="Przed detekcją (max m=14) dla n=30")
    
    # n=40 -> max m=19 (Tu był błąd w oryginalnej liście)
    generate_evolution(n=40, m_list=[14, 15, 16, 17, 18, 19], filename="trig_detekcja_n40.png", title="Zbliżanie do pełnej detekcji przy n=40 (max m=19)")
    
    # n=50 -> max m=24
    generate_evolution(n=50, m_list=[15, 18, 20, 22, 23, 24], filename="trig_detekcja_n50.png", title="Pełna stabilność ekstremów dla n=50")
    
    # n=60 -> max m=29 (Poprawione)
    generate_evolution(n=60, m_list=[20, 22, 24, 26, 28, 29], filename="trig_detekcja_n60.png", title="Wysoka precyzja lokalizacji przy n=60 (max m=29)")
    
    # n=80 -> max m=39 (Poprawione)
    generate_evolution(n=80, m_list=[25, 30, 35, 37, 38, 39], filename="trig_detekcja_n80.png", title="Bardzo dobra lokalizacja przy n=80 (max m=39)")
    
    # n=90 -> max m=44 (Poprawione)
    generate_evolution(n=90, m_list=[30, 35, 40, 42, 43, 44], filename="trig_detekcja_n90.png", title="Praktyczna granica dla n=90 (max m=44)")
    
    # n=100 -> max m=49
    generate_evolution(n=100, m_list=[20, 30, 40, 45, 48, 49], filename="trig_granica_n100.png", title="Granica numeryczna dla gęstej siatki n=100")
    plot_heatmaps()
    plot_best_fit()
    plot_fixed_m_trig()
    print("Zrobione!")