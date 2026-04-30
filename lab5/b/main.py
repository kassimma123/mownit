import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 12

# dane do zadania
m_param = 5.0
k_param = 0.5
INTERVAL = (-5.0, 5.0)
L = INTERVAL[1] - INTERVAL[0] # Długość przedziału = 10
X_DENSE = np.linspace(INTERVAL[0], INTERVAL[1], 1000)
SAVE_DIR = "wykresy_trygonometryczna"

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

def f(x):
    # funkcja z zadania
    return x**2 - m_param * np.cos((np.pi * x) / k_param)

def get_nodes(n, include_endpoint=False):
    x = np.linspace(INTERVAL[0], INTERVAL[1], n, endpoint=True)
    if not include_endpoint:
        x = x[:-1]
    return x, f(x)

# algorytm aproksymacji (z uzyciem macierzy Grama)
def trig_approximation(x_nodes, y_nodes, m):
    n = len(x_nodes)
    if n < 2 * m + 1:
        return None, None

    # transformacja na przedzial [-pi, pi] zeby uniknac zjawiska Gibbsa na brzegach
    v_nodes = 2 * np.pi * (x_nodes - INTERVAL[0]) / L - np.pi

    B = [np.ones_like(v_nodes) / np.sqrt(2)]
    for k in range(1, m + 1):
        B.append(np.cos(k * v_nodes))
        B.append(np.sin(k * v_nodes))
    
    B = np.column_stack(B)
    
    G = B.T @ B
    b = B.T @ y_nodes

    try:
        # rozwiazanie G * c = b
        coeffs = np.linalg.solve(G, b)
    except np.linalg.LinAlgError:
        return None, None

    def approx_func(x):
        v = 2 * np.pi * (x - INTERVAL[0]) / L - np.pi
        res = coeffs[0] / np.sqrt(2)
        idx = 1
        for k in range(1, m + 1):
            res += coeffs[idx] * np.cos(k * v) + coeffs[idx+1] * np.sin(k * v)
            idx += 2
        return res

    return approx_func, np.linalg.cond(G)

# funkcje do rysowania wykresow do sprawozdania

def plot_1_interpolation_comparison():
    # rysuje wykres 1 (przypadek n = 2m+1)
    n, m = 22, 10
    
    x_inc, y_inc = get_nodes(n, include_endpoint=True)
    x_exc, y_exc = get_nodes(n, include_endpoint=False)
    
    approx_inc, cond_inc = trig_approximation(x_inc, y_inc, m)
    approx_exc, cond_exc = trig_approximation(x_exc, y_exc, m)
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    # Wykres 1: Z węzłem końcowym
    ax[0].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał f(x)')
    ax[0].scatter(x_inc, y_inc, c='red', s=40, zorder=5, label=f'Węzły (n={n})')
    if approx_inc: ax[0].plot(X_DENSE, approx_inc(X_DENSE), 'r', label="Aproksymacja")
    ax[0].set_title(f"Z WĘZŁEM x=5\nMacierz układu równań bliska osobliwości (cond={cond_inc:.1e})", color='darkred')
    ax[0].legend(loc='lower center')
    ax[0].set_ylim(-10, 35)
    
    # Wykres 2: BEZ węzła końcowego
    ax[1].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał f(x)')
    ax[1].scatter(x_exc, y_exc, c='blue', s=40, zorder=5, label=f'Węzły (n={n})')
    if approx_exc: ax[1].plot(X_DENSE, approx_exc(X_DENSE), 'b', label="Aproksymacja")
    ax[1].set_title(f"BEZ WĘZŁA x=5\nIdealna dyskretna ortogonalność (cond={cond_exc:.1f})", color='darkblue')
    ax[1].legend(loc='lower center')
    ax[1].set_ylim(-10, 35)
    
    plt.suptitle("Etap 1: Warunek rozwiązalności układu równań (Interpolacja n=21, m=10)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "1_interpolacja_porownanie.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_2_error_comparison():
    # wykres z bledami na osi y w skali logarytmicznej
    m = 5
    n_vals = np.arange(12, 40)
    err_inc, err_exc = [], []
    
    for n in n_vals:
        x_in, y_in = get_nodes(n, include_endpoint=True)
        app_in, _ = trig_approximation(x_in, y_in, m)
        err_inc.append(np.max(np.abs(f(X_DENSE) - app_in(X_DENSE))) if app_in else np.nan)
            
        x_ex, y_ex = get_nodes(n, include_endpoint=False)
        app_ex, _ = trig_approximation(x_ex, y_ex, m)
        err_exc.append(np.max(np.abs(f(X_DENSE) - app_ex(X_DENSE))) if app_ex else np.nan)
            
    plt.figure(figsize=(10, 6))
    plt.plot(n_vals, err_inc, 'r-o', label='Z WĘZŁEM (endpoint=True)')
    plt.plot(n_vals, err_exc, 'b-o', label='BEZ WĘZŁA (endpoint=False)')
    
    plt.yscale('log')
    plt.xlabel('Liczba węzłów (n)')
    plt.ylabel('Błąd maksymalny (skala logarytmiczna)')
    plt.title(f'Etap 2: Błąd aproksymacji (m={m})\nDla zbyt małego stopnia bazy oba modele dają ogromny błąd, a wariant bez węzła jest nieznacznie gorszy.', fontsize=13)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig(os.path.join(SAVE_DIR, "2_blad_porownanie.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_2b_error_comparison_proper_m():
    # wykres pokazujacy udana probe przy prawidlowym m=12
    m = 12
    n_vals = np.arange(26, 60)
    err_inc, err_exc = [], []
    
    for n in n_vals:
        x_in, y_in = get_nodes(n, include_endpoint=True)
        app_in, _ = trig_approximation(x_in, y_in, m)
        err_inc.append(np.max(np.abs(f(X_DENSE) - app_in(X_DENSE))) if app_in else np.nan)
            
        x_ex, y_ex = get_nodes(n, include_endpoint=False)
        app_ex, _ = trig_approximation(x_ex, y_ex, m)
        err_exc.append(np.max(np.abs(f(X_DENSE) - app_ex(X_DENSE))) if app_ex else np.nan)
            
    plt.figure(figsize=(10, 6))
    plt.plot(n_vals, err_inc, 'r-o', label='Z WĘZŁEM (endpoint=True)')
    plt.plot(n_vals, err_exc, 'b-o', label='BEZ WĘZŁA (endpoint=False)')
    
    plt.yscale('log')
    plt.xlabel('Liczba węzłów (n)')
    plt.ylabel('Błąd maksymalny (skala logarytmiczna)')
    plt.title(f'Etap 2b: Błąd aproksymacji po pokonaniu underfittingu (m={m})\nBłąd drastycznie spada, a ortogonalna siatka (bez węzła) staje się docelowo lepsza.', fontsize=13)
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig(os.path.join(SAVE_DIR, "2b_blad_porownanie_m12.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_3_extrema_threshold():
    # pokazujemy jak dla m=10 nagle model lapie wszystkie ekstrema
    n = 31
    x_nodes, y_nodes = get_nodes(n, include_endpoint=False)
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    m1 = 9
    app1, _ = trig_approximation(x_nodes, y_nodes, m1)
    ax[0].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał z falą cos(2πx)')
    ax[0].scatter(x_nodes, y_nodes, c='black', s=20, label=f'Węzły (n={n})')
    if app1: ax[0].plot(X_DENSE, app1(X_DENSE), 'orange', linewidth=2, label=f'Aproksymacja m={m1}')
    ax[0].set_title(f"Zbyt niski stopień m={m1}\n(Brak detekcji wysokiej częstotliwości)")
    ax[0].legend(loc='lower center')
    ax[0].set_ylim(-10, 35)
    
    m2 = 10
    app2, _ = trig_approximation(x_nodes, y_nodes, m2)
    ax[1].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał z falą cos(2πx)')
    ax[1].scatter(x_nodes, y_nodes, c='black', s=20, label=f'Węzły (n={n})')
    if app2: ax[1].plot(X_DENSE, app2(X_DENSE), 'green', linewidth=2, label=f'Aproksymacja m={m2}')
    ax[1].set_title(f"Trafienie w dziesiątą harmonikę m={m2}\n(Idealne pokrycie oscylacji)")
    ax[1].legend(loc='lower center')
    ax[1].set_ylim(-10, 35)
    
    plt.suptitle("Etap 3a: Moment detekcji ekstremów (Granica to m=10 ze względu na częstotliwość fali)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "3a_detekcja_prog.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_4_extrema_with_without():
    # jak zachowuje sie przy optymalnym m=10 dla z wezlem i bez
    m = 10 
    n = 31 # Bierzemy n=31, gdzie różnica błędu brzegowego jest największa i wyraźnie widoczna
    
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    # Lewa strona - Z WĘZŁEM
    x_in, y_in = get_nodes(n, include_endpoint=True)
    app_in, _ = trig_approximation(x_in, y_in, m)
    ax[0].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał')
    ax[0].scatter(x_in, y_in, c='red', s=40, label=f'Węzły (n={len(x_in)})', zorder=5)
    if app_in: 
        err_in = np.max(np.abs(f(X_DENSE) - app_in(X_DENSE)))
        ax[0].plot(X_DENSE, app_in(X_DENSE), 'r', linewidth=2, label=f'Aproksymacja m={m}')
        ax[0].set_title(f"Z WĘZŁEM x=5\nZaburzenia widoczne na krańcach (Max błąd: {err_in:.2f})", color='darkred')
    
    ax[0].legend(loc='lower left')
    # Zbliżenie na prawy kraniec, żeby pokazać "telepanie"
    ax[0].set_xlim(3.0, 5.2)
    ax[0].set_ylim(8, 28)
    
    # Prawa strona - BEZ WĘZŁA
    x_ex, y_ex = get_nodes(n, include_endpoint=False)
    app_ex, _ = trig_approximation(x_ex, y_ex, m)
    ax[1].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5, label='Oryginał')
    ax[1].scatter(x_ex, y_ex, c='blue', s=40, label=f'Węzły (n={len(x_ex)})', zorder=5)
    if app_ex: 
        err_ex = np.max(np.abs(f(X_DENSE) - app_ex(X_DENSE)))
        ax[1].plot(X_DENSE, app_ex(X_DENSE), 'b', linewidth=2, label=f'Aproksymacja m={m}')
        ax[1].set_title(f"BEZ WĘZŁA x=5\nKrzywa stabilniejsza na krańcach (Max błąd: {err_ex:.2f})", color='darkblue')
    
    ax[1].legend(loc='lower left')
    # Zbliżenie na prawy kraniec
    ax[1].set_xlim(3.0, 5.2)
    ax[1].set_ylim(8, 28)
    
    plt.suptitle("Etap 3b: Zbliżenie na prawy kraniec (Warianty nadokreślone n=31)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "3b_detekcja_porownanie.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_5_n_impact_comparison():
    # rosnace n przy stalym m=10
    m = 10
    n_list = [22, 30, 50, 80]
    
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']
    
    # Lewy wykres: Z węzłem
    ax[0].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.4, label='Oryginał', linewidth=3)
    for i, n in enumerate(n_list):
        x_in, y_in = get_nodes(n, include_endpoint=True)
        app, _ = trig_approximation(x_in, y_in, m)
        if app: ax[0].plot(X_DENSE, app(X_DENSE), color=colors[i], label=f'n={n}', alpha=0.8, linewidth=1.5)
    ax[0].set_title("Z WĘZŁEM\nDla n > 2m+1 metoda najmniejszych kwadratów świetnie sobie radzi", color='darkred')
    ax[0].legend(loc='lower center', ncol=2)
    ax[0].set_ylim(-10, 35)

    # Prawy wykres: BEZ węzła
    ax[1].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.4, label='Oryginał', linewidth=3)
    for i, n in enumerate(n_list):
        x_ex, y_ex = get_nodes(n, include_endpoint=False)
        app, _ = trig_approximation(x_ex, y_ex, m)
        if app: ax[1].plot(X_DENSE, app(X_DENSE), color=colors[i], label=f'n={n}', alpha=0.8, linewidth=1.5)
    ax[1].set_title("BEZ WĘZŁA\nRównież stabilna od samego początku, różnice są marginalne", color='darkblue')
    ax[1].legend(loc='lower center', ncol=2)
    ax[1].set_ylim(-10, 35)
    
    plt.suptitle(f"Etap 4: Wpływ zagęszczania siatki n przy stałym stopniu m={m}", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "4_wplyw_n_porownanie.png"), bbox_inches='tight', dpi=150)
    plt.close()

def plot_6_heatmaps_comparison():
    # mapy ciepla dla n i m
    # n bierzemy co 5 zeby bylo szybciej
    n_vals = list(range(12, 105, 5))
    m_vals = list(range(2, 41, 2))
    
    err_inc = np.full((len(n_vals), len(m_vals)), np.nan)
    err_exc = np.full((len(n_vals), len(m_vals)), np.nan)

    for i, n in enumerate(n_vals):
        x_in, y_in = get_nodes(n, include_endpoint=True)
        x_ex, y_ex = get_nodes(n, include_endpoint=False)
        
        for j, m in enumerate(m_vals):
            if len(x_in) >= 2*m + 1:
                app_in, _ = trig_approximation(x_in, y_in, m)
                if app_in:
                    err = np.max(np.abs(f(X_DENSE) - np.nan_to_num(app_in(X_DENSE), nan=1e10)))
                    err_inc[i, j] = err if err < 1e5 else np.nan
                    
            if len(x_ex) >= 2*m + 1:
                app_ex, _ = trig_approximation(x_ex, y_ex, m)
                if app_ex:
                    err = np.max(np.abs(f(X_DENSE) - np.nan_to_num(app_ex(X_DENSE), nan=1e10)))
                    err_exc[i, j] = err if err < 1e5 else np.nan

    fig, axes = plt.subplots(1, 2, figsize=(22, 10))
    
    # Rozszerzona skala od 1e-4 do 1e2 by uwydatnić subtelne różnice
    norm = plt.matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e2)
    
    sns.heatmap(err_inc, xticklabels=m_vals, yticklabels=n_vals, annot=True, 
                fmt=".1e", annot_kws={"size": 6}, cmap="turbo", norm=norm, ax=axes[0])
    axes[0].set_title("Błąd Maksymalny - Z WĘZŁEM\nZnacznie gorsza precyzja, strefy dużego błędu", color='darkred', fontsize=14)
    axes[0].set_xlabel("Stopień szeregu (m)", fontsize=12)
    axes[0].set_ylabel("Liczba węzłów (n)", fontsize=12)

    sns.heatmap(err_exc, xticklabels=m_vals, yticklabels=n_vals, annot=True, 
                fmt=".1e", annot_kws={"size": 6}, cmap="turbo", norm=norm, ax=axes[1])
    axes[1].set_title("Błąd Maksymalny - BEZ WĘZŁA\nStabilny spadek błędu do 10^-3, duża precyzja", color='darkblue', fontsize=14)
    axes[1].set_xlabel("Stopień szeregu (m)", fontsize=12)
    axes[1].set_ylabel("Liczba węzłów (n)", fontsize=12)

    plt.suptitle("Etap 5: Gęsta Mapa Ciepła Błędu Maksymalnego (Z Węzłem vs Bez Węzła)", fontsize=18, y=1.03)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, "5_mapy_ciepla_porownanie.png"), bbox_inches='tight', dpi=150)
    plt.close()

if __name__ == "__main__":
    print("Generowanie pełnego porównania z/bez węzła...")
    plot_1_interpolation_comparison()
    plot_2_error_comparison()
    plot_2b_error_comparison_proper_m()
    plot_3_extrema_threshold()
    plot_4_extrema_with_without()
    plot_5_n_impact_comparison()
    plot_6_heatmaps_comparison()
    print("Zakończono. Pliki w folderze", SAVE_DIR)