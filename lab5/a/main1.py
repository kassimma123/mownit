import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
import os

# Ukrywamy ostrzeżenia o źle uwarunkowanych macierzach
# Robimy to celowo, bo w ramach laboratorium właśnie badamy to złe uwarunkowanie
warnings.filterwarnings('ignore')

# Ustawienie globalnego stylu, aby wykresy wyglądały profesjonalnie ("sprawozdaniowo")
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13

# ==========================================
# 1. PARAMETRY I FUNKCJA BAZOWA
# ==========================================
# Parametry funkcji zadane w poleceniu do laboratorium
m_param = 5.0
k_param = 0.5
INTERVAL = (-5.0, 5.0)

# Gęsta siatka punktów do rysowania płynnych wykresów "prawdziwej" funkcji
X_DENSE = np.linspace(INTERVAL[0], INTERVAL[1], 1000)
SAVE_DIR = "wykresy_sprawozdanie"

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

def f(x):
    """
    Nasza oryginalna, badana funkcja f(x).
    """
    return x**2 - m_param * np.cos((np.pi * x) / k_param)

def get_nodes(n):
    """
    Funkcja generująca n równoodległych węzłów aproksymacji na zadanym przedziale.
    Zwraca tablice (wektory) współrzędnych x oraz wartości y=f(x) dla tych węzłów.
    """
    x = np.linspace(INTERVAL[0], INTERVAL[1], n)
    return x, f(x)

# ==========================================
# 2. SILNIK APROKSYMACJI (Metoda Najmniejszych Kwadratów)
# ==========================================
def poly_approximation(x_nodes, y_nodes, degree):
    """
    Główna funkcja wykonująca aproksymację średniokwadratową wielomianami algebraicznymi.
    Zwraca wyliczoną funkcję wielomianową (callable), którą można potem narysować.
    """
    m = degree
    
    # Warunek konieczny: liczba węzłów (n) musi być ściśle większa od stopnia wielomianu (m)
    # W przeciwnym razie układ jest niedookreślony (mamy za mało punktów do narysowania tak "powyginanej" krzywej)
    if len(x_nodes) <= m:
        return None 

    # --- KROK 1: Budowa układu równań normalnych (Macierz Grama i Wektor wyrazów wolnych) ---
    # W aproksymacji średniokwadratowej nie szukamy wielomianu przechodzącego przez punkty, 
    # ale takiego, który minimalizuje sumę kwadratów błędów.
    # W tym celu budujemy układ postaci: A * c = b
    # gdzie:
    # A - macierz układu (tzw. Macierz Grama), zbudowana z sum potęg zmiennej x (od x^0 do x^(2m))
    # c - wektor niewiadomych współczynników naszego wielomianu 
    # b - wektor wyrazów wolnych, zależący od wartości x oraz f(x) (iloczyny f(x_i) * x_i^j)

    # Najpierw liczymy wszystkie potrzebne sumy potęg iksów, żeby nie liczyć ich w pętli wielokrotnie (optymalizacja)
    x_powers = [np.sum(x_nodes**k) for k in range(2 * m + 1)]
    
    # Inicjalizacja pustej macierzy Grama (rozmiar m+1 na m+1)
    gram_matrix = np.zeros((m + 1, m + 1))
    for j in range(m + 1):
        for k in range(m + 1):
            gram_matrix[j, k] = x_powers[j + k]

    # Budowa wektora wyrazów wolnych
    moment_vector = np.zeros(m + 1)
    for j in range(m + 1):
        moment_vector[j] = np.sum(y_nodes * (x_nodes**j))

    # --- KROK 2: Rozwiązywanie układu równań ---
    try:
        # DO ROZWIĄZANIA UKŁADU UŻYWAMY FUNKCJI BIBLIOTECZNEJ: numpy.linalg.solve()
        # Funkcja ta pod spodem wykorzystuje zaawansowany rozkład LU (z częściowym wyborem elementu głównego),
        # co gwarantuje w miarę dobrą stabilność numeryczną dla dobrze uwarunkowanych macierzy.
        coeffs = np.linalg.solve(gram_matrix, moment_vector)
        
    except np.linalg.LinAlgError:
        # Jeżeli stopień wielomianu (m) jest bardzo duży, wartości w macierzy Grama (typu x^30) 
        # są tak ogromne, że macierz staje się osobliwa (wyznacznik bliski zeru). 
        # Układ jest źle uwarunkowany i numpy nie potrafi go rozwiązać metodą LU.
        # W takim przypadku wstawiamy wektor zer.
        coeffs = np.zeros(m + 1) 

    # --- KROK 3: Zwrócenie wyniku jako funkcji ---
    # Definiujemy lokalnie nową funkcję, która po prostu wylicza wartość W(x) wstawiając policzone współczynniki c_k
    def approx_func(x):
        return sum(coeffs[i] * x**i for i in range(m + 1))
        
    return approx_func

# ==========================================
# 3. METRYKI BŁĘDÓW 
# ==========================================
def calculate_mse(func1, func2, x_vals):
    """
    Oblicza błąd średniokwadratowy (MSE) pomiędzy dwiema funkcjami.
    """
    # nan_to_num zabezpiecza przed wybuchem błędu do nieskończoności 
    y_approx = np.nan_to_num(func2(x_vals), nan=1e20, posinf=1e20, neginf=-1e20)
    return np.mean((func1(x_vals) - y_approx)**2)

def calculate_max_error(func1, func2, x_vals):
    """
    Oblicza maksymalne odchylenie (błąd maksymalny E_max) pomiędzy dwiema funkcjami.
    """
    y_approx = np.nan_to_num(func2(x_vals), nan=1e20, posinf=1e20, neginf=-1e20)
    return np.max(np.abs(func1(x_vals) - y_approx))


# ==========================================
# 4. NOWE, KLAROWNE NARZĘDZIA DO RYSOWANIA PANELI
# ==========================================
def draw_single_isolated_panel(ax, n, m, color):
    """Rysuje jeden, pojedynczy wykres dla konkretnego n i konkretnego m"""
    x_nodes, y_nodes = get_nodes(n)
    
    ax.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginał', alpha=0.4, linewidth=1.5)
    ax.scatter(x_nodes, y_nodes, color='black', s=20, label=f'Węzły (n={n})', zorder=5)
    
    pf = poly_approximation(x_nodes, y_nodes, m)
    if pf:
        y_approx = pf(X_DENSE)
        # Obcinamy wykresy w osi Y, żeby jeden "wystrzał" wynikający ze złego uwarunkowania 
        # nie spłaszczył nam całej osi i nie zepsuł rysunku.
        y_approx = np.clip(y_approx, -20, 50) 
        ax.plot(X_DENSE, y_approx, color=color, label=f'Aproks. (m={m})', linewidth=2)
    else:
         ax.text(0, 15, "Zbyt mało węzłów!", ha='center', color='red', fontsize=12)
        
    ax.set_ylim(-10, 40)
    ax.set_title(f"Stopień m={m}")
    ax.legend(loc='lower center', fontsize=9, ncol=3)

def generate_isolated_grid(n, m_list, filename):
    """Generuje całą siatkę (np. 2x3 lub 1x3) z izolowanymi wykresami"""
    num_plots = len(m_list)
    cols = 3
    rows = (num_plots + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
    
    # Bezpieczeństwo iteracji po osiach
    if rows == 1: axes = axes.reshape(1, -1)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
    
    for idx, m in enumerate(m_list):
        r = idx // cols
        c = idx % cols
        draw_single_isolated_panel(axes[r, c], n, m, colors[idx % len(colors)])
        
    # Usunięcie pustych, niewykorzystanych kratek w siatce
    for idx in range(len(m_list), rows * cols):
        r = idx // cols
        c = idx % cols
        fig.delaxes(axes[r, c])
        
    plt.suptitle(f"Analiza: Stała liczba węzłów n={n}", fontsize=18, y=1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(SAVE_DIR, filename), bbox_inches='tight', dpi=150)
    plt.close()


print(f"Zapisywanie wykresów do folderu: {SAVE_DIR} ...")

# ==========================================
# 5. GENEROWANIE WYKRESÓW - ZADANIA Z LABORATORIUM
# ==========================================

# --- ZADANIE 1: Zmiana stopnia m dla stałego n ---
generate_isolated_grid(n=10, m_list=[2, 3, 4, 5, 6, 7], filename="fig_n10_stale.png")
generate_isolated_grid(n=20, m_list=[2, 5, 8, 10, 15, 18], filename="fig_n20_stale.png")
generate_isolated_grid(n=50, m_list=[5, 10, 15, 20, 30, 40], filename="fig_n50_stale.png")
generate_isolated_grid(n=100, m_list=[10, 20, 30, 40, 60, 80], filename="fig_n100_stale.png")


# --- ZADANIE 2: Zmiana liczby węzłów n dla stałego m ---
def plot_panel_fixed_m(ax, m, n_list, title_prefix=""):
    ax.plot(X_DENSE, f(X_DENSE), 'k--', label='Oryginalna', alpha=0.6)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    for idx, n in enumerate(n_list):
        x_nodes, y_nodes = get_nodes(n)
        if idx == 0: ax.scatter(x_nodes, y_nodes, color='black', s=10, alpha=0.5, label='Węzły')
            
        pf = poly_approximation(x_nodes, y_nodes, m)
        if pf:
            y_approx = np.clip(pf(X_DENSE), -20, 50)
            ax.plot(X_DENSE, y_approx, color=colors[idx % len(colors)], label=f'Aproks. n={n}')
            
    ax.set_ylim(-10, 40)
    ax.set_title(f"{title_prefix} m={m}, n={n_list}")
    ax.legend(loc='lower center', fontsize=8, ncol=2)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
plot_panel_fixed_m(axes[0,0], 5, [10, 20, 50])
plot_panel_fixed_m(axes[0,1], 10, [15, 30, 60])
plot_panel_fixed_m(axes[1,0], 20, [25, 50, 80])
plot_panel_fixed_m(axes[1,1], 40, [45, 70, 100])
plt.suptitle("Aproksymacja dla stałego stopnia wielomianu (m)", fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(SAVE_DIR, "fig_stale_m.png"), bbox_inches='tight', dpi=150)
plt.close()


# --- ZADANIE 3: Wykres spadku błędu od zagęszczania węzłów ---
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

run_error_vs_nodes_plot()


# --- ZADANIE 4: Generowanie Map Ciepła i Tabeli ---
def generate_heatmaps_and_table():
    n_values = list(range(10, 105, 10))
    m_values = list(range(2, 41, 4))
    
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
                    
                    # Zabezpieczenie wizualne - ucinamy ogromne błędy żeby nie psuły kolorów mapy
                    mse_results[i, j] = max(1e-4, min(mse, 1e12))
                    max_results[i, j] = max(1e-4, min(max_err, 1e12))
                    
                    if n_val in [20, 50, 100] and m_val in [4, 10, 18, 30]:
                        table_data.append({"n": n_val, "m": m_val, "MSE": mse, "E_max": max_err})

    # Mapa Ciepła MSE
    plt.figure(figsize=(12, 8))
    sns.heatmap(mse_results, xticklabels=m_values, yticklabels=n_values, 
                annot=False, cmap="viridis_r", 
                norm=plt.matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e8))
    plt.title("Błąd Średniokwadratowy (MSE) - 'Mozaika' Złego Uwarunkowania", fontsize=14)
    plt.xlabel("Stopień wielomianu (m)")
    plt.ylabel("Liczba węzłów (n)")
    plt.savefig(os.path.join(SAVE_DIR, "heatmap_mse.png"), bbox_inches='tight', dpi=150)
    plt.close()
    
    # Mapa Ciepła MAX ERROR
    plt.figure(figsize=(12, 8))
    sns.heatmap(max_results, xticklabels=m_values, yticklabels=n_values, 
                annot=False, cmap="plasma", 
                norm=plt.matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e8))
    plt.title("Błąd Maksymalny (E_max) - Aproksymacja Wielomianowa", fontsize=14)
    plt.xlabel("Stopień wielomianu (m)")
    plt.ylabel("Liczba węzłów (n)")
    plt.savefig(os.path.join(SAVE_DIR, "heatmap_max.png"), bbox_inches='tight', dpi=150)
    plt.close()

    # Zapis wybranych danych do CSV, dla dowodu w sprawozdaniu
    df = pd.DataFrame(table_data)
    df['MSE'] = df['MSE'].map('{:.3e}'.format)
    df['E_max'] = df['E_max'].map('{:.3e}'.format)
    df.to_csv(os.path.join(SAVE_DIR, "tabela_bledow.csv"), index=False)

generate_heatmaps_and_table()


# --- DODATEK: Rysowanie wariantu optymalnego ---
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
        plt.plot(X_DENSE, y_approx, 'b-', label=f'Najlepsze Dopasowanie (m={m})')
        plt.title(f"Aproksymacja optymalna: n={n}, m={m} (MSE = {mse:.3f})")
    
    plt.legend()
    plt.ylim(-10, 40)
    plt.savefig(os.path.join(SAVE_DIR, "najlepsze_dopasowanie.png"), dpi=150)
    plt.close()

best_fit_plot()

print("Gotowe! Wszystkie pliki zostały zapisane w folderze:", SAVE_DIR)