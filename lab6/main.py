import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from matplotlib.colors import ListedColormap

# Wyłączenie ostrzeżeń o dzieleniu przez zero w numpy
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================
# 1. DEFINICJE FUNKCJI I USTAWIENIA
# ==========================================
def f(x):
    return x**10 - (1 - x)**15

def df(x):
    return 10 * x**9 + 15 * (1 - x)**14

TRUE_ROOT = 0.4301597090019467340886000
INTERVAL_A = -1.2
INTERVAL_B = 0.8

RHOS = [10**(-i) for i in range(2, 10)]
START_POINTS = np.round(np.linspace(INTERVAL_A, INTERVAL_B, 21), 2)
MAX_ITER = 100

output_dir = "output_results"
plot_dir = "wykresy_analiza"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# Kody statusów dla mapy ciepła:
STATUS_SUCCESS = 0
STATUS_FALSE_CONV_INC = 1
STATUS_STUCK_RES = 2
STATUS_THROWN_OUT = 3
STATUS_DIV_ZERO = 4

# ==========================================
# 2. IMPLEMENTACJA METOD
# ==========================================
def run_newton_experiments():
    results = []
    for pt in START_POINTS:
        for rho in RHOS:
            for criterion in ['inc', 'res']:
                x = pt
                status = STATUS_SUCCESS
                iters = MAX_ITER
                
                for i in range(1, MAX_ITER + 1):
                    dfx = df(x)
                    
                    if abs(dfx) < 1e-12:
                        status = STATUS_THROWN_OUT
                        iters = i; break
                    
                    x_new = x - f(x) / dfx
                    
                    if x_new < INTERVAL_A - 10 or x_new > INTERVAL_B + 10:
                        status = STATUS_THROWN_OUT
                        iters = i; break
                    
                    if criterion == 'inc' and abs(x_new - x) <= rho:
                        if abs(x_new - TRUE_ROOT) > 0.05: # Fałszywa zbieżność (daleko od celu)
                            status = STATUS_FALSE_CONV_INC
                        iters = i; x = x_new; break
                    
                    elif criterion == 'res' and abs(f(x_new)) <= rho:
                        if abs(x_new - TRUE_ROOT) > 0.05: # Ugrzęźnięcie w płaskim
                            status = STATUS_STUCK_RES
                        iters = i; x = x_new; break

                    x = x_new
                
                results.append({
                    'start_point': pt, 'rho': rho, 'criterion': criterion,
                    'iterations': iters, 'error': abs(x - TRUE_ROOT), 'status': status
                })
    return pd.DataFrame(results)

def run_secant_experiments():
    results = []
    for fixed_pt, var_points, mode_name in [
        (INTERVAL_A, START_POINTS[1:], "Usztywnione A (-1.2)"),
        (INTERVAL_B, START_POINTS[:-1], "Usztywnione B (0.8)")
    ]:
        for var_pt in var_points:
            for rho in RHOS:
                for criterion in ['inc', 'res']:
                    x0, x1 = fixed_pt, var_pt
                    status = STATUS_SUCCESS
                    iters = MAX_ITER
                    
                    for i in range(1, MAX_ITER + 1):
                        f0, f1 = f(x0), f(x1)
                        if abs(f1 - f0) < 1e-15:
                            status = STATUS_DIV_ZERO
                            iters = i; break
                            
                        x_new = x1 - f1 * (x1 - x0) / (f1 - f0)
                        
                        if x_new < INTERVAL_A - 10 or x_new > INTERVAL_B + 10:
                            status = STATUS_THROWN_OUT
                            iters = i; break
                        
                        if criterion == 'inc' and abs(x_new - x1) <= rho:
                            if abs(x_new - TRUE_ROOT) > 0.05:
                                status = STATUS_FALSE_CONV_INC
                            x1 = x_new; iters = i; break
                            
                        elif criterion == 'res' and abs(f(x_new)) <= rho:
                            if abs(x_new - TRUE_ROOT) > 0.05:
                                status = STATUS_STUCK_RES
                            x1 = x_new; iters = i; break
                        
                        x0, x1 = x1, x_new
                    
                    results.append({
                        'mode': mode_name, 'var_pt': var_pt, 'rho': rho, 'criterion': criterion,
                        'iterations': iters, 'error': abs(x1 - TRUE_ROOT), 'status': status
                    })
    return pd.DataFrame(results)

# ==========================================
# 3. GENEROWANIE DANYCH
# ==========================================
print("Obliczanie danych (Newton i Sieczne)...")
df_newton = run_newton_experiments()
df_secant = run_secant_experiments()

df_newton.to_csv(f'{output_dir}/Newton_results.csv', index=False)
df_secant.to_csv(f'{output_dir}/Secant_results.csv', index=False)
# ==========================================
# 4. WIZUALIZACJE (POTĘŻNE I CZYTELNE)
# ==========================================
print("Generowanie potężnych, czytelnych wykresów...")
sns.set_theme(style="whitegrid")

# --- WIZUALIZACJA 1: Wykres badanej funkcji ---
plt.figure(figsize=(10, 5))
x_vals = np.linspace(INTERVAL_A, INTERVAL_B, 500)
plt.plot(x_vals, f(x_vals), label='f(x) = x^10 - (1-x)^15', color='blue', lw=2)
plt.axhline(0, color='black', lw=1)
plt.axvline(TRUE_ROOT, color='red', linestyle='--', label=f'True Root ({TRUE_ROOT:.4f})')
plt.title("Kształt funkcji i miejsce zerowe", fontsize=14)
plt.legend()
plt.tight_layout()
plt.savefig(f'{plot_dir}/01_Funkcja.png', dpi=150)
plt.close()

selected_pts = [-1.2, 0.0, 0.8]

def plot_lines(df_data, group_col, title_prefix, filename):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for i, crit in enumerate(['inc', 'res']):
        ax = axes[i]
        d_crit = df_data[df_data['criterion'] == crit]
        for pt in selected_pts:
            d_pt = d_crit[np.isclose(d_crit[group_col], pt)]
            if not d_pt.empty:
                ax.plot(d_pt['rho'], d_pt['error'], marker='o', lw=2, label=f'x0 = {pt}')
        
        ax.set_xscale('log'); ax.set_yscale('log'); ax.invert_xaxis()
        ax.set_xlabel('Dokładność (rho)', fontsize=12)
        ax.set_ylabel('Błąd bezwzględny', fontsize=12)
        ax.set_title(f'{title_prefix} - Kryterium {"Przyrostowe" if crit=="inc" else "Rezydualne"}', fontsize=14)
        ax.legend()
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/{filename}.png', dpi=150)
    plt.close()

# Wykresy błędów
plot_lines(df_newton, 'start_point', 'Błąd Newton', '02_Newton_Errors')
plot_lines(df_secant[df_secant['mode']=="Usztywnione A (-1.2)"], 'var_pt', 'Błąd Sieczne (Sztywne A)', '03_Secant_Errors_A')

# --- WIZUALIZACJA 4, 5, 6: CZYTELNE MAPY ALERTÓW Z TEKSTEM ---
# Słowniki tłumaczące kody liczbowe na krótki, dosadny tekst w komórkach
status_text_map = {
    0: "OK", 
    1: "Fałsz\n(Inc)", 
    2: "Ugrzązł\n(Res)", 
    3: "Wyrzucony", 
    4: "Dziel/0"
}

# Paleta kolorów dopasowana do czytelności tekstu (jasne tła)
cmap_alerts = ListedColormap(['#a8e6cf', '#ffd3b6', '#ffaaa5', '#ff8b94', '#dcedc1'])

def plot_annotated_heatmap(df_data, y_col, title, filename):
    fig, axes = plt.subplots(1, 2, figsize=(22, 12)) # Ogromny rozmiar dla czytelności
    
    for i, crit in enumerate(['inc', 'res']):
        d_crit = df_data[df_data['criterion'] == crit]
        pivot = d_crit.pivot(index=y_col, columns='rho', values='status')
        
        # Tworzenie macierzy z tekstami zamiast liczb
        annot_matrix = pivot.map(lambda v: status_text_map.get(v, "BŁĄD"))
        
        # Rysowanie z parametrem annot=annot_matrix
        ax = sns.heatmap(pivot, cmap=cmap_alerts, vmin=0, vmax=4, ax=axes[i], 
                         annot=annot_matrix, fmt="", cbar=False, 
                         linewidths=1, linecolor='black', 
                         annot_kws={"size": 10, "weight": "bold", "color": "black"})
        
        axes[i].set_title(f'Kryterium: {"Przyrostowe" if crit=="inc" else "Rezydualne"}', fontsize=16, pad=15)
        axes[i].invert_xaxis() # Rho malejące
        
        # Obrót etykiet na osi Y, by były w poziomie
        axes[i].tick_params(axis='y', rotation=0, labelsize=12)
        axes[i].tick_params(axis='x', labelsize=12)
        axes[i].set_xlabel('Dokładność (rho)', fontsize=14)
        axes[i].set_ylabel('Punkt startowy (x0)', fontsize=14)
    
    plt.suptitle(title, fontsize=22, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Marginesy na tytuł główny
    plt.savefig(f'{plot_dir}/{filename}.png', dpi=200) # Wysokie DPI
    plt.close()

# Generowanie czytelnych tabel-map
plot_annotated_heatmap(df_newton, 'start_point', 'STATUS ALGORYTMU - Metoda Newtona', '04_Newton_Tabela_Zdarzen')
plot_annotated_heatmap(df_secant[df_secant['mode']=="Usztywnione A (-1.2)"], 'var_pt', 'STATUS ALGORYTMU - Sieczne (Usztywnione A = -1.2)', '05_Secant_Tabela_Zdarzen_A')
plot_annotated_heatmap(df_secant[df_secant['mode']=="Usztywnione B (0.8)"], 'var_pt', 'STATUS ALGORYTMU - Sieczne (Usztywnione B = 0.8)', '06_Secant_Tabela_Zdarzen_B')

print(f"\nGotowe! Pliki zapisano w '{plot_dir}'.")
print("Sprawdź pliki 04, 05 i 06 - to teraz wielkie tabele ze statusem wypisanym w każdym okienku.")