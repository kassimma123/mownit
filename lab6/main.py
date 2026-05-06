import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from matplotlib.colors import ListedColormap

# Wywalamy ostrzeżenia o zerze
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Nasza funkcja - x^10 - (1-x)^15. 
# Matematycznie: to jest wielomian 15 stopnia, będzie się działo przy pochodnych.
def f(x):
    return x**10 - (1 - x)**15

# Pochodna f'(x) = 10x^9 + 15(1-x)^14. 
# Potrzebna do Newtona, żeby wiedzieć w którą stronę i jak szybko "spadać" do zera.
def df(x):
    return 10 * x**9 + 15 * (1 - x)**14

# Prawdziwy pierwiastek (wyliczony wcześniej z wolframa)
TRUE_ROOT = 0.4301597090019467340886000
INTERVAL_A = -1.2
INTERVAL_B = 0.8

# Rho to nasze upragnione przybliżenie (epsilon/dokładność). Testujemy różne rzędy wielkości.
RHOS = [10**(-i) for i in range(2, 10)]
# Punkty startowe - sprawdzamy czy metoda Newtona ruszy z każdego miejsca
START_POINTS = np.round(np.linspace(INTERVAL_A, INTERVAL_B, 21), 2)
MAX_ITER = 100

output_dir = "output_results"
plot_dir = "wykresy_analiza"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# Statusy, żeby na mapie ciepła było widać co poszło nie tak
STATUS_SUCCESS = 0          # Wszystko pykło
STATUS_FALSE_CONV_INC = 1   # Niby się nie rusza (x_new - x bliskie 0), ale jesteśmy daleko od pierwiastka! (pułapka przyrostowa)
STATUS_STUCK_RES = 2        # Niby f(x) bliskie 0, ale to tylko płaski teren, a nie pierwiastek (pułapka rezydualna)
STATUS_THROWN_OUT = 3       # Metoda wystrzeliła nas w kosmos poza zakres
STATUS_DIV_ZERO = 4         # Pochodna = 0, czyli styczna jest pozioma i nie przetnie osi X.

# Metoda Newtona (Metoda Stycznych)
# Matematycznie: x_{n+1} = x_n - f(x_n)/f'(x_n)
# Przybliżamy funkcję linią styczną i patrzymy gdzie ta linia przecina zero.
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
                    
                    # Jak pochodna jest prawie zerem, to Newton głupieje (dzielenie przez zero)
                    if abs(dfx) < 1e-12:
                        status = STATUS_THROWN_OUT
                        iters = i; break
                    
                    x_new = x - f(x) / dfx
                    
                    # Jak nas wywali za daleko, to nie ma sensu dalej mielić
                    if x_new < INTERVAL_A - 10 or x_new > INTERVAL_B + 10:
                        status = STATUS_THROWN_OUT
                        iters = i; break
                    
                    # Kryterium przyrostowe: patrzymy czy kolejne X-sy są już blisko siebie
                    if criterion == 'inc' and abs(x_new - x) <= rho:
                        # Sprawdzamy czy to prawdziwy sukces czy "fałszywy alarm"
                        if abs(x_new - TRUE_ROOT) > 0.05: 
                            status = STATUS_FALSE_CONV_INC
                        iters = i; x = x_new; break
                    
                    # Kryterium rezydualne: patrzymy czy wartość funkcji jest już bliska 0
                    elif criterion == 'res' and abs(f(x_new)) <= rho:
                        if abs(x_new - TRUE_ROOT) > 0.05: 
                            status = STATUS_STUCK_RES
                        iters = i; x = x_new; break

                    x = x_new
                
                results.append({
                    'start_point': pt, 'rho': rho, 'criterion': criterion,
                    'iterations': iters, 'error': abs(x - TRUE_ROOT), 'status': status
                })
    return pd.DataFrame(results)

# Metoda Siecznych
# Matematycznie: zamiast pochodnej używamy ilorazu różnicowego z dwóch punktów.
# x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
def run_secant_experiments():
    results = []
    # Tu robimy trik: jeden punkt sztywny (brzeg przedziału), a drugi "lata"
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
                        # Jak mianownik zerowy, to sieczna jest pozioma.
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
                        
                        # Przesuwamy punkty: x0 bierze miejsce x1, x1 bierze nowy wynik
                        x0, x1 = x1, x_new
                    
                    results.append({
                        'mode': mode_name, 'var_pt': var_pt, 'rho': rho, 'criterion': criterion,
                        'iterations': iters, 'error': abs(x1 - TRUE_ROOT), 'status': status
                    })
    return pd.DataFrame(results)

# Odpalamy te wszystkie obliczenia
print("Liczymy Newtona i Sieczne, chwilę to zajmie...")
df_newton = run_newton_experiments()
df_secant = run_secant_experiments()

# Zapisujemy do CSV
df_newton.to_csv(f'{output_dir}/Newton_results.csv', index=False)
df_secant.to_csv(f'{output_dir}/Secant_results.csv', index=False)

# Rysowanie wykresów
print("Generujemy wykresy...")
sns.set_theme(style="whitegrid")

# Wykres funkcji 
plt.figure(figsize=(10, 5))
x_vals = np.linspace(INTERVAL_A, INTERVAL_B, 500)
plt.plot(x_vals, f(x_vals), label='nasza funkcja f(x)', color='blue', lw=2)
plt.axhline(0, color='black', lw=1)
plt.axvline(TRUE_ROOT, color='red', linestyle='--', label=f'Tu jest zero ({TRUE_ROOT:.4f})')
plt.title("Podgląd funkcji - tu szukamy zera", fontsize=14)
plt.legend()
plt.tight_layout()
plt.savefig(f'{plot_dir}/01_Funkcja.png', dpi=150)
plt.close()

# Funkcja pomocnicza do rysowania błędów w skali logarytmicznej
# Log-log to podstawa, żeby zobaczyć rząd zbieżności
def plot_lines(df_data, group_col, title_prefix, filename):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for i, crit in enumerate(['inc', 'res']):
        ax = axes[i]
        d_crit = df_data[df_data['criterion'] == crit]
        for pt in [-1.2, 0.0, 0.8]: # Wybieramy kilka punktów startowych do porównania
            d_pt = d_crit[np.isclose(d_crit[group_col], pt)]
            if not d_pt.empty:
                ax.plot(d_pt['rho'], d_pt['error'], marker='o', lw=2, label=f'start x0 = {pt}')
        
        ax.set_xscale('log'); ax.set_yscale('log'); ax.invert_xaxis()
        ax.set_xlabel('Zadana dokładność (rho)', fontsize=12)
        ax.set_ylabel('Błąd (jak daleko od celu)', fontsize=12)
        ax.set_title(f'{title_prefix} - Kryterium {"Przyrostowe" if crit=="inc" else "Rezydualne"}', fontsize=14)
        ax.legend()
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/{filename}.png', dpi=150)
    plt.close()

plot_lines(df_newton, 'start_point', 'Błąd Newtona', '02_Newton_Errors')
plot_lines(df_secant[df_secant['mode']=="Usztywnione A (-1.2)"], 'var_pt', 'Błąd Siecznych', '03_Secant_Errors_A')

# Mapy statusów - najciekawsza część sprawozdania
status_text_map = {
    0: "OK", 
    1: "FAŁSZ\n(INC)", # Przyrostowe oszukuje!
    2: "STUCK\n(RES)", # Rezydualne ugrzęzło!
    3: "OUT",          # Wystrzeliło w kosmos
    4: "DIV/0"         # Zero w mianowniku
}

cmap_alerts = ListedColormap(['#a8e6cf', '#ffd3b6', '#ffaaa5', '#ff8b94', '#dcedc1'])

def plot_annotated_heatmap(df_data, y_col, title, filename):
    fig, axes = plt.subplots(1, 2, figsize=(22, 12)) 
    
    for i, crit in enumerate(['inc', 'res']):
        d_crit = df_data[df_data['criterion'] == crit]
        pivot = d_crit.pivot(index=y_col, columns='rho', values='status')
        annot_matrix = pivot.map(lambda v: status_text_map.get(v, "??"))
        
        ax = sns.heatmap(pivot, cmap=cmap_alerts, vmin=0, vmax=4, ax=axes[i], 
                         annot=annot_matrix, fmt="", cbar=False, 
                         linewidths=1, linecolor='black', 
                         annot_kws={"size": 10, "weight": "bold", "color": "black"})
        
        axes[i].set_title(f'Kryterium: {"Przyrostowe" if crit=="inc" else "Rezydualne"}', fontsize=16, pad=15)
        axes[i].invert_xaxis()
        axes[i].tick_params(axis='y', rotation=0, labelsize=12)
        axes[i].set_xlabel('Rho (dokładność)', fontsize=14)
        axes[i].set_ylabel('Punkt startowy', fontsize=14)
    
    plt.suptitle(title, fontsize=22, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{plot_dir}/{filename}.png', dpi=200)
    plt.close()

# Generujemy te wielkie tabele, z których wszystko widać jak na dłoni
plot_annotated_heatmap(df_newton, 'start_point', 'Gdzie Newton daje radę, a gdzie pada?', '04_Newton_Tabela_Zdarzen')
plot_annotated_heatmap(df_secant[df_secant['mode']=="Usztywnione A (-1.2)"], 'var_pt', 'Sieczne - jak punkty startowe wpływają na wynik?', '05_Secant_Tabela_Zdarzen_A')

print(f"\nRobota skończona. Wykresy są w folderze '{plot_dir}'.")
