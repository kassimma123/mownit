import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d, CubicSpline

# 1. Definicja funkcji i jej pochodnej
def f(x):
    return x**2 - 5 * np.cos(2 * np.pi * x)

def df(x):
    return 2*x + 10 * np.pi * np.sin(2 * np.pi * x)

def generuj_i_zapisz_splajny(n_nodes):
    a, b = -5, 5
    x_fine = np.linspace(a, b, 1000)
    y_true = f(x_fine)
    
    x_w = np.linspace(a, b, n_nodes)
    y_w = f(x_w)
    
    # --- OBLICZENIA ---
    # Stopień 3: Naturalny (2-ga pochodna na końcach = 0)
    s3_nat = CubicSpline(x_w, y_w, bc_type='natural')
    
    # Stopień 3: Ustalony (pochodne na końcach zgodne z funkcją)
    s3_cla = CubicSpline(x_w, y_w, bc_type=((1, df(a)), (1, df(b))))
    
    # Stopień 2: Standardowy (kwadratowy)
    s2_std = interp1d(x_w, y_w, kind='quadratic', fill_value="extrapolate")
    
    # Stopień 2: Inny wariant (np. slinear dla porównania lub inna metoda)
    # W scipy interp1d 'quadratic' ma jedno podejście, ale możemy zasymulować
    # różnicę używając splajnu 2 stopnia przez CubicSpline z redukcją stopnia (teoretycznie)
    # Dla potrzeb zadania skupmy się na tych trzech kluczowych różnicach.

    # --- WIZUALIZACJA ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    fig.suptitle(f'Porównanie metod interpolacji: n = {n_nodes}', fontsize=16)

    # Definicja danych do pętli wykresów
    wykresy = [
        (s3_nat(x_fine), 'S3: Naturalny', 'royalblue', '--'),
        (s3_cla(x_fine), 'S3: Ustalony (Clamped)', 'crimson', '-'),
        (s2_std(x_fine), 'S2: Kwadratowy', 'forestgreen', '-.')
    ]

    for ax, (y_spline, tytul, kolor, styl) in zip(axes, wykresy):
        ax.plot(x_fine, y_true, color='gray', alpha=0.3, label='Oryginał f(x)')
        ax.plot(x_fine, y_spline, color=kolor, linestyle=styl, label=tytul, lw=2)
        ax.scatter(x_w, y_w, color='black', s=10, zorder=5, label='Węzły')
        
        ax.set_title(tytul)
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend(loc='lower center', fontsize='small')
        
        # Obliczanie błędu maksymalnego dla opisu
        err = np.max(np.abs(y_true - y_spline))
        ax.text(0.05, 0.95, f'Max Err: {err:.2e}', transform=ax.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # --- ZAPIS DO PLIKU ---
    filename = f"porownanie_n{n_nodes}.png"
    plt.savefig(filename, dpi=150)
    print(f"Zapisano wykres: {filename}")
    plt.show()

# Lista n do testów
nodes_list = [3, 5, 6, 7, 8, 9]

for n in nodes_list:
    generuj_i_zapisz_splajny(n)