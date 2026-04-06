import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d, CubicSpline

# 1. Funkcja i węzły
def funkcja_f10(x):
    m = 5; k = 0.5
    return x**2 - m * np.cos((np.pi * x) / k)

def pochodna_f10(x):
    m = 5; k = 0.5
    return 2*x + m * (np.pi / k) * np.sin((np.pi * x) / k)

def badaj_splajny(n_lista):
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)
    
    wyniki = []

    for n in n_lista:
        x_w = np.linspace(a, b, n)
        y_w = funkcja_f10(x_w)
        
        # --- SPLINE 2 STOPNIA ---
        # Używamy interp1d z kind='quadratic'
        spline2 = interp1d(x_w, y_w, kind='quadratic', fill_value="extrapolate")
        y_spline2 = spline2(x_geste)
        
        # --- SPLINE 3 STOPNIA (Warunek Naturalny) ---
        # bc_type='natural' -> druga pochodna na brzegach = 0
        spline3_nat = CubicSpline(x_w, y_w, bc_type='natural')
        y_spline3_nat = spline3_nat(x_geste)
        
        # --- SPLINE 3 STOPNIA (Warunek Ustalony / Clamped) ---
        # bc_type=((1, pochodna_lewa), (1, pochodna_prawa))
        poch_lewa = pochodna_f10(a)
        poch_prawa = pochodna_f10(b)
        spline3_cla = CubicSpline(x_w, y_w, bc_type=((1, poch_lewa), (1, poch_prawa)))
        y_spline3_cla = spline3_cla(x_geste)
        
        # Błędy MSE
        mse_2 = np.mean((y_prawdziwe - y_spline2)**2)
        mse_3_nat = np.mean((y_prawdziwe - y_spline3_nat)**2)
        mse_3_cla = np.mean((y_prawdziwe - y_spline3_cla)**2)
        
        # Błędy Max
        max_2 = np.max(np.abs(y_prawdziwe - y_spline2))
        max_3_nat = np.max(np.abs(y_prawdziwe - y_spline3_nat))
        max_3_cla = np.max(np.abs(y_prawdziwe - y_spline3_cla))
        
        wyniki.append({
            "n": n,
            "Max Błąd (2 st.)": f"{max_2:.2e}",
            "MSE (2 st.)": f"{mse_2:.2e}",
            "Max Błąd (3 st. Natural)": f"{max_3_nat:.2e}",
            "MSE (3 st. Natural)": f"{mse_3_nat:.2e}",
            "Max Błąd (3 st. Clamped)": f"{max_3_cla:.2e}",
            "MSE (3 st. Clamped)": f"{mse_3_cla:.2e}"
        })
        
        # Generowanie wykresów dla ciekawych n
        if n in [10, 20]:
            plt.figure(figsize=(10, 5))
            plt.plot(x_geste, y_prawdziwe, 'k--', label='Oryginał', alpha=0.6)
            plt.plot(x_geste, y_spline2, 'g-', label='Spline 2 st.')
            plt.plot(x_geste, y_spline3_nat, 'b-', label='Spline 3 st. (Natural)')
            plt.plot(x_geste, y_spline3_cla, 'r-.', label='Spline 3 st. (Clamped)')
            plt.plot(x_w, y_w, 'ko', label='Węzły')
            plt.title(f"Funkcje Sklejane n={n}")
            plt.ylim(min(y_prawdziwe) - 5, max(y_prawdziwe) + 15)
            plt.legend()
            plt.grid(True)
            plt.savefig(f"spline_porownanie_n{n}.png", bbox_inches='tight')
            plt.close()

    df = pd.DataFrame(wyniki)
    print(df.to_string(index=False))

n_do_testow = [5, 8, 10, 15, 20, 30, 45]
badaj_splajny(n_do_testow)