import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 1. Definicja badanej funkcji f10(x)
def funkcja_f10(x):
    m = 5
    k = 0.5
    return x**2 - m * np.cos((np.pi * x) / k)

# 2. Generowanie węzłów równomiernych
def wezly_rownomierne(a, b, n):
    return np.linspace(a, b, n)

# 3. Generowanie węzłów Czebyszewa
def wezly_czebyszewa(a, b, n):
    k = np.arange(n)
    # przeliczanie zera Czebyszewa na przedział [a, b]
    punkty = 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k + 1) * np.pi / (2 * n))
    return punkty

def interpolacja_lagrangea(x_wezly, y_wezly, x_punkty):
    n = len(x_wezly)
    y_wynik = np.zeros_like(x_punkty, dtype=float)

    for i in range(n):
        L_i = np.ones_like(x_punkty, dtype=float)

        for j in range(n):
            if i != j:
                L_i *= (x_punkty - x_wezly[j]) / (x_wezly[i] - x_wezly[j])
        
        y_wynik += y_wezly[i] * L_i
    return y_wynik

def ilorazy_roznicowe(x_wezly, y_wezly):
    n = len(x_wezly)
    tabela = np.zeros([n, n])
    tabela[:, 0] = y_wezly
    for j in range(1, n):
        for i in range(n-j):
            tabela[i][j] = (tabela[i+1][j-1] - tabela[i][j-1]) / (x_wezly[i+j] - x_wezly[i])

    return tabela[0, :]

def interpolacja_newtona(x_wezly, y_wezly, x_punkty):
    wspolczynniki = ilorazy_roznicowe(x_wezly, y_wezly)
    n = len(x_wezly)

    y_wynik = np.full_like(x_punkty, wspolczynniki[0], dtype=float)

    for i in range(1, n):
        kolejna_warstwa = np.full_like(x_punkty, wspolczynniki[i], dtype=float)

        for j in range(i):
            kolejna_warstwa *= (x_punkty - x_wezly[j])

        y_wynik += kolejna_warstwa
    return y_wynik

def generuj_wyniki_do_sprawozdania(n_lista):
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)
    
    wyniki = []

    for n in n_lista:
        # 1. Węzły
        xr = wezly_rownomierne(a, b, n)
        xc = wezly_czebyszewa(a, b, n)
        
        # 2. Wielomiany
        y_lag_r = interpolacja_lagrangea(xr, funkcja_f10(xr), x_geste)
        y_lag_c = interpolacja_lagrangea(xc, funkcja_f10(xc), x_geste)
        y_newt_r = interpolacja_newtona(xr, funkcja_f10(xr), x_geste)
        y_newt_c = interpolacja_newtona(xc, funkcja_f10(xc), x_geste) 
        
        # 3. Pomiary Błędów Równomierne (Max i MSE)
        err_max_r = np.max(np.abs(y_prawdziwe - y_lag_r))
        err_mse_r = np.mean((y_prawdziwe - y_lag_r)**2)
        
        # 4. Pomiary Błędów Czebyszew (Max i MSE)
        err_max_c = np.max(np.abs(y_prawdziwe - y_lag_c))
        err_mse_c = np.mean((y_prawdziwe - y_lag_c)**2)
        
        # 5. Różnica metody Lagrange vs Newton (dla obu rodzajów węzłów)
        roznica_metod_r = np.max(np.abs(y_lag_r - y_newt_r))
        roznica_metod_c = np.max(np.abs(y_lag_c - y_newt_c)) 
        
        # 6. Dodanie do tabeli
        wyniki.append({
            "n": n,
            "Max Błąd (Równ)": f"{err_max_r:.2e}",
            "MSE (Równ)": f"{err_mse_r:.2e}",
            "Max Błąd (Czeb)": f"{err_max_c:.2e}",
            "Różnica Lag-Newt (Równ)": f"{roznica_metod_r:.2e}",
            "Różnica Lag-Newt (Czeb)": f"{roznica_metod_c:.2e}" 
        })
        
        # 7. Zapisywanie wykresów tylko dla ciekawych 'n' 
        if n in [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 45]:
            plt.figure(figsize=(8, 4))
            plt.plot(x_geste, y_prawdziwe, 'k--', label='Funkcja f10(x)')
            plt.plot(x_geste, y_lag_r, 'b-', label=f'Lagrange (Równomierne n={n})')
            plt.plot(x_geste, y_lag_c, 'g-', label=f'Lagrange (Czebyszew n={n})')
            plt.plot(xr, funkcja_f10(xr), 'bo', label='Węzły Równ.')
            
            # Ogranicznik
            plt.ylim(min(y_prawdziwe) - 10, max(y_prawdziwe) + 15)
            
            plt.title(f"Interpolacja n={n} | Max błąd (Równ): {err_max_r:.1f}")
            plt.legend()
            plt.grid(True)
            plt.savefig(f"wykres_interpolacja_n{n}.png", bbox_inches='tight')
            plt.close()

    df = pd.DataFrame(wyniki)
    print(df.to_string(index=False))

# ODPALENIE TESTÓW
n_do_testow = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 45]
generuj_wyniki_do_sprawozdania(n_do_testow)