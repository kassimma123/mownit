import numpy as np
import matplotlib.pyplot as plt

def funkcja_f10(x):
    m = 5; k = 0.5
    return x**2 - m * np.cos((np.pi * x) / k)

def pochodna_f10(x):
    m = 5; k = 0.5
    return 2*x + m * (np.pi / k) * np.sin((np.pi * x) / k)

def wezly_rownomierne(a, b, n):
    return np.linspace(a, b, n)

def wezly_czebyszewa(a, b, n):
    k = np.arange(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k + 1) * np.pi / (2 * n))

# --- NEWTON ---
def ilorazy_roznicowe_newton(x_wezly, y_wezly):
    n = len(x_wezly)
    tabela = np.zeros([n, n])
    tabela[:, 0] = y_wezly
    for j in range(1, n):
        for i in range(n-j):
            tabela[i][j] = (tabela[i+1][j-1] - tabela[i][j-1]) / (x_wezly[i+j] - x_wezly[i])
    return tabela[0, :]

def interpolacja_newtona(x_wezly, y_wezly, x_punkty):
    wspolczynniki = ilorazy_roznicowe_newton(x_wezly, y_wezly)
    n = len(x_wezly)
    y_wynik = np.full_like(x_punkty, wspolczynniki[0], dtype=float)
    for i in range(1, n):
        kolejna_warstwa = np.full_like(x_punkty, wspolczynniki[i], dtype=float)
        for j in range(i):
            kolejna_warstwa *= (x_punkty - x_wezly[j])
        y_wynik += kolejna_warstwa
    return y_wynik

# --- HERMITE ---
def ilorazy_roznicowe_hermite(x_wezly, y_wezly, yp_wezly):
    n = len(x_wezly)
    z = np.zeros(2*n)
    Q = np.zeros((2*n, 2*n))
    for i in range(n):
        z[2*i] = x_wezly[i]
        z[2*i+1] = x_wezly[i]
        Q[2*i][0] = y_wezly[i]
        Q[2*i+1][0] = y_wezly[i]
        Q[2*i+1][1] = yp_wezly[i] 
        if i != 0:
            Q[2*i][1] = (Q[2*i][0] - Q[2*i-1][0]) / (z[2*i] - z[2*i-1])
    for j in range(2, 2*n):
        for i in range(j, 2*n):
            Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (z[i] - z[i-j])
    return z, np.diag(Q)

def interpolacja_hermite(x_wezly, y_wezly, yp_wezly, x_punkty):
    z, wspolczynniki = ilorazy_roznicowe_hermite(x_wezly, y_wezly, yp_wezly)
    n = len(z)
    y_wynik = np.full_like(x_punkty, wspolczynniki[0], dtype=float)
    iloczyn = np.ones_like(x_punkty, dtype=float)
    for i in range(1, n):
        iloczyn *= (x_punkty - z[i-1])
        y_wynik += wspolczynniki[i] * iloczyn
    return y_wynik

# --- RYSOWANIE WYKRESÓW ---
def generuj_wykresy():
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)

    for n in [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]: 
        # 1. RÓWNOMIERNE 
        xr = wezly_rownomierne(a, b, n)
        yn_r = interpolacja_newtona(xr, funkcja_f10(xr), x_geste)
        yh_r = interpolacja_hermite(xr, funkcja_f10(xr), pochodna_f10(xr), x_geste)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.plot(x_geste, y_prawdziwe, 'k--', label='Oryginał')
        ax1.plot(x_geste, yn_r, 'b-', label='Newton')
        ax1.plot(xr, funkcja_f10(xr), 'ro'); ax1.set_title(f'Newton (Równomierne n={n})')
        ax1.set_ylim(-10, 40); ax1.grid(True); ax1.legend()
        
        ax2.plot(x_geste, y_prawdziwe, 'k--', label='Oryginał')
        ax2.plot(x_geste, yh_r, 'g-', label='Hermite')
        ax2.plot(xr, funkcja_f10(xr), 'ro'); ax2.set_title(f'Hermite (Równomierne n={n})')
        ax2.set_ylim(-10, 40); ax2.grid(True); ax2.legend()
        plt.savefig(f"newt_vs_herm_rown_n{n}.png", bbox_inches='tight')
        plt.close()

        # 2. CZEBYSZEWA 
        xc = wezly_czebyszewa(a, b, n)
        yn_c = interpolacja_newtona(xc, funkcja_f10(xc), x_geste)
        yh_c = interpolacja_hermite(xc, funkcja_f10(xc), pochodna_f10(xc), x_geste)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        ax1.plot(x_geste, y_prawdziwe, 'k--', label='Oryginał')
        ax1.plot(x_geste, yn_c, 'b-', label='Newton')
        ax1.plot(xc, funkcja_f10(xc), 'ro'); ax1.set_title(f'Newton (Czebyszew n={n})')
        ax1.set_ylim(-10, 40); ax1.grid(True); ax1.legend()
        
        ax2.plot(x_geste, y_prawdziwe, 'k--', label='Oryginał')
        ax2.plot(x_geste, yh_c, 'g-', label='Hermite')
        ax2.plot(xc, funkcja_f10(xc), 'ro'); ax2.set_title(f'Hermite (Czebyszew n={n})')
        ax2.set_ylim(-10, 40); ax2.grid(True); ax2.legend()
        plt.savefig(f"newt_vs_herm_czeb_n{n}.png", bbox_inches='tight')
        plt.close()

# --- GENEROWANIE TABELI WYNIKÓW ---
def generuj_tabele():
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)
    
    n_lista = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]
    
    print("GOTOWY KOD DO TYPSTA (skopiuj i wklej wewnątrz `table()`):")
    print("-" * 100)
    
    for n in n_lista:
        # Pomiary dla węzłów równomiernych
        xr = wezly_rownomierne(a, b, n)
        yr = funkcja_f10(xr)
        ypr = pochodna_f10(xr)
        
        yn_r = interpolacja_newtona(xr, yr, x_geste)
        yh_r = interpolacja_hermite(xr, yr, ypr, x_geste)
        
        err_max_n_r = np.max(np.abs(y_prawdziwe - yn_r))
        err_max_h_r = np.max(np.abs(y_prawdziwe - yh_r))
        mse_h_r = np.mean((y_prawdziwe - yh_r)**2)
        
        # Pomiary dla węzłów Czebyszewa
        xc = wezly_czebyszewa(a, b, n)
        yc = funkcja_f10(xc)
        ypc = pochodna_f10(xc)
        
        yn_c = interpolacja_newtona(xc, yc, x_geste)
        yh_c = interpolacja_hermite(xc, yc, ypc, x_geste)
        
        err_max_n_c = np.max(np.abs(y_prawdziwe - yn_c))
        err_max_h_c = np.max(np.abs(y_prawdziwe - yh_c))
        mse_h_c = np.mean((y_prawdziwe - yh_c)**2)
        
        # Wypisywanie z odpowiednim formatowaniem
        print(f"  [{n}],  [{err_max_n_r:.2e}], [{err_max_h_r:.2e}], [{mse_h_r:.2e}], [{err_max_n_c:.2e}], [{err_max_h_c:.2e}], [{mse_h_c:.2e}],")

    print("-" * 100)


# Uruchomienie skryptu
generuj_tabele()
generuj_wykresy()