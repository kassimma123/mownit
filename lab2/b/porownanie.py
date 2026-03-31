import numpy as np
import matplotlib.pyplot as plt

#  1. Funkcja, pochodna, węzły
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

#  2. ALGORYTMY INTERPOLACJI 
#  NEWTON 
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

#  HERMITE 
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

#  3. RYSOWANIE PORÓWNANIA 
def porownaj_wykresy(n):
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)
    
    # Przygotowanie okna z dwoma wykresami (lewy i prawy)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    #  WYKRES 1: Węzły Równomierne 
    x_row = wezly_rownomierne(a, b, n)
    y_row = funkcja_f10(x_row)
    yp_row = pochodna_f10(x_row)
    
    y_newt_row = interpolacja_newtona(x_row, y_row, x_geste)
    y_herm_row = interpolacja_hermite(x_row, y_row, yp_row, x_geste)
    
    ax1.plot(x_geste, y_prawdziwe, 'k--', linewidth=2, label='Funkcja oryginalna')
    ax1.plot(x_geste, y_newt_row, 'b-', label=f'Newton (Stopień {n-1})')
    ax1.plot(x_geste, y_herm_row, 'r-', label=f'Hermite (Stopień {2*n-1})')
    ax1.plot(x_row, y_row, 'ko', markersize=6, label='Węzły')
    
    ax1.set_title(f"Węzły Równomierne (n={n})")
    ax1.set_ylim(min(y_prawdziwe) - 10, max(y_prawdziwe) + 20)
    ax1.grid(True)
    ax1.legend()

    #  WYKRES 2: Węzły Czebyszewa 
    x_czeb = wezly_czebyszewa(a, b, n)
    y_czeb = funkcja_f10(x_czeb)
    yp_czeb = pochodna_f10(x_czeb)
    
    y_newt_czeb = interpolacja_newtona(x_czeb, y_czeb, x_geste)
    y_herm_czeb = interpolacja_hermite(x_czeb, y_czeb, yp_czeb, x_geste)
    
    ax2.plot(x_geste, y_prawdziwe, 'k--', linewidth=2, label='Funkcja oryginalna')
    ax2.plot(x_geste, y_newt_czeb, 'b-', label=f'Newton (Stopień {n-1})')
    ax2.plot(x_geste, y_herm_czeb, 'r-', label=f'Hermite (Stopień {2*n-1})')
    ax2.plot(x_czeb, y_czeb, 'ko', markersize=6, label='Węzły')
    
    ax2.set_title(f"Węzły Czebyszewa (n={n})")
    ax2.set_ylim(min(y_prawdziwe) - 10, max(y_prawdziwe) + 20)
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.show()

porownaj_wykresy(n=20)