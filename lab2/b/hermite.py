import numpy as np
import matplotlib.pyplot as plt

#  1. FUNKCJA I POCHODNA 
def funkcja_f10(x):
    m = 5; k = 0.5
    return x**2 - m * np.cos((np.pi * x) / k)

def pochodna_f10(x):
    m = 5; k = 0.5
    return 2*x + m * (np.pi / k) * np.sin((np.pi * x) / k)

#  2. GENEROWANIE WĘZŁÓW 
def wezly_rownomierne(a, b, n):
    return np.linspace(a, b, n)

def wezly_czebyszewa(a, b, n):
    k = np.arange(n)
    return 0.5 * (a + b) + 0.5 * (b - a) * np.cos((2 * k + 1) * np.pi / (2 * n))

#  3. INTERPOLACJA HERMITE'A 
def ilorazy_roznicowe_hermite(x_wezly, y_wezly, yp_wezly):
    """Tworzy tabelę ilorazów różnicowych dla węzłów podwójnych."""
    n = len(x_wezly)
    z = np.zeros(2*n)
    Q = np.zeros((2*n, 2*n))
    
    for i in range(n):
        z[2*i] = x_wezly[i]
        z[2*i+1] = x_wezly[i]
        
        Q[2*i][0] = y_wezly[i]
        Q[2*i+1][0] = y_wezly[i]
        Q[2*i+1][1] = yp_wezly[i] # Pochodna w miejsce pierwszego ilorazu
        
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

#  4. POMIARY 
def badaj_hermite(n_lista):
    a, b = -5, 5
    x_geste = np.linspace(a, b, 1000)
    y_prawdziwe = funkcja_f10(x_geste)
    
    print(f"{'n':<5} | {'Błąd Równomierne':<20} | {'Błąd Czebyszewa':<20}")
    print("-" * 52)
    
    for n in n_lista:
        # Pomiary dla węzłów równomiernych
        x_row = wezly_rownomierne(a, b, n)
        y_herm_row = interpolacja_hermite(x_row, funkcja_f10(x_row), pochodna_f10(x_row), x_geste)
        blad_row = np.max(np.abs(y_prawdziwe - y_herm_row))
        
        # Pomiary dla węzłów Czebyszewa
        x_czeb = wezly_czebyszewa(a, b, n)
        y_herm_czeb = interpolacja_hermite(x_czeb, funkcja_f10(x_czeb), pochodna_f10(x_czeb), x_geste)
        blad_czeb = np.max(np.abs(y_prawdziwe - y_herm_czeb))
        
        print(f"{n:<5} | {blad_row:<20.2e} | {blad_czeb:<20.2e}")
        
        # Zapisywanie wykresu dla wybranych n 
        if n in [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]:
            plt.figure(figsize=(10,5))
            plt.plot(x_geste, y_prawdziwe, 'k--', label='Funkcja oryginalna')
            plt.plot(x_geste, y_herm_row, 'b-', label='Hermite (równomierne)')
            plt.plot(x_row, funkcja_f10(x_row), 'ro')
            plt.ylim(min(y_prawdziwe)-5, max(y_prawdziwe)+5)
            plt.title(f"Hermite, n={n}, Węzły równomierne, Błąd: {blad_row:.2e}")
            plt.legend()
            plt.savefig(f"hermite_row_{n}.png")
            plt.close()

n_do_sprawdzenia = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30]
badaj_hermite(n_do_sprawdzenia)