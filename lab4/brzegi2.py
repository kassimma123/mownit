import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

# 1. Definicja funkcji bazowej i jej pochodnych
def f(x):
    return x**2 - 5 * np.cos(2 * np.pi * x)

def df(x):
    return 2*x + 10 * np.pi * np.sin(2 * np.pi * x)

def d2f(x):
    return 2 + 20 * np.pi**2 * np.cos(2 * np.pi * x)

# 2. Własna implementacja splajnu 2. stopnia z warunkami brzegowymi
def fit_quadratic_spline(xs, ys, bc_type='free', d2f_a=None):
    n = len(xs)
    h = np.diff(xs) 
    
    A = np.zeros((n, n))
    B = np.zeros(n)
    
    if bc_type == 'free':
        A[0, 0] = 1.0
        B[0] = 0.0
    elif bc_type == 'clamped':
        A[0, 0] = -1.0
        A[0, 1] = 1.0
        B[0] = d2f_a * h[0]  # POPRAWKA: Usunięto 0.5 wg poprawnego wzoru!
        
    for i in range(1, n):
        A[i, i-1] = 1.0
        A[i, i] = 1.0
        B[i] = 2 * (ys[i] - ys[i-1]) / h[i-1]
        
    b_coeffs = np.linalg.solve(A, B)
    
    def eval_spline(x_eval):
        idx = np.searchsorted(xs, x_eval) - 1
        idx = np.clip(idx, 0, n - 2) 
        dx = x_eval - xs[idx]
        a_coeffs = (b_coeffs[1:] - b_coeffs[:-1]) / (2 * h)
        c_coeffs = ys[:-1]
        return a_coeffs[idx] * dx**2 + b_coeffs[idx] * dx + c_coeffs[idx]
        
    return eval_spline

# 3. Główna funkcja eksperymentu
def generuj_i_zapisz_splajny(n_nodes):
    a, b = -5, 5
    x_fine = np.linspace(a, b, 1000)
    y_true = f(x_fine)
    
    x_w = np.linspace(a, b, n_nodes)
    y_w = f(x_w)
    
    s3_nat = CubicSpline(x_w, y_w, bc_type='natural')
    y_s3_nat = s3_nat(x_fine)
    
    s3_cla = CubicSpline(x_w, y_w, bc_type=((1, df(a)), (1, df(b))))
    y_s3_cla = s3_cla(x_fine)
    
    s2_free = fit_quadratic_spline(x_w, y_w, bc_type='free')
    y_s2_free = s2_free(x_fine)
    
    s2_cla = fit_quadratic_spline(x_w, y_w, bc_type='clamped', d2f_a=d2f(a))
    y_s2_cla = s2_cla(x_fine)

    # POPRAWKA: Dodano np.sqrt() by faktycznie liczyć RMSE, a nie MSE!
    wyniki = {
        'n': n_nodes,
        'S3_Nat_Max': np.max(np.abs(y_true - y_s3_nat)),
        'S3_Nat_RMSE': np.sqrt(np.mean((y_true - y_s3_nat)**2)),
        'S3_Cla_Max': np.max(np.abs(y_true - y_s3_cla)),
        'S3_Cla_RMSE': np.sqrt(np.mean((y_true - y_s3_cla)**2)),
        'S2_Free_Max': np.max(np.abs(y_true - y_s2_free)),
        'S2_Free_RMSE': np.sqrt(np.mean((y_true - y_s2_free)**2)),
        'S2_Cla_Max': np.max(np.abs(y_true - y_s2_cla)),
        'S2_Cla_RMSE': np.sqrt(np.mean((y_true - y_s2_cla)**2))
    }
    
    return wyniki

# --- URUCHOMIENIE ---
nodes_list = [5, 7, 8, 10, 12, 14, 15, 18, 25, 30, 60, 100]
tabela_wynikow = []

for n in nodes_list:
    tabela_wynikow.append(generuj_i_zapisz_splajny(n))

df_wyniki = pd.DataFrame(tabela_wynikow)
pd.options.display.float_format = '{:.2e}'.format
print("\n--- NOWA TABELA BŁĘDÓW (Z PRAWDZIWYM RMSE) ---")
print(df_wyniki.to_string(index=False))