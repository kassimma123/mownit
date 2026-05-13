import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import time
import tracemalloc
import random
import gc

# --- 1. ALGORYTMY ROZWIĄZYWANIA UKŁADÓW RÓWNAŃ ---

def gauss(A, b, precision):
    """Metoda eliminacji Gaussa z częściowym wyborem elementu głównego."""
    A = A.astype(precision)
    b = b.astype(precision)
    n = A.shape[0]
    
    for i in range(n):
        # Wybór elementu głównego
        _max_i = np.argmax(np.abs(A[i:, i])) + i
        A[[i, _max_i], :] = A[[_max_i, i], :]
        b[i], b[_max_i] = b[_max_i], b[i]
        
        if A[i, i] == 0:
            continue
            
        b[i] /= A[i, i]
        A[i, :] /= A[i, i]
        
        for j in range(i + 1, n):
            b[j] -= b[i] * A[j, i]
            A[j, :] -= A[i, :] * A[j, i]
            
    # Postępowanie odwrotne (Back substitution)
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            b[j] -= b[i] * A[j, i]
            A[j, i] = 0
            
    return b

def thomas_memory(A_banded, b, precision):
    """
    Zoptymalizowany algorytm Thomasa.
    A_banded to macierz n x 3, gdzie:
    kolumna 0 to poddiagonala, 1 to diagonala, 2 to naddiagonala.
    """
    A = A_banded.astype(precision)
    b = b.astype(precision)
    n = A.shape[0]
    
    # Eliminacja w przód
    if A[0, 1] != 0:
        b[0] /= A[0, 1]
        A[0, 2] /= A[0, 1]
        A[0, 1] = 1.0
    
    for i in range(1, n):
        factor = A[i, 1] - A[i, 0] * A[i-1, 2]
        if factor == 0:
            factor = 1e-15 # zapobieganie dzieleniu przez 0
        if i < n - 1:
            A[i, 2] /= factor
        b[i] = (b[i] - A[i, 0] * b[i-1]) / factor
        
    # Podstawianie wstecz
    for i in range(n - 2, -1, -1):
        b[i] -= A[i, 2] * b[i+1]
        
    return b

# --- 2. GENERATORY MACIERZY ---

def get_X_vector(n, precision):
    """Generuje wektor x_zad złożony z losowych wartości 1 i -1."""
    random.seed(42) # Dla powtarzalności wyników
    return np.array([1 if random.randint(0, 1) == 1 else -1 for _ in range(n)], dtype=precision)

def task_1_matrix(n, X, precision):
    """Zadanie 1: Macierz bliska macierzy Hilberta (źle uwarunkowana)."""
    A = np.zeros((n, n), dtype=precision)
    for i in range(n):
        for j in range(n):
            if i == 0: A[i, j] = 1.0
            else: A[i, j] = 1.0 / (i + 1 + j + 1 - 1)
    b = A @ X
    return A, b 

def task_2_matrix(n, X, precision):
    """Zadanie 2: Dobrze uwarunkowana macierz."""
    A = np.zeros((n, n), dtype=precision)
    for i in range(n):
        for j in range(n):
            if j >= i: A[i, j] = (2.0 * (i + 1)) / (j + 1)
            else: A[i, j] = (2.0 * (j + 1)) / (i + 1)
    b = A @ X
    return A, b 

def task_3_matrix_full(n, X, precision, k=5, m=4):
    """Zadanie 3: Macierz trójdiagonalna (wersja pełna n x n dla Gaussa)."""
    A = np.zeros((n, n), dtype=precision)
    for i in range(n):
        math_i = i + 1
        A[i, i] = -m * math_i - k
        if i < n - 1: 
            A[i, i + 1] = math_i
        if i > 0: 
            A[i, i - 1] = m / math_i
    b = A @ X
    return A, b

def task_3_matrix_banded(n, X, precision, k=5, m=4):
    """Zadanie 3: Macierz trójdiagonalna (wersja pasmowa n x 3 dla Thomasa)."""
    A_full, b = task_3_matrix_full(n, X, precision, k, m)
    A_banded = np.zeros((n, 3), dtype=precision)
    for i in range(n):
        math_i = i + 1
        A_banded[i, 1] = -m * math_i - k
        if i > 0: 
            A_banded[i, 0] = m / math_i
        if i < n - 1: 
            A_banded[i, 2] = math_i
    return A_banded, b

# --- 3. INFRASTRUKTURA BADAWCZA ---

def maximum_error(X_zad, X_obl):
    """Zwraca maksymalny błąd bezwzględny (norma nieskończoność)."""
    return np.max(np.abs(X_zad - X_obl))

def run_benchmark(task_name, matrix_func, method, n_range, precisions, track_cond=True):
    """Przeprowadza eksperyment dla zadanych parametrów i zapisuje do CSV."""
    os.makedirs('data', exist_ok=True)
    
    for precision in precisions:
        prec_str = 'float32' if precision == np.float32 else 'float64'
        results = []
        
        for n in n_range:
            X_zad = get_X_vector(n, precision)
            A, b = matrix_func(n, X_zad, precision)
            
            cond_val = np.nan
            if track_cond:
                # Uwarunkowanie = ||A|| * ||A^-1||
                try:
                    cond_val = np.linalg.cond(A, p=1)
                except np.linalg.LinAlgError:
                    cond_val = float('inf')
            
            tracemalloc.start()
            tracemalloc.reset_peak()
            
            gc.collect()
            start_time = time.perf_counter()
            X_obl = method(A.copy(), b.copy(), precision)
            solve_time = time.perf_counter() - start_time
            
            _, peak_memory = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            memory_kb = peak_memory / 1024.0
            
            error = maximum_error(X_zad, X_obl)
            
            results.append({
                'n': n,
                'error_max': error,
                'cond': cond_val,
                'time_s': solve_time,
                'memory_kb': memory_kb
            })
            
        df = pd.DataFrame(results)
        df.to_csv(f"data/{task_name}_{prec_str}.csv", index=False, float_format='%.6e')
        print(f"Zakończono: {task_name} | {prec_str} | Max n={n_range[-1]}")

# --- 4. WYKONANIE I WIZUALIZACJA ---

def generate_plots():
    """Odczytuje wyniki z CSV i generuje wykresy dla wszystkich zadań."""
    
    # --- ZADANIE 1 ---
    if os.path.exists("data/Zadanie_1_Gauss_float32.csv") and os.path.exists("data/Zadanie_1_Gauss_float64.csv"):
        df_z1_f32 = pd.read_csv("data/Zadanie_1_Gauss_float32.csv")
        df_z1_f64 = pd.read_csv("data/Zadanie_1_Gauss_float64.csv")
        
        # Wykres błędów Z1
        plt.figure(figsize=(8,4.5))
        plt.plot(df_z1_f32['n'], df_z1_f32['error_max'], 'o-', color='red', label='float32')
        plt.plot(df_z1_f64['n'], df_z1_f64['error_max'], 's-', color='blue', label='float64')
        plt.yscale('log')
        plt.axhline(1.0, color='gray', linestyle='--', label='Błąd krytyczny = 1.0')
        plt.xticks(np.arange(2, 21, 1))
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('Błąd maksymalny ||x_zad - x_obl||')
        plt.title('Zadanie 1: Załamanie precyzji (Macierz źle uwarunkowana)')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad1_error.png', dpi=150)
        plt.close()
        
        # Wykres uwarunkowania Z1
        plt.figure(figsize=(8,4.5))
        plt.plot(df_z1_f64['n'], df_z1_f64['cond'], 'D-', color='purple', label='Wskaźnik uwarunkowania')
        plt.yscale('log')
        plt.xticks(np.arange(2, 21, 1))
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('cond(A) = ||A|| * ||A⁻¹||')
        plt.title('Zadanie 1: Złe uwarunkowanie macierzy')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad1_cond.png', dpi=150)
        plt.close()

    # --- ZADANIE 2 ---
    if os.path.exists("data/Zadanie_2_Gauss_float32.csv") and os.path.exists("data/Zadanie_2_Gauss_float64.csv"):
        df_z2_f32 = pd.read_csv("data/Zadanie_2_Gauss_float32.csv")
        df_z2_f64 = pd.read_csv("data/Zadanie_2_Gauss_float64.csv")
        
        # Wykres błędów Z2
        plt.figure(figsize=(8,4.5))
        plt.plot(df_z2_f32['n'], df_z2_f32['error_max'], 'o-', color='red', label='float32')
        plt.plot(df_z2_f64['n'], df_z2_f64['error_max'], 's-', color='blue', label='float64')
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('Błąd maksymalny')
        plt.title('Zadanie 2: Stabilne rozwiązanie (Dobre uwarunkowanie)')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad2_error.png', dpi=150)
        plt.close()

        # Wykres uwarunkowania Z2
        plt.figure(figsize=(8,4.5))
        plt.plot(df_z2_f64['n'], df_z2_f64['cond'], 'D-', color='green', label='Wskaźnik uwarunkowania')
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('cond(A)')
        plt.title('Zadanie 2: Wskaźnik uwarunkowania')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad2_cond.png', dpi=150)
        plt.close()

    # --- ZADANIE 3 ---
    if os.path.exists("data/Zadanie_3_Gauss_float64.csv") and os.path.exists("data/Zadanie_3_Thomas_float64.csv"):
        df_g = pd.read_csv("data/Zadanie_3_Gauss_float64.csv")
        df_t = pd.read_csv("data/Zadanie_3_Thomas_float64.csv")
        
        # Wykres czasu Z3
        plt.figure(figsize=(8,4))
        plt.plot(df_g['n'], df_g['time_s'], 'o-', color='purple', label='Gauss O(n³)')
        plt.plot(df_t['n'], df_t['time_s'], 's-', color='green', label='Thomas O(n)')
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('Czas [s] (skala LOG)')
        plt.title('Zadanie 3: Przepaść wydajnościowa - Czas Obliczeń')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad3_time_log.png', dpi=150)
        plt.close()
        
        # Wykres pamięci Z3
        plt.figure(figsize=(8,4))
        plt.plot(df_g['n'], df_g['memory_kb'], 'o-', color='purple', label='Gauss O(n²)')
        plt.plot(df_t['n'], df_t['memory_kb'], 's-', color='green', label='Thomas O(n)')
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('Pamięć [KB] (skala LOG)')
        plt.title('Zadanie 3: Różnica w zużyciu pamięci RAM')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad3_mem_log.png', dpi=150)
        plt.close()
        
        # Wykres błędów Z3
        plt.figure(figsize=(8,4))
        plt.plot(df_g['n'], df_g['error_max'], 'o-', color='purple', label='Błąd - Gauss', alpha=0.7)
        plt.plot(df_t['n'], df_t['error_max'], 's-', color='green', label='Błąd - Thomas', alpha=0.7)
        plt.yscale('log')
        plt.xlabel('Rozmiar macierzy (n)')
        plt.ylabel('Błąd maksymalny')
        plt.title('Zadanie 3: Porównanie dokładności (Gauss vs Thomas)')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig('zad3_error.png', dpi=150)
        plt.close()

def main():
    precisions = [np.float32, np.float64]
    
    print("Rozpoczęto Zadanie 1...")
    run_benchmark("Zadanie_1_Gauss", task_1_matrix, gauss, range(2, 21), precisions, track_cond=True)
    
    print("Rozpoczęto Zadanie 2...")
    run_benchmark("Zadanie_2_Gauss", task_2_matrix, gauss, range(2, 151, 2), precisions, track_cond=True)
    
    print("Rozpoczęto Zadanie 3...")
    n_large = [10, 50, 100, 200, 300, 400, 500]
    run_benchmark("Zadanie_3_Gauss", task_3_matrix_full, gauss, n_large, [np.float64], track_cond=False)
    run_benchmark("Zadanie_3_Thomas", task_3_matrix_banded, thomas_memory, n_large, [np.float64], track_cond=False)
    
    print("\nGenerowanie wykresów...")
    generate_plots()
    
    print("Gotowe. Wyniki zapisane w folderze /data/, a wykresy w głównym katalogu.")

if __name__ == "__main__":
    main()