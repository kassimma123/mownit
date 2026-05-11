import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
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
    b[0] /= A[0, 1]
    A[0, 2] /= A[0, 1]
    A[0, 1] = 1.0
    
    for i in range(1, n):
        factor = A[i, 1] - A[i, 0] * A[i-1, 2]
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
            else: A[i, j] = 1.0 / (i + 1 + j + 1 - 1) # i, j są indeksowane od 0
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
        A[i, i] = -m * (i + 1) - k
        if i > 0: A[i, i - 1] = (i + 1) / m
        if i < n - 1: A[i, i + 1] = (i + 1) / m # Zgodnie z symetrią podaną w treści
    b = A @ X
    return A, b

def task_3_matrix_banded(n, X, precision, k=5, m=4):
    """Zadanie 3: Macierz trójdiagonalna (wersja pasmowa n x 3 dla Thomasa)."""
    A_full, b = task_3_matrix_full(n, X, precision, k, m)
    A_banded = np.zeros((n, 3), dtype=precision)
    for i in range(n):
        A_banded[i, 1] = -m * (i + 1) - k
        if i > 0: A_banded[i, 0] = (i + 1) / m
        if i < n - 1: A_banded[i, 2] = (i + 1) / m
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
            
            # Pamięć start
            tracemalloc.start()
            tracemalloc.reset_peak()
            
            # Czas start
            gc.collect()
            start_time = time.perf_counter()
            X_obl = method(A.copy(), b.copy(), precision)
            solve_time = time.perf_counter() - start_time
            
            # Pamięć stop
            _, peak_memory = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            memory_kb = peak_memory / 1024.0
            
            error = maximum_error(X_zad, X_obl)
            
            cond_val = np.nan
            if track_cond:
                # Uwarunkowanie liczymy tylko dla pełnych macierzy
                cond_val = np.linalg.cond(A, p=1)
                
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

def main():
    precisions = [np.float32, np.float64]
    
    # Doświadczenie 1 (Macierz źle uwarunkowana, awaria float32 w okolicach n=7, float64 n=13)
    run_benchmark("Zadanie_1_Gauss", task_1_matrix, gauss, range(2, 21), precisions, track_cond=True)
    
    # Doświadczenie 2 (Macierz dobrze uwarunkowana, większe rozmiary)
    run_benchmark("Zadanie_2_Gauss", task_2_matrix, gauss, range(2, 151, 5), precisions, track_cond=True)
    
    # Doświadczenie 3 (Porównanie złożoności: Thomas vs Gauss)
    n_large = [10, 50, 100, 200, 300, 400, 500]
    run_benchmark("Zadanie_3_Gauss", task_3_matrix_full, gauss, n_large, [np.float64], track_cond=False)
    run_benchmark("Zadanie_3_Thomas", task_3_matrix_banded, thomas_memory, n_large, [np.float64], track_cond=False)
    
    print("\nGotowe. Wyniki zapisane w folderze /data/ do tabel w raporcie.")

if __name__ == "__main__":
    main()