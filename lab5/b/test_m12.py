import numpy as np
from main import get_nodes, trig_approximation, f, X_DENSE

m = 12
n_vals = np.arange(25, 60)
err_inc, err_exc = [], []

for n in n_vals:
    x_in, y_in = get_nodes(n, True)
    app_in, _ = trig_approximation(x_in, y_in, m)
    err_inc.append(np.max(np.abs(f(X_DENSE) - app_in(X_DENSE))) if app_in else np.nan)
    
    x_ex, y_ex = get_nodes(n, False)
    app_ex, _ = trig_approximation(x_ex, y_ex, m)
    err_ex.append(np.max(np.abs(f(X_DENSE) - app_ex(X_DENSE))) if app_ex else np.nan)

for i in range(5):
    print(f"n={n_vals[i]} | IN: {err_inc[i]:.2f} | EX: {err_ex[i]:.2f}")
print("...")
for i in range(-5, 0):
    print(f"n={n_vals[i]} | IN: {err_inc[i]:.2f} | EX: {err_ex[i]:.2f}")
