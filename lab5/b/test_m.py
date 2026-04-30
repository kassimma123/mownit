import numpy as np
from main import get_nodes, trig_approximation, f, X_DENSE

m = 10
n = 30
x_in, y_in = get_nodes(n, include_endpoint=True)
app_in, cond_in = trig_approximation(x_in, y_in, m)
err_in = np.max(np.abs(f(X_DENSE) - app_in(X_DENSE)))

x_ex, y_ex = get_nodes(n, include_endpoint=False)
app_ex, cond_ex = trig_approximation(x_ex, y_ex, m)
err_ex = np.max(np.abs(f(X_DENSE) - app_ex(X_DENSE)))

print(f"For m={m}, n={n}:")
print(f"With endpoint: error={err_in}, cond={cond_in}")
print(f"Without endpoint: error={err_ex}, cond={cond_ex}")
