import numpy as np
import matplotlib.pyplot as plt
from main import get_nodes, trig_approximation, f, X_DENSE

m = 10
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

for i, n in enumerate([21, 22, 26]):
    x_in, y_in = get_nodes(n, True)
    app_in, c_in = trig_approximation(x_in, y_in, m)
    
    ax[i].plot(X_DENSE, f(X_DENSE), 'k--', alpha=0.5)
    if app_in: 
        ax[i].plot(X_DENSE, app_in(X_DENSE), 'r')
    ax[i].set_title(f"n={n}, cond={c_in:.1e}")
    ax[i].set_ylim(-10, 35)

plt.savefig("test_wobble.png")
