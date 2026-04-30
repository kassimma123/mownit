import re

with open('main.py', 'r') as f:
    content = f.read()

new_get_nodes = """def get_nodes(n, include_endpoint=False):
    x = np.linspace(INTERVAL[0], INTERVAL[1], n, endpoint=True)
    if not include_endpoint:
        x = x[:-1]
    return x, f(x)"""

content = re.sub(r"def get_nodes\(n, include_endpoint=False\):\n    x = np.linspace\(INTERVAL\[0\], INTERVAL\[1\], n, endpoint=include_endpoint\)\n    return x, f\(x\)", new_get_nodes, content)

with open('main.py', 'w') as f:
    f.write(content)
