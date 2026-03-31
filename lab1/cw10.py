
def get_wilkindon_coeffs(n):
    coeffs = [1.0] 
    for r in range(1, n+1):
        new_coeffs = [0.0] * (len(coeffs) + 1)
        for i in range(len(coeffs)):
            new_coeffs[i] += coeffs[i] * (-r)
            new_coeffs[i+1] += coeffs[i]
        coeffs = new_coeffs
    return coeffs

def horner(coeffs, x):
    result = 0.0
    for i in range(len(coeffs) -1, -1, -1):
        result = result * x + coeffs[i]
    return result

N = 20
a = get_wilkindon_coeffs(N)

print(f"Współczynniki wielomianu W(x) dla n={N}:")
for i, val in enumerate(a):
    print(f"a_{i} = {val:.0f}")

print("\n--- TEST SCHEMATU HORNERA ---")
test_points = [1.0, 2.0, 19.0, 20.0]

for x in test_points:
    value_horner = horner(a, x)
    print(f"Dla x = {x}:")
    print(f"  Wartość (Horner): {value_horner}")
    print(f"  Błąd bezwzględny: {abs(value_horner)}")