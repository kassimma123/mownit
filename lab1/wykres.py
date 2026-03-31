import numpy as np
import matplotlib.pyplot as plt

#Dokładne pierwiastki
roots_exact = np.arange(1, 21)

#Wyliczamy współczynniki w 64-bitach
coeffs_64 = np.poly(roots_exact)

#ZMUSZAMY PYTHONA DO BŁĘDU (Symulacja typu float 32-bit z C++)
# Obcinamy precyzję współczynników, tak jakbyśmy zapisali je w małej pamięci
coeffs_32 = np.float32(coeffs_64)

#Obliczamy pierwiastki numerycznie z obciętych współczynników
roots_computed = np.roots(coeffs_32)

#Rysowanie wykresu
plt.figure(figsize=(10, 6))
plt.scatter(roots_exact, np.zeros(20), c='green', s=100, label='Dokładne (1 do 20)', zorder=5)
plt.scatter(roots_computed.real, roots_computed.imag, c='red', marker='x', s=100, label='Wyliczone numerycznie (float32)', zorder=5)

plt.axhline(0, color='black', lw=0.5)
plt.title('Pierwiastki Wielomianu Wilkinsona (Zjawisko ucieczki - float32)', fontsize=14)
plt.xlabel('Część Rzeczywista (Re)', fontsize=12)
plt.ylabel('Część Urojona (Im)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)

# Zapis do pliku
plt.savefig('wilkinson_plot.png', dpi=300, bbox_inches='tight')
print("Wykres zapisany jako wilkinson_plot.png - teraz musi być widać ucieczkę!")