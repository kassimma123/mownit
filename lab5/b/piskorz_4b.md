#set document(title: "Sprawozdanie - Aproksymacja Trygonometryczna", author: "Katarzyna Piskorz")
#set page(paper: "a4", margin: 2cm)
#set text(font: "Linux Libertine", size: 10pt, lang: "pl")
#set heading(numbering: "1.")

#align(center)[
  #text(16pt, weight: "bold")[Metody Obliczeniowe w Nauce i Technice]\
  #v(0.3em)
  #text(13pt)[Aproksymacja Średniokwadratowa Trygonometryczna]
  #v(0.5em)
  #grid(
    columns: (1fr, 1fr),
    align: (left, right),
    [Autorka: *Katarzyna Piskorz*],
    [Data: 28.04.2026]
  )
]

= Analiza przypadku i cel badania

Celem doświadczenia jest wyznaczenie przybliżenia zadanej funkcji za pomocą aproksymacji średniokwadratowej przy użyciu wielomianów trygonometrycznych oraz analiza błędów tego przybliżenia[cite: 4, 286]. 

Rozpatrywana jest funkcja bazowa[cite: 6, 276]:
$ f(x) = x^2 - m cos((pi x) / k) $
gdzie $m = 5.0$ oraz $k = 0.5$, określona na przedziale $x in [-5, 5]$[cite: 7, 277]. Składowa trygonometryczna funkcji wprowadza wysokoczęstotliwościowe oscylacje, które stanowią wyzwanie dla stabilności numerycznej modelu[cite: 389].

= Metodologia i Aparat Matematyczny

Aproksymacja realizowana jest metodą wyznaczania współczynników kombinacji liniowej funkcji bazy trygonometrycznej[cite: 287].

*Szczegółowy opis aparatu matematycznego:*
+ *Transformacja dziedziny (Mapowanie):* Aby uniknąć zjawiska Gibbsa oraz błędu wynikającego z braku okresowości funkcji bazowej na brzegach przedziału, zastosowano liniowe odwzorowanie przedziału $[-5, 5]$ na półokres $[0, pi]$ wzorem:
  $ v_i = pi (x_i + 5) / 10 $
  , gdzie:
  - $x_i$: oryginalna współrzędna węzła w dziedzinie $[-5, 5]$.
  - $v_i$: przekształcona współrzędna w dziedzinie $[0, pi]$. 
  *Znaczenie:* Takie "rozciągnięcie" pozwala funkcjom trygonometrycznym na naturalne dopasowanie się do trendu kwadratowego funkcji $f(x)$ bez wymuszania sztucznej periodyzacji na krańcach[cite: 289].

+ *Baza modelu ($B$):* Przyjęto układ $2m+1$ funkcji trygonometrycznych stanowiących bazę aproksymacji[cite: 285]:
  $ B(v) = [1, cos(v), sin(v), cos(2v), sin(2v), ..., cos(m v), sin(m v)] $
  *Znaczenie:* Każdy element wektora $B(v)$ to kolejna harmoniczna. Stopień $m$ określa najwyższą częstotliwość dostępną dla modelu. Aproksymacja jest sumą wałoną tych funkcji: $W(v) = sum c_j B_j (v)$.

+ *Budowa układu równań:* Współczynniki $c$ wyznaczono z układu równań normalnych (Macierz Grama)[cite: 10, 12]:
  $ G c = b, quad "gdzie" quad G = B^T B, quad b = B^T y $
  - $G$ – macierz Grama (zbudowana z iloczynów skalarnych funkcji bazy)[cite: 17].
  - $b$ – wektor wyrazów wolnych[cite: 18].
  - $c$ – wektor szukanych współczynników wagowych harmonicznych[cite: 16].

+ *Rozwiązanie i warunki:* Układ rozwiązano funkcją `numpy.linalg.solve()`, wykorzystującą *stabilny rozkład LU*[cite: 22, 23]. Warunkiem koniecznym jest $n >= 2m + 1$[cite: 24, 288].

= Wyznaczenie dokładności aproksymacji

Pomiar przeprowadzono na gęstej siatce $l = 1000$ punktów pomiarowych[cite: 25, 294].

*Błąd średniokwadratowy (MSE):* [cite: 27, 296, 298]
$ E_"MSE" = sqrt(sum_(i=1)^l (W(v_i) - f(x_i))^2) / l $

*Błąd maksymalny:* [cite: 28, 304, 306]
$ E_"max" = max_(i in {1..l}) |W(v_i) - f(x_i)| $

, gdzie:
- $l$ – liczba punktów pomiaru błędu [cite: 29, 307]
- $W$ – funkcja wielomianu trygonometrycznego aproksymującego [cite: 32, 301]
- $f$ – badana funkcja bazowa [cite: 33, 300]
- $x_i$ – punkt, w którym dokonywany jest pomiar [cite: 302, 311]

= Analiza węzła końcowego (Parametr endpoint)

#figure(
  image("trig_wezel_koncowy.png", width: 90%),
  caption: [Różnica w aproksymacji z włączonym (lewo) i wyłączonym (prawo) węzłem brzegowym $x=5$.]
)
*Wniosek:* Zastosowana metoda macierzowa #underline[wymaga dołączenia ostatniego węzła]. Odpięcie prawego brzegu przerywa ciągłość informacyjną, co skutkuje wyraźnym odchyleniem wielomianu od funkcji na skraju przedziału.

= Wizualizacja dla stałej liczby węzłów $n$ (Wykrywanie ekstremów)

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("trig_detekcja_n40.png", width: 100%),
    image("trig_detekcja_n50.png", width: 100%)
  ),
  caption: [Proces osiągania rezonansu częstotliwościowego: brak detekcji ($n=40$) vs pełna detekcja ($n=50, m=20$).]
)

*Dowód wizualny:*
- *Dla $n=40$:* Warunek oznaczoności ogranicza stopień do $m <= 19$. Baza nie osiąga częstotliwości składowej $cos(2pi x)$, co skutkuje #underline[wysokim błędem niedopasowania (bias)].
- *Dla $n=50$:* Model zyskuje dostęp do harmonicznej $m=20$. Następuje #underline[idealny rezonans]. Składowa $cos(2pi x)$ (która w skali $v$ odpowiada 20-stej harmonicznej) jest idealnie odwzorowana, a błąd skokowo spada[cite: 379, 380].

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("trig_detekcja_n80.png", width: 100%),
    image("trig_granica_n100.png", width: 100%)
  ),
  caption: [Stabilność modelu vs wybuch błędów maszyny cyfrowej na granicy oznaczoności układu.]
)

*Fakt:* - Dla $m >= 20$ wykresy niemal się pokrywają – #underline[model matematyczny osiągnął pełną zbieżność].
- Dla $n=100, m=49$: Na brzegach pojawiają się wibracje. To dowód na #underline[złe uwarunkowanie macierzy Grama] przy zbliżaniu się do granicy oznaczoności układu ($n approx 2m+1$)[cite: 407, 408].

= Analiza uwarunkowania układu (Mapy Ciepła)

#figure(
  image("heatmap_trig_mse.png", width: 95%),
  caption: [Mapa Ciepła MSE dla bazy trygonometrycznej (Puste pola to $n < 2m+1$).]
)

*Analiza zjawisk:*
- Mapa przyjmuje kształt trójkątny ze względu na matematyczne ograniczenie algorytmu[cite: 368].
- #underline[Dolina błędów]: Obszar optymalny ($E_"MSE" approx 10^(-4)$) jest rozległy, co dowodzi wysokiej stabilności bazy trygonometrycznej.
- #underline[Ściana błędów]: Gwałtowny skok błędu następuje wyłącznie na krawędzi matematycznej rozwiązywalności układu.

= Wybór wariantu optymalnego 

#figure(
  image("trig_najlepsze.png", width: 85%),
  caption: [Najbardziej stabilny i dokładny model: $n=100, m=20$.]
)

*Naukowe uzasadnienie wyboru:*
+ *Zgodność widmowa:* Stopień $m=20$ to najniższa wartość w pełni pokrywająca częstotliwość składowej bazowej, co eliminuje błąd modelu.
+ *Nadokreśloność:* Siatka $n=100$ zapewnia ogromny margines bezpieczeństwa ($100 limits(>>) 41$), gwarantując #underline[doskonałe uwarunkowanie macierzy Grama] i odporność na błędy zmiennoprzecinkowe[cite: 395].