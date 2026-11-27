# Zakres prototypu dynamicznego quasi-3D

Stan na 2026-11-26, po audytach 744–788 (prototyp prowadzenia 2D) i 835–866 (postać spinu w pomiarze p814).
Szkic zakresu przyjęty przez autora; pierwszy etap (E1) to audyt 868.

Cel: obserwable wszystkich poziomów n — wybór torusa, rozszczepienia, docelowo zanik i kaskada — mają wynikać z dynamiki pola
i źródła, a nie z wpisanej energii Sommerfelda. n = 1: stabilność stanu 1S (E1), rozszczepienie nadsubtelne wodoru F = 0 / F = 1
i p-Ps / o-Ps (E4–E5), anihilacja (E5–E6); n = 2: pierwszy poziom z degeneracją l — główny test struktury subtelnej (E3); n = 3:
skalowanie 1/n³ tym samym przepisem (bramka po E3; pomiar do dobrania z cytowaniem). Każdy kanał osobno: wodór F = 0 i F = 1,
p-Ps i o-Ps; porównanie z pomiarami.

Zastrzeżenie (uwaga autora po 867): przy danym widmie pola (baza modów wiązanych, Langer, faza kaustyk −π/2 wpisane) dynamika ma
wyprowadzić wybór stanu i rozszczepienia, a nie samą dyskretność. Pytanie o przyczynę kwantowania: test przesunięcia fazy kaustyk
(czy wybrany torus idzie za fazą źródła, czy zostaje na modach), docelowo pole bez bazy modów (poza E1–E7).
Prototyp nie zastępuje produkcji i nie walczy o 0,015 MHz, dopóki nie zgadzają się rzędy i znaki.

---

## 1. Punkt wyjścia

Prototyp prowadzenia (guide770–784, offline): pole ruchu względnego rozpisane na dokładne mody wiązane 2D −1/r + β/r²; cząstka na
torusie z wag modów; źródło z zegarem fazy S/ħ i skokiem −π/2 w punktach zwrotnych (licznik kaustyk, 769).

| | stan | audyty |
|---|---|---|
| działa | wybór torusa EBK (własny mod wygrywa 269/297); mechanizm skoku = pętla lokalna (widmo napędu trafia zwycięzcę w 85 %); przy fizycznej dobroci jednomodowe 19/20 | 778–780 |
| nie działa | kaskada promienista — cztery próby; najlepiej stała norma pola, 2/5; pole napędzane przez źródło broni modu wzbudzonego | 781–784 |
| sztuczne | β/r² 10³–10⁴ razy za duże i o odwrotnym znaku (w 2D fizycznego znaku nie da się sprawdzić — m = 0 spada na centrum); dobroć ~400 razy za niska; brak spinu | 780 |

Z serii p814: z Diracem w α⁶ zgadza się tylko postać ZDC (L_T + dipol własny + czynnik BMT + pęd kanoniczny); przyczynowa postać O
jest ~55σ od Diraca w tempie apsyd (866); szum trybu O w Ps ~0,0002 MHz (865); start pierwszego rzędu wystarcza (864, 865).

## 2. Etapy

Każdy etap = osobny audyt z preregiem i bramką; następny dopiero po przejściu poprzedniego.

| etap | co dodaje | bramka |
|---|---|---|
| E1 kinematyka quasi-3D | orbita w płaszczyźnie ⊥ L, \|L\| = (l + ½)ħ (Langer); mody radialne R_{n_r,l} 3D Coulomba z (l + ½)²; energia zależy tylko od n = n_r + l + 1 | poziomy = Bohr; wybór torusa EBK co najmniej jak w 2D; stany l = 0 nie spadają na centrum |
| E2 fizyczna poprawka relatywistyczna | zamiast β/r² poprawka o fizycznym znaku i wielkości (kinematyka γ, Darwin; bez spinu) | tylko wodór: precesja apsyd zgodna z Sommerfeldem bez spinu (znak, 1 %) — w Ps (m₁ = m₂) nie ma jednociałowego Sommerfelda; wybór torusa przetrwa fizyczną małość rozszczepienia |
| E3 struktura subtelna przez więz | rozszczepienie z precesji przez zasadę I (701, 711, 767), bez spinu | przedział 2P₃/₂ − 2P₁/₂ = 10 969,045(15) MHz (Hagley i Pipkin, PRL 72, 1172 (1994)) = 45,364 µeV = 2,7331 cm bez energii Sommerfelda, najpierw do 10⁻³ (ryzyko 807) |
| E4 spin jawny | precesja Thomasa–BMT, siły spinowe jednej postaci (decyzja a) | wodór F = 0 / F = 1: Dirac w α⁴ i znak α⁶; Ps p-Ps / o-Ps: ν₀ (Gurung, PRL 125, 073002 (2020)), ν₁ (Hagena, PRL 71, 2887 (1993)), ν₂ (Sheldon, PRL 131, 043001 (2023)) co do rzędu |
| E5 pozytonium | równe masy, dwa spiny, odrzut | dwuciałowe α⁴ (Breit + anihilacja wirtualna ³S₁): ν₀, ν₁, ν₂ (o-Ps, tryplet) i rozszczepienie 1S p-Ps / o-Ps 203 389,10(74) MHz (Ritter, PRA 30, 1331 (1984)); kanał p-Ps n = 2 — pomiar do sprawdzenia; szum trybu O (865) jako dolna granica |
| E6 zanik i kaskada | ziarno zaniku z fluktuacji próżni / Weisskopf–Wigner ze stałą normą (784) / skok jako pomiar (755) | 2p wodoru 1,600(4) ns (Bickel i Goodman 1966) w skali prototypu co do rzędu; zejście przez torusy EBK z Δl = ±1 |
| E7 decyzja o produkcji | co przenieść do silnika, z wyłącznikiem | reguły 848 (bajt w bajt, walidacja 72/72) |

## 3. Decyzje

- (a) Postać spinu w E4: O z poprawkami 841 — **rozstrzygnięte przez 867** (Δ_O / Δ_Z = 1,0000 / 0,9998 w wodorze F = 0 / F = 1,
  1,002 w o-Ps; odchyłka od Diraca +5,3 / −6,9σ₂P; reszta = brakujący pęd kanoniczny Tulczyjewa, osobny audyt).
- (b) E6 nie blokuje E1–E5 (osobny tor badawczy).
- (c) Prototyp offline (osobny plik w `audits/data/section-NNN`, jak guide784); silnik dopiero w E7.
- (d) Kolejność: E1 najpierw (najtańszy; od razu odpowiada, czy l = 0 z Langerem przestaje spadać na centrum).
- (e) Etykiety kanałów od E4: wodór n ≥ 2 jako (j, F) (2P₁/₂ F = 0/1, 2P₃/₂ F = 1/2 — w stanach P F nie jest spinem pary);
  pozytonium jako (S, J) osobno p-Ps i o-Ps (S dobre dzięki parzystości C). Odpowiedniość F = 0 ↔ p-Ps, F = 1 ↔ o-Ps tylko w stanach S.

## 4. Ryzyka

- E3 może powtórzyć 807: żaden przepis precesja → przedział nie daje członu α² Diraca; wtedy struktura subtelna zostaje importem i
  trzeba to zapisać.
- Koszt modów (n_r, l) w quasi-3D; 7,8 GB RAM — jeden stan = jeden proces, `ulimit`; audyty nie równolegle.
- Fizyczna dobroć wydłuża czas ustalania (otwarte; do zmierzenia w E2).

## 5. Dziennik etapów

| etap | audyt | wynik |
|---|---|---|
| E1 | 868 | **niezaliczony**: Bohr dokładnie i bez spadku na centrum (H868k, H868l), licznik θ konieczny (bez niego 0/20), ale wybór torusa EBK 10/40 wobec 24/40 w 2D (sukces 37 %); POST HOC atraktory mieszane blisko modów, l = 0 nieodwiedzane; E2 wstrzymany do rozbicia przyczyny |
