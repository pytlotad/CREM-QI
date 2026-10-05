# Semi-klasyczny model pozytonium (CREM)

Ten plik opisuje **model fizyczny** w jednym miejscu: z czego się składa, co
jest w nim klasyczne, co jest importem kwantowym i jak wypada wobec pomiarów.
Szczegóły pomiarów, historia zmian i wycofane twierdzenia są w `README.md` i w
pliku audytu `audits/2026-09-10-para-ortho-after-dipole-fix.txt` (z danymi w
`audits/data/section-NNN/`); tutaj są tylko odsyłacze do nich. Stan opisu:
audyt 323 (2026-10-04).

## 1. Czym jest model

Model jest **klasyczną elektrodynamiką dwóch punktowych cząstek**, elektronu i
pozytonu, z których każda niesie ładunek i moment magnetyczny. Do tej
dynamiki dołożono **nazwaną, policzalną listę importów kwantowych**: miejsc,
w których wynik kwantowy wchodzi jako założenie. Model jest semi-klasyczny
dokładnie w tym sensie: trajektoria jest klasyczna, a dyskretność (kwant
emisji, moment pędu fotonu, przygotowanie stanu, reguła anihilacji) jest
wpisana i oznaczona.

Cel jest jeden: sprawdzić, **jak daleko klasyczna elektrodynamika sięga w
odtwarzaniu wyników doświadczalnych pozytonium**, zanim zabraknie jej opisu.
Model ocenia się wobec **pomiarów**, nie wobec innych modeli (sekcja 7).
Żaden człon nie jest dopasowywany do danych; jedynymi wejściami są stałe
CODATA i zmierzone momenty magnetyczne.

## 2. Składniki i skale

Każda cząstka ma masę \(m\), ładunek \(\pm e\) i moment magnetyczny
\(\mu=(g/2)\mu_B\) ze zmierzonym \(g\). Para ma masę zredukowaną
\(\mu_{\rm red}=m_e/2\). Wszystkie długości modelu wynikają z jednej skali,
\(\lambda_C\), i jednej liczby bezwymiarowej, \(\alpha\) (relacje są
zablokowane asercjami w `modules/particle_species.hpp`):

| skala | wyrażenie | wartość dla e⁺e⁻ |
|---|---|---|
| sprzężenie Coulomba | \(k=\alpha\hbar c\) | — |
| promień Bohra pary | \(a_{\rm pary}=2\lambda_C/\alpha\) | \(2a_0=105{,}8\) pm |
| promień „bariery” | \(r^*=\mu/(ec)=(g/2)\lambda_C/2\) | \(193{,}3\) fm |
| rozmycie pola dipola | \(\varepsilon=0{,}96682\,r^*\) | \(186{,}9\) fm |

\(r^*\) jest dokładnie promieniem, na którym energia dipol–dipol równa się
kulombowskiej: ich stosunek to \((r^*/r)^2\). Długość rozmycia \(\varepsilon\)
nie jest parametrem dopasowania: moment magnetyczny traktuje się jak pętlę
prądu o najmniejszym promieniu, przy którym ładunek \(e\) niesie \(\mu\) z
prędkością nie większą niż \(c\) (\(R=2r^*\)), a profil Plummera ma tę samą
wagę kontaktową co jednorodny dysk tej pętli (audyt 86).

## 3. Dynamika klasyczna

Całkowane są położenia, pędy i wektory momentów magnetycznych obu cząstek
(`modules/electrodynamics.hpp`, `modules/crem_engine.hpp`).

1. **Kinematyka relatywistyczna:** \(\mathbf p=\gamma m\mathbf v\); warunek
   \(v<c\) wynika z tej relacji, nie z obcinania.
2. **Oddziaływanie ładunków:** pełne, retardowane pola Liénarda–Wiecherta
   drugiej cząstki i siła Lorentza
   \(\mathbf F_i=q_i[\mathbf E_j+\mathbf v_i\times\mathbf B_j]_{t_{\rm ret}}\),
   dokładne we wszystkich rzędach \(v/c\). Przybliżenie Darwina występuje
   tylko w diagnostyce energii.
3. **Oddziaływanie momentów:** pole dipola rozmyte profilem Plummera
   (\(\varepsilon\) wyżej), z członem magnetyzacji (kontaktowym), tak żeby
   \(\nabla\cdot\mathbf B=0\). Daleko daje to zwykłą energię
   \(U_{dd}\propto[\boldsymbol\mu_1\cdot\boldsymbol\mu_2-3(\boldsymbol\mu_1\cdot\hat r)(\boldsymbol\mu_2\cdot\hat r)]/r^3\),
   której znak dla orbit kołowych zależy od \(P_2(\cos\theta)\). Przy
   \(r\lesssim r^*\) dominuje człon kontaktowy, który **przyciąga momenty
   równoległe w każdym ułożeniu** (audyt 312).
4. **Ruch momentów:** precesja BMT w lokalnym polu. Część orbitalna precesji
   jest identyczna dla obu cząstek, a część od dipola partnera wynosi
   \(\pm1{,}38\cdot10^{11}\) rad/s, ze znakiem zależnym od kanału (audyt 313).
   W estymatorze sekularnym wektory \(\mathbf L,\mathbf S_1,\mathbf S_2\) są
   transportowane po orbicie uśrednionej (`modules/secular_spin_orbit.hpp`).
5. **Promieniowanie** (sekcja 4).

**Konwencja kanałów.** W p-Ps spiny są antyrównoległe, więc momenty
magnetyczne elektronu i pozytonu są **równoległe**
(\(|\boldsymbol\mu_1+\boldsymbol\mu_2|=2\mu\)). W o-Ps momenty są
antyrównoległe. Kanały różnią się na starcie tylko znakiem
\(\boldsymbol\mu_2\). Domyślnie kąt wzajemny jest losowany swobodnie (audyt
91), a kanał jest klasyfikacją wylosowanego układu. `--spin-quantization`
narzuca \(\cos=\pm1\).

**Przygotowanie stanu związanego.** Start na \(a_n=n^2a_{\rm pary}\), domyślnie
\(n=1\), ostro kołowy: \(L=n\hbar\) z konstrukcji (`--level`).
`CREM_INITIAL_ANGULAR_MOMENTUM=J₀` daje ten sam poziom energii przy
\(L=J_0\,n\hbar\), czyli orbitę o \(e^2=1-J_0^2\) (audyt 301).
`--microcanonical-start` (audyt 326) losuje \(J_0^2\) jednostajnie z \([0,1]\):
to klasyczny zespół mikrokanoniczny przy energii \(E_n\), czyli odpowiednik
całej powłoki \(n\), a nie jej kołowego członu. Jako jedyny klasyczny stan
\(n=1\) ma gęstość w zerze rzędu \(|\psi_{1s}(0)|^2\) (audyt 324). Strumień
losowy jest ten sam co przy okręgu, więc oba starty da się porównywać parami.

## 4. Promieniowanie: dwa tryby

**Ciągły (domyślny):** zredukowana siła Landaua–Lifshitza na każdą cząstkę
(`--radiation-reaction individual`). Orbita traci energię płynnie po
obwiedni Larmora; nie ma fotonów ani kwantu.

**Fotonowy (`--radiation-reaction stochastic`):** ta sama moc E1 trafia do
hazardu i jest wypłacana dyskretnymi fotonami, które zachowują pęd
(estymator sekularny `estimateCremCollapse`, `modules/crem_collapse.hpp`).
Ten tryb był używany we wszystkich badaniach kaskady i para/orto. Jego
reguły, w kolejności działania:

- **Hazard:** \(\lambda=S(e)\,P_{E1}(a,e)/\hbar\omega_{\rm orb}\), gdzie
  \(S(e)\) uwzględnia, że orbita eliptyczna rozdziela moc między harmoniczne
  (tablica interpolowana w \(e^2\), bo przy małym \(e\) mamy \(S\approx1-2e^2\); audyt 340).
  Domyślnie foton pada, gdy skumulowany hazard osiągnie 1
  (deterministycznie); `--emission poisson` losuje próg z
  \(\mathrm{Exp}(1)\).
- **Energia:** \(E_\gamma=k\,\hbar\omega_{\rm orb}\). Harmoniczna \(k\) jest
  losowana z widma dipola orbity Keplera, \(d(t)=e\,\mathbf r(t)\), które
  dla \(e>0\) ma składowe \(k\omega\) (tablica z DFT,
  `eccentricOrbitHarmonicNumber`). Na okręgu \(k=1\). Źródłem harmonicznych
  jest wyłącznie kształt orbity, a nie momenty magnetyczne; fotony M1 mają
  \(k=1\).
- **Kierunek i skrętność:** kierunek pochodzi z wzoru kątowego pola dalekiego,
  \(\propto1+\cos^2\theta\) względem \(\mathbf L\). Skrętność \(h=\pm1\) jest
  losowana z \(P(+)=(1+c)^2/(2(1+c^2))\).
- **Moment pędu:** foton zabiera dokładnie \(\hbar\) wzdłuż kierunku lotu,
  więc \(\mathbf L'=\mathbf L-h\hbar\hat d\) (domyślne od 2026-09-26;
  `CREM_NO_SPIN_MAGNITUDE=1` wyłącza).
- **Sufit energii:** orbita z \(L'\) musi istnieć (\(e^2\ge0\)), więc
  \(E'\ge-B_1(\hbar/L')^2\). Foton większy od tej różnicy jest przycinany, a
  gdy sufit wychodzi \(\le0\), losowanie jest odrzucane.
- **Między fotonami:** dynamika zachowawcza (audyt 288). Całe promieniowanie
  wychodzi fotonami.
- **Stan końcowy:** gdy \(|L/\hbar-1|\ge n\), żaden foton nie przechodzi
  sufitu. Kanał E1 jest zamknięty i przebieg kończy się z przyczyną
  `EmissionChannelClosed` (audyt 308). Dla orbity kołowej zachodzi to przy
  \(n<1/2\).

**Czego tryb fotonowy nie ma.** Nie ma dyskretnych poziomów: kwant
\(\hbar\omega\) własnej orbity z koła \(n\) prowadzi do \(n\to\sqrt{n^3/(n+2)}\).
Nie ma też przejścia \(2\to1\) z koła \(n=2\): foton zabiera najwyżej
\(\hbar\), a \(L=2\hbar\) nie zejdzie do \(L=\hbar\) inaczej niż przy
idealnym ustawieniu, czyli z prawdopodobieństwem zero. Zmierzone: \(0\) z
\(27\) pierwszych fotonów (audyt 321).

## 5. Importy kwantowe

Każdy import ma przełącznik albo jest nazwany jako nieusuwalny. Ich
zmierzone ceny są w README („Sześć importów kwantowych i ich cena”).

| import | co wpisuje | domyślnie |
|---|---|---|
| promień i moment startowy | \(a_n=n^2a_{\rm pary}\), \(L=n\hbar\) (obraz Bohra) | **tak**, bez przełącznika dla promienia; `--microcanonical-start` zastępuje \(L=n\hbar\) zespołem \(L^2\) jednostajnym |
| kwant emisji | \(E_\gamma=k\hbar\omega_{\rm orb}\) | tylko w trybie fotonowym |
| \(\hbar\) fotonu | \(\mathbf L'=\mathbf L-h\hbar\hat d\) | tak, w trybie fotonowym |
| kwantyzacja spinu | \(\cos(\boldsymbol\mu_1,\boldsymbol\mu_2)=\pm1\) | nie (`--spin-quantization`) |
| podłoga stanu podstawowego | emisja gaśnie przy \(n=1\); **L podnoszone do \(\hbar\)**, orbita kołowa | nie (`--ground-state-floor`) |
| drabina Bohra | \(E_\gamma=E(n)-E(n-k)\) | nie (`--bohr-photon-energy`; poprawiona w audycie 322) |

**Jedno założenie skalowe.** Skala całego modelu bierze się z \(L=\hbar\):
przy ustalonym \(L\) minimum \(L^2/(2\mu r^2)-k/r\) leży dokładnie w
\(a_{\rm pary}\). To jest argument Bohra. Stan 1s ma \(L=0\), a jego skalę
kwantowo daje energia lokalizacji, której model nie ma. Pozostałe importy
to reguły dyskretności i skali nie ustalają (README, „Sześć importów to jedno
założenie skalowe i pięć reguł”).

**Konsekwencja tego samego założenia:** klasyczny okrąg \(n\) ma
\(L=n\hbar\), o jeden kwant więcej niż kwantowe \(l\le n-1\). Stąd zamknięcie
kanału E1 przy \(n\approx1/2\) zamiast w \(n=1\) (audyt 308) i brak przejścia
\(2\to1\) z okręgu (audyt 321).

## 6. Anihilacja

**Klasyczna dynamika nie zawiera anihilacji**: nie ma w niej ani kanału
kontaktowego, ani tempa. Czasy podawane przez tryby 1 i 2 to czasy kaskady
albo kolapsu, a nie czasy życia. Generator rozpadu w idealnej próżni jest
osobną, jawnie kwantową receptą.

**Eksperyment 6** (`--mode statistical --phenomenon 6`, audyty 316–320)
daje modelowi zdarzenie anihilacji, za cenę trzech importów:

1. anihilacja przy pierwszym wejściu w \(r\le r^*\);
2. tłumienie \(3\gamma\) czynnikiem Ore–Powella
   \(\varepsilon_{\rm OP}=4(\pi^2-9)\alpha/(9\pi)=1/1113{,}9\) (teoria wiodącego
   rzędu QED);
3. rozkład płaszczyzny \(3\gamma\) względem spinu,
   \(dN/d\cos\theta_n\propto1-\tfrac13\cos^2\theta_n\) dla \(m=\pm1\),
   policzony z drzewowej amplitudy QED (audyt 320).

Z modelu pochodzą: orbita (start bliski 1s, \(L\approx0{,}07\hbar\), dynamika
zachowawcza), siła kontaktowa i momenty w chwili kontaktu. Waga
\(w_{2\gamma}=(|\boldsymbol\mu_1+\boldsymbol\mu_2|/2\mu)^2\) wybiera kanał z
\(P(2\gamma)=w/(w+(1-w)\varepsilon_{\rm OP})\). Fotony powstają w układzie
spoczynkowym pary: \(2\gamma\) wzdłuż izotropowej osi, \(3\gamma\) z dokładnym
rozkładem energii Ore–Powella i płaszczyzną związaną z
\(\mathbf S_1+\mathbf S_2\). Bezwzględne czasy nie mają sensu fizycznego;
obserwablą jest stosunek.

## 7. Wobec pomiarów

Pomiary odniesienia: \(\tau_{\rm para}=125{,}14\) ps (Al-Ramadhan i Gidley,
PRL 72, 1632 (1994)); \(\tau_{\rm orto}=142{,}04\) ns (Vallery i in., PRL 90,
203402 (2003)); odstęp 1S–2S \(5{,}10179\) eV (Fee i in., PRL 70, 1397
(1993)); \(m_ec^2=510{,}99895\) keV (CODATA 2018).

| obserwabla | model | pomiar | ocena |
|---|---|---|---|
| tempo w stanie podstawowym | \(P/\hbar\omega\) daje \(186{,}74\) ps | \(\tau_{\rm para}=125{,}14\) ps | ta sama potęga \(\alpha^5\), czynnik \(1{,}49\); współczynnik \(O(1)\) niewyjaśniony, to nie jest wyprowadzenie anihilacji |
| skala \(\alpha^5\) bez kwantu | \(P/E_{\rm wiąz}\) | — | przeżywa usunięcie \(\hbar\omega\) (czynnik \(0{,}746\)) |
| rozpraszanie | kształt Rutherforda, jedna normalizacja | — | klasyczne; danych niskoenergetycznych e⁺e⁻ projekt nie używa |
| \(\tau_{\rm orto}/\tau_{\rm para}\), kaskada | \(1{,}000006\) | \(1135{,}0\) | **nie odtworzone**: mechanizm klasyczny jest o \(\sim7\) rzędów za słaby |
| \(\tau_{\rm orto}/\tau_{\rm para}\), eksperyment 6, \(n=1\) | \(27\,215\pm42\%\) = geometria \(24{,}4\) × reguła \(1113\) | \(1135{,}0\) | reguła wyboru jest importem; geometria, jedyne własne przewidywanie, powinna wynosić \(\approx1\) |
| to samo przy \(n=2,3\) | \(196\,755\), \(86\,011\) | \(1135{,}0\) | stosunek zależy od \(n\) (\(3{,}2\sigma\)), geometria \(177\), \(77\) |
| odstęp \(n=2\to1\) | brak przejścia z okręgu | \(5{,}10179\) eV | drabina Bohra (import) daje \(5{,}10214\) eV, \(6{,}8\cdot10^{-5}\) od pomiaru |
| rozszczepienie nadsubtelne 1s | koło \(L=\hbar\): \(6{,}7\cdot10^{-7}\) pomiaru; zespół mikrokanoniczny \(n=1\): \(290{,}5\) GHz | \(203{,}3942\) GHz (Ishida i in., PLB 734, 338 (2014)) | znak poprawny (orto wyżej), człon kontaktowy; wynik zależy od tego, który stan klasyczny gra rolę \(1s\) — przy zespole mikrokanonicznym czynnik \(1{,}43\), ale przez kompensację braków, nie z wyprowadzenia (audyt 324) |

**Wynik negatywny jest wynikiem.** Różnica para/orto w przyrodzie pochodzi z
reguły wyboru \(2\gamma/3\gamma\), czyli z zachowania parzystości
ładunkowej w procesie, którego klasyczna dynamika nie zawiera. Jedyne miejsce,
w którym kanały rozdzielają się w modelu z właściwym znakiem, to człon
kontaktowy przy \(r\lesssim r^*\) na orbitach niemal radialnych (audyty
312–317). Jego siła jest klasycznym przesadzeniem, bo kwantowo
\(|\psi(0)|^2\) w wiodącym rzędzie nie zależy od spinu.

## 8. Granice dziedziny

- **Bariera Comptona** \(r^*\): poniżej niej klasyczna elektrodynamika
  punktowa przestaje obowiązywać, a przebieg się kończy.
- **Granica retardacji:** okres / czas przelotu światła \(\le150\). To
  margines przybliżenia Darwina i całkowania, a nie fizyka.
- **Poniżej \(n=1\)** nie istnieje żaden poziom. Bez podłogi i bez fotonu
  niosącego \(\hbar\) model przechodzi klasyczną katastrofę promienistą
  (\(\approx199\) ps), której prawdziwe pozytonium nie przechodzi.
- **Nie ma w modelu:** funkcji falowej, zasady Pauliego, poprawek
  radiacyjnych, rozpraszania Bhabhy, odpowiedzi detektora ani zjawisk w
  materiale (pełna lista w README, „Ograniczenia”).

## 9. Gdzie szukać

| temat | miejsce |
|---|---|
| prawa sił, pola, rozmycie dipola | `modules/electrodynamics.hpp` |
| integrator, historia retardowana | `modules/crem_engine.hpp` |
| przygotowanie stanu, pętla trajektorii | `modules/crem_trajectory.hpp` |
| tryb fotonowy, sufit, zamknięcie kanału | `modules/crem_collapse.hpp` |
| transport \(\mathbf L,\mathbf S_1,\mathbf S_2\) | `modules/secular_spin_orbit.hpp` |
| eksperyment 6, fotony anihilacji | `modules/contact_annihilation.hpp` |
| wszystkie przełączniki i ich uzasadnienia | `modules/configuration_panel.hpp` |
| pomiary, wyceny importów, historia | `README.md` |
| rozumowanie krok po kroku, korekty | plik audytu i `audits/data/` |
