# Semi-klasyczny model pozytonium (CREM)

Ten plik opisuje **model fizyczny** w jednym miejscu: z czego się składa, co
jest w nim klasyczne, co jest importem kwantowym i jak wypada wobec pomiarów.
Szczegóły pomiarów, historia zmian i wycofane twierdzenia są w `README.md` i w
pliku audytu `audits/2026-09-10-para-ortho-after-dipole-fix.txt` (z danymi w
`audits/data/section-NNN/`); tutaj są tylko odsyłacze do nich. Stan opisu:
audyt 353 (2026-10-06). Każda zmiana modelu trafia tutaj w tym samym
commicie co audyt.

## 1. Czym jest model

Model jest **klasyczną elektrodynamiką dwóch punktowych cząstek**, elektronu i
pozytonu, z których każda niesie ładunek i moment magnetyczny. Do tej
dynamiki dołożono **nazwaną, policzalną listę importów kwantowych**: miejsc,
w których wynik kwantowy wchodzi jako założenie. Model jest semi-klasyczny
dokładnie w tym sensie: trajektoria jest klasyczna, a dyskretność (kwant
działania fotonu, moment pędu fotonu, przygotowanie stanu, reguła anihilacji)
jest wpisana i oznaczona.

Cel jest jeden: sprawdzić, **jak daleko klasyczna elektrodynamika sięga w
odtwarzaniu wyników doświadczalnych pozytonium**, zanim zabraknie jej opisu.
Model ocenia się wobec **pomiarów**, nie wobec innych modeli (sekcja 8).
Żaden człon nie jest dopasowywany do danych; jedynymi wejściami są stałe
CODATA i zmierzone momenty magnetyczne. p-Ps i o-Ps są zawsze liczone i
raportowane osobno.

## 2. Składniki i skale

Każda cząstka ma masę \(m\), ładunek \(\pm e\) i moment magnetyczny
\(\mu=(g/2)\mu_B\) ze zmierzonym \(g\). Para ma masę zredukowaną
\(\mu_{\rm red}=m_e/2\). Wszystkie długości modelu wynikają z jednej skali,
\(\lambda_C\), i jednej liczby bezwymiarowej, \(\alpha\) (relacje są
zablokowane asercjami w `modules/particle_species.hpp`):

| skala | wyrażenie | wartość dla e⁺e⁻ |
|---|---|---|
| sprzężenie Coulomba | \(k=\alpha\hbar c\) | — |
| promień Bohra pary | \(a=a_{\rm pary}=2\lambda_C/\alpha\) | \(2a_0=105{,}8\) pm |
| energia wiązania \(n=1\) | \(B=k/(2a)\) | \(6{,}80285\) eV |
| promień „bariery” | \(r^*=\mu/(ec)=(g/2)\lambda_C/2\) | \(193{,}3\) fm |
| rozmycie pola dipola | \(\varepsilon=0{,}96682\,r^*\) | \(186{,}9\) fm |

\(r^*\) jest dokładnie promieniem, na którym energia dipol–dipol równa się
kulombowskiej: ich stosunek to \((r^*/r)^2\). Długość rozmycia \(\varepsilon\)
nie jest parametrem dopasowania: moment magnetyczny traktuje się jak pętlę
prądu o najmniejszym promieniu, przy którym ładunek \(e\) niesie \(\mu\) z
prędkością nie większą niż \(c\) (\(R=2r^*\)), a profil Plummera ma tę samą
wagę kontaktową co jednorodny dysk tej pętli (audyt 86).

**Zmienne orbity.** Orbita Keplera ma energię \(E=-B/n^2\) i moment pędu
\(L\); \(n=\sqrt{B/|E|}\) jest ciągłe. Główne działanie Delaunaya to
\(J=J_r+L=n\hbar\), a \(e^2=1-(L/n\hbar)^2\). Energia zależy **tylko od
\(J\)**: \(E=-B(\hbar/J)^2\), a częstość orbity to \(\omega=\partial E/\partial J
=2B/(\hbar n^3)\).

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
5. **Transport sekularny** (estymator `estimateCremCollapse`, między
   fotonami): wektory \(\mathbf L,\mathbf S_1,\mathbf S_2\) są przenoszone po
   orbicie uśrednionej z zachowaniem \(\mathbf J=\mathbf L+\mathbf S_1+\mathbf S_2\)
   (`modules/secular_spin_orbit.hpp`), z dwiema regułami bilansu:
   - **\(E_{\rm orb}+U\) stałe** (audyt 334): zmiana uśrednionej energii
     dipolowej \(U\) między checkpointami jest odejmowana od energii orbity
     (`CREM_NO_SPIN_ENERGY_EXCHANGE=1` wyłącza);
   - **\(|\mathbf L|\le L_{\rm kol}(E)=n\hbar\)** (audyt 337): nadmiar momentu
     pędu przechodzi do spinów (`CREM_NO_L_CAP=1` wyłącza).
   Sprzężenie \(\mathbf L\cdot\mathbf S\) obraca \(\mathbf L\), ale nie zmienia
   \(|\mathbf L|\); zmienia je tylko część tensorowa spin–spin:
   \(d|L|/dt=-\partial U/\partial g\) (\(g\) – argument perycentrum).
   Zmierzone przy \(n=1\): \(|L|\) drga z amplitudą \(<0{,}07\hbar\), a dryf
   wynosi \(\sim10^{-3}\hbar\)/ns (audyt 351) — **\(L\) nie miesza się** w
   czasie życia.
6. **Promieniowanie** (sekcja 4).

**Konwencja kanałów.** W p-Ps spiny są antyrównoległe, więc momenty
magnetyczne elektronu i pozytonu są **równoległe**
(\(|\boldsymbol\mu_1+\boldsymbol\mu_2|=2\mu\)). W o-Ps momenty są
antyrównoległe. Kanały różnią się na starcie tylko znakiem
\(\boldsymbol\mu_2\). Domyślnie kąt wzajemny jest losowany swobodnie (audyt
91), a kanał jest klasyfikacją wylosowanego układu. `--spin-quantization`
narzuca \(\cos=\pm1\). Waga kanału \(2\gamma\) to
\(w=(|\boldsymbol\mu_1+\boldsymbol\mu_2|/2\mu)^2\). W o-Ps (momenty
antyrównoległe, układ stabilny) \(w=0\) do \(10^{-8}\); w p-Ps (układ
niestabilny) \(w\) **oscyluje** między \(\approx0{,}62\) a \(1\) z okresem
\(\sim10\) ps, ze średnią czasową \(\approx0{,}87\) (audyty 344, 351). W QM
singlet ma \(w=1\) stale.

**Przygotowanie stanu związanego.** Start na \(a_n=n^2a\), domyślnie
\(n=1\), ostro kołowy: \(L=n\hbar\) z konstrukcji (`--level`).
`CREM_INITIAL_ANGULAR_MOMENTUM=J₀` daje ten sam poziom energii przy
\(L=J_0\,n\hbar\), czyli orbitę o \(e^2=1-J_0^2\) (audyt 301).
`--microcanonical-start` (audyt 326) losuje \(J_0^2\) jednostajnie z \([0,1]\):
to klasyczny zespół mikrokanoniczny przy energii \(E_n\), czyli odpowiednik
całej powłoki \(n\). W zmiennej \(L\) gęstość wynosi
\(p(L)=2L/(n\hbar)^2\) na \([0,n\hbar]\). Strumień losowy jest ten sam co
przy okręgu, więc oba starty da się porównywać parami.

**Start Langera** (audyt 353; przez `CREM_INITIAL_ANGULAR_MOMENTUM` z
\(J_0=(l+\tfrac12)/n\)): \(L=(l+\tfrac12)\hbar\), \(J_r=(n_r+\tfrac12)\hbar\),
\(n=n_r+l+1\) — kwantowanie półklasyczne WKB z indeksem Maslowa (import).
Reguły fotonu zmieniają \(L\) o \(\pm\hbar\) albo odbijają je
(\(L\to\hbar-L\)), więc część ułamkowa \(f=(L/\hbar)\bmod1\) przechodzi
\(f\to f\) lub \(f\to1-f\): sieć \(L\in\{\tfrac12,\tfrac32,\dots\}\hbar\) jest
**zachowywana**, a \(f=\tfrac12\) jest punktem stałym odbicia. Każda para
startująca z siatki Langera dochodzi do \(n=1\) z \(L=\hbar/2\)
(zmierzone: \(|L-\tfrac12|\le8\cdot10^{-3}\)). Dynamika tego \(L\) **nie
wybiera** (brak atraktora, audyt 351) — wybiera je przygotowanie.

## 4. Promieniowanie: dwa tryby

**Ciągły:** zredukowana siła Landaua–Lifshitza na każdą cząstkę
(`--radiation-reaction individual`); domyślny w trybie wizualnym i w
wiązkach. Orbita traci energię płynnie po obwiedni Larmora; nie ma fotonów
ani kwantu, a orbita spada aż do bariery (klasyczna katastrofa promienista).

**Fotonowy (`--radiation-reaction stochastic`; od audytu 341 domyślny w
statystycznych eksperymentach 1 i 2, czyli dla stanów związanych):** moc
E1 trafia do hazardu i jest wypłacana dyskretnymi fotonami, które zachowują
pęd (estymator sekularny `estimateCremCollapse`, `modules/crem_collapse.hpp`).
Reguły, w kolejności działania:

- **Hazard** (wyłącznie fotonowy od audytu 328; strata mierzonej orbity też
  trafia do hazardu, `CREM_CONTINUOUS_ORBIT_CREDIT=1` przywraca kredyt
  ciągły): \(\lambda=S(e)\,P_{E1}(a,e)/E_{\rm ref}\), gdzie \(E_{\rm ref}\) to
  energia fotonu z reguły niżej, a \(S(e)\) uwzględnia rozkład mocy orbity
  eliptycznej na harmoniczne (tablica interpolowana w \(e^2\), przy małym
  \(e\): \(S\approx1-2e^2\); audyt 340). Domyślnie foton pada, gdy skumulowany
  hazard osiągnie 1 (deterministycznie); `--emission poisson` losuje próg z
  \(\mathrm{Exp}(1)\). Moc \(P_{E1}\) liczy się po węzłach orbity z siłami
  opóźnionymi; równanie Keplera rozwiązuje Newton zabezpieczony przedziałem
  \(E-M\in[-e,e]\) (audyt 348; goły Newton rozbiegał się przy \(e\ge0{,}985\)).
- **Energia — reguła działania (domyślna od audytu 352;
  `CREM_NO_ACTION_PHOTON=1` przywraca \(k\hbar\omega\)).** Foton zabiera
  \(\hbar\) głównego działania \(J\). Ponieważ \(E\) zależy tylko od \(J\),
  \[\Delta E_{\rm orb}=\int_{J-\hbar}^{J}\omega\,dJ=E(n)-E(n-1)\]
  **przy dowolnym mimośrodzie**. Harmoniczna \(k\), losowana z widma dipola
  orbity Keplera (`eccentricOrbitHarmonicNumber`), wybiera
  \(\Delta n=k\), z podłogą \(n'=1\): \(\Delta=B[1/\max(1,n-k)^2-1/n^2]\).
  Odrzut: orbita traci dokładnie \(\Delta\), więc foton dostaje
  \(E_\gamma=\Delta-\Delta^2/(2W)\), gdzie \(W=Mc^2+\mu_{\rm red}E\).
  Klasyczne \(E_\gamma=\hbar\omega\) jest granicą tej reguły dla małego skoku.
- **Kierunek i skrętność:** kierunek z wzoru kątowego pola dalekiego,
  \(\propto1+\cos^2\theta\) względem \(\mathbf L\); skrętność \(h=\pm1\) z
  \(P(+)=(1+c)^2/(2(1+c^2))\).
- **Moment pędu (\(\Delta l=-1\), od audytu 330):** foton zabiera \(\hbar\)
  wzdłuż osi orbity, \(\mathbf L'=\mathbf L-\hbar\hat{\mathbf L}\), czyli
  \(|L'|=\big||L|-\hbar\big|\) (`CREM_DIRECTIONAL_PHOTON_SPIN=1` przywraca
  \(\hbar h\) wzdłuż kierunku lotu). Przy \(|L|\ge\hbar\) daje to
  \(\Delta J_r=0\) (działanie radialne zachowane, jak w kwantowej kaskadzie
  yrast). Przy \(|L|<\hbar\) reguła **odbija** \(L\) do \(\hbar-|L|\) — tego
  QM nie ma (z \(l=0\) dozwolone jest tylko \(\Delta l=+1\)); wariant testowy
  `CREM_DL_PLUS_BELOW_HBAR=1` daje wtedy \(L'=L+\hbar\) (audyt 348–349).
- **Sufit energii:** orbita z \(L'\) musi istnieć (\(e^2\ge0\)), czyli
  \(E'\ge-B(\hbar/L')^2\). Przy regule działania foton, którego sufit nie
  mieści, jest **odrzucany**, a nie przycinany (przycięcie złamałoby
  \(\Delta J=\hbar\)); przy \(n\le1\) każdy foton jest odrzucany.
- **Między fotonami:** dynamika zachowawcza i transport sekularny (sekcja 3).
- **Stan końcowy (`EmissionChannelClosed`):** przy regule działania kanał
  jest zamknięty, gdy \(n\le1\) albo gdy najmniejszy krok (\(k=1\)) nie
  zmieści \(L'\): \(L'>\max(1,n-1)\,\hbar\). Bez reguły działania warunek
  brzmi \(|L/\hbar-1|\ge n\) (audyt 308).

**Co daje reguła działania** (audyt 349, start mikrokanoniczny \(n=2\),
oba kanały identycznie):
- każdy foton \(2\to1\) ma \(E_\gamma=5{,}10212\) eV, a każda kaskada ląduje
  dokładnie na \(n=1\) z \(L'=|L_0-\hbar|\);
- **poziomy \(n<1\) nie istnieją**: z \(n=1\) nie da się zabrać \(\hbar\)
  działania, więc nie ma rodzin kołowych \(n\approx0{,}44\) ani eliptycznych
  z audytów 343–345, które były artefaktami reguły \(k\hbar\omega\);
- **start przy `--level 1` nie emituje nic** (stan końcowy od razu);
  kaskadę daje dopiero \(n\ge2\). Kaskada \(2\to1\) w programie głównym
  (start mikrokanoniczny, 4 pary) trwa średnio \(9{,}2\pm2{,}7\) ns w obu
  kanałach i potrzebuje `--crem-wallclock-budget-s` rzędu 1200 s — przy
  domyślnych 90 s wszystkie przebiegi są cenzurowane (audyt 352);
- dyskretność jest **siecią niezmienniczą, nie atraktorem**: całkowite \(n\)
  zostaje całkowite, ale model nie wybiera całkowitych \(n\) sam — robi to
  przygotowanie stanu;
- z \(\Delta l=+1\) przy \(|L|<\hbar\) stany \(L_0<\hbar\) przy \(n=2\) są
  zamknięte (odpowiednik metastabilnego \(2s\), bo E1 \(2s\to1s\) jest
  zabronione) i \(J_r\) jest zachowane do \(5\cdot10^{-5}\).

## 5. Importy kwantowe

Każdy import ma przełącznik albo jest nazwany jako nieusuwalny. Ich
historyczne ceny są w README („Sześć importów kwantowych i ich cena”).

| import | co wpisuje | domyślnie |
|---|---|---|
| promień i moment startowy | \(a_n=n^2a\), \(L=n\hbar\) (obraz Bohra) | **tak**; `--microcanonical-start` zastępuje \(L=n\hbar\) zespołem \(L^2\) jednostajnym |
| \(\hbar\) działania na foton | \(\Delta J=-\hbar\), więc \(\Delta E=E(n)-E(n-k)\) i \(n=1\) jako stan końcowy | **tak**, w trybie fotonowym (`CREM_NO_ACTION_PHOTON=1`: \(E_\gamma=k\hbar\omega\)) |
| \(\hbar\) momentu pędu fotonu | \(\mathbf L'=\mathbf L-\hbar\hat{\mathbf L}\) (\(\Delta l=-1\)) | **tak**, w trybie fotonowym |
| kwantyzacja spinu | \(\cos(\boldsymbol\mu_1,\boldsymbol\mu_2)=\pm1\) | nie (`--spin-quantization`) |
| tempo anihilacji | \(\sigma v=4\pi r_e^2c\), \(\varepsilon_{\rm OP}\), \(w\) (sekcja 6) | tak, w trybie fotonowym |
| podłoga / drabina Bohra | dawniej dwa importy (`--ground-state-floor`, `--bohr-photon-energy`) | zastąpione regułą działania; flagi zostały dla przebiegów z `CREM_NO_ACTION_PHOTON=1` |

**Jedno założenie skalowe.** Skala całego modelu bierze się z \(\hbar\):
przy ustalonym \(L=\hbar\) minimum \(L^2/(2\mu r^2)-k/r\) leży dokładnie w
\(a\) (argument Bohra), a reguła działania przenosi to samo \(\hbar\) na
skoki energii. Stan 1s ma kwantowo \(l=0\), a jego skalę daje energia
lokalizacji, której model nie ma.

## 6. Anihilacja

**Klasyczna dynamika nie zawiera anihilacji.** Model ma dwie jawnie
kwantowe recepty.

### 6.1. Tempo w estymatorze (audyt 342; tryb fotonowy)

\[\Gamma=\sigma v\,[\,w+(1-w)\,\varepsilon_{\rm OP}\,]\,\langle n\rangle_{\rm orb},\qquad
E[T]=\int_0^\infty S(t)\,dt,\quad S=e^{-\int\Gamma dt}\]

- \(\sigma v=4\pi r_e^2c\) (import; z \(|\psi_{1s}(0)|^2=1/(\pi a^3)\) daje
  \(T_0=1/(\sigma v|\psi_{1s}(0)|^2)=124{,}49\) ps);
- \(\varepsilon_{\rm OP}=4(\pi^2-9)\alpha/(9\pi)=1/1113{,}9\) — tłumienie
  \(3\gamma\) (Ore–Powell, wiodący rząd QED);
- \(\langle n\rangle_{\rm orb}\) — gęstość kontaktowa uśredniona po orbicie z
  jądrem Plummera \(3\varepsilon^2/(4\pi\rho^5)\).

**Ta recepta nie pasuje do stanu \(n=1\)** (audyty 343–347, 349d): orbita
\(n=1\) z \(L\sim\hbar\) prawie nie dochodzi do zera, więc
\(\langle n\rangle_{\rm orb}\) wynosi \(10^{-6}\)–\(10^{-3}|\psi_{1s}(0)|^2\), a
czasy życia wychodzą \(10^4\)–\(10^5\) razy za długie. \(|\psi_{1s}(0)|^2\)
odtwarza tylko wąskie \(L^*=0{,}178\hbar\), które jednak przesuwa się z
\(\varepsilon\) jak \(\varepsilon^{0{,}28}\), czyli jest dopasowaniem, a nie
przewidywaniem (audyt 347). W języku reguły działania kontakt to \(J=0\)
(„\(n=0\)”), a anihilacja to ostatni krok \(1\to0\), którego tempa estymator
jeszcze nie ma.

### 6.2. Krok \(1\to0\) z twierdzenia Quigga–Rosnera (audyt 350; analiza, nie w silniku)

**Kinematyka:** fotony anihilacyjne (\(2\gamma\)/\(3\gamma\)) zabierają cały
\(\mathbf J=\mathbf L+\mathbf S_1+\mathbf S_2\), więc reguła \(\Delta l=-1\)
tego kroku nie ogranicza; kanał wybiera \(w\).

**Tempo:** twierdzenie Quigga–Rosnera (QM, ścisłe dla stanów s):
\(|\psi(0)|^2=(\mu_{\rm red}/2\pi\hbar^2)\langle dV/dr\rangle\). Klasyczna
orbita ma \(\langle k/r^2\rangle=k\hbar/(n^3a^2L)\), więc
\[n_{\rm QR}=\frac{1}{2\pi a^3n^3}\,\frac{\hbar}{L},\qquad
\tau(L)=\frac{1}{\Gamma}=T_0\,\frac{2n^3L/\hbar}{w+(1-w)\varepsilon_{\rm OP}}.\]
Przy \(n=1\) i \(L=\hbar/2\) (Langer) \(n_{\rm QR}=|\psi_{1s}(0)|^2\). Wynik
nie zależy od \(\varepsilon\); importem jest tożsamość QR, nie
\(|\psi(0)|^2\). Przy \(n=2\) czynnik \(n^3=8\) daje
\(\tau_{2s}=8\tau_{1s}\) dla tego samego \(L\), tak jak
\(|\psi_{2s}(0)|^2=|\psi_{1s}(0)|^2/8\).

**Średnie po zespole.** Start mikrokanoniczny przy \(n=2\): w jednostkach
\(\hbar\) \(L_0=2j_0\), \(j_0^2\sim U(0,1)\), czyli \(p(L_0)=L_0/2\) na
\([0,2]\). Średni czas życia mieszaniny to
\(\langle\tau\rangle=\int p(L)\,\tau(L)\,dL\propto\langle L'\rangle\).

- **A′** (odbicie, \(L'=|L_0-1|\), wszystkie pary na \(n=1\)):
  \[\langle L'\rangle=\int_0^2|x-1|\,\frac{x}{2}\,dx
  =\int_0^1\frac{(1-x)x}{2}dx+\int_1^2\frac{(x-1)x}{2}dx=\frac1{12}+\frac5{12}=\frac12,\]
  więc \(\langle\tau\rangle_p=T_0=124{,}5\) ps (przy \(w=1\)) i
  \(\langle\tau\rangle_o=T_0\cdot1113{,}9=138{,}7\) ns.
- **C′** (\(\Delta l=+1\) przy \(L_0<1\)): zespół dzieli się na dwie
  populacje. \(P(L_0<1)=\int_0^1x/2\,dx=1/4\) zostaje przy \(n=2\) (\(2s\));
  \(3/4\) schodzi na \(n=1\) z \(L'=L_0-1\):
  \[\langle L'\rangle=\frac{\int_1^2(x-1)x\,dx}{\int_1^2x\,dx}=\frac{5/6}{3/2}=\frac59
  \;\Rightarrow\;138{,}3\ \text{ps},\ 154\ \text{ns}.\]
  Populacja \(2s\) ma \(\langle L_0\rangle=\int_0^1x^2dx/\int_0^1x\,dx=2/3\),
  więc \(\langle\tau\rangle_{2s}=T_0\cdot16\cdot\tfrac23=10{,}7\,T_0\).

**Wagi.** A′ i C′ to **alternatywne reguły**, nie populacje — nie uśrednia
się ich. Wewnątrz wariantu wagi zespołu startowego są już w
\(\langle L'\rangle\); w C′ jeden czas życia całego zespołu wymagałby wag
\(3/4\) (\(n=1\)) i \(1/4\) (\(n=2\)). Pomiary \(125{,}14\) ps i \(142{,}04\)
ns dotyczą stanu podstawowego, więc porównuje się populację \(n=1\). p-Ps i
o-Ps się nie waży: to osobne składowe widma czasów życia.

**Zbieżność, nie wyprowadzenie.** \(\langle L'\rangle=1/2\) w A′ jest
wartością Langera, ale wynika z dwóch arbitralnych wyborów (zespół
mikrokanoniczny przy \(n=2\), reguła odbicia). Mieszanina **nie jest
jednowykładnicza**: gęstość \(L'\) przy zerze jest skończona, więc
\(\langle1/L'\rangle\) (tempo początkowe) rozbiega się logarytmicznie, a ogon
\(S(t)\) ma czas 1,5–1,9× średniej. Zmierzone rozpady są jednowykładnicze,
a \(L\) nie miesza się w czasie życia (audyt 351).

**Zmierzone na trajektoriach** (stany końcowe 349, bez bariery):

| wariant | p-Ps \(\langle\tau\rangle\) | o-Ps \(\langle\tau\rangle\) | \(\tau_o/\tau_p\) |
|---|---|---|---|
| A′ | 160,8 ps (1,28×; \(\langle L\rangle=0{,}55\), \(w\approx0{,}92\)) | 150,0 ns (1,06×) | 933 |
| C′ | 187,5 ps (1,50×) | 177,6 ns (1,25×) | 947 |

**Ze startem Langera** (audyt 353; \(2p\), \(3p\), \(3d\) i \(2s\) w A′;
26 par na kanał na \(n=1\), wszystkie z \(L=\hbar/2\)):

| kanał | \(\langle\tau\rangle\) | mediana | \(\langle1/\Gamma\rangle\langle\Gamma\rangle\) | pomiar |
|---|---|---|---|---|
| p-Ps | 136,3 ps | 127,65 ps | 1,0225 (oscylacja \(w\)) | 125,14 ps |
| o-Ps | 138,65 ns | 138,67 ns | 1,00002 (jednowykładniczo) | 142,04 ns |

o-Ps różni się od pomiaru o \(-2{,}4\%\), czyli tyle, ile
\(1/\varepsilon_{\rm OP}=1113{,}9\) od \(1135\). Rozrzut p-Ps pochodzi
wyłącznie z oscylacji \(w\) (min 0,59, średnia 0,934), nie z \(L\). Stan
\(2s\) w C′ (zamknięty przy \(n=2\)) ma \(\tau_{\rm QR}=995{,}95\) ps (p-Ps)
i \(1{,}109\) µs (o-Ps), czyli \(8\times\) więcej, jak
\(|\psi_{2s}(0)|^2=|\psi_{1s}(0)|^2/8\).

### 6.3. Eksperyment 6 (zdarzenie anihilacji)

`--mode statistical --phenomenon 6` (audyty 316–320): anihilacja przy
pierwszym wejściu w \(r\le r^*\); \(P(2\gamma)=w/(w+(1-w)\varepsilon_{\rm OP})\);
fotony w układzie spoczynkowym pary — \(2\gamma\) wzdłuż izotropowej osi,
\(3\gamma\) z dokładnym rozkładem energii Ore–Powella i płaszczyzną
\(dN/d\cos\theta_n\propto1-\tfrac13\cos^2\theta_n\) względem
\(\mathbf S_1+\mathbf S_2\) (drzewowa amplituda QED). Bezwzględne czasy nie
mają tu sensu fizycznego; obserwablą jest stosunek.

## 7. Przełączniki testowe i pomiarowe (domyślnie wyłączone)

| przełącznik | działanie | audyt |
|---|---|---|
| `CREM_NO_ACTION_PHOTON=1` | energia fotonu \(k\hbar\omega\) zamiast reguły działania | 352 |
| `CREM_CLOSE_BELOW_HBAR=1` | kanał zamknięty przy \(|L|<\hbar\) | 346 |
| `CREM_DL_PLUS_BELOW_HBAR=1` | \(\Delta l=+1\) przy \(|L|<\hbar\) | 348–349 |
| `CREM_HOLD_AFTER_CLOSURE=<s>` | po zamknięciu trzymaj parę bez emisji, drukuj \(L\), \(n\), \(w\) | 351 |
| `CREM_NO_L_CAP=1`, `CREM_NO_SPIN_ENERGY_EXCHANGE=1` | wyłączają reguły bilansu z sekcji 3 | 337, 334 |
| `CREM_DIRECTIONAL_PHOTON_SPIN=1` | \(\hbar h\) wzdłuż kierunku fotonu zamiast osi | 330 |
| `CREM_CONTINUOUS_ORBIT_CREDIT=1` | kredyt ciągły straty mierzonej orbity | 328 |
| `CREM_EMISSION_REACH=1` | wydruk każdego fotonu (\(E\), \(L\), \(e^2\), \(k\), \(n\)) | — |

Pełna lista z uzasadnieniami: `modules/configuration_panel.hpp`.

## 8. Wobec pomiarów

Pomiary odniesienia: \(\tau_{\rm para}=125{,}14\) ps (Al-Ramadhan i Gidley,
PRL 72, 1632 (1994)); \(\tau_{\rm orto}=142{,}04\) ns (Vallery i in., PRL 90,
203402 (2003)) — średnie czasy życia, nie półokresy; odstęp 1S–2S
\(5{,}10179\) eV (Fee i in., PRL 70, 1397 (1993)); HFS \(203{,}3942\) GHz
(Ishida i in., PLB 734, 338 (2014)); \(m_ec^2=510{,}99895\) keV (CODATA 2018).

| obserwabla | model | pomiar | ocena |
|---|---|---|---|
| linia \(2\to1\) | \(5{,}10212\) eV w 100% kaskad (reguła działania) | \(5{,}10179\) eV | \(+0{,}33\) meV (\(6{,}5\cdot10^{-5}\)) — rząd \(\alpha^2B\), poprawki subtelne i QED, których drabina \(-B/n^2\) nie ma (audyt 349) |
| stan końcowy kaskady | dokładnie \(n=1\) | \(1s\) | tak; przy starcie mikrokanonicznym \(L=|L_0-\hbar|\) rozrzucone, przy starcie Langera \(L=\hbar/2\) w każdej parze |
| linie \(3\to1\), \(3\to2\) | \(6{,}04698\), \(0{,}94484\) eV (kaskada yrast \(3d\to2p\to1s\), audyt 353) | — | pomiaru nie ma w repozytorium |
| \(\tau\) z kroku \(1\to0\) (QR, start Langera) | p-Ps 127,65 ps (mediana), o-Ps 138,65 ns, o-Ps jednowykładniczo | 125,14 ps / 142,04 ns | \(+2{,}0\%\) / \(-2{,}4\%\); p-Ps średnia \(+8{,}9\%\) przez oscylację \(w\) |
| \(2s\) metastabilne | w C′ brak fotonu z \(2s\) (E1 zabronione) | QM: tak | zgodne; w A′ \(2s\) schodzi przez odbicie (niezgodne z QM) |
| \(\tau\) z recepty 6.1 przy \(n=1\) | \(10^4\)–\(10^5\times\) za długie | 125,14 ps / 142,04 ns | recepta nie pasuje do \(n=1\) (sekcja 6.1) |
| \(\tau\) z kroku \(1\to0\) (QR, A′) | p-Ps 160,8 ps, o-Ps 150,0 ns | 125,14 ps / 142,04 ns | właściwa skala bez \(\varepsilon\); rozpad niejednowykładniczy; p-Ps wydłużone przez \(w\approx0{,}9\) |
| \(\tau_{\rm orto}/\tau_{\rm para}\) | \(1/\varepsilon_{\rm OP}=1113{,}9\) przy \(w_p=1\); 933–953 z oscylacją \(w\) | \(1135{,}0\) | reguła wyboru \(\varepsilon_{\rm OP}\) jest importem; własny wkład modelu to \(w\) |
| \(\tau_{\rm orto}/\tau_{\rm para}\), eksperyment 6, \(n=1\) | \(27\,215\pm42\%\) = geometria \(24{,}4\) × reguła \(1113\) | \(1135{,}0\) | geometria powinna wynosić \(\approx1\) |
| rozszczepienie nadsubtelne 1s | koło \(L=\hbar\): \(6{,}7\cdot10^{-7}\) pomiaru; zespół mikrokanoniczny: \(290{,}5\) GHz | \(203{,}3942\) GHz | znak poprawny, człon kontaktowy; zależy od tego, który stan klasyczny gra \(1s\) (audyt 324) |
| rozpraszanie | kształt Rutherforda, jedna normalizacja | — | klasyczne |

**Wynik negatywny jest wynikiem.** Różnica para/orto w przyrodzie pochodzi z
reguły wyboru \(2\gamma/3\gamma\) (parzystość ładunkowa), której klasyczna
dynamika nie zawiera; w modelu wchodzi przez \(\varepsilon_{\rm OP}\) i
\(w\). Skala bezwzględna wymaga gęstości w zerze, której klasyczna orbita
nie ma; tożsamość QR daje ją z \(\langle dV/dr\rangle\), ale jedno \(L\) na
parę daje rozrzut temp.

## 9. Granice dziedziny i otwarte problemy

- **Bariera Comptona** \(r^*\): poniżej niej klasyczna elektrodynamika
  punktowa przestaje obowiązywać, a przebieg się kończy.
- **Granica retardacji:** okres / czas przelotu światła \(\le150\) — margines
  przybliżenia Darwina i całkowania, nie fizyka.
- **Poniżej \(n=1\):** w trybie fotonowym z regułą działania nie ma
  poziomów; w trybie ciągłym orbita przechodzi klasyczną katastrofę
  promienistą (\(\approx199\) ps), której prawdziwe pozytonium nie przechodzi.
- **Otwarte:** (1) jedno \(L=\hbar/2\) na \(n=1\) daje tylko przygotowanie
  Langera (import WKB) — dynamika go nie wybiera; (2) oscylacja \(w\) w p-Ps (QM: \(w=1\));
  (3) tempo kroku \(1\to0\) w samym silniku; (4) reguła przy \(|L|<\hbar\)
  (odbicie vs \(\Delta l=+1\)); (5) czas kaskady \(2\to1\) wobec zmierzonego
  czasu życia \(2P\) (brak źródła pomiarowego w repozytorium).
- **Nie ma w modelu:** funkcji falowej, zasady Pauliego, poprawek
  radiacyjnych, rozpraszania Bhabhy, odpowiedzi detektora ani zjawisk w
  materiale (pełna lista w README, „Ograniczenia”).

## 10. Gdzie szukać

| temat | miejsce |
|---|---|
| prawa sił, pola, rozmycie dipola | `modules/electrodynamics.hpp` |
| integrator, historia retardowana | `modules/crem_engine.hpp` |
| przygotowanie stanu, pętla trajektorii | `modules/crem_trajectory.hpp` |
| tryb fotonowy, reguła działania, sufit, zamknięcie, tempo anihilacji | `modules/crem_collapse.hpp` |
| transport \(\mathbf L,\mathbf S_1,\mathbf S_2\) | `modules/secular_spin_orbit.hpp` |
| eksperyment 6, fotony anihilacji | `modules/contact_annihilation.hpp` |
| wszystkie przełączniki i ich uzasadnienia | `modules/configuration_panel.hpp` |
| pomiary, wyceny importów, historia | `README.md` |
| rozumowanie krok po kroku, korekty | plik audytu i `audits/data/` |
