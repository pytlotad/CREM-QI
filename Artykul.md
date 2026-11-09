# Klasyczna elektrodynamika pozytonium i wodoru z czterema zasadami kwantowymi

**Tadeusz Sławomir Pytlos**

*Oprogramowanie i pełny zapis audytów: CREM-QI, <https://github.com/pytlotad/CREM-QI>
(model fizyczny: `Model.md`; audyty: `audits/2026-09-10-para-ortho-after-dipole-fix.txt`, dane:
`audits/data/section-NNN/`). Numery w nawiasach kwadratowych [NNN] odsyłają do audytów.*

---

## Streszczenie

Badamy, jak daleko relatywistyczna klasyczna elektrodynamika dwóch naładowanych cząstek
z momentami magnetycznymi odtwarza pomiary pozytonium (Ps) i wodoru, jeśli treść kwantową
ograniczyć do jawnie zadeklarowanego minimum. Dynamika obejmuje retardowane pola
Liénarda–Wiecherta, oddziaływanie dipol–dipol z członem kontaktowym, kowariantne pędy
spinowe, precesję BMT i promieniowanie rozłożone na mody klasycznego pola orbity. Treść
kwantowa to cztery zasady w kolejności postulatów mechaniki kwantowej: **(I) wartości**
(kwantowanie EBK z fazą Maslowa), **(II) elementy macierzowe** (korespondencja Borna w
punkcie środkowym z jednym członem porządku), **(III) prawdopodobieństwa** (reguła Borna)
i **(IV) punkt $r=0$** (gęstość kontaktu oraz liczby QED z literatury). Żaden parametr nie
jest dopasowany do danych; wejściami są $\alpha$, masy, zmierzone momenty magnetyczne i
czynniki $g$. Model daje czasy życia p-Ps i o-Ps do $1{,}5\cdot10^{-4}$ i $6\cdot10^{-5}$,
rozszczepienie nadsubtelne 1s Ps do $6\cdot10^{-5}$, strukturę $n=2$ Ps do 5 MHz (od importu przesunięcia Lamba Ps [720]), przesunięcie Lamba wodoru do 0,5 MHz i 1S–2S wodoru do 4 MHz (import [721]),
strukturę subtelną i nadsubtelną $2P$ wodoru oraz czynnik $g$ elektronu związanego
w granicach niepewności pomiarów. Granice są nazwane: tempo E1 stanów $p$ (2p wodoru
$+4{,}7\,\%$, $\sim19\sigma$; od [734] tempo E1 z fali EBK: $1{,}608$ ns, $+2{,}1\sigma$), przesunięcie Lamba tylko jako import liczb QED i skład stanów własnych
w polu pośrednim, którego żadna z sześciu sprawdzonych dróg klasycznych nie odtwarza.

---

## 1. Wstęp

Klasyczny atom wodoru jest niestabilny. Promieniuje i kolapsuje w czasie rzędu
$10^{-11}$ s, nie ma dyskretnego widma ani reguł wyboru. Znane jest też, gdzie klasyce
brakuje treści: w fazie de Broglie'a, która daje warunek kwantowania, w
nieprzemienności, która daje elementy macierzowe, i w regule Borna. Mniej znane jest,
**ile z obserwowanych liczb** daje sama klasyczna elektrodynamika, gdy te braki uzupełni
się możliwie oszczędnie i każde uzupełnienie nazwie.

Tę ilość mierzymy. Model łączy pełną relatywistyczną dynamikę dwóch ciał z momentami
magnetycznymi z czterema zasadami kwantowymi. Każdą sprawdzamy pod kątem tego, czy da się
ją wyprowadzić z dynamiki. Wynikiem są zarówno liczby zgodne z pomiarem, jak i wyniki
negatywne: miejsca, w których klasyka zawodzi w sposób mierzalny, z oszacowaną skalą.
Porównujemy się z pomiarami. Wyjątkiem są tempa E1 wodoru: wzór Gordona jest tam ścisły
do $\alpha^2$, a stare pomiary metodą folii odbiegają od niego o $\sim9\,\%$ [616].

---

## 2. Dynamika klasyczna

### 2.1. Skale

Cząstka $i$ ma masę $m_i$, ładunek $q_i$ i moment magnetyczny
$\boldsymbol\mu_i=(g_i/2)\,(q_i/2m_i)\,\mathbf S_i$ o długości ze zmierzonego $g_i$.
Dla pary o masie zredukowanej $\mu$:

$$
k=\alpha\hbar c,\qquad a=\frac{\hbar^2}{\mu k},\qquad B=\frac{k}{2a},\qquad
\omega_n=\frac{2B}{\hbar n^3},\qquad E=-\frac{B\hbar^2}{J^2},
$$

gdzie $J=J_r+L$ jest głównym działaniem Delaunaya, a mimośród spełnia
$e^2=1-(L/J)^2$. Moment magnetyczny jest pętlą prądu o najmniejszym promieniu
$R=\hbar/(m c)$, przy którym ładunek $e$ niesie $\mu$ z prędkością nie większą niż $c$.
Pole dipola jest rozmyte profilem Plummera o długości $\varepsilon=0{,}967\,r^*$, gdzie
$r^*=\mu/(ec)$ (dla e⁺e⁻ $r^*=193$ fm, $a=105{,}8$ pm). Długość $\varepsilon$ nie jest
dopasowana: profil ma tę samą wagę kontaktową co jednorodny dysk pętli.

### 2.2. Ruch ładunków

Pędy są relatywistyczne, $\mathbf p_i=\gamma_i m_i\mathbf v_i$. Pole cząstki $j$ to pole
Liénarda–Wiecherta w chwili opóźnionej:

$$
\mathbf E_j=\frac{q_j}{4\pi\varepsilon_0}\left[\frac{(\hat{\mathbf n}-\boldsymbol\beta)(1-\beta^2)}{\kappa^3R^2}
+\frac{\hat{\mathbf n}\times\big((\hat{\mathbf n}-\boldsymbol\beta)\times\dot{\boldsymbol\beta}\big)}{c\,\kappa^3R}\right]_{\rm ret},
\qquad \mathbf B_j=\frac{\hat{\mathbf n}\times\mathbf E_j}{c},\qquad \kappa=1-\hat{\mathbf n}\cdot\boldsymbol\beta ,
$$

a siła to $\mathbf F_i=q_i(\mathbf E_j+\mathbf v_i\times\mathbf B_j)$. Wszystko jest
dokładne we wszystkich rzędach $v/c$.

### 2.3. Ruch i siły momentów

Moment partnera daje rozmyte pole dipola z członem magnetyzacji, tak że
$\nabla\cdot\mathbf B=0$; daleko

$$
U_{dd}=\frac{\mu_0}{4\pi r^3}\Big[\boldsymbol\mu_1\cdot\boldsymbol\mu_2-3(\boldsymbol\mu_1\cdot\hat{\mathbf r})(\boldsymbol\mu_2\cdot\hat{\mathbf r})\Big]
-\frac{2\mu_0}{3}\,\boldsymbol\mu_1\cdot\boldsymbol\mu_2\,\delta^3(\mathbf r).
$$

Siła na moment zawiera gradient energii oraz tempo pędu ukrytego. W postaci kowariantnej
[570] pędy spinowe to $p'^{\mu}_{\rm hid}=\mathbf m_0\times\mathbf E'/c^2$ oraz człon
typu Tulczyjewa $p'_T=q(\mathbf E'\times\mathbf S_0)/(mc^2)$, gdzie
$\mathbf E'=\gamma(\mathbf E+\mathbf v\times\mathbf B)-\tfrac{\gamma^2}{\gamma+1}(\mathbf v\cdot\mathbf E)\mathbf v/c^2$
jest polem w układzie spoczynkowym. Pędy te są transformowane jak czterowektory, a siła
to $-dp^\mu/d\tau$ rzutowane prostopadle do $u$. Sektor dipolowy jest kowariantny do
$4\cdot10^{-5}$, a moment siły zgodny z precesją do $6\cdot10^{-4}$.

Spiny precesują zgodnie z równaniem BMT:

$$
\frac{d\mathbf S}{dt}=\frac{q}{m}\,\mathbf S\times\left[\Big(\frac g2-1+\frac1\gamma\Big)\mathbf B
-\Big(\frac g2-1\Big)\frac{\gamma}{\gamma+1}(\boldsymbol\beta\cdot\mathbf B)\boldsymbol\beta
-\Big(\frac g2-\frac{\gamma}{\gamma+1}\Big)\frac{\boldsymbol\beta\times\mathbf E}{c}\right].
$$

W estymatorze sekularnym (torze produkcyjnym) ruch jest uśredniany po orbicie. Precesja
zgadza się wtedy człon po członie z klasycznym hamiltonianem Breita–Pauliego rzędu
$1/c^2$ [467], na przykład dla spinu–orbity

$$
\boldsymbol\Omega_{SO,i}=\frac{k}{c^2}\left[\frac{g_i-1}{2m_i^2}+\frac{g_i}{2m_1m_2}\right]\langle r^{-3}\rangle\,\mathbf L .
$$

Na ścieżce mechanicznej siły ładunek–moment z sumą chwilową wynikają z dwóch lagranżjanów liniowych
w prędkości [706]:

$$
L_{q\mu}=q(\mathbf v_q-\mathbf v_\mu)\cdot\mathbf A_\mu,\qquad
L_T=-\boldsymbol\omega_T\cdot\mathbf S_p=\frac{k\,q_oq_p}{2m_pc^2}\,
\mathbf v_p\cdot\frac{\mathbf d\times\mathbf S_p}{d^3},
$$

gdzie $\mathbf A_\mu=(\mu_0/4\pi)\,\boldsymbol\mu\times\mathbf d/d^3$, a $L_T$ to odrzut precesji
Thomasa z wagami Barkera–O'Connella. Oba lagranżjany wchodzą do pędu kanonicznego, a nie do energii.
Energia ładunek–moment $U_{qd}=-\mathbf p\cdot\mathbf E$ jest jednorodna stopnia 1 w prędkościach i
do energii zachowanej sumy chwilowej nie należy [145, 708]. W sumie opóźnionej pędy spinowe są
czterowektorami [570]. Ta suma liczy prędkość mechaniczną bez nich, więc $U_{qd}$ wraca w jej energii
zachowanej [708]. Obie postacie dają tę samą orbitę wiekową: precesja peryapsu wodoru i antywodoru
zgadza się do $4\cdot10^{-4}$ [707]. Poziom $n=\sqrt{Ry/|E|}$, od którego zależą tempa i
okno drabiny ścieżki mechanicznej, silnik czyta z energii zachowanej całkowanej sumy: z $U_{qd}$ o wadze 0
przy sumie chwilowej i 1 przy opóźnionej, jak w księdze energii. Przebiegi regresji są przy tym identyczne
bajt w bajt, a walidacja daje 72/72 [719].

### 2.4. Promieniowanie

Dipol pary $\mathbf d=q_{\rm ef}\mathbf r$ rozkłada się na harmoniczne orbity Keplera:

$$
x+iy=\sum_{s\ge1}\Big(c_{s,+}\,e^{-is\Omega t}+c_{s,-}\,e^{+is\Omega t}\Big),\qquad
P_{s,\pm}=\frac{q_{\rm ef}^2(s\Omega)^4|c_{s,\pm}|^2}{6\pi\varepsilon_0c^3}.
$$

Mod współbieżny ($+$) niesie $L_z/E=+1/(s\Omega)$, a przeciwbieżny ($-$) $-1/(s\Omega)$.
Uśredniony po orbicie moment siły reakcji promieniowania spełnia w każdym punkcie
$\mathbf r\times\dot{\mathbf a}=-\mathbf L/r^3$ [597]. Moment pędu nigdy więc nie zmienia
kierunku, a mod współbieżny niesie 95–98 % momentu siły [598]. Klasycznie orbita
promieniuje w sposób ciągły i kolapsuje. Dyskretność wprowadzają dopiero zasady z
sekcji 3.

### 2.5. Tryb fotonowy

W trybie produkcyjnym ciągła reakcja Landaua–Lifszyca jest zastąpiona procesem
punktowym. Orbita porusza się zachowawczo, a foton jest emitowany w chwili, gdy
skumulowany hazard przekroczy próg $\mathrm{Exp}(1)$:

$$
\lambda(t)=\sum_{s,\pm}\frac{P_{s,\pm}}{\hbar\,s\Omega},\qquad
\int_{t_0}^{t_\gamma}\lambda\,dt=\xi,\quad \xi\sim\mathrm{Exp}(1).
$$

Mod $(s,m)$ jest losowany z wagą proporcjonalną do liczby jego kwantów. Emisja zabiera
$\Delta n=\min(s,\lfloor n-1\rfloor)$ kwantów głównego działania $J$ i $\hbar$ momentu
pędu wzdłuż osi orbity. Ponieważ $E$ zależy tylko od $J$, energia orbity zmienia się o
całkę wzdłuż krzywej działania, przy dowolnym mimośrodzie (reguła działania):

$$
\Delta=\int_{J-\Delta n\,\hbar}^{J}\omega(J')\,dJ'=B\Big[\frac{1}{(n-\Delta n)^2}-\frac1{n^2}\Big],\qquad
\mathbf L'=\mathbf L-\hbar\hat{\mathbf L}.
$$

Przy $|\mathbf L|\ge\hbar$ działanie radialne jest zachowane, $\Delta J_r=0$, jak w
kwantowej kaskadzie yrast. Orbita traci dokładnie $\Delta$, a odrzut pary o energii
$W=Mc^2+\mu E$ zmniejsza energię fotonu:

$$
E_\gamma=\Delta-\frac{\Delta^2}{2W},\qquad \mathbf P_{\rm odrz}=-\frac{E_\gamma}{c}\,\hat{\mathbf k}.
$$

Klasyczne $E_\gamma=\hbar\omega$ jest granicą tej reguły dla małego skoku. Kierunek
fotonu względem $\hat{\mathbf L}$ ma rozkład pola dalekiego obracającego się dipola:

$$
\frac{dP}{d\Omega}\propto1+\cos^2\theta,\qquad
P(h=+1)=\frac{(1+\cos\theta)^2}{2(1+\cos^2\theta)} .
$$

Orbita końcowa musi istnieć ($e^2\ge0$, czyli $E'\ge-B(\hbar/L')^2$). Foton, który nie
spełnia tego warunku, jest odrzucany, a nie przycinany, bo przycięcie złamałoby
całkowite $\Delta J$. Kanał emisji jest zamknięty, gdy $n<2$ albo gdy $L'>(n-1)\hbar$.
Przejścia $\Delta l=+1$ i wychwyt działają tak samo, z modem przeciwbieżnym.

Foton spełnia algebrę spinu 1 [M-WB] w dwóch miejscach. Jego polaryzacja
$\boldsymbol\varepsilon$, poprzeczna część $(\mathbf e_1+im\,\mathbf e_2)/\sqrt2$, jest stanem spinu 1 [644]. Zabierany
moment pędu $m\hbar$ na osi $\hat{\mathbf L}$ (w polu na osi $\hat{\mathbf B}$, z $m=\Delta m$ [504]) to rzut $m_S\hbar$.
$\hbar$ jest rzutem, nie długością. Odpowiednikiem $S^2=2\hbar^2$ w języku modelu jest długość reguły Q
$\tfrac32\hbar$. Ta długość działa przy wyborze kanału (trójkąt z $L_{ph}=\tfrac32\hbar$ dla E1 i $\tfrac52\hbar$ dla E2,
zero różnic od QM [598]), ale nie w rozgałęzieniu $J$. Tam foton jako wektor na stożku trójkąta dopuszcza $\Delta J=2$
i odbiega od 6j o 0,15–0,26, a sam rzut wzdłuż $\hat{\mathbf L}$ daje 0,002 (wodór) i 0,028 (Ps) [611, 725]. Rozgałęzienie nadsubtelne $F'$ stanu 1s wodoru po fotonie $2p\to1s$ też nie wynika z geometrii klasycznych spinów (próg kąta daje 2/64 zamiast ~1/4 [739]); od [740] $F'$ jest losowane z wagą 6j, jak $J'$ trypletu Ps [631], . Od [741] basen $F$ stanu 2p jest etykietą trzymaną przez księgę poziomu (klasyczna precesja przenosi kierunki spinów przez granicę basenu), i w 64 kaskadach wychodzi 0/8, 5/19, 8/16, 0/21 wobec 0, 1/3, 2/3, 0, a każda linia jest składową $F$.

Dynamika z §2.2–2.5 jest niezmiennicza względem C i P z dokładnością do bitu. Obraz stanu
przygotowanego pod C (antycząstki, $\boldsymbol\mu\to-\boldsymbol\mu$, $\mathbf B\to-\mathbf B$) i pod P
($\mathbf r,\mathbf v\to-\mathbf r,-\mathbf v$) daje identyczny przebieg, także dla antywodoru. Losowania
kierunku fotonu i stanu własnego pola są symetryczne w rozkładzie, choć nie punktowo [623].

---

## 3. Cztery zasady kwantowe

Treść kwantowa modelu to cztery zasady. Każdy element jest zapisany jednym wierszem,
wraz z audytem, który go ustalił albo sprawdził.

### Zasada I — wartości (kwantowanie EBK z fazą Maslowa)

Treść zasady: fala niesie fazę $S/\hbar$, $S=\oint\mathbf p\cdot d\mathbf q$ (z potencjałem), i skok
$-\pi/2$ w każdym punkcie zwrotnym. Suma fal po kolejnych obiegach, $\sum_N e^{iN\Phi}$, znosi się,
chyba że $\Phi=2\pi k$: wartości działania to maksima tej interferencji (Keller; Gutzwiller).
Importem jest sama faza i skok Maslowa; liczba całkowita i połówka są ich wnioskiem (§5.1) [627].
Audyty 723–730 policzyły to na orbitach i spinach: kaustyki $r$ dają $J_r=(n_r+\tfrac12)\hbar$, kaustyki $\vartheta$
nachylonej płaszczyzny — długość Langera $L=(l+\tfrac12)\hbar$, precesja węzła — $L_z=m\hbar$, faza spinu — trójkąt $J$,
a w reprezentacji Schwingera (spin = para oscylatorów) ten sam skok daje długość spinu $S+\tfrac12$ i parzystość trójkąta
(singlet $\tfrac12$, tryplet $\tfrac32$) [723, 726–730]. Fazę zasady I da się wyprowadzić. Równanie pola ruchu względnego samo daje skok Maslowa: jego widmo (różnice skończone, bez wzoru na
$E_n$) przy klasycznym działaniu radialnym daje $J_r=n_r+\tfrac12$ do $4\cdot10^{-8}$, a w 3D wymaga długości Langera $L=l+\tfrac12$
(z $\sqrt{l(l+1)}$ odchylenie 0,035–0,174); tam, gdzie EBK nie jest ścisłe, reszta to $10^{-3}$ i maleje z $n_r$ [743, 748]. Przez retardowane
pole pamięci ta faza przechodzi dynamicznie do źródła: poślizg fazy źródła wynosi $-3{,}103$ rad na okres radialny (do $-\pi$ przy
silniejszym dostrajaniu) [744, 748]; jedynym wyjątkiem jest orbita radialna $L=0$, na której mody $\pm m$ tworzą wzdłuż toru falę stojącą i wygrywają z własnym modem [749]. Gdy orbita traci energię klasycznie (reakcja promienista), pole samo przechodzi między modami, ale nie przenosi orbity na nowy torus EBK — „lądowanie” po skoku zostaje importem [750]. Sprzężone z polem EM w przybliżeniu dipolowym pole promieniuje tylko dudnieniem dwóch modów: częstość linii to różnica ich energii (10,19872 eV dla 2p → 1s wodoru, bez importu częstości), tempo to reguła Fermiego ($\tau(2p)$ 1,59620 ns wobec pomiaru 1,600(4) ns), a ziarno emisji daje źródło na orbicie; zanik jest jednak logistyczny, nie wykładniczy, więc statystyka skoków (zasada III) zostaje [751]; w silniku energia fotonu jest odtąd zapisana jako ta różnica energii modów (zmiana opisowa, liczby bez zmian) [752]. Przeniesienie losowości do klasycznego pola próżni (widmo $\hbar\omega/2$, każdy bieg deterministyczny) nie działa: zanik z 2p zaczyna się od połowy tempa, stan 1s jest wzbudzany do tej samej równowagi (~0,66), a chwile emisji nie są wykładnicze — klasyczny szum działa jak promieniowanie cieplne, nie jak próżnia QED [753]. Z polem jako operatorem (trajektorie bohmowskie konfiguracji pola) statystyka zespołu jest kwantowa, a stan podstawowy stabilny, ale pojedynczy bieg nie zawiera skoku: prawdopodobieństwo emisji rośnie w nim płynnie, a chwila emisji pojawia się dopiero przy detektorze [754]. Z idealnym detektorem pojedynczy bieg ma skok (gałęzie się rozdzielają), statystyka rejestracji jest kwantowa, a chwila skoku jest funkcją początkowej konfiguracji wskaźnika i pola; ogon chwil jest bez pamięci — zasada III jest w tym obrazie skutkiem pomiaru, a nie dynamiki CREM [755]. W silniku foton zabiera całkowite kwanty działania i zachowuje części ułamkowe $n$ i $L$, więc sam nie wybiera torusa; robi to wyłącznie zatrzask fali przed fotonem — zasady I i III są rozdzielnymi importami [756]. Zatrzask nie wynika też z pola z przejęciem fazy: źródło kopiujące fazę pola pamięci nie czuje niedopasowania fazy, więc praca pola nad orbitą ma stały znak i nie ma punktu stałego przy EBK — przejęcie fazy Maslowa i wybór torusa wykluczają się [757]. Skan siły przejmowania fazy potwierdza to dla każdej wartości: poślizg przechodzi do $-\pi$ ostro, a zmiana znaku pracy przy EBK występuje tylko przed tym przejściem [758]. Źródło o dynamice solitonu (faza zmieniana przez amplitudę) tego nie zmienia: zmiana znaku pracy przy EBK pojawia się tylko na gałęzi bez Maslowa [759]. Ta gałąź nie jest torusem z inną „połówką”: poślizg zmienia się na niej ciągle z $J_r$, a zera pracy nie spełniają samouzgodnienia z fazą modu i nie powtarzają się z okresem 1 [760]. Precesja spinu (AC + Thomas) na torusie EBK daje strukturę subtelną 2P dokładnie (wodór 10 968,69 wobec 10 969,045 MHz, stała LS o-Ps do 5·10⁻⁵), ale tylko z $\langle r^{-3}\rangle$ fali — orbita EBK daje 8/9; basen $F$ utrzymuje się w precesji dynamicznie, a faza spinowa na okres ($2{,}4\cdot10^{-5}\pi$) jest za mała, by nieść fazę Maslowa orbity [761]. Z falą zasady I w rzędzie $\hbar^2$ (iteracja Schwarza) $\langle r^{-3}\rangle$ zgadza się z QM do $10^{-5}$–$10^{-6}$, a precesja daje 2P$_{3/2}$ − 2P$_{1/2}$ wodoru 10 967,97 MHz wobec zmierzonych 10 969,045 MHz i stałą LS o-Ps do $10^{-4}$ — część spin–orbita struktury subtelnej wynika z zasady I i BMT, bez Sommerfelda [762]. Wstawiona do silnika jako opcja, część LS z precesji daje 10 967,97 MHz — o 0,74 MHz dalej od pomiaru niż Sommerfeld, z powodu resztkowego błędu $\langle r^{-3}\rangle$ fali; domyślna ścieżka zostaje bez zmian [763]. Druga iteracja Schwarza (rząd $\hbar^4$), policzona z wyciętymi oknami wokół punktów zwrotnych (metoda dobrana po fakcie), daje $\langle r^{-3}\rangle$ do $\sim2\cdot10^{-6}$ i 2P$_{3/2}$ − 2P$_{1/2}$ 10 968,72 MHz — równe Sommerfeldowi do 0,01 MHz; reszta wobec pomiaru (−0,33 MHz) to QED [764]. Powtórzona z prerejestracją dla wszystkich 28 stanów n = 2–8 metoda jest stabilna i daje rozszczepienia równe Sommerfeldowi do 0,013 MHz, ale kontrola metody dla 5p–8p przekracza ustalony próg, więc tablica nie została podstawiona do silnika [765]. Importem zasady I zostają warunek na orbitę (wybór torusa), którego pole nie daje [743–747],
i wybór polaryzacji Schwingera dla spinów [731].
W zmiennych Kustaanheimo–Stiefela ruch Keplera to cztery oscylatory z więzem, $\sum_i(n_i+\tfrac12)=2n$,
więc połówka Langera $(l+\tfrac12)\hbar$ jest tą samą fazą Maslowa ½ co połówka radialna, a nie
poprawką doraźną [630].

- Pętle działania: $\oint p\,dq=(k+\tfrac12)h$, co daje całkowite $n$ i długości Langera $|\mathbf L|=(l+\tfrac12)\hbar$ [462].
- Długości momentów: $|\mathbf S|=(S+\tfrac12)\hbar$, $|\mathbf J|=(J+\tfrac12)\hbar$, $|\mathbf F|=(F+\tfrac12)\hbar$, $k=j+\tfrac12$ [462, 515, 520].
- Mod pola: odstęp $\hbar\omega_{\rm mod}$ na foton, czyli $E/\omega=\hbar$ dla modu $(s,m)$. Mod zabiera $s\hbar$ działania i $m\hbar$ momentu pędu [408, 592]. Z EBK oscylatora pola bierzemy **tylko odstęp**; człon $\tfrac12\hbar\omega$ (pole punktu zerowego) jest opcją, odrzuconą dynamicznie [387, 465].
- Trójkąt długości dla fotonu $j$: $|L_l-L_{\rm ph}|\le L_{l'}\le L_l+L_{\rm ph}$, $L_{\rm ph}=(j+\tfrac12)\hbar$; potrzebny dla E2, gdzie zabrania $s\to s$ [598].

### Zasada II — elementy macierzowe (korespondencja Borna)

- Wielkość przejścia liczy się w punkcie środkowym, $L_m=l_>$. Pochodna po działaniu staje się różnicą między sąsiednimi punktami środkowymi [519, 594]:
  $\langle r^{-3}\rangle_B=\dfrac{\mu k}{L\hbar}\big[\langle r^{-2}\rangle(L-\tfrac\hbar2)-\langle r^{-2}\rangle(L+\tfrac\hbar2)\big]=\dfrac{1}{a^3n^3\,l(l+\frac12)(l+1)}$, czyli QM do $2\cdot10^{-15}$.
- Kwadraty stają się iloczynami wartości na obu końcach: $L^2\to(L_l-\tfrac12)(L_l+\tfrac12)=l(l+1)\hbar^2$, a $|c|^2\to c_ic_f$ (średnia geometryczna wag Zeemana; $\sqrt{P_T}$ Zeemana pary) [518, 539, 594].
- Człon porządku $(\mathbf S\cdot\mathbf L)^2\to(\mathbf S\cdot\mathbf L)(\mathbf S\cdot\mathbf L+\tfrac12)$ w $\langle(\mathbf S\cdot\hat{\mathbf r})^2\rangle$ to jedyna część spoza zasady punktu środkowego [521, 594]. Daje $\langle S_{12}\rangle(^3P_0)=-4$ dokładnie. Suma Heisenberga–Kramersa odtwarza go tylko asymptotycznie, z resztą $\sim\hbar^2/4L^2$ [595].
- Energia poziomu ponad Bohrem to klasyczny Breit–Pauli policzony regułami Q i B [532, 590]. Dla trypletu P Ps $E_J=a\langle\mathbf L\cdot\mathbf S\rangle_J+b\langle S_{12}\rangle_J$, gdzie $a=(2g-1)k\hbar^2\langle r^{-3}\rangle/(2m^2c^2)$ i $b=(\mu_0/4\pi)\mu^2\langle r^{-3}\rangle$; dla par o różnych masach $F(L_l)+(g-1)[F(k\hbar)-F(L_l)]$. Wynik równa się wartości oczekiwanej QM operatorów BP do $3\cdot10^{-4}$ [590].
- Tempo $\Delta l=-1$ to liczba kwantów modów współbieżnych na orbicie Langera, $\Gamma=\sum_sP_{s,+}/(\hbar s\Omega)$ [409, 432].
- Tempo $\Delta l=+1$ (stany $s$, kanały w górę z $l\ge1$, wychwyt) liczy się z modu przeciwbieżnego w punkcie środkowym $n_m=(n+n')/2$, $L_m=l_>$, z wagą [560, 579, 584, 585]:
  $A=w\,\tfrac23\,\omega^3|c_s^{(m)}|^2/c^3$, przy czym $w=2l_>/(2l+1)=L_m/L_l$ to kątowy czynnik QM względem granicy klasycznej $\tfrac12$.
- Podział tempa $\Delta l=-1$ na kanały $n'$ pochodzi z punktu środkowego w średniej geometrycznej, $n_m=\sqrt{nn'}$; rozgałęzienia p, d, f wychodzą w 0,016 od QM [599, 613].
- E2 stanów $s$ to mod przeciwbieżny $m=-2$ kwadrupola w punkcie środkowym, z $P=\omega^6q_Q^2|C_s|^2/(40c^5)$ [583]. Od [736] tempa E2 wszystkich stanów z elementu $r^2$ fali EBK: $A=kq_Q^2\omega^5G(R_2a^2)^2/(90\hbar c^5)$, $G=6(2l'+1)\left(\begin{smallmatrix}l&2&l'\\0&0&0\end{smallmatrix}\right)^2$.

### Zasada III — prawdopodobieństwa (reguła Borna)

- Wybór stanu przy powstaniu: $p_k=|\langle k|\psi_{\rm koh}(\hat{\mathbf L},\hat{\mathbf s}_1,\hat{\mathbf s}_2)\rangle|^2$, czyli rzut stanu koherentnego klasycznej konfiguracji na stany stacjonarne QM. Konfiguracja jest następnie obracana do basenu wylosowanej etykiety [547, 591, 619]. Rzut jednego momentu na oś jest szczególnym przypadkiem tej reguły.
- Wybór $J'$ trypletu przy fotonie E1 ($l,l'\ge1$): $P(J')\propto(2J'+1)\{l\,J\,1;J'\,l'\,1\}^2$, ten sam mechanizm co przy powstaniu [631].
- Kanał wolnych spinów: $P_S=(1-\hat{\mathbf s}_1\cdot\hat{\mathbf s}_2)/4$ [421].
- Stany własne pola: hamiltonian $2P$ (12×12: $\sum_JE_JP_J$ + Zeeman) i $1S$ (Breit–Rabi); po fotonie E1 stan spinowy to $\mathrm{Tr}_L|k\rangle\langle k|$ [547, 601].
- Wybór modu i kanału fotonu jest proporcjonalny do liczby kwantów modu.
- Chwila zdarzenia bez pamięci przy stałym tempie: $P(t)=e^{-\Gamma t}$, z progiem $\mathrm{Exp}(1)$ dla fotonu, anihilacji i wychwytu [561].

### Zasada IV — punkt $r=0$

- Gęstość kontaktu stanu $s$ to granica $l\to0$ członu odśrodkowego reguły B w klasycznym bilansie sił $\langle V'\rangle=\langle L^2/\mu r^3\rangle$ [593]:
  $n_{QR}=\dfrac{1}{2\pi}\lim_{L\to\hbar/2}\dfrac{(L^2-\hbar^2/4)\langle r^{-3}\rangle_B}{\hbar^2}=\dfrac{\hbar}{2\pi a^3n^3L}=|\psi_{ns}(0)|^2$. Importem zostaje liczba $1/(2\pi)$ i umieszczenie tej granicy w $r=0$.
- $L$ w gęstości kontaktu jest $L$ kanonicznym, niosącym etykietę: $\mathbf L_{\rm kan}=\mathbf L_{\rm mech}+\Delta\mathbf J_{\rm ukr}$, a zachowane całkowite to $\mathbf L_{\rm kan}+\mathbf S$. Wtedy 1S osiągnięty kaskadą jest tym samym 1S co przygotowany: $\tau_{\rm o\text{-}Ps}$ ze startów 2p i 3d równa się startowi 1s do $2\cdot10^{-10}$ dla każdego $J'$ (sprawdzenie walidacji) [636, 637].
- Kontakt Fermiego w HFS przy $g=2$ jest dokładnie granicą $l\to0$ wzorów reguły B dla $l\ge1$ (wodór $16/3$, Ps $8/3$) [622]. Gęstość jest potrzebna tylko dla części anomalnej.
- Liczby QED z literatury dla temp anihilacji i HFS, do rzędu $\alpha^2$ i jednej pętli [575, 587, 620]:
  $\Gamma_p=4\pi r_e^2c\,n_{QR}\,f_p$, $\Gamma_o=\varepsilon_{\rm OP}\Gamma_p\,f_o/f_p$, $\varepsilon_{\rm OP}=\dfrac{4(\pi^2-9)\alpha}{9\pi}$,
  $f_p=1-\big(5-\tfrac{\pi^2}{4}\big)\tfrac\alpha\pi+2\alpha^2\ln\tfrac1\alpha+5{,}1243\big(\tfrac\alpha\pi\big)^2$,
  $f_o=1-10{,}286606\,\tfrac\alpha\pi-\tfrac13\alpha^2\ln\tfrac1\alpha+45{,}06\big(\tfrac\alpha\pi\big)^2$.
- Wirtualna anihilacja trypletu: $\Delta E=\hbar\Gamma_p/(2\alpha)=\dfrac{2\pi\alpha\hbar^3}{m^2c}\,n_{QR}$ [427].
- Jedna pętla HFS 1s Ps: $-\tfrac\alpha\pi\big(\tfrac89+\tfrac12\ln2\big)\alpha^4mc^2$ (Karplus i Klein), minus część anomalna $\alpha^4mc^2[(g/2)^2-1]/3$, którą model ma już w kontakcie magnetycznym [575].
- Prawo płaszczyzny $3\gamma$: $W(\hat{\mathbf n}|\boldsymbol\xi)\propto|\boldsymbol\xi|^2+|\boldsymbol\xi\cdot\hat{\mathbf n}|^2$ dla każdej polaryzacji trypletu. Z diagramów drzewowych do $5\cdot10^{-15}$; z QED zostaje tylko liczba 2 [587, 588].

### Czego zasady nie obejmują

**Wybór kanału E1 przy emisji jest klasyczny** [596–598] (rys. 1b). Warunek „$l'$ mieści się w $n'$” to istnienie elipsy Keplera, $L'\le n'\hbar$. Przy wartościach z zasady I daje on $l'\le n'-1$ w 17 545 kanałach bez wyjątku.

$\Delta l=\pm1$ to zawartość modów dipola orbity płaskiej: wokół $\hat{\mathbf L}$ występują tylko $m=\pm1$.

Przypadek $l=0$ wynika z tego, że reakcja promieniowania nie odwraca orbity. Przy długościach Langera odwrócenie jest równoważne złej parzystości $\Delta l$.

**Wejścia mierzone, a nie importy:** $\alpha$, masy, momenty magnetyczne, czynniki $g$ (od [724] CODATA 2022: $\alpha^{-1}=137{,}035999177$; wcześniej $\varepsilon_0$ z 2018 i obcięte $\hbar$ dawały $137{,}035999000$).

---

## 4. Zgodność z pomiarami

Wartości pochodzą z domyślnego silnika. W kolumnie „ocena” podajemy odchylenie od pomiaru,
a dla temp E1 wodoru od dokładnego QM.

### 4.1. Pozytonium

| obserwabla | model | pomiar | ocena |
|---|---|---|---|
| $\tau_{\rm p\text{-}Ps}$ | 125,161 ps | 125,14 ps [AG94] | $+0{,}015\,\%$ |
| $\tau_{\rm o\text{-}Ps}$ (start 1s) | 142,046 ns | 142,04 ns [V03] | $+0{,}006\,\%$ |
| $\tau_o/\tau_p$ (tempa) | 1134,9 | 1135,0 | $-0{,}01\,\%$; w rzędzie wiodącym $1/\varepsilon_{\rm OP}=1113{,}97$ |
| rozkład czasu rozpadu | wykładniczy (KS $p=0{,}90/0{,}99$) | wykładniczy | zgodny |
| HFS 1s | 203,381 GHz | 203,3942 GHz [I14] | $-6\cdot10^{-5}$; w rzędzie wiodącym $+0{,}62\,\%$ |
| 1S–2S | reszta $+11{,}68$ MHz [724] ($+15{,}62$ przed stałymi CODATA 2022 [720]; przed importem #24 $+1737{,}0$); rząd $\alpha^4$ = teoria do 0,03 MHz [629] | 1 233 607 216,4 MHz [F93] | $1{,}3\cdot10^{-8}$; reszta ≈ brakujące $\alpha^6$–$\alpha^7$ z [K04] |
| $2^3P_1-2^3P_0$ / $2^3P_2-2^3P_1$ | 5,4852 / 4,3856 GHz | 5,488 / 4,388 GHz [H93] | do $5\cdot10^{-4}$ |
| $\nu_J=E(2^3S_1)-E(2^3P_J)$ (o-Ps) | 18 496,12 / 13 010,90 / 8 625,26 MHz [720] | 18 501,02(61) [G20] / 13 012,42(1,67) [H93] / 8 627,94(95) MHz [S23] | $-4{,}90/-1{,}52/-2{,}68$ MHz; silnik − teoria [K04] $-2{,}1$ do $-1{,}5$ MHz ($\alpha^6$) |
| $2^3S_1-2^1P_1$ (p-Ps) | 11 183,90 MHz [720] | 11 180 MHz [L94] | $+3{,}90$ MHz |
| linia o-Ps $2p\to1s$ wobec 1S–2S | $-0{,}081/-0{,}058/-0{,}040$ meV ($^3P_{0,1,2}$) | — | w całości $\nu_J$ ($2P$ leży poniżej $2^3S_1$); reszta $-4{,}4$ µeV |
| gaszenie magnetyczne (niskie pole) | $\lambda_p/\lambda_o=1134{,}9$ | $1302\pm15\,\%$ [HMW55] | $0{,}86\sigma$ |
| liniowy Zeeman singletu | brak (bit w bit przy $\pm B$) | brak (QM) | zgodny [589] |
| płaszczyzna $3\gamma$: $\langle P_2(\hat{\mathbf n}\cdot\hat{\mathbf S})\rangle$ | $-0{,}041$ (1000 rozpadów) | QED $-0{,}05$ | w niepewności próby |
| korelacja CPT $\hat{\mathbf S}\cdot(\mathbf k_1\times\mathbf k_2)$ | $-0{,}0023\pm0{,}0175$ | QED 0 | zgodny |
| korelacja CP $(\hat{\mathbf S}\cdot\hat{\mathbf k}_1)(\hat{\mathbf S}\cdot\hat{\mathbf n}_{12})$ | $(+0{,}1\pm3{,}3)\cdot10^{-5}$ (generator, $6\cdot10^7$ losowań) | $0{,}0013\pm0{,}0022$ [Y10] | zgodny; zero z konstrukcji [623] |

### 4.2. Wodór

| obserwabla | model | odniesienie | ocena |
|---|---|---|---|
| $\tau(2p)$ | 1,6082 ns (od [734]; wcześniej 1,6747) | 1,600(4) ns [BG66]; QM 1,595 ns | $+0{,}5\,\%$, $+2{,}1\sigma$ (wcześniej $+4{,}7\,\%$, $\sim19\sigma$ — granica reguły tempa) |
| $\tau(3p)$ | 5,5848 ns | QM 5,269 ns; pomiar 5,58(13) ns [BG66] | $1{,}060\times$ QM |
| $\tau(np)$, $n=4..8$ | $1{,}009\times$ QM od [734] (wcześniej $1{,}061$–$1{,}062$) | QM (Gordon) | fala EBK [733]; wcześniej stały nadmiar 6 % [612] |
| $\tau(nd)$ / $\tau(nf)$, $n\le8$ | $1{,}02$–$1{,}03$ / $1{,}01$–$1{,}02\times$ QM | QM | — |
| $\tau(ns)$, $n=3..6$ | $0{,}977/1{,}016/1{,}066/1{,}123\times$ QM | QM | dryf z $n$ [560] |
| rozgałęzienie $3p\to2s$ / $4p\to2s$ | $0{,}139\pm0{,}019$ / $0{,}097\pm0{,}024$ (silnik) | QM 0,118 / 0,119 | zgodne [614] |
| rozgałęzienie w $J$ przy E1 ($j'$) | $=6j$ do $10^{-3}$ | QM | dokładne dla spinu ½ [610] |
| rozgałęzienie $F'$ przy $2p\to1s$ wodoru | 0/8, 5/19, 8/16, 0/21 (losowanie z wagą 6j [740], etykieta z księgi [741]) | QM 0, 1/3, 2/3, 0 | −0,65σ, −1,41σ; waga jest importem |
| $2P_{3/2}-2P_{1/2}$ | 10 968,705 MHz | 10 969,045(15) MHz [HP94, LP81] | $+0{,}34$ MHz = Lamb $P$ i odrzut |
| $2S_{1/2}-2P_{1/2}$ (Lamb) | 1058,29 MHz (import [721]) | 1057,8298(32) MHz [B19] | $+0{,}46$ MHz (+1,9 neV), $\alpha(Z\alpha)^6$ pominięte |
| 1S–2S (środek ciężkości) | reszta $-3{,}94$ MHz (import [721], stałe [724]) | 2 466 061 413,187 MHz [P11] | $1{,}6\cdot10^{-9}$; $\alpha(Z\alpha)^6$ pominięte (przed [724] $+3{,}94$, w tym $+7{,}9$ ze stałych) |
| HFS $2P_{1/2}$ | 59,169 MHz | 59,22(14) MHz [LJP75] | $0{,}37\sigma$ |
| HFS 1s | 1,420485 GHz | 1,420406 GHz [H70] | $+5{,}6\cdot10^{-5}$ |
| $g_J/g_s-1$, H 1s | $-17{,}70\cdot10^{-6}$ | $-17{,}709(13)\cdot10^{-6}$ [TR77] | w niepewności |

![Rys. 1. Kaskada E1 wodoru i płaszczyzna działań](rysunki/rys1_kaskada.svg)

**Rys. 1.** (a) Kaskada E1 wodoru dla $n\le4$, ze schematyczną skalą energii. Grubość
strzałki jest proporcjonalna do udziału kanału w QM (wzór Gordona). Niebieskie strzałki
to $\Delta l=-1$ (mod współbieżny, podział z punktu środkowego $\sqrt{nn'}$), pomarańczowe
to $\Delta l=+1$ (mod przeciwbieżny w punkcie środkowym). Liczby przy strzałkach podają
udział QM i udział silnika [613]; liczby przy poziomach to stosunek $\tau_{\rm silnik}/\tau_{\rm QM}$
[560, 612]. (b) Ta sama kaskada na płaszczyźnie działań $(|\mathbf L|,J)$. Szare odcinki
to przedziały $0\le L\le n\hbar$, w których istnieje elipsa Keplera. Punkty pełne to
długości Langera zasady I; punkty puste ($l=n$) leżą poza przedziałem. Foton E1 zabiera
$\hbar$ z $|\mathbf L|$ i $\Delta n$ kwantów $J$ (§2.5). Kanał $3d\to1p$ (linia czerwona)
wymagałby $L'=1{,}5\hbar$ przy $n'=1$. Takiej orbity nie ma, więc zakaz $l'\le n'-1$
wynika z klasyki (§3) [596]. Skrypt rysunku: `rysunki/rys1_kaskada.py`.

### 4.3. Rząd QED według obserwabli

Rząd QED w modelu zależy od obserwabli [620]:
- tempa anihilacji: do $\alpha^2$, łącznie z $\alpha^2\ln\alpha$;
- HFS 1s Ps: jedna pętla (człon $\alpha^6\ln\alpha$, $+0{,}019$ GHz, pominięty);
- poziomy: rząd wiodący $\alpha^4$; w Ps od [720] także blok $m\alpha^5$ niezależny od spinu (przesunięcie Lamba Ps, import #24: Pineda–Soto [PS98], logarytmy Bethego [JM05]); w wodorze od [721] przesunięcie Lamba jako import #25 (teoria poziomów H CODATA 2006 [MTN08], bez członu momentu anomalnego, który silnik już ma; $R_p$ z [A13]).

Stąd reszty: 1S–2S Ps $+11{,}7$ MHz (od [724]; $+15{,}6$ przed stałymi CODATA 2022; przed [720] $+1{,}74$ GHz — dokładnie brakujące $\alpha^5$ bez części HFS, $\alpha^6$ i $\alpha^7$ z [K04] [629]), $\nu_J$ $-1{,}5$ do $-4{,}9$ MHz (przed [720] $-0{,}27$ GHz), wodór $2P$ $+0{,}34$ MHz, Lamb wodoru $+0{,}46$ MHz (przed [721] $-1041$ MHz), 1S–2S wodoru $-3{,}94$ MHz (od [724]; $+3{,}94$ przed stałymi CODATA 2022; przed [721] $+7{,}13$ GHz — nie $+6{,}82$ GHz, jak podawał [429]: pomiar [P11] to środek ciężkości nadsubtelnej),
$3D_{5/2}-3S_{1/2}$ $-0{,}23$ MHz (przed [721] $+311$ MHz).

---

## 5. Czego klasyczna dynamika nie daje (wyniki negatywne)

Każda zasada była sprawdzana pod kątem tego, czy da się ją wyprowadzić z dynamiki. Poniżej
zestawienie prób rozstrzygających.

### 5.1. Warunek Bohra a faza zegara wewnętrznego

Model ma jedną fazę, której nie wprowadza żadna z zasad: fazę pętli momentu
magnetycznego. Pętla ma promień $R=\hbar/(mc)$, a ładunek obiega ją z prędkością $c$,
więc jej częstość spoczynkowa to $\Omega_0=c/R=mc^2/\hbar$, czyli częstość de Broglie'a.
Sprawdzamy, czy ten zegar wybiera całkowite $n$ [380, 603, 604].

**Faza czasu własnego.** Wzdłuż toru cząstki $i$ faza w laboratorium to
$\varphi_i(t)=(m_ic^2/\hbar)\int_0^t dt'/\gamma_i(t')$. W ciągu okresu $P$ zegar
orbitujący traci względem nieruchomego zegara o tej samej częstości

$$
\Delta\varphi_i=\frac{m_ic^2}{\hbar}\oint\Big(1-\frac{1}{\gamma_i}\Big)dt
=\frac{1}{\hbar}\oint T_i\,dt\;\big[1+O(v^2/c^2)\big],\qquad T_i=\tfrac12m_iv_i^2 .
\tag{5.1}
$$

**Krok 1: podział między cząstki.** W układzie środka masy
$\mathbf v_1=(m_2/M)\,\mathbf v$ i $\mathbf v_2=-(m_1/M)\,\mathbf v$, więc

$$
T_i=\frac{m_j}{M}\,T,\qquad T=\tfrac12\mu v^2,\qquad \sum_iT_i=T .
\tag{5.2}
$$

**Krok 2: działanie.** Energia kinetyczna jest jednorodna stopnia 2 w prędkościach,
więc $\mathbf p\cdot\mathbf v=2T$ przy dowolnym potencjale. Z definicji zmiennej
działania $\oint\mathbf p\cdot d\mathbf q=2\pi J$, stąd

$$
\oint T\,dt=\tfrac12\oint\mathbf p\cdot d\mathbf q=\pi J
\quad\Longrightarrow\quad
N_i\equiv\frac{\Delta\varphi_i}{2\pi}=\frac{m_j}{M}\,\frac{J}{2\hbar}.
\tag{5.3}
$$

Zegar liczy więc **główne działanie** $J=J_r+L$, a nie $L$. W rzędzie wiodącym nie zależy
od mimośrodu. W pozytonium $m_j/M=\tfrac12$ i $N_i=J/(4\hbar)$. W wodorze zegar
elektronu daje $N_e=0{,}49973\,J/\hbar$ ($=\tfrac12m_p/M$), a protonu $0{,}00027\,J/\hbar$.
Kwadratura na orbitach Keplera ($n=1$–$4$, $e\le0{,}95$) daje $0{,}250000$–$0{,}250005$ na
zegar w Ps [380, 603].

**Krok 3: następny rząd.** Z rozwinięcia $1-1/\gamma=\tfrac12\beta^2+\tfrac18\beta^4+\dots$
oraz średnich Keplera $\langle v^2\rangle=k/(\mu a_n)$,
$\langle v^4\rangle=(k/\mu a_n)^2\,(4/\sqrt{1-e^2}-3)$, przy $1/\sqrt{1-e^2}=J/L$:

$$
N_i=\frac{m_j}{M}\,\frac{J}{2\hbar}\,(1+\delta),\qquad
\delta=\frac{\langle v_i^4\rangle}{4c^2\langle v_i^2\rangle}=\frac{\alpha^2}{16n^2}\Big(\frac{4J}{L}-3\Big)\quad(\mathrm{Ps}).
\tag{5.4}
$$

Wzór (5.4) wyjaśnia zmierzone odchyłki od (5.3), sięgające $3\cdot10^{-5}$, z resztą
$\le5\cdot10^{-9}$, a po dodaniu wyrazu $\beta^6$ z resztą $\le2\cdot10^{-9}$ [380].
Czynnik $4J/L-3$ ma postać relatywistycznej poprawki kinetycznej Sommerfelda
$n/l-\tfrac34$, bo faza zegara jest działaniem swobodnej cząstki
$-mc^2\!\int d\tau/\hbar$.

**Krok 4: blokada fazy.** Warunek $N_i\in\mathbb Z$ daje $J=4N\hbar$ w pozytonium i
$J\simeq2N\hbar$ w wodorze, a nie $J=n\hbar$. Dla $n=1,2,3$ zegar Ps przesuwa się o
$\tfrac14,\tfrac12,\tfrac34$ cyklu. Dwa zegary pary nie dają dudnień: w układzie środka
masy $|\mathbf v_1|=|\mathbf v_2|$, więc $d(\varphi_1-\varphi_2)/dt=0$ dokładnie. Faza
opóźnienia $\Omega_0r/c$ między nimi nie ma struktury przy całkowitym $L/\hbar$ [380].

**Krok 5: co dokładnie brakuje.** Weźmy fazę lagranżjanu obu cząstek, liczoną względem
odniesienia o częstości poziomu $E/\hbar$. Z tożsamości Legendre'a $\mathcal L+H=\mathbf p\cdot\mathbf v$:

$$
\Phi=\frac1\hbar\oint(\mathcal L+E)\,dt=\frac1\hbar\oint\mathbf p\cdot\mathbf v\,dt=\frac{2\pi J}{\hbar},
\tag{5.5}
$$

co numerycznie zachodzi do $10^{-12}$ dla $n=1$–$3$ i $e\le0{,}9$ [603]. Odejmując (5.3)
zsumowane po cząstkach,

$$
\Phi-\sum_i\Delta\varphi_i=\frac1\hbar\oint T\,dt=\frac1\hbar\oint\big(E-V\big)\,dt .
\tag{5.6}
$$

Brakująca połowa jest więc dokładnie sprzężeniem fazy z potencjałem, $-q\phi/\hbar$
(faza cechowania), wraz z odniesieniem o częstości $E/\hbar$ zamiast $mc^2/\hbar$.
Dopiero wtedy jednoznaczność fazy, $\Phi=2\pi n$, daje $J=n\hbar$. Faza Maslowa $\tfrac12$
(zasada I) jest jeszcze osobnym efektem falowym. Jednoznaczna jest tylko suma po obu
cząstkach, a podział potencjału między nie jest nieokreślony. To zgadza się z mechaniką
kwantową, w której fala żyje we współrzędnej względnej.

Zegar nie jest przy tym fragmentem brakującej fali. Wzdłuż linii świata swobodnej
cząstki faza fali de Broglie'a $(\mathbf p\cdot\mathbf x-Et)/\hbar$ równa się dokładnie
$-mc^2\tau/\hbar$. Zegar jest więc tą samą falą, odczytaną wzdłuż toru. Brakuje pola
rozciągniętego w przestrzeni, o liczbie falowej $p/\hbar$, którego jednoznaczność wokół
orbity daje warunek Bohra.

**Krok 6: czy klasyczna pętla czuje potencjał?** Rozważmy elektron jako naładowaną pętlę
$R=\hbar/mc$ na orbicie wodoru [604]. Jednorodna część pola Coulomba znosi się z siłą
bezwładności, bo środek pętli spada swobodnie. Zostaje część pływowa:

$$
\Big\langle\frac{\delta\omega}{\omega}\Big\rangle\sim\frac{\alpha^4}{n^6\beta_{\rm wewn}^2}\;(+),
\qquad\text{wymagane}\qquad \frac{\langle V\rangle}{mc^2}=-\frac{\alpha^2}{n^2}
=-5{,}3\cdot10^{-5}\;(n=1).
$$

Zmierzone przesunięcie wynosi $10^{-11}$–$1{,}6\cdot10^{-7}$ na orbitę. Ma zły znak,
skaluje się jak $n^{-6}$ zamiast $n^{-2}$ i maleje z prędkością wewnętrzną.

**Wniosek.** Zegar wybiera właściwą zmienną: $J$, a nie $L$, niezależnie od $e$. Uzasadnia
to, że reguła działania (§2.5) odejmuje kwanty $J$. Wartości całkowitej i fazy Maslowa
zegar nie daje; zostają one treścią zasady I. Pole punktu zerowego (SED) też nie daje
stanu stacjonarnego $n=1$ [387].

Brak z kroku 5 można zapisać bez fali rozciągniętej w przestrzeni [627]. Suma po kolejnych
obiegach $\sum_N e^{iN\Phi}$ jest niezerowa tylko przy $\Phi=2\pi k$, więc jednoznaczność fali
to warunek interferencji fazy niesionej wzdłuż toru. Wynik zależy od sumowanej fazy. Zegar
pętli ($\Phi=(m_j/M)\pi J/\hbar$) daje $J=4N\hbar$. Faza (5.5) daje $J=k\hbar$, a ze skokiem
$-\pi/2$ w punktach zwrotnych — wartości reguły Q. Interferencja zastępuje więc jednoznaczność
fali, ale potencjał w fazie i odniesienie $E/\hbar$ zostają jednym importem: fazą $S/\hbar$.
Suma wybiera wartości, a nie prowadzi do nich orbity; w modelu to wystarcza, bo zasada I działa
jako reguła wyboru.

Odniesienie $E/\hbar$ ma przy tym sam zegar czasu własnego. Z wiriału dla $1/r$,
$\langle\sum_i\mathbf p_i\cdot\mathbf v_i\rangle=-\langle V\rangle$, wynika

$$
\sum_i\Big\langle\frac{m_ic^2}{\gamma_i}-m_ic^2\Big\rangle=E
$$

dokładnie i relatywistycznie, z udziałem cząstki $(m_j/M)E$ (do $7{,}6\cdot10^{-12}$ na orbitach
z $e\le0{,}9$) [712]. Brakuje więc tylko pamięci. Składowa pola o pierwotnej częstości, propagująca
swobodnie (także progowa Kleina–Gordona), ma pamięć $\sim k^{-3/2}$ na obieg i rezonanse przy $J=2nh$.
Składowa sprzężona z potencjałem partnera, napędzana tym zegarem, ma mody związane, które rosną
rezonansowo dokładnie przy orbitach o energii poziomów. Pasma boczne $E+j\hbar\Omega=E_m$ odtwarzają
przy tym zasadę korespondencji [712]. Jest to jednak naładowana fala materii, a nie część pola EM.
Dla Ps działa ona tylko we współrzędnej względnej, a Maslowa nie daje.

Człon Darwina stanów S pozostaje importem także po próbie wyprowadzenia go z geometrii modelu.
Ładunek rozmyty szybkim ruchem daje kontakt $\tfrac16\langle\delta^2\rangle\nabla^2V$, a Darwin
Diraca odpowiada $\langle\delta^2\rangle=\tfrac34\lambda_C^2$. Pętla momentu $R=(g/2)\lambda_C$ daje
$\tfrac43(g/2)^2$ tej wartości, czyli 2S$_{1/2}-$2P$_{1/2}=+7142{,}69$ MHz, promień $r^*=R/2$ daje
$-14\,584{,}93$ MHz, a profil Plummera przeniesiony na ładunek $+69\,051{,}84$ MHz, wobec zmierzonych
$1057{,}845(9)$ MHz [Lundeen, Pipkin 1981] [713]. Izotropowe drganie 3D o amplitudzie $\hbar/2mc$
na oś odtwarza Darwina, z resztą $-137{,}66$ MHz o znaku przeciwnym do Lamba. Na orbicie Langera
pierścień nie jest osiągany (peryapsis $/R=17$–186), kontakt znika, a człon tensorowy
$(R^2/4r^3)(3(\hat s\cdot\hat r)^2-1)$ zmienia 2P$_{3/2}-$2P$_{1/2}$ o $+2434{,}82$ MHz (22 %) [713].
Dla spinu ½ ten tensor znika tożsamościowo, więc Darwin wymaga drgania 3D, gęstości w zerze i
spinu kwantowego. W Ps nadmiar pętli jest równy w obu kanałach: $+1796{,}60$ MHz dla 2$^1S_0$ i 2$^3S_1$.
Struktura $n=2$ Ps zawiera wirtualną anihilację i kontakt $C_n=\alpha^4mc^2/(4n^3)$ od [427–431]. Bieżący silnik
daje $\nu_J=18\,230{,}08$, $12\,744{,}86$ i $8\,359{,}22$ MHz wobec zmierzonych $18\,501{,}02(61)$ [Gurung 2020],
$13\,012{,}42(1{,}67)$ [Hagena 1993] i $8\,627{,}94(95)$ MHz [Sheldon 2023]. Wspólna reszta $\approx-268$ MHz to blok
$m\alpha^5$ niezależny od spinu [640, 714]. Od [720] jest on importem #24 (zasada IV): równanie (3.8) Pinedy–Soto [PS98]
z logarytmami Bethego [JM05], dodane do energii poziomu Ps. Poziomy $n=2$ przesuwają się dokładnie o import ($2^3S_1$ $+263{,}13$, $2P$
$-2{,}91$ MHz). o-Ps: $\nu_J=18\,496{,}12$ / $13\,010{,}90$ / $8\,625{,}26$ MHz, reszty $-4{,}90$ / $-1{,}52$ / $-2{,}68$ MHz.
p-Ps: $2^3S_1-2^1P_1$ $+3{,}90$ MHz. 1S–2S: reszta $+15{,}62$ MHz. Linie $2p\to1s$ przesuwają się o $-8{,}219$ µeV; o-Ps
$J=2$ to 5,1017412 eV = 243,0233 nm, p-Ps 5,1025717 eV = 242,9838 nm. Czasy życia i fotony anihilacji się nie zmieniają,
a reszta wobec teorii ($-1{,}5$ do $-2{,}1$ MHz) to rząd $\alpha^6$.

Zegar drgania Diraca, o częstości $2(mc^2+E)/\hbar$, jako źródło pola swobodnego o progu $2mc^2$ ma największe
rezonanse pamięci dokładnie przy $J=nh$, bez sprzężenia pola z potencjałem [715]. Rachunek jest jednak tożsamy
z zegarem Comptona przy $\hbar\to\hbar/2$, co jest czynnikiem ½ z [603]. Kontrast rezonansów spada od 5,7 do 1,7
wraz ze wzrostem $L$, a $L$ nie jest wybierane. Opóźnienie zegara $2E$ dają tylko zegary obu cząstek razem.
Zegar samego elektronu przesuwa 1S–2S wodoru o $-2{,}71$ THz wobec pomiaru [Parthey 2011], a w Ps daje czynnik 4 [715].
Spin nie przejmuje roli pamięci. Kąt jego składowej w płaszczyźnie orbity, mierzony od peryapsis, zmienia się
na obieg o $(\pi\alpha^2/L^2)[(g-2)+3(g-1)S_\parallel/L]$, ponieważ człon spin–orbita przesuwa także apsydy. Silnik
odtwarza ten wzór dla wodoru do 3 % [716]. Jedyny zatrzask, $S_\parallel\approx0$, daje $|L+s|=1{,}803$, więc nie
wybiera $j$. W p-Ps klasyczny moment dipol–dipol tworzy spin całkowity, rosnący o $1{,}8\cdot10^{-5}$ na okres [716].
Relatywistycznie warunek rezonansu wspólnego zegara, $2|E|=j\hbar\Omega_r$, wyznacza
$N=J_r+\sqrt{L^2-\alpha^2\hbar^2}$, czyli energię $E=-\alpha^2/(2n^2)-\tfrac38\alpha^4/n^4$ niezależną od $L$ [717].
Struktura subtelna wodoru albo znika, albo, przy jawnych członach Darwina i spin–orbita, przesunięcie Lamba wychodzi
$+29\,198$ MHz, 27,6 raza więcej niż zmierzone. W Ps wszystkie przedziały $n=2$ odbiegają o $+14{,}3$ GHz. Reguła
EBK przy Langerze [628] pozostaje jedynym warunkiem zgodnym z rzędem $\alpha^4$.
Izotropowe drganie Diraca (amplituda $\hbar/2mc$ na oś, częstość $2mc^2/\hbar$) przyjmujemy więc jako fizyczną
treść dwóch importów — kontaktu Darwina [713] i reguły działania w rzędzie Bohra [715] — bez zmiany ich wartości
liczbowych, które w rzędzie $\alpha^4$ zostają przy EBK [717, 718].
Dlaczego zegar trafia w $J=nh$, choć nie jest sprzężony z potencjałem: w granicy nierelatywistycznej jego warunek to
$2\langle T\rangle=j\hbar\Omega_r$ — energetyczny warunek Bohra z 1913 r. Zegar zna potencjał przez wirial. W Keplerze
$\langle T\rangle=|E|$, a $E(J)$ jest jednorodne stopnia $-2$, więc $J\Omega=2|E|$ i warunek jest dokładnie $J=j\hbar$
(numerycznie do $10^{-13}$) [722]. Poza Keplerem ta równoważność znika: oscylator izotropowy dostaje $E=2j\hbar\omega$ bez
energii punktu zerowego $\tfrac32\hbar\omega$, a zaburzenie $\beta/r^2$ przesuwa poziomy $(1+L/\Lambda)\approx2$ razy mocniej niż QM [722].
Zegar ten już niesie fazę działania: $2T=\mathbf p\cdot\dot{\mathbf q}$, więc jego faza to $-S_0/\hbar$, a
$2\langle T\rangle=J_r\Omega_r+L\Omega_\phi$ [723]. Odczytany w czasie (źródło punktowe i pole) daje jeden warunek. Odczytany
w przestrzeni (ta sama faza w tym samym punkcie na kolejnych obiegach, skok $-\pi/2$ w punkcie zwrotnym) jest jednoznaczny
dokładnie w torusach EBK, także na orbitach niezamkniętych (Coulomb + $\beta/r^2$: ścisłe widmo 2D QM), a w rezonansach
odczytu czasowego dokładnie antysymetryczny. Przy skończonym czasie obserwacji filtr przepuszcza jednak wąskie pasy
przez punkty EBK, które znikają dopiero przy porównaniu kilku czasów [723]. Zasada I wymaga więc fali w przestrzeni
konfiguracyjnej i skoku Maslowa, a nie tylko zegara.
Fizyczny nośnik tej fali — pole ruchu względnego pary ze źródłem na cząstce o fazie cechowania $S_0-Et$ — daje Maslowa w swoich modach ($J_r=n_r+\tfrac12$, w 3D Langer), ale nie w fazie źródła: na własnym torusie EBK źródło i mod różnią się dokładnie o $\pi$ na okres radialny (2D), więc źródło rezonuje na krzywych jednego warunku, a orbita z siłą zwrotną nie zatrzaskuje się (0/8) [743]. Gdy źródło przejmuje fazę retardowanego pola pamięci, skok przechodzi do niego dynamicznie (poślizg $-3{,}103$ rad na okres radialny w torusie EBK, niezależnie od opóźnienia), ale źródło śledzi dominujący mod na każdej orbicie, więc torus nadal nie jest wybrany (0/8) [744]. Siła z gradientu natężenia pola pamięci (jak u kropli chodzących po fali) też nie zatrzaskuje orbity (0/12); pole wielomodowe wykonuje pracę [745]. Wyróżnia ona torus EBK (0,5; 1) jako najwolniejszy z 49 startów, ale nie jako atraktor i nie w drugim punkcie EBK, także przy sile przeskalowanej do tej samej względnej siły pola [746, 747]. Nośnik daje więc fazę Maslowa, a nie warunek na orbitę. Skok $-\pi/2$ jako warunek na orbitę zostaje importem.
Na orbicie osadzonej w 3D ten sam odczyt daje połówkę Langera. Współrzędna $\vartheta$ nachylonej płaszczyzny ma
dwa punkty zwrotne na obrót, więc skok Maslowa dodaje $-\pi$ na obrót, a jednoznaczność wymaga $J_r=(n_r+\tfrac12)\hbar$
i $L=(l+\tfrac12)\hbar$. Na Coulombie + $\beta/r^2$ to ścisłe widmo 3D QM: 18/18 punktów Langera spójnych przy
każdej długości przebiegu, 18/18 punktów 2D antysymetrycznych, 0/200 losowych torusów fałszywie spójnych [726].
Z precesją węzła (liniowy Zeeman orbitalny, $H=H_0+\omega L_z$) trajektoria obiega trzeci cykl torusa, który nie ma
kaustyki, więc ten sam odczyt daje $L_z=m\hbar$ całkowite: 16/16 torusów spójnych, $L_z=m+\tfrac12$ antysymetryczne,
0/100 losowych fałszywie spójnych; zarejestrowany test przeszedł dopiero po poprawce przekroju (orbity biegunowe) [727].
Cały warunek EBK 3D wynika więc z orbit, z jedynym importem: skokiem $-\pi/2$. Mechanizm działa dla wodoru; para e⁺e⁻
nie ma orbitalnego Zeemana (p-Ps), a w o-Ps węzeł precesuje wokół $\mathbf J$.
Tam sama faza orbitalna kwantowałaby $\mathbf L\cdot\hat{\mathbf J}$ — źle. Z fazą geometryczną spinu jako bąka warunek
jednoznaczności to $|\mathbf J|-|\mathbf S|\in\mathbb Z$, więc przy długości spinu reguły Q $|\mathbf S|=\tfrac32$ wychodzi dokładnie
trójkąt silnika $|\mathbf J|=J+\tfrac12$ (6/6 punktów 2³P_J spójnych, 0/100 losowych), a przy $|\mathbf S|=1$ — $|\mathbf J|$ całkowite [728].
Połówka $|\mathbf J|$ pochodzi z długości spinu $S+\tfrac12$. Ta z kolei wynika z warunku Diraca ($2|\mathbf S|$ całkowite,
inaczej faza geometryczna nie jest jednoznaczna) i skoku Maslowa na czapie sfery spinu: $S_z=|\mathbf S|-\tfrac12-n$, $2|\mathbf S|$ stanów,
przy liczbie stanów elektronu 2 [729]. Sprzężenie dwóch spinów Ps tym mechanizmem nie wychodzi: suma faz geometrycznych
wybiera $|\mathbf S|=1,2$ zamiast singletu $\tfrac12$ i trypletu $\tfrac32$; potrzebna jest parzystość trójkąta
$|\mathbf J_1|+|\mathbf J_2|+|\mathbf J|-\tfrac32\in\mathbb Z$ [LY95] [729]. Ta jednak też jest skokiem Maslowa: w reprezentacji Schwingera
(spin = para oscylatorów, zwykła przestrzeń konfiguracyjna) cykl precesji dwóch spinów ma indeks $\mu=6$, długość spinu
$\mu=4$, rzut $S_z$ $\mu=4$; EBK daje singlet $\tfrac12$ i tryplet $\tfrac32$, a dla spinów 10 + 10 widmo $\mathbf s_1\cdot\mathbf s_2$
ze stałym błędem $-A/8$ [730]. Cała reguła Q wynika więc z jednoznaczności fali i jednego skoku $-\pi/2$, dla spinów
liczonego w polaryzacji Schwingera.

Zastosowana do hamiltonianu silnika bez spinu (kinetyka relatywistyczna, Coulomb, człon Darwina;
pęd kanoniczny $\gamma m\mathbf v+\partial E_D/\partial\mathbf v$), ta interferencja daje strukturę subtelną
bez spinu dokładnie w rzędzie $\alpha^4$ [628]:

$$
\Delta E_{nl}=m\alpha^4c^2\Big[\frac{11}{64n^4}-\frac{1}{2n^3(2l+1)}\Big]\quad(\mathrm{Ps}),
$$

przy każdym $l$, także $l=0$, a dla wodoru $l\ge1$ z resztą ścisłego wzoru Sommerfelda. Reszty skalują
się jak $\alpha^6$. Reguła Q nie jest więc od dynamiki niezależna: przy jednej zaimportowanej fazie poziomy
wynikają z klasycznego hamiltonianu modelu.

Część spinowa wychodzi tak samo [701–708]. Na ścieżce mechanicznej, okres radialny po okresie,
liczymy działanie kanoniczne $J+\oint(L_{q\mu}+L_T)\,dt$ (§2.3) i energię zachowaną. Przesunięcie
poziomu wynosi $\Delta E=-(2Ry/n^3)(n_{\rm dz}-n_E)$, a spin–orbitę wyodrębnia różnica wodór − antywodór
przy tych samych spinach:

$$
\frac{\Delta E_{\rm H}-\Delta E_{\bar{\rm H}}}{\Delta E_{2P_{3/2}-2P_{1/2}}\cos(\hat{\mathbf L},\hat{\mathbf s})}
=0{,}8912\quad\text{wobec}\quad \tfrac89=0{,}8889 .
$$

Połowę Thomasa daje $L_T$, który odejmuje dokładnie połowę $L_{q\mu}$. Dla pozytonium różnica
o-Ps − p-Ps daje $0{,}958$–$1{,}019$ (przy 8192 krokach na orbitę $0{,}9993$–$0{,}9994$ [711]) wartości $\tfrac89(a\,\mathbf L\cdot\mathbf S+b\,S_{12})$ z multipletu
$2^3P_J$ [H93] po odjęciu klasycznej energii dipol–dipol p-Ps. Czynnik $8/9$ to stosunek klasycznego
$\langle r^{-3}\rangle$ na orbicie Langera do kwantowego. Jest to efekt kwantowy, w modelu nazwany regułą
Borna $9/8$ [699, 700].

Część spinowa struktury subtelnej jest więc w klasycznej elektrodynamice modelu. Brakuje tylko kwantowego
$\langle r^{-3}\rangle$. Energia z $U_{qd}$ daje natomiast $g/(g-1)\cdot\tfrac89=1{,}775$ [701, 704], bo przy
sumie chwilowej nie jest zachowana. Żadna energia fotonu modelu nie jest jednak z niej liczona [708]. Na sumie
opóźnionej, gdzie $U_{qd}$ należy do energii zachowanej, wartość $\tfrac89$ daje działanie $J+\oint L_T\,dt$
(0,8876) [709]. Relacja Hamiltona $dJ_r=T\,dE-\Delta\varphi\,dL$ na rodzinie orbit o stałym $L$ tego wyboru nie
rozstrzyga. Całki spinowe na okres zależą tam tylko od $L$, bo $\oint r^{-3}dt=2\pi m/(Lp)$, więc spełnia ją
każdy kandydat [709]. Rozstrzygnąć może dopiero $\partial J_r/\partial L|_E=-\Delta\varphi$, czyli precesja peryapsu.
Na rodzinie zmieniającej $L$, dla pozytonium bez pola zewnętrznego, para sumy chwilowej
ma część spinową reszty $3\cdot10^{-9}$, a pozostali kandydaci $5\cdot10^{-6}$ [711]. Ta para daje o − p
$=0{,}9993$ i $0{,}9994$ wartości $\tfrac89(a\,\mathbf L\cdot\mathbf S+b\,S_{12})-U_{dd}(p)$ [H93], na przykład $-1819{,}36$ MHz wobec
$-1820{,}57$ MHz [711]. Na sumie opóźnionej żaden z czterech kandydatów nie tworzy pary z energią z $U_{qd}$:
$J+\oint L_T\,dt$ daje w wodorze 0,8876, ale dla Ps 0,33 [711]. Para tej sumy pozostaje otwarta.

### 5.2. Pozostałe zasady

**Element macierzowy (II).**
- Żadna jedna zasada punktu środkowego nie daje temp E1 wszystkich stanów w 10 % [599, 600, 615]. Sprawdzone były cztery średnie orbity, każda z dwiema częstościami, średnia Kramersa wzdłuż drogi i iloczyn $c_ic_f$.
- Element macierzowy samej fali EBK (Langer, jednorodne Airy) daje $\tau(2p)=1{,}653$ ns wobec 1,600(4) ns zmierzonych i 1,596 ns QM; elementy radialne są 1,7–5,3 % za małe [732]. Ten błąd pochodził jednak ze zszycia gałęzi Airy'ego: fala EBK jednorodna przy obu punktach zwrotnych (odwzorowanie na oscylator) daje elementy 0,1–0,5 % od QM i $\tau(2p)=1{,}608$ ns, 2,1σ od pomiaru; czasy życia $np$ stale ~0,9 % nad QM [733]. Nadmiar silnika (5–6 %) jest więc własnością reguły B, nie semiklasyki [733].
- Człon porządku jest dokładny tylko w granicy $l\to\infty$ [595].
- Rozgałęzienie „wbrew” $\Delta l$ dla spinu 1 maleje jak $l^{-4}$, ale nie znika [611]. Dwa kroki przez stan pośredni (zmiana $L$, potem przestawienie sprzężenia ze spinem), każdy w punkcie środkowym reguły B, dają $p=1/(4l_>^4)$: asymptotykę dokładnie, ale nie wartość przy skończonym $l$ [625].

**Reguła Borna (III).**
- Klasyczny układ $2P$ Ps w polu jest chaotyczny przy 0,01–0,3 T: od 53 do 87 % trajektorii ma dodatni wykładnik Lapunowa [602]. Torusów EBK tam nie ma.
- Bilans energii fotonu wyznacza wagi tylko w blokach $M$ o co najwyżej dwóch stanach (1S) [605].
- Średnie czasowe w klasycznej dynamice stanów koherentnych mylą wagi o 0,10–0,23 [609].
- Ergodyczność stanów własnych nie odtwarza składu stanów nawet w słabym polu [607, 608].
- Reguła Q wystarcza tylko przy $B\lesssim0{,}003$ T (baseny $J$) i $B\gtrsim1$ T (rzuty iloczynowe) [606].

**Punkt $r=0$ (IV).** Mod pola orbity przy $\hbar\omega=mc^2$ ma amplitudę $\sim10^{-804}$
(harmoniczna 37 557 orbity 1s), więc anihilacja nie jest modem klasycznego pola [587].
Logarytm Bethego z widma orbity wychodzi $2{,}62/2{,}49$ wobec $2{,}984/2{,}812$ (1s/2s), z
niedoborem $\approx0{,}34$ niezależnym od $n$. Wersja z różnicami reguły B nie ma
kontinuum, które dominuje [621].

---

## 6. Granice dziedziny

1. **Bariera Comptona** $r^*$: poniżej niej punktowa elektrodynamika przestaje obowiązywać i przebieg się kończy. Momenty są porównywalne z polem Coulomba dopiero dla $r\lesssim r^*\approx2\cdot10^{-3}a$.
2. **Retardacja:** stosunek okresu do czasu przelotu światła musi być $\ge150$; to margines całkowania, a nie fizyka.
3. **Tempo E1 stanów $p$:** 66–95 % tempa $\Delta l=-1$ pochodzi z harmonicznych poza siatką ($s\ge n$). Całość zgadza się dzięki kompensacji, a kanały daje podział z punktu środkowego [612, 613]. Zostawał nadmiar 6 %, w tym 2p $+4{,}7\,\%$ ($\sim19\sigma$); od [734] tempo E1 bierze się z elementu macierzowego fali EBK (zasada I) — 2p $1{,}608$ ns ($+2{,}1\sigma$), stany $np$ ~1 % nad QM [733, 734]. Ta sama fala daje E2 (elementy $r^2$ 0,04–0,57 % od QM; klasyczne E2 silnika w 3d ~29 % za małe), a fala spinu — M1 obrotu pojedynczego spinu dokładnie; M1 między sprzężonymi stanami spinu (21 cm, ²⁰⁹Bi⁸²⁺) wymaga współczynników 6j [735]. Od [736] silnik liczy E2 z elementu $r^2$ tej fali (3d: udział 9,87·10⁻⁶ wobec QM 9,95·10⁻⁶; fałszywy kanał E2 ze stanów $p$ usunięty, $\tau(2p)$ $1{,}60826$ ns), a M1 przewrócenia spinu jako tempo przejścia tryplet → singlet z elementu fali spinu i wagi 6j (Ps 1³S₁ → 1¹S₀ $3{,}38\cdot10^{-8}$ s⁻¹, 21 cm z momentem protonu); regresja bez zmian [736]. Resztę ~1 % zamyka rząd $\hbar^2$ tej samej fali (człon Schwarza odwzorowania, bez parametru): elementy E1 do ~10⁻⁴ QM, $\tau(np)$ 0,02–0,03 % pod QM, 2p $1{,}5957$ ns ($-1{,}1\sigma$) — na razie offline [737]. Od [738] silnik liczy elementy E1 i E2 z tej fali: $\tau(2p)$ $1{,}5957$ ns ($-1{,}07\sigma$), Ps 2P $3{,}190$ ns; zmiana czasu emisji ujawniła, że basen spinów 1s po emisji zależy od fazy precesji (linia $2P_{1/2}$ jednej kaskady −4,6 µeV poza składowe $F$) [738]. Próg kąta spinów dla basenu $F$ stanu 1s daje energię składowej $F$, ale nie rozgałęzienie: kąt przy emisji wyznacza basen $(j,F)$ stanu 2p, więc $F'=0$ wychodzi 0 / 0,125 / 0 zamiast 6j 1/3 / 2/3 / 0 — reguła zostaje opcją [739]. Stany $s$ dryfowały przy $n\ge7$ o 18–25 % [560]. Poza oknem fali EBK ($|{\rm poziom}-n|\ge0{,}05$) zostaje reguła punktu środkowego z $\omega_m$; w zbadanych sondach i scenariuszach ta ścieżka nie decyduje o żadnym fotonie, a wymiana $\omega_m$ na częstość fotonu pogarsza czasy życia stanów s [742].
4. **Przesunięcie Lamba tylko jako import** (§4.3): Ps [720], wodór [721] (Lamb $+0{,}46$ MHz, 1S–2S $-3{,}94$ MHz od pomiaru, reszta = $\alpha(Z\alpha)^6$; stałe silnika CODATA 2022 od [724]). Logarytm Bethego z widma orbity daje skalę z dokładnością do 12 %, ale nie dokładną wartość [621].
5. **Pole pośrednie:** skład stanów własnych przy $\ge3$ stanach na blok $M$ (Ps $2P$, $\sim$0,01–0,3 T) jest importem III i nie redukuje się do klasyki [602, 605–609].
6. **Tryplet $T_0$ nie ma klasycznego odpowiednika:** antyrównoległe spiny to w modelu singlet. Prawo płaszczyzny $3\gamma$ dla dowolnej polaryzacji jest gotowe [588].
7. **Zakres wyboru Borna przy powstaniu:** wodór w polu i Ps $l\ge2$ w polu używają jeszcze basenów reguły Q [619].
8. **Bilans pędu w sumie retardowanej z momentami:** reszta $\sim\beta^2$ (2–6 % siły dipolowej przy $\beta=0{,}1$–$0{,}2$) wynika z różnicy między siłą na dipol punktowy a siłą na rozłożone źródła [571–573]. W produkcji to rząd $10^{-10}$ siły Coulomba, a estymator sekularny jej nie dziedziczy. To znana niespójność reprezentacji bez skutku obserwowalnego.
9. **Poza modelem:** funkcja falowa, zasada Pauliego, poprawki radiacyjne poza #11, rozpraszanie Bhabhy, ośrodek.

---

## 7. Wnioski

Relatywistyczna klasyczna elektrodynamika dwóch ciał z momentami magnetycznymi zachowuje
wartości kwantowe przez całą kaskadę i czas życia. Dotyczy to siatki długości $L$,
całkowitego $n$ i całkowitego $J$. Daje ona też wybór kanałów E1 przy emisji i cały
hamiltonian Breita–Pauliego rzędu $1/c^2$.

Czterech rzeczy nie daje: widma (I), elementów macierzowych przy skończonym $l$ (II),
reguły Borna w obszarach chaotycznych (III) i QED w punkcie kontaktu (IV). Te braki są
zmierzone, a nie założone. Przy czterech zasadach i wyłącznie zmierzonych stałych model
odtwarza pomiary pozytonium do $10^{-4}$, strukturę $2P$ wodoru w granicach niepewności,
a granice dziedziny są nazwane i mają oszacowaną skalę. Strukturę subtelną, bez spinu i ze spinem,
łącznie z połową Thomasa, daje klasyczny hamiltonian przy jednej zaimportowanej fazie $S/\hbar$. Zostaje
tylko kwantowy czynnik $\langle r^{-3}\rangle$ [628, 706].

---

## Literatura

- [AG94] A. H. Al-Ramadhan, D. W. Gidley, Phys. Rev. Lett. 72, 1632 (1994).
- [V03] R. S. Vallery, P. W. Zitzewitz, D. W. Gidley, Phys. Rev. Lett. 90, 203402 (2003).
- [I14] A. Ishida i in., Phys. Lett. B 734, 338 (2014).
- [F93] M. S. Fee i in., Phys. Rev. Lett. 70, 1397 (1993).
- [H93] R. Hagena i in., Phys. Rev. Lett. 71, 2887 (1993).
- [L94] R. Ley i in., Hyperfine Interact. 89, 327 (1994).
- [K04] S. G. Karshenboim, Int. J. Mod. Phys. A 19, 3879 (2004) (teoria Ps, tab. 4 i 5).
- [G20] L. Gurung i in., Phys. Rev. Lett. 125, 073002 (2020).
- [S23] R. E. Sheldon i in., Phys. Rev. Lett. 131, 043001 (2023).
- [PS98] A. Pineda, J. Soto, Phys. Rev. D 59, 016005 (1998); T. Fulton, P. C. Martin, Phys. Rev. 95, 811 (1954).
- [MTN08] P. J. Mohr, B. N. Taylor, D. B. Newell, Rev. Mod. Phys. 80, 633 (2008) (CODATA 2006, teoria poziomów wodoru).
- [A13] A. Antognini i in., Science 339, 417 (2013) (promień protonu z wodoru mionowego).
- [B19] N. Bezginov i in., Science 365, 1007 (2019).
- [P11] C. G. Parthey i in., Phys. Rev. Lett. 107, 203001 (2011).
- [LY95] R. G. Littlejohn, K.-L. Yu, J. Phys. Chem. 99, 2938 (1995) (semiklasyczne sprzężenie momentów pędu, parzystość trójkąta).
- [M-WB] M. Makowiecki, *Mechanika kwantowa*, rozdz. 33 „Wstęp do teorii promieniowania kwantów pola elektromagnetycznego”, pl.wikibooks.org (CC BY-SA, wersja 558852).
- [JM05] U. D. Jentschura, P. J. Mohr, Phys. Rev. A 72, 012110 (2005) (logarytmy Bethego $n\le200$).
- [Y10] T. Yamazaki, T. Namba, S. Asai, T. Kobayashi, Phys. Rev. Lett. 104, 083401 (2010).
- [HMW55] V. W. Hughes, S. Marder, C. S. Wu, Phys. Rev. 98, 1840 (1955).
- [BG66] W. S. Bickel, A. S. Goodman, Phys. Rev. 148, 1 (1966).
- [HP94] E. W. Hagley, F. M. Pipkin, Phys. Rev. Lett. 72, 1172 (1994); [LP81] S. R. Lundeen, F. M. Pipkin, Phys. Rev. Lett. 46, 232 (1981).
- [LJP75] S. R. Lundeen, P. E. Jessop, F. M. Pipkin, Phys. Rev. Lett. 34, 377 (1975).
- [H70] H. Hellwig i in., IEEE Trans. Instrum. Meas. IM-19, 200 (1970).
- [TR77] J. S. Tiedeman, H. G. Robinson, Phys. Rev. Lett. 39, 602 (1977).
- R. Karplus, A. Klein, Phys. Rev. 87, 848 (1952); W. E. Caswell, G. P. Lepage, J. Sapirstein, Phys. Rev. Lett. 38, 488 (1977); G. S. Adkins i in., Phys. Rev. A 68, 032512 (2003).
- P. Kustaanheimo, E. Stiefel, J. reine angew. Math. 218, 204 (1965) (regularyzacja KS).
- J. B. Keller, Ann. Phys. 4, 180 (1958); M. C. Gutzwiller, J. Math. Phys. 12, 343 (1971) (kwantowanie EBK, sumy po orbitach).
- G. W. F. Drake, R. A. Swainson, Phys. Rev. A 41, 1243 (1990) (logarytmy Bethego).
- C. Quigg, J. L. Rosner, Phys. Rep. 56, 167 (1979) (tożsamość kontaktu).
