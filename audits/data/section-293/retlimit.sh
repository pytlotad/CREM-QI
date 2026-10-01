#!/bin/bash
# Audit 293: skan progu retardacyjnego przeciw wyprowadzeniu z 292d.
#
# 292d: dla orbity kolowej T/t_light = 2 pi c sqrt(a/k), wiec
#   a(X) = k (X/(2 pi c))^2        -- KWADRATOWO w progu
#   L(X) = mu sqrt(k a) = mu k X/(2 pi c)  -- LINIOWO w progu
# Kalibracja na progu domyslnym: a(150)=16,617 r*, L(150)=0,17421 hbar.
#
# ZAREJESTROWANE PRZED URUCHOMIENIEM:
#   R293a: r_p/r* ma isc jak X^2, czyli r_p(X)/r_p(150) = (X/150)^2.
#   R293b: L ma isc jak X, czyli L(X) = 0,17421 * X/150.
#   R293c: kontakt (L <= 0,042737 hbar) ma byc pierwszy raz osiagniety
#          przy X = 150*0,042737/0,17421 = 36,80, czyli
#          contactPassageAtBarrierSeconds przestanie byc NaN miedzy
#          X = 40 i X = 36.  To jest przewidywanie ILOSCIOWE i
#          falsyfikowalne co do jednej cyfry.
#   R293d: przyczyna stopu ma przejsc z RetardationLimit na
#          ComptonBarrier gdy a(X) spadnie pod 1 r*, czyli przy
#          X = 150/sqrt(16,617) = 36,80 -- TA SAMA liczba, bo dla
#          orbity kolowej oba warunki pokrywaja sie dokladnie tam.
#          To wyjasnia obserwacje 292d o T/t=36,1 i jest jej testem.
#
# Model CIAGLY, bo tylko on konczy w budzecie (292c).  Wyniki sa
# probami ksiegowosci modelu, nie twierdzeniami fizycznymi -- tak je
# kwalifikuje komentarz przy CREM_RETARDATION_LIMIT.
S="$1"
for X in 150 100 70 50 40 36 30; do
    python3 -c "print('### prog=%s  przewidywane r_p=%.3f r*  L=%.5f hbar'
        %($X,16.617*($X/150.0)**2,0.17421*$X/150.0))"
    CREM_CONTINUOUS=1 CREM_RETARDATION_LIMIT=$X timeout 900 "$S/stopcause" 2 300 1 \
        | grep -v "^# nan\|^# ziaren\|^# bariera\|^ziarno"
done
