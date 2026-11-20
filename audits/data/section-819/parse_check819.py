# Audit 819 gate: every quantity of H819s must be found in the OLD files (run before reading new).  Usage: python3 parse_check819.py DIR PREFIX
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import stat819 as S
d, pre = sys.argv[1], sys.argv[2]
miss = S.check({t: f'{d}/{pre}{t}.out' for t in S.PAT})
for t in S.PAT: print(t, {k: v[1] for k, v in S.parse(f'{d}/{pre}{t}.out', t).items()})
print('PARSER OK' if not miss else f'PARSER BRAK {miss}')
