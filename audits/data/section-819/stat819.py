# Audit 819: parser of the regression outputs (written before the runs; checked on the 818 old files for the format only).
import re
F = re.S
PAT = {  # set -> list of (key, regex, kind); kind 'pm' = value +/- sigma, 'eq' = must equal old, 'zero' = must be 0 in old and new
    '1':   [('collapse', r'mean of completed runs (\S+) \+/- (\S+) ', 'pm'), ('completion', r'Completion fraction: (\S+)%', 'eq'),
            ('tau_lab', r'annihilation time on the lab clock: tau \(MLE\) (\S+) \+/- (\S+)', 'pm')],
    '2':   [('collapse', r'mean of completed runs (\S+) \+/- (\S+) ', 'pm'), ('completion', r'Completion fraction: (\S+)%', 'eq'),
            ('tau_lab', r'annihilation time on the lab clock: tau \(MLE\) (\S+) \+/- (\S+)', 'pm')],
    'L2':  [('cascade', r'Photon-cascade time.*?mean\s+(\S+) \+/- (\S+) \w+ over \d+ completed runs', 'pm')],
    '5Ps': [('sigma_capture', r'sigma_capture = (\S+) \+/- (\S+) barn', 'pm'),
            ('sigma_over_K', r'<sigma/sigma_Kramers\(K\)> event by event \(audit 551\) = (\S+) \+/- (\S+)', 'pm')],
    '6':   [('tau_para', r'para\s*: \d+ decays \(\d+ failed\).*?tau \(MLE = mean t\) (\S+) \+/- (\S+)', 'pm'),
            ('tau_ortho', r'ortho\s*: \d+ decays \(\d+ failed\).*?tau \(MLE = mean t\) (\S+) \+/- (\S+)', 'pm'),
            ('failed_para', r'para\s*: \d+ decays \((\d+) failed\)', 'zero'), ('failed_ortho', r'ortho\s*: \d+ decays \((\d+) failed\)', 'zero')]}
def parse(path, t):
    s = open(path).read(); out = {}
    for k, rx, kind in PAT[t]:
        m = re.search(rx, s, F)
        if m: out[k] = (kind, tuple(float(g) for g in m.groups()))
    return out
def check(paths):  # paths: set -> file; returns list of missing (set, key)
    return [(t, k) for t, p in paths.items() for k, _, _ in PAT[t] if k not in parse(p, t)]
