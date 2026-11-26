"""Run every ancillary script of the article and compare each script's summary with expected.json.

Usage: python run_all.py [--skip-long | --only-long] [script names or tags...]
Each entry of expected.json: {"script": ..., "args": [...], "cwd": "kod", "long": false,
                              "expect": [regex, ...], "forbid": [regex, ...],
                              "agree": [{"regex": ..., "tol": ..., "min_count": n}, ...], "save_stdout": file}
Every "expect" regex must match the output and no "forbid" regex may match it. For each "agree" entry, every match
of the regex must have |a - b| <= tol (named groups a and b) or |d| <= tol (named group d), and there must be at
least min_count matches; this checks whole printed tables of exact against numerical values.
Exit code 1 if any script fails, times out, or misses an expected line. Logs go to ci_logs/.
"""
import json, re, subprocess, sys, time, os
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# a failed check is printed as 'FEL ...' or 'FAIL ...' by the scripts; a crash as a traceback
DEFAULT_FORBID = [r'^\s*(FEL|FAIL)\b', r'^Traceback']
def main(argv):
    skip_long = '--skip-long' in argv; only_long = '--only-long' in argv
    names = [a for a in argv if not a.startswith('--')]
    cases = json.loads((ROOT / 'expected.json').read_text(encoding='utf-8'))
    (ROOT / 'ci_logs').mkdir(exist_ok=True)
    bad = []
    for c in cases:
        tag = c.get('tag', Path(c['script']).stem)
        if names and c['script'] not in names and tag not in names: continue
        if skip_long and c.get('long'): continue
        if only_long and not c.get('long'): continue
        t = time.time()
        try:
            p = subprocess.run([sys.executable, '-X', 'utf8', c['script'], *c.get('args', [])],
                               cwd=ROOT / c.get('cwd', 'kod'), capture_output=True, text=True,
                               encoding='utf-8', timeout=c.get('timeout', 3600))
            out = p.stdout + p.stderr; rc = p.returncode
        except subprocess.TimeoutExpired as e:
            out = f'TIMEOUT {e}'; rc = -1
        (ROOT / 'ci_logs' / f'{tag}.log').write_text(out, encoding='utf-8')
        if c.get('save_stdout'):   # some later scripts read this output as a file
            (ROOT / c.get('cwd', 'kod') / c['save_stdout']).write_text(p.stdout if rc != -1 else out, encoding='utf-8')
        missing = [r for r in c.get('expect', []) if not re.search(r, out, re.M)]
        missing += ['forbidden: ' + r for r in DEFAULT_FORBID + c.get('forbid', []) if re.search(r, out, re.M)]
        for ag in c.get('agree', []):
            ms = list(re.finditer(ag['regex'], out, re.M)); tol = Decimal(str(ag['tol']))
            bad_m = [m.group(0) for m in ms if (abs(Decimal(m.group('d'))) if 'd' in m.groupdict()
                     else abs(Decimal(m.group('a')) - Decimal(m.group('b')))) > tol]
            if len(ms) < ag.get('min_count', 1): missing.append(f"agree: {len(ms)} matches of {ag['regex']!r}")
            missing += [f'agree: tolerance {tol} exceeded: {x}' for x in bad_m]
        ok = rc == 0 and not missing
        print(f"{'OK ' if ok else 'FAIL'} {tag:40s} {time.time()-t:7.1f}s" + (f'  rc={rc} missing={missing}' if not ok else ''))
        if not ok: bad.append(tag)
    print(f'{len(bad)} failure(s)' if bad else 'ALL OK')
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
