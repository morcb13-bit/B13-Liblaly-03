"""第14章　ずらすと窓になる ── 検算スクリプト

使い方:  python3 ch14_window.py <b13_chain_units.py のあるフォルダ>
  平面（五角形 6030 枚）を、番地 d だけずらした写しと重ね、両方に五角形がある場所（重なり）を数える。
  m(d) = d だけずらしても五角形が五角形に重なる枚数。ずらし d は、4 成分がどれも −4〜4 の番地 6560 本（0 を除く）。
  番地の組：a0 − a1 + a2 − a3 を 5 で割った余り（第6章）。残渣：ζ → ζ³ と読み替えた番地（第4章）。
  長さ² は p + qφ、残渣の長さ² は (p + q) − qφ。判定はすべて整数。

事前登録
  検定1 二つの組：五角形の番地の組は 2 種類だけで、枚数は等しい。ずらしの番地の組が 1 か 4（組をまたいで
        どちらにも行けない）なら、重なりは必ず 0。外れれば NG
  検定2 窓の開き：組を変えないずらし（番地の組 0）では、重なりの枚数は残渣の長さ² だけで決まり（同じ残渣なら同じ枚数）、
        残渣が長いほど少なくなり、ある長さを越えると 0 になる（窓が閉じる）。外れれば NG
  検定3 短いずらしほど閉じる：ずらしの長さ² × 残渣の長さ² は整数（ノルム N）になる。組を変えないずらしのうち、
        実際の長さが最も短いものは重ならず、重なりの最も大きいずらしは実際の長さが長い。外れれば NG
"""
import sys, itertools
from collections import Counter, defaultdict
from functools import cmp_to_key
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

T = (-2, 2, 0, 3); S0 = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S0, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = set(B.fits(U))

def group(d): return (d[0] - d[1] + d[2] - d[3]) % 5
def residue2(d):
    n = B.norm2(d); return (n[0] + n[1], -n[1])
def qcmp(a, b): return B.phi_sign((a[0] - b[0], a[1] - b[1]))
def norm(x): return x[0] * x[0] + x[0] * x[1] - x[1] * x[1]       # (p+qφ)((p+q)−qφ)
def overlap(d): return sum(1 for c in cells if B.zadd(c, d) in cells)

SHIFTS = [d for d in itertools.product(range(-4, 5), repeat=4) if d != (0, 0, 0, 0)]

if __name__ == "__main__":
    print("0. 平面：五角形", len(cells), "枚。ずらし", len(SHIFTS), "本")
    M = {d: overlap(d) for d in SHIFTS}

    print("検定1 二つの組")
    G = Counter(group(c) for c in cells)
    byg = defaultdict(Counter)
    for d, m in M.items(): byg[group(d)][m > 0] += 1
    print("   五角形の番地の組:", dict(sorted(G.items())))
    for g in range(5): print("   ずらしの組 %d：重なる %4d 本、重ならない %4d 本" % (g, byg[g][True], byg[g][False]))
    ok1 = len(G) == 2 and len(set(G.values())) == 1 and byg[1][True] == 0 and byg[4][True] == 0
    print("   →", "OK" if ok1 else "NG")

    print("検定2 窓の開き（組を変えないずらし）")
    I0 = [d for d in SHIFTS if group(d) == 0]
    byr = defaultdict(set)
    for d in I0: byr[residue2(d)].add(M[d])
    order = sorted(byr, key=cmp_to_key(qcmp))
    single = all(len(v) == 1 for v in byr.values())
    seq = [next(iter(byr[r])) for r in order]
    opened = [r for r in order if next(iter(byr[r])) > 0]
    closed = [r for r in order if next(iter(byr[r])) == 0]
    for r in order[:16]:
        print("   残渣² %3d%+4dφ ：重なり %s 枚" % (r[0], r[1], sorted(byr[r])))
    dec = all(seq[i] > seq[i + 1] for i in range(len(opened) - 1)) and all(x == 0 for x in seq[len(opened):])
    print("   同じ残渣なら同じ枚数:", single, "  開いている残渣の最大:", opened[-1], "  閉じている残渣の最小:", closed[0])
    print("   残渣が長いほど少なく、越えると 0:", dec)
    ok2 = single and dec and qcmp(opened[-1], closed[0]) < 0
    print("   →", "OK" if ok2 else "NG")

    print("検定3 短いずらしほど閉じる")
    prod_ok = all(norm(B.norm2(d)) >= 1 for d in SHIFTS)
    by_len = sorted(I0, key=cmp_to_key(lambda a, b: qcmp(B.norm2(a), B.norm2(b))))
    shortest = B.norm2(by_len[0])
    short_all = [d for d in I0 if B.norm2(d) == shortest]
    best = max(I0, key=lambda d: M[d])
    first_open = next(d for d in by_len if M[d] > 0)
    print("   長さ² × 残渣² ＝ N（整数）が 1 以上:", prod_ok)
    print("   最も短いずらし：長さ² %d%+dφ（%d 本）、残渣² %d%+dφ、N = %d、重なり %s"
          % (shortest[0], shortest[1], len(short_all), *residue2(by_len[0]), norm(shortest), sorted({M[d] for d in short_all})))
    print("   最初に重なるずらし：長さ² %d%+dφ、重なり %d 枚" % (*B.norm2(first_open), M[first_open]))
    print("   重なりの最大：", best, "長さ² %d%+dφ、残渣² %d%+dφ、N = %d、重なり %d / %d 枚"
          % (*B.norm2(best), *residue2(best), norm(B.norm2(best)), M[best], len(cells)))
    ok3 = prod_ok and all(M[d] == 0 for d in short_all) and qcmp(B.norm2(best), B.norm2(first_open)) > 0
    print("   →", "OK" if ok3 else "NG")
