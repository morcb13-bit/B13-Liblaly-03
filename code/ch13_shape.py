"""第13章　形を見分ける ── 検算スクリプト

使い方:  python3 ch13_shape.py <b13_chain_units.py のあるフォルダ>
  網膜は 910環の図（五角形 6030 枚）。構造の点は、五角形の中心・五芒星・円環の中心。
  判定はすべて整数（4成分の番地。長さ² と内積は p + qφ で比べ、比は掛けて比べる）。

事前登録
  検定1 正多角形の角：五芒星 g を中心に、同じ円の上にある構造の点の組 (a, b) で、g から見た角が
        360°/n（n = 3, 4, 5, 6, 10）にちょうどなるものを数える（2Re(conj(a−g)(b−g)) = |a−g|²·2cos(360°/n) を整数で比べる）。
        n = 5, 10 は見つかり、n = 3, 4, 6 は一つも見つからない。外れれば NG
  検定2 三角形を探す：真ん中から |c|² ≤ R の五角形の中心から 3 つを選ぶ全部の組で、三辺の長さ² がすべて等しい組
        （正三角形）を数え、最も歪みの小さい組（(最長² − 最短²) ÷ 最長² が最小）を探す。R = 20, 40, 80, 160。
        正三角形は 0 組。歪みの最小は R を大きくするほど小さくなる。外れれば NG
  検定3 検図：図に描いた五角形の組を、網膜で調べる。
        隣（距離² 1 + φ）は頂点 2 つを共有（辺で接する）、距離² 2 + φ は頂点 1 つを共有、
        どの頂点もほかの五角形の内側（辺の上を除く）に入らない、を満たせば「設計どおり」。
        この本の平面と半径 φ² の円環は設計どおり。半径 φ（重なる）・半径 φ³（離れる）・一枚だけ 36° 回した円環は、
        どれも設計どおりにならない。外れれば NG
"""
import sys, itertools
from collections import Counter
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
cells = set(B.fits(U)); Uset = set(U)
stars = {}
for p in cells:
    for j in range(10):
        g = B.zsub(p, B.zrot(B.PHI, j))
        if g in cells or g in stars: continue
        if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)): stars[g] = j % 2
STEP = [B.zrot(B.PHI, j) for j in range(10)]
VO = {c: 1 - [j for j in range(10) if B.zadd(c, STEP[j]) in cells][0] % 2 for c in cells}   # 頂点の向き
POINTS = cells | set(stars) | Uset

def qmul(a, b): return (a[0] * b[0] + a[1] * b[1], a[0] * b[1] + a[1] * b[0] + a[1] * b[1])   # (p+qφ)(r+sφ)
def qsub(a, b): return (a[0] - b[0], a[1] - b[1])
def qcmp(a, b): return B.phi_sign(qsub(a, b))
def im_sign(w): return B.phi_sign((w[1], w[2] + w[3]))
def re2(u, v):
    """2Re(conj(u)·v) を p + qφ で。"""
    z = B.zmul(B.zconj(u), v); w = B.zadd(z, B.zconj(z))
    return (w[0] - w[2], w[2])
TWO_COS = {3: (-1, 0), 4: (0, 0), 5: (-1, 1), 6: (1, 0), 10: (0, 1)}   # 2cos(360°/n)

def verts(c, o): return [B.zadd(c, B.zt(o + 2 * i)) for i in range(5)]
def strictly_inside(c, o, p):
    V = verts(c, o)
    return all(im_sign(B.zmul(B.zconj(B.zsub(V[(i + 1) % 5], V[i])), B.zsub(p, V[i]))) > 0 for i in range(5))
def inspect(pent):
    """図に描いた五角形（中心, 頂点の向き）の組を調べる。(接し方の内訳, 内側に入った頂点の数) を返す。"""
    shared, inside = Counter(), 0
    for (a, oa), (b, ob) in itertools.combinations(pent, 2):
        n = B.norm2(B.zsub(a, b))
        if qcmp(n, (4, 0)) >= 0: continue                                        # 距離 2 以上は重ならない
        va, vb = verts(a, oa), verts(b, ob)
        shared[(n, len(set(va) & set(vb)))] += 1
        inside += sum(strictly_inside(b, ob, v) for v in va) + sum(strictly_inside(a, oa, v) for v in vb)
    return shared, inside
def as_designed(res):
    shared, inside = res
    ok = {((1, 1), 2), ((2, 1), 1)}
    return inside == 0 and len(shared) > 0 and set(shared) <= ok

if __name__ == "__main__":
    sys.path.insert(0, ".")
    import ch04_views as V4
    far = [g for g in stars if V4.to_outer(g) >= 12]
    print("0. 網膜：五角形", len(cells), "枚、構造の点", len(POINTS), "個、外周から 12 以上内側の五芒星", len(far), "個")

    print("検定1 正多角形の角（五芒星を中心に、同じ円の上の構造の点の組）")
    found = Counter()
    pts = list(POINTS)
    for g in far:
        ring = {}
        for p in pts:
            if p == g: continue
            r = B.norm2(B.zsub(p, g))
            if qcmp(r, (40, 0)) <= 0: ring.setdefault(r, []).append(B.zsub(p, g))
        for r, vs in ring.items():
            for u in vs:
                for v in vs:
                    if u == v: continue
                    w = re2(u, v)
                    for n, c in TWO_COS.items():
                        if w == qmul(r, c): found[n] += 1
    for n in (3, 4, 5, 6, 10): print("   360°/%-2d の角をなす組: %d" % (n, found[n]))
    ok1 = found[5] > 0 and found[10] > 0 and found[3] == found[4] == found[6] == 0
    print("   →", "OK" if ok1 else "NG")

    print("検定2 三角形を探す")
    prev, ok2 = None, True
    for R in (20, 40, 80, 160):
        P = [c for c in cells if qcmp(B.norm2(c), (R, 0)) <= 0]
        best, exact = None, 0
        for a, b, c in itertools.combinations(P, 3):
            s = sorted([B.norm2(B.zsub(a, b)), B.norm2(B.zsub(b, c)), B.norm2(B.zsub(c, a))], key=cmp_to_key(qcmp))
            if s[0] == s[2]: exact += 1
            d = qsub(s[2], s[0])
            if best is None or qcmp(qmul(d, best[1]), qmul(best[0], s[2])) < 0: best = (d, s[2], (a, b, c))
        d, mx, tri = best
        print("   R = %3d：五角形 %3d 枚、正三角形 %d 組、最も歪みの小さい組の (最長² − 最短², 最長²) = (%d+%dφ, %d+%dφ)"
              % (R, len(P), exact, d[0], d[1], mx[0], mx[1]), tri)
        ok2 &= exact == 0
        if prev is not None: ok2 &= qcmp(qmul(d, prev[1]), qmul(prev[0], mx)) < 0     # 歪みが小さくなる
        prev = (d, mx)
    print("   →", "OK" if ok2 else "NG")

    print("検定3 検図")
    def ring_of(R, flip=None):
        return [(B.zrot(R, k), (k % 2) ^ (1 if k == flip else 0)) for k in range(10)]
    PHI3 = B.zmul(B.PHI2, B.PHI)
    cases = [("この本の平面（6030 枚）", [(c, VO[c]) for c in cells], True),
             ("半径 φ² の円環", ring_of(B.PHI2), True),
             ("半径 φ の円環", ring_of(B.PHI), False),
             ("半径 φ³ の円環", ring_of(PHI3), False),
             ("一枚だけ 36° 回した円環", ring_of(B.PHI2, flip=3), False)]
    ok3 = True
    for name, pent, expect in cases:
        res = inspect(pent)
        got = as_designed(res)
        print("   %s：接し方 %s、内側に入った頂点 %d → %s"
              % (name, {"距離²%d+%dφ・共有頂点%d" % (n[0], n[1], s): v for (n, s), v in res[0].items()} or "なし", res[1],
                 "設計どおり" if got else "設計と違う"))
        ok3 &= got == expect
    print("   →", "OK" if ok3 else "NG")
