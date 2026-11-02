"""第17章　初期値の違いと折り返し（1, 1 から生まれる組み合わせ）── 検算スクリプト

使い方:  python3 ch17_residue.py <b13_chain_units.py のあるフォルダ>
  最初の 1 ＝ 頂点の五芒星（真ん中 (0,0,0,0)）、二番目の 1 ＝ 段 1 の円環 (−2, 2, 0, 3)（第1章）。
  円環をつなぐ一歩（連続：距離² 7 + 11φ、2つ飛ばし：3 + 4φ）を足すことだけで、910 環の平面の円環をたどる。
  残渣 ＝ 番地の ζ を ζ³ に読み替えた点。判定はすべて整数。

事前登録
  F1 足し算だけで届く：段 1 の円環から一歩を足していくだけで、910 環すべてに届く。足した回数の内訳を記録する
  F2 二度と同じにならない：910 環の残渣は 4 つの整数として全部違う。対照：残渣の長さ² だけ比べると同じ値が出る
  F3 平面では広がり、残渣は小さいまま：残渣の長さ² の最大の 1000 倍 ＜ 平面の長さ² の最大（p + qφ の整数比較）
  F4 平面では戻らず、番地では折り返す：二番目の 1 に φ を n 回掛けた点（n = 1〜60）は平面で二度と同じ場所に来ない。
     13 の番地で読むとちょうど 28 回ごとに元に戻る。負の対照：13 で割らなければ戻らない
"""
import sys
from collections import Counter, deque
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

T0 = (-2, 2, 0, 3); S0 = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T0 if k % 2 == 0 else S0, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
US = set(U)
ONE2 = T0                                          # 二番目の 1（段 1 の円環）
Z3 = [B.zt(3 * k) for k in range(4)]               # ζ^0, ζ^3, ζ^6, ζ^9
def sigma(a):
    r = (0, 0, 0, 0)
    for k in range(4):
        if a[k]: r = B.zadd(r, B.zmul((a[k], 0, 0, 0), Z3[k]))
    return r
def qlt(a, b): return B.phi_sign((a[0] - b[0], a[1] - b[1])) < 0

STEPS = sorted({B.zsub(b, a) for a in U for b in U if B.norm2(B.zsub(b, a)) in ((7, 11), (3, 4))})

if __name__ == "__main__":
    print("最初の 1：頂点の五芒星 (0, 0, 0, 0)   二番目の 1：段 1 の円環", ONE2, " 距離²", B.norm2(ONE2))
    print("円環をつなぐ一歩：", len(STEPS), "本（連続", sum(1 for s in STEPS if B.norm2(s) == (7, 11)), "・2つ飛ばし", sum(1 for s in STEPS if B.norm2(s) == (3, 4)), "）")
    print("F1 足し算だけで届く")
    gen = {ONE2: 0}; q = deque([ONE2])
    while q:
        c = q.popleft()
        for s in STEPS:
            x = B.zadd(c, s)
            if x in US and x not in gen: gen[x] = gen[c] + 1; q.append(x)
    cnt = Counter(gen.values())
    print("   届いた円環:", len(gen), "/", len(U), "  足した回数ごとの数:", [cnt[i] for i in range(max(cnt) + 1)])
    ok1 = len(gen) == len(U)
    print("   →", "OK" if ok1 else "NG")
    print("F2 二度と同じにならない")
    res = [sigma(c) for c in U]
    ok2a = len(set(res)) == len(U)
    lens = Counter(B.norm2(r) for r in res)
    print("   残渣（4 つの整数）の種類:", len(set(res)), "/", len(U), "  対照：残渣の長さ² の種類", len(lens), "（同じ長さ² が最多", max(lens.values()), "個）")
    ok2 = ok2a and len(lens) < len(U)
    print("   →", "OK" if ok2 else "NG")
    print("F3 平面では広がり、残渣は小さいまま")
    mr = max((B.norm2(c) for c in U), key=lambda v: (v[0] + v[1] * 1.618033988749895))   # 候補を選ぶだけ（判定は下の整数比較）
    ms = max((B.norm2(r) for r in res), key=lambda v: (v[0] + v[1] * 1.618033988749895))
    ok_max = all(not qlt(mr, B.norm2(c)) for c in U) and all(not qlt(ms, B.norm2(r)) for r in res)
    ok3 = ok_max and qlt((1000 * ms[0], 1000 * ms[1]), mr)
    print("   平面の長さ² の最大:", mr, "  残渣の長さ² の最大:", ms, "  （最大であることの整数確認:", ok_max, "）  1000 倍しても平面より小さい:", ok3)
    print("   →", "OK" if ok3 else "NG")
    print("F4 平面では戻らず、番地では折り返す")
    pts = [ONE2]; m13 = lambda z: tuple(x % 13 for x in z)
    for n in range(60): pts.append(B.zmul(pts[-1], B.PHI))
    distinct = len(set(pts)) == len(pts)
    adr = [m13(p) for p in pts]
    per = next(k for k in range(1, 61) if adr[k] == adr[0])
    okper = all(adr[i] == adr[i + per] for i in range(len(adr) - per))
    print("   φ を 0〜60 回掛けた 61 点：平面で全部違う場所:", distinct, "  13 の番地で初めて元に戻る回数:", per, "  以後も", per, "回ごと:", okper)
    print("   平面の長さ²（n = 0, 1, 2, 3, 28）:", [B.norm2(pts[n]) for n in (0, 1, 2, 3, 28)])
    print("   残渣の長さ²（n = 0, 1, 2, 3, 28）:", [B.norm2(sigma(pts[n])) for n in (0, 1, 2, 3, 28)])
    ok4 = distinct and per == 28 and okper
    print("   →", "OK" if ok4 else "NG")
