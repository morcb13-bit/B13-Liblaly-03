"""第5章（もう一枚の模様）・第6章（向きの番地の組）── 検算スクリプト（ch00_ring.py・ch04_join.py・ch05_grow.py と同じフォルダで走る）

使い方:  python3 ch06_name.py <b13_chain_units.py のあるフォルダ> rings_integer.json

事前登録
  検定1 円環の中心を連続のつなぎ方で結んだ網の面：正五角形とひし形だけで、ほかの形が出れば NG
  検定2 正五角形の面の重心は五芒星の位置と一致（5倍した整数4成分で比べる）。一つでも外れれば NG
  検定3 ひし形の面の短い対角線は2つ飛ばしの組。長い対角線になっていれば NG
  検定4 大きさの比：大きい図の辺² ＝ 小さい図の辺² × φ⁶（距離の比 φ³）。合わなければ NG
  検定5 指数 a0 − a1 + a2 − a3（5で割った余り）：五角形は二つの値だけで、向きの番地の偶奇と一対一。
        円環の中心は一つの値だけ。五芒星は残りの二つの値だけ（5 つの値が 2・1・2 に分かれる）。外れれば NG
  負の対照：指数の式の符号を一つ変える（a0 + a1 + a2 − a3）と、偶奇との一対一は崩れる
"""
import sys, json
from collections import Counter, defaultdict
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B
from ch00_ring import add, sub, turn, dist2, ZERO, ring
from ch05_grow import Plane

DIR = {}
for k in range(10):
    DIR[B.CONT[k]] = (9 + 2 * k) % 20          # 18° 目盛りの向きの番号
    DIR[B.SKIP[k]] = (7 + 2 * k) % 20

def faces(points, vecs):
    """点 points を、差が vecs に入る組で結んだ網の面（向きの番号の順で回る）。"""
    P = set(points)
    nb = {p: sorted([v for v in vecs if add(p, v) in P], key=lambda v: DIR[v]) for p in points}
    seen, out = set(), []
    for p in points:
        for v in nb[p]:
            if (p, v) in seen: continue
            cyc, a, e = [], p, v
            while (a, e) not in seen:
                seen.add((a, e)); cyc.append(a)
                b = add(a, e); back = tuple(-x for x in e)
                lst = nb[b]; i = lst.index(back)
                e = lst[(i - 1) % len(lst)]                 # 戻ってきた向きのすぐ手前（時計回り側）へ曲がる
                a = b
            out.append(cyc)
    return out

def idx(a): return (a[0] - a[1] + a[2] - a[3]) % 5
def idx_bad(a): return (a[0] + a[1] + a[2] - a[3]) % 5

if __name__ == "__main__":
    rec = [tuple(r) for r in json.load(open(sys.argv[2]))["rings"]]
    P = Plane(rec)
    stars = P.pentagrams()
    print("0. 60環：五角形", len(P.cells), "枚、五芒星", len(stars), "個")

    print("検定1 円環の中心を連続で結んだ網の面")
    F = faces(rec, set(B.CONT))
    F.sort(key=len)
    outer = F[-1]; inner = F[:-1]
    print("   面の頂点の数ごとの個数:", dict(Counter(len(f) for f in inner)), "  外側の輪の頂点", len(outer))
    print("   →", "OK" if set(len(f) for f in inner) <= {4, 5} else "NG")

    print("検定2 正五角形の面の重心と五芒星")
    pent = [f for f in inner if len(f) == 5]
    cen5 = {tuple(sum(p[i] for p in f) for i in range(4)) for f in pent}
    star5 = {tuple(5 * x for x in g) for g in stars}
    print("   正五角形の面", len(pent), "個、重心×5 と 五芒星×5 の一致:", cen5 == star5)
    print("   正五角形の頂点（円環の中心）から重心（五芒星）までの距離²:",
          {B.norm2(sub(tuple(5 * x for x in f[0]), c)) for f, c in
           [(f, tuple(sum(p[i] for p in f) for i in range(4))) for f in pent]}, "（5倍した長さの2乗。25φ⁶ = (125, 200)）")
    print("   →", "OK" if cen5 == star5 else "NG")

    print("検定3 ひし形の面の対角線")
    rh = [f for f in inner if len(f) == 4]
    short = [(f[0], f[2]) if sub(f[2], f[0]) in set(B.SKIP) else (f[1], f[3]) for f in rh]
    ok3 = all(sub(b, a) in set(B.SKIP) for a, b in short)
    longd = [(f[1], f[3]) if sub(f[2], f[0]) in set(B.SKIP) else (f[0], f[2]) for f in rh]
    print("   ひし形の面", len(rh), "個。短い対角線がすべて2つ飛ばし:", ok3,
          "  長い対角線の距離²:", {B.norm2(sub(b, a)) for a, b in longd})
    print("   →", "OK" if ok3 else "NG")

    print("検定4 大きさの比")
    s2 = B.norm2(sub(add(ZERO, turn((1, 0, 0, 0), 2)), (1, 0, 0, 0)))    # 五角形の辺²（頂点 1 と ζ² の距離²）
    big = B.N_CONT                                                        # 網の辺²（連続）
    phi6 = (5, 8)
    prod = (s2[0] * phi6[0] + s2[1] * phi6[1], s2[0] * phi6[1] + s2[1] * phi6[0] + s2[1] * phi6[1])
    print("   五角形の辺² =", s2, "  ×φ⁶ =", prod, "  網の辺² =", big, "  一致:", prod == big)
    print("   五角形の外接半径² = 1、網の正五角形の外接半径² = 5 + 8φ ＝ φ⁶")
    print("   →", "OK" if prod == big else "NG")

    print("検定5 指数 a0 − a1 + a2 − a3（5で割った余り）")
    c = Counter((idx(p), a % 2) for p, a in P.cells.items())
    print("   五角形（指数, 向きの番地の偶奇）:", dict(c))
    print("   円環の中心:", dict(Counter(idx(r) for r in rec)), "  五芒星（指数, 向き k）:", dict(Counter((idx(g), k) for g, k in stars.items())))
    def split(cells, rings, sts):
        pv = {a % 2: {idx(p) for p, b in cells.items() if b % 2 == a % 2} for a in (0, 1)}
        one = all(len(v) == 1 for v in pv.values()) and pv[0] != pv[1]
        rv = {idx(r) for r in rings}; sv = {idx(g) for g in sts}
        return one and len(rv) == 1 and len(sv) == 2 and (pv[0] | pv[1] | rv | sv) == {0, 1, 2, 3, 4}
    print("   →", "OK" if split(P.cells, rec, stars) else "NG")
    cb = Counter((idx_bad(p), a % 2) for p, a in P.cells.items())
    print("   負の対照（a0 + a1 + a2 − a3）:", dict(cb))

    print("6. 910環の図でも同じか")
    T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
    rows, place, offs = B.build_stack()
    U, seen = [], set()
    for k in range(10):
        t = B.zrot(T if k % 2 == 0 else S, k)
        for cc in sum(place, []):
            c2 = B.zadd(B.zrot(cc, k), t)
            if c2 not in seen: seen.add(c2); U.append(c2)
    cells = B.fits(U)
    st = {}
    for p in cells:
        for j in range(10):
            g = B.zsub(p, B.zrot(B.PHI, j))
            if g in cells or g in st: continue
            if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)): st[g] = j % 2
    print("   五角形（指数, 偶奇）:", dict(Counter((idx(p), a % 2) for p, a in cells.items())),
          "  円環の中心:", dict(Counter(idx(r) for r in U)), "  五芒星:", dict(Counter(idx(g) for g in st)))
    print("   →", "OK" if split(cells, U, st) else "NG")
    print("   60環は 910環を 36° 回したものの一部：36° 回すと指数はそろって一つずれる（60環で 1,3 / 2 / 0,4 → 910環で 2,4 / 3 / 0,1）")
