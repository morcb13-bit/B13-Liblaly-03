"""第2章　2つを足す ── 検算スクリプト

使い方:  python3 ch02_add.py <b13_chain_units.py のあるフォルダ>
  判定はすべて整数（4成分の番地、長さ² は p + qφ）。

事前登録
  検定F1 頂点の五芒星の 1 と段 1 の 1 から、隣り合う 2 つを足していくと 1, 1, 2, 3, 5, 8, 13 で、どれも 13 段の中の段になり、
         次の 21 は 13 を越える。段の円環の数も同じ足し算になる（段 a ＋ 段 b の円環 ＝ 段 a + b の円環）。外れれば NG
  検定T4 五芒星から測る：φ^n を p + qφ で書くと (F(n−1), F(n))（n = 0..13）。真ん中の五芒星の位置には五角形がない。
         そこから φ^n（n = 1..6）だけ離れた 5 点は、向きをそろえるか 36° 回すかのどちらかで、すべて構造の点
         （五角形の中心・五芒星・円環の中心）に着地する。外れれば NG
  検定F2 一歩のノルム p² + pq − q²：頂点の五芒星から段 1 の円環への一歩（距離² 5 + 8φ と 34 + 55φ）は 1、
         円環どうしの連続（7 + 11φ）と 2つ飛ばし（3 + 4φ）は 5。外れれば NG
"""
import sys
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, ".")
import b13_chain_units as B
import ch04_side as C5

T = (-2, 2, 0, 3); S0 = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen, tri = [], set(), []
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S0, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2); tri.append(k)
cells = set(B.fits(U))
stars = set()
for p in cells:
    for j in range(10):
        g = B.zsub(p, B.zrot(B.PHI, j))
        if g in cells or g in stars: continue
        if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)): stars.add(g)
POINTS = cells | stars | set(U)

def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a

if __name__ == "__main__":
    print("検定F1 2 つを足す")
    seq = [1, 1]
    while seq[-1] + seq[-2] <= 13: seq.append(seq[-1] + seq[-2])
    nxt = seq[-1] + seq[-2]
    sizes = [len(r) for r in place]
    cnt = {n: sizes[n - 1] for n in range(1, 14)}
    adds = [(seq[i], seq[i + 1], seq[i + 2]) for i in range(len(seq) - 2)]
    ring = lambda i, v: 1 if i == 0 else cnt[v]       # 最初の 1 は頂点の五芒星（1 個）
    lines = ["%s %d ＋ 円環 %d ＝ 円環 %d" % ("五芒星" if i == 0 else "円環", ring(i, a), cnt[b], cnt[c]) for i, (a, b, c) in enumerate(adds)]
    for s in lines: print("   " + s)
    print("   足していった列:", seq, "  次:", nxt)
    okf = seq == [1, 1, 2, 3, 5, 8, 13] and nxt == 21 and all(ring(i, a) + cnt[b] == cnt[c] for i, (a, b, c) in enumerate(adds))
    print("   →", "OK" if okf else "NG")
    print("検定T4 五芒星から測る（φⁿ ＝ (F(n−1), F(n))）")
    p, ok3 = (1, 0), True                                        # φ⁰ = 1
    pows = []
    for n in range(14):
        pows.append(p); ok3 &= p == (fib(n - 1) if n else 1, fib(n))
        p = (p[1], p[0] + p[1])                                  # ×φ：(p, q) → (q, p + q)
    print("   φ^n = (F(n−1), F(n))  n = 0..13:", pows)
    O = (0, 0, 0, 0)
    print("   真ん中 (0,0,0,0)：五芒星", O in stars, "  五角形", O in cells)
    ok3 &= O in stars and O not in cells
    R = B.ONE
    for n in range(1, 7):
        R = B.zmul(R, B.PHI)
        land = [off for off in (0, 1) if all(B.zrot(R, 1 + off + 2 * i) in POINTS for i in range(5))]
        print("   φ^%d = %d + %dφ の 5 点：着地する向き %s" % (n, pows[n][0], pows[n][1], ["そろえる", "36°回す"][land[0]] if land else "なし"))
        ok3 &= bool(land)
    print("   →", "OK" if ok3 else "NG")

    print("検定F2 一歩のノルム")
    nrm = lambda d: d[0] * d[0] + d[0] * d[1] - d[1] * d[1]
    steps = [("五芒星→段 1 の円環（φ³）", (5, 8)), ("五芒星→段 1 の円環（φ⁵）", (34, 55)), ("円環の連続", (7, 11)), ("円環の 2つ飛ばし", (3, 4))]
    for s, d in steps: print("   %-18s 距離² %s  ノルム %d" % (s, d, nrm(d)))
    okn = [nrm(d) for _, d in steps] == [1, 1, 5, 5]
    print("   →", "OK" if okn else "NG")
