"""第1章　整数を並べる ── 検算スクリプト

使い方:  python3 ch01_numbers.py <b13_chain_units.py のあるフォルダ>
  910 環の平面（円環を 1 から 13 まで積んだ三角形を 10 方向に組んだもの）の上で調べる。
  判定はすべて整数（4成分の番地、長さ² は p + qφ）。

事前登録
  検定T0 扇の頂点：どの扇でも、段 1 の円環から扇の軸を段 2 と反対向きに距離 φ³ 進んだ所に五芒星がある
         （扇の頂点の五芒星＝最初の1、段 1 の円環＝二番目の1）。頂点の五芒星は、5 枚の扇では真ん中の (0,0,0,0)、
         残る 5 枚では真ん中から距離² 13 + 21φ。段 1 の円環は 10 個とも違い、距離² は 5 + 8φ が 5 個、34 + 55φ が 5 個（ノルム 1）。
         段の数 1〜13 のうちフィボナッチ数は 1, 2, 3, 5, 8, 13。外れれば NG
  検定T1 段 n の円環はちょうど n 個（n = 1..13）、三角形の円環は 91 個。10 方向に組むと重ならずに 910 個。外れれば NG
  検定T5 左右対称：頂点 (0,0,0,0) を通る扇の軸で鏡に写す（z → ζ^m · z̄）と、ある m ただ一つで、
         三角形の円環の集まりも、各段の円環の集まりも、平面の五角形の集まりも、それぞれ自分に重なる。
         ほかの m（負の対照）では三角形が自分に重ならない。外れれば NG
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
    print("検定T0 扇の頂点は最初の1")
    from collections import Counter
    T, S0 = (-2, 2, 0, 3), (-5, 5, 0, 8)
    apex, tips = [], []
    def im_sign(w): return B.phi_sign((w[1], w[2] + w[3]))
    for k in range(10):
        t = B.zrot(T if k % 2 == 0 else S0, k)
        a = B.zadd(B.zrot(place[0][0], k), t)
        r1, r2 = [B.zadd(B.zrot(c, k), t) for c in place[1]]
        d = B.zsub(B.zadd(a, a), B.zadd(r1, r2))                 # 段 2 の真ん中から段 1 へ向かう向き（2 倍）
        ps = [p for p in stars if B.norm2(B.zsub(p, a)) == (5, 8) and im_sign(B.zmul(B.zconj(d), B.zsub(p, a))) == 0
              and B.phi_sign((lambda z: (z[0] - z[2], z[2]))(B.zadd(B.zmul(B.zconj(d), B.zsub(p, a)), B.zconj(B.zmul(B.zconj(d), B.zsub(p, a)))))) > 0]
        apex.append(a); tips.append(ps[0] if len(ps) == 1 else None)
    print("   扇の頂点の五芒星:", "すべての扇で 1 個" if all(tips) else "見つからない扇がある",
          "  真ん中からの距離²:", dict(Counter(B.norm2(p) for p in tips if p)))
    dist = [B.norm2(a) for a in apex]
    cnt = Counter(dist)
    nrm = {d: d[0] * d[0] + d[0] * d[1] - d[1] * d[1] for d in cnt}
    fibs = [n for n in range(1, 14) if n in {fib(k) for k in range(1, 9)}]
    print("   段 1 の円環（二番目の1）:", len(set(apex)), "個  真ん中からの距離²:", dict(cnt), "  ノルム:", nrm)
    print("   段の数 1〜13 のうちフィボナッチ数:", fibs)
    ok0 = all(tips) and Counter(B.norm2(p) for p in tips) == Counter({(0, 0): 5, (13, 21): 5}) and len(set(apex)) == 10 and cnt == Counter({(5, 8): 5, (34, 55): 5}) and set(nrm.values()) == {1} and fibs == [1, 2, 3, 5, 8, 13]
    print("   →", "OK" if ok0 else "NG")
    print("検定T1 三角形は数そのもの")
    sizes = [len(r) for r in place]
    print("   段ごとの円環の数:", sizes, " 合計:", sum(sizes), "  10 方向に組んだ円環:", len(U))
    ok1 = sizes == list(range(1, 14)) and sum(sizes) == 91 and len(U) == 910
    print("   →", "OK" if ok1 else "NG")

    print("検定T5 左右対称（扇の軸で鏡に写す）")
    tri0 = [[B.zadd(c, T) for c in r] for r in place]
    S = set(sum(tri0, []))
    good = [m for m in range(10) if {B.zrot(B.zconj(c), m) for c in S} == S]
    m = good[0] if good else None
    rowok = m is not None and all({B.zrot(B.zconj(c), m) for c in r} == set(r) for r in tri0)
    cellok = m is not None and {B.zrot(B.zconj(c), m) for c in cells} == cells
    print("   三角形が自分に重なる m:", good, "  各段も重なる:", rowok, "  平面の五角形 %d 枚も重なる:" % len(cells), cellok,
          "  頂点は動かない:", m is not None and B.zrot(B.zconj((0, 0, 0, 0)), m) == (0, 0, 0, 0))
    ok5 = len(good) == 1 and rowok and cellok
    print("   →", "OK" if ok5 else "NG")
