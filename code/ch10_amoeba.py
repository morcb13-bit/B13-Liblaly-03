"""第10章　アメーバは最短路を知っている ── 検算スクリプト

使い方:  python3 ch10_amoeba.py <b13_chain_units.py のあるフォルダ>
  910環の図の五角形 6030 枚を地図とし、辺を共有する五角形（中心の距離 φ）だけを隣とする。
  餌の匂いとアメーバの体は、どちらも一刻ごとに隣へ広がる。広がるとき、量は隣へそのまま足すだけ。
  判定はすべて整数（4成分の番地、距離² は p + qφ、匂いの量は整数の足し算）。

事前登録
  検定1 隣：どの五角形も隣は 2, 3, 4 枚のどれか。地図は一つにつながり、
        五角形は2色に塗り分けられる（隣どうしは必ず別の色）。外れれば NG
  検定2 匂いの量＝道の数：餌から匂いを足し算で広げると、ある五角形に初めて届いた刻の量は、
        餌からそこまでの最短の道の本数（層ごとに数えた別の計算）にすべて一致する。
        負の対照：2刻あとの量は、遠回りの道が加わって一致しない
  検定3 二つの波が出会う：餌から匂い、アメーバから体を同時に広げる。届いた刻の和が最短距離 d に
        等しい五角形（帯）について、体から k 刻の帯の五角形の「体の量 × 匂いの量」を足すと、
        k = 0..d のどの段でも最短の道の総数に等しい。負の対照：和が d + 2 の五角形では等しくならない
  検定4 しっぽが縮む：帯の中で、体の量 × 匂いの量が最も大きい隣へ一歩ずつ進むと、d 歩ちょうどで餌に着く。
        負の対照：匂いの量（足し続けた量）だけを見て濃い隣へ進むと、d 歩で着かない組が出る
  検定5 障害物：五角形を外して壁を2枚置く（壁は、向き 0°・72°・108° への射影 p + qφ の比較だけで選ぶ）。
        地図はつながったまま。壁のない地図で残った道は壁を通る。壁のある地図で検定2〜4 をやり直すと
        すべて成り立ち、残った道は壁を一枚も通らず、最短距離は壁のないときより長い。外れれば NG
  試す組：乱数の種 10 で選んだ 200 組（アメーバ, 餌）
"""
import sys, random
from collections import Counter, deque
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

# 910環の図の五角形（第4章・第9章と同じ組み方）
T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = set(B.fits(U))
STEP = [B.zrot(B.PHI, j) for j in range(10)]            # 中心の距離 φ：辺を共有する隣
nb = {c: [B.zadd(c, s) for s in STEP if B.zadd(c, s) in cells] for c in cells}

def spread(src):
    """src から一刻ずつ隣へ足し算で広げる。初めて届いた刻と、その刻の量を返す。"""
    t_first, amount = {src: 0}, {src: 1}
    W, t = {src: 1}, 0
    hist = [dict(W)]
    while len(t_first) < len(cells):
        W2 = Counter()
        for u, w in W.items():
            for v in nb[u]: W2[v] += w
        W, t = W2, t + 1
        hist.append(dict(W))
        for v, w in W.items():
            if v not in t_first: t_first[v] = t; amount[v] = w
    return t_first, amount, hist

# 壁：向き 36k° への射影の2倍 2Re(c ζ^(−k)) を p + qφ で書き、整数の比較だけで選ぶ
COS2 = [(2, 0), (0, 1), (-1, 1), (1, -1), (0, -1), (-2, 0), (0, -1), (1, -1), (-1, 1), (0, 1)]
def proj2(c, k):
    p = q = 0
    for i in range(4):
        a, b = COS2[(i - k) % 10]; p += c[i] * a; q += c[i] * b
    return (p, q)
def lt(u, v): return B.phi_sign((u[0] - v[0], u[1] - v[1])) < 0
def is_wall(c):
    x = proj2(c, 0)
    w1 = lt((-44, 0), x) and lt(x, (-38, 0)) and lt(proj2(c, 3), (60, 0))     # 左の壁：下から立ち、上が開く
    w2 = lt((38, 0), x) and lt(x, (44, 0)) and lt((-60, 0), proj2(c, 2))     # 右の壁：上から下り、下が開く
    return w1 or w2

def layer_count(src):
    """別の計算：層（隣の段数）を幅優先で求め、一つ内側の層の本数を足して最短の道の本数を出す。"""
    d = {src: 0}; q = deque([src]); order = []
    while q:
        u = q.popleft(); order.append(u)
        for v in nb[u]:
            if v not in d: d[v] = d[u] + 1; q.append(v)
    n = {src: 1}
    for u in order[1:]:
        n[u] = sum(n[v] for v in nb[u] if d[v] == d[u] - 1)
    return d, n

if __name__ == "__main__":
    print("0. 地図：五角形", len(cells), "枚。隣の距離²", B.norm2(STEP[0]))

    print("検定1 隣")
    deg = Counter(len(v) for v in nb.values())
    src0 = min(cells)
    d0, _ = layer_count(src0)
    same = sum(1 for u in cells for v in nb[u] if d0[u] % 2 == d0[v] % 2)
    print("   隣の枚数の内訳:", dict(sorted(deg.items())), "  つながり:", len(d0), "/", len(cells),
          "  同じ色どうしの隣:", same)
    print("   →", "OK" if set(deg) <= {2, 3, 4} and len(d0) == len(cells) and same == 0 else "NG")

    random.seed(10)
    pool = sorted(cells)
    pairs = [tuple(random.sample(pool, 2)) for _ in range(200)]

    print("検定2 匂いの量＝道の数")
    ok2, neg2, color2 = True, 0, True
    for A, F in pairs[:20]:
        tF, aF, hist = spread(F)
        dF, nF = layer_count(F)
        ok2 &= all(tF[v] == dF[v] and aF[v] == nF[v] for v in cells)
        color2 &= all(all(t % 2 == dF[v] % 2 for v in hist[t]) for t in range(len(hist)))
        neg2 += sum(1 for v in cells if tF[v] + 2 < len(hist) and hist[tF[v] + 2].get(v, 0) != nF[v])
    print("   餌 20 か所 × 五角形 6030 枚で、初めて届いた刻と量が層の計算と一致:", ok2)
    print("   どの刻も、匂いは片方の色の五角形にだけある:", color2)
    print("   負の対照：2刻あとの量が道の本数と違う五角形の数:", neg2)
    print("   →", "OK" if ok2 and color2 and neg2 > 0 else "NG")

    print("検定3 二つの波が出会う")
    ok3, neg3, ds = True, 0, []
    for A, F in pairs:
        dA, nA = layer_count(A); dF, nF = layer_count(F)
        d = dA[F]; total = nA[F]; ds.append(d)
        assert nF[A] == total
        band = Counter(); off = Counter()
        for v in cells:
            if dA[v] + dF[v] == d: band[dA[v]] += nA[v] * nF[v]
            elif dA[v] + dF[v] == d + 2: off[dA[v]] += nA[v] * nF[v]
        ok3 &= all(band[k] == total for k in range(d + 1))
        neg3 += sum(1 for k in off if off[k] != total)
    print("   200 組の最短距離 d：", min(ds), "〜", max(ds), "  帯のどの段でも 体×匂い の和＝道の総数:", ok3)
    print("   負の対照：和が d+2 の五角形で、段の和が道の総数と違う段の数:", neg3)
    print("   →", "OK" if ok3 and neg3 > 0 else "NG")

    print("検定4 しっぽが縮む")
    def shrink(A, F, dA, nA, dF, nF):
        p, path = A, [A]
        while p != F and len(path) <= dA[F] + 1:
            band_nb = [v for v in nb[p] if dA[v] == dA[p] + 1 and dA[v] + dF[v] == dA[F]]
            p = max(band_nb, key=lambda v: nA[v] * nF[v]); path.append(p)
        return path
    def greedy(A, F, hist, dF):
        p, steps = A, 0; t0 = dF[A]
        while p != F and steps < 3 * dF[A]:
            Sx = hist[min(t0 + steps, len(hist) - 1)]
            p = max(nb[p], key=lambda v: Sx.get(v, 0)); steps += 1
        return p == F and steps == dF[A]
    ok4 = True
    for A, F in pairs:
        dA, nA = layer_count(A); dF, nF = layer_count(F)
        path = shrink(A, F, dA, nA, dF, nF)
        ok4 &= path[-1] == F and len(path) - 1 == dA[F]
    neg_ok = 0
    for A, F in pairs[:20]:
        tF, aF, hist = spread(F)
        cum = []; acc = Counter()
        for W in hist:
            for v, w in W.items(): acc[v] += w
            cum.append(dict(acc))
        neg_ok += greedy(A, F, cum, tF)
    print("   200 組すべてで d 歩ちょうどで餌に着く:", ok4)
    A, F = (-43, -2, -28, 24), (43, 3, 28, -24)            # 図の組
    dA, nA = layer_count(A); dF, nF = layer_count(F)
    path = shrink(A, F, dA, nA, dF, nF)
    print("   図の組：アメーバ", A, "→ 餌", F, " 最短距離", dA[F], " 最短の道の本数", nA[F],
          " 帯の五角形", sum(1 for v in cells if dA[v] + dF[v] == dA[F]), "枚  縮んだ道", len(path) - 1, "歩")
    print("   負の対照：足し続けた匂いの濃い隣へ進むと d 歩で着く組:", neg_ok, "/ 20")
    print("   →", "OK" if ok4 and neg_ok < 20 else "NG")

    print("検定5 障害物")
    A, F = (-43, -2, -28, 24), (43, 3, 28, -24)
    dA, nA = layer_count(A); dF, nF = layer_count(F)
    free_path = shrink(A, F, dA, nA, dF, nF)
    d_free = dA[F]
    WALL = {c for c in cells if is_wall(c)}
    cells = cells - WALL                                       # 壁の五角形を地図から外す
    nb = {c: [B.zadd(c, s) for s in STEP if B.zadd(c, s) in cells] for c in cells}
    hits = sum(1 for c in free_path if c in WALL)
    d0, _ = layer_count(A)
    print("   壁の五角形", len(WALL), "枚  残りの地図", len(cells), "枚  つながり:", len(d0), "/", len(cells))
    print("   壁のない地図で残った道が通る壁の五角形:", hits, "枚")
    pool = sorted(cells); random.seed(10)
    wpairs = [tuple(random.sample(pool, 2)) for _ in range(200)]
    ok5 = len(d0) == len(cells) and hits > 0
    for A2, F2 in wpairs[:20]:
        tF, aF, hist = spread(F2); dF2, nF2 = layer_count(F2)
        ok5 &= all(tF[v] == dF2[v] and aF[v] == nF2[v] for v in cells)
    for A2, F2 in wpairs:
        dA2, nA2 = layer_count(A2); dF2, nF2 = layer_count(F2); d = dA2[F2]
        band = Counter()
        for v in cells:
            if dA2[v] + dF2[v] == d: band[dA2[v]] += nA2[v] * nF2[v]
        ok5 &= all(band[k] == nA2[F2] for k in range(d + 1))
        p = shrink(A2, F2, dA2, nA2, dF2, nF2)
        ok5 &= p[-1] == F2 and len(p) - 1 == d
    dA, nA = layer_count(A); dF, nF = layer_count(F)
    wall_path = shrink(A, F, dA, nA, dF, nF)
    ok5 &= dA[F] > d_free and not any(c in WALL for c in wall_path)
    print("   壁のある地図の 200 組で、検定2（20 か所）・検定3・検定4 がすべて成り立つ:", ok5)
    print("   図の組：最短距離", d_free, "→", dA[F], "  最短の道の本数", nA[F],
          "  帯の五角形", sum(1 for v in cells if dA[v] + dF[v] == dA[F]), "枚  残った道", len(wall_path) - 1,
          "歩  壁を通る数", sum(1 for c in wall_path if c in WALL))
    print("   →", "OK" if ok5 else "NG")
