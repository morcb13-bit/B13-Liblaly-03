"""第9章　五芒星をニューロンに ── 検算スクリプト

使い方:  python3 ch09_neuron.py <b13_chain_units.py のあるフォルダ>
  910環の図の五芒星をニューロン、距離 φ⁴ でとなり合う五芒星どうしをシナプスとみなす。
  五芒星の5つの先端を平衡5進数の桁 −2, −1, 0, +1, +2 に対応させ、信号の行き先と記憶をこの桁で書く。
  判定はすべて整数（Z[ζ₁₀] の4成分、距離² は p + qφ）。

事前登録
  検定1 シナプスは先端の向きだけ：五芒星（向き k）から距離 φ⁴ にある五芒星は、必ず k+1+2i（i = 0..4）の向き、
        すなわち五芒星の先端の向きにある。五角形の向き（k+2i）には一つもない。外れれば NG
  検定2 シナプスの数：真ん中の五芒星は 5 本、外周から離れた五芒星はすべて 3 本。ほかが出れば NG
  検定3 隣だけで伝える：真ん中から信号を隣へ隣へと渡すと、416 個すべてに届く。n 歩目に初めて届く数を数える
  検定4 記憶：信号が最初に届いたとき、来た先端の桁（平衡5進数）を五芒星に書き残す。
        どの五芒星からでも、書き残した桁をたどると、届いた歩数ちょうどで真ん中に戻る。戻らなければ NG
  検定5 番地をヘッダに：真ん中から目的の五芒星までの道を桁の列（ヘッダ）で書き、真ん中から桁どおりに進むと
        目的の五芒星に着く（416 個すべて）。桁を一つ変えると別の五芒星に着く（負の対照）
"""
import sys, math
from collections import Counter, deque
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = B.fits(U)
stars = {}
for p in cells:
    for j in range(10):
        g = B.zsub(p, B.zrot(B.PHI, j))
        if g in cells or g in stars: continue
        if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)): stars[g] = j % 2
SS = set(stars)
PHI4 = B.zmul(B.PHI2, B.PHI2)                       # φ⁴ = (5, 0, 3, −3)

def tip_digit(k, j):
    """向き j（36° 単位）の先端の、平衡5進数の桁。先端は k+1+2i、桁は i − 2。"""
    i = ((j - k - 1) % 10) // 2
    return i - 2
def tip_dir(k, d):
    return (k + 1 + 2 * (d + 2)) % 10

def synapses(g):
    return [j for j in range(10) if B.zadd(g, B.zrot(PHI4, j)) in SS]

def bal5(n, w):
    out = []
    for _ in range(w):
        r = n % 5
        if r > 2: r -= 5
        out.append(r); n = (n - r) // 5
    return out

if __name__ == "__main__":
    O = (0, 0, 0, 0)
    print("0. ニューロン（五芒星）", len(stars), "個。シナプスの長さ² φ⁸ =", B.norm2(PHI4))

    print("検定1 シナプスの向き")
    cnt = Counter()
    for g, k in stars.items():
        for j in synapses(g): cnt["先端" if (j - k) % 2 == 1 else "五角形の向き"] += 1
    print("   向きの内訳（両端から数える）:", dict(cnt), "→", "OK" if set(cnt) == {"先端"} else "NG")

    print("検定2 シナプスの数")
    # 外周の索引（観察者側）：第4章と同じく、図の外周の輪からの距離で選ぶ
    import ch04_views as V4
    def far(g, R): return V4.to_outer(g) >= R
    deg = Counter((g == O, len(synapses(g))) for g in stars if far(g, 8))
    print("   外周から 8 以上離れた五芒星（真ん中か, シナプスの数）:", dict(deg))
    print("   すべての五芒星:", dict(Counter(len(synapses(g)) for g in stars)))
    print("   →", "OK" if deg == Counter({(True, 5): 1, (False, 3): sum(v for (c, n), v in deg.items() if not c)}) and
          all(n == 3 for (c, n) in deg if not c) else "NG")

    print("検定3 隣だけで伝える（真ん中から）")
    hop, back = {O: 0}, {}
    q = deque([O])
    while q:
        u = q.popleft()
        for j in synapses(u):
            v = B.zadd(u, B.zrot(PHI4, j))
            if v not in hop:
                hop[v] = hop[u] + 1
                back[v] = tip_digit(stars[v], (j + 5) % 10)     # v から見て、来た向き（逆向き）の先端の桁
                q.append(v)
    layers = Counter(hop.values())
    print("   届いた数:", len(hop), "/", len(stars), "  n 歩目に初めて届く数:", [layers[i] for i in range(max(layers) + 1)])
    print("   →", "OK" if len(hop) == len(stars) else "NG")

    print("検定4 記憶：来た先端の桁を書き残し、たどって戻る")
    ok4 = True
    for g in stars:
        if g == O: continue
        cur, steps = g, 0
        while cur != O and steps < 100:
            j = tip_dir(stars[cur], back[cur])
            cur = B.zadd(cur, B.zrot(PHI4, j)); steps += 1
        ok4 &= cur == O and steps == hop[g]
    print("   書き残した桁の内訳:", dict(Counter(back.values())))
    print("   どの五芒星からも、届いた歩数で真ん中に戻る:", ok4, "→", "OK" if ok4 else "NG")

    print("検定5 番地をヘッダに")
    def header(g):
        path, cur = [], g
        while cur != O:
            j = tip_dir(stars[cur], back[cur]); nxt = B.zadd(cur, B.zrot(PHI4, j))
            path.append(tip_digit(stars[nxt], (j + 5) % 10)); cur = nxt
        return path[::-1]
    def follow(h):
        cur = O
        for d in h:
            j = tip_dir(stars[cur], d); nxt = B.zadd(cur, B.zrot(PHI4, j))
            if nxt not in SS: return None
            cur = nxt
        return cur
    ok5 = all(follow(header(g)) == g for g in stars)
    far_g = max(stars, key=lambda g: hop[g])
    h = header(far_g)
    print("   例：", far_g, "（", hop[far_g], "歩）のヘッダ", h)
    alt = list(h); alt[-1] = next(d for d in (-2, -1, 0, 1, 2) if d != h[-1] and follow(h[:-1] + [d]) is not None)
    print("   最後の桁を", h[-1], "→", alt[-1], "に変えると", follow(alt), "に着く（別の五芒星）")
    print("   416 個すべてでヘッダどおりに着く:", ok5, "→", "OK" if ok5 and follow(alt) != far_g else "NG")
    print("   平衡5進数の例：", 13, "=", bal5(13, 3), "（下の桁から）、", 2026, "=", bal5(2026, 6))
