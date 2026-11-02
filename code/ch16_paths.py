"""第16章　すべての経路 ── 検算スクリプト

使い方:  python3 ch16_paths.py <b13_chain_units.py のあるフォルダ>
  地図・規則・出発は第15章（ch15_wave.py）と同じ。一本の経路の重みは、通った五角形ごとの係数の積：
  隣 2 枚：抜ける 12・戻る 0／隣 3 枚：曲がる 8・戻る −4／隣 4 枚：曲がる 6・戻る −6

経路の数え始め
  五芒星から出発の 5 枚に入る形にする：出発の五角形（隣 2 枚）から両隣へ、1 歩目の重みは 12（戻りはない）。
  波の値は 1 歩目以降、第15章と同じ（第15章の出発は、この 1 歩目と同じ値を作る）。

事前登録
  G1 経路の和：1〜10 歩のすべての辺で、波の値 ＝ そこに着く全経路の重みの和（整数で一致）。
     負の対照：2〜10 歩の各歩で、重みの符号を外す（絶対値の積）と一致しない
  G2 打ち消し：経路の本数・重み 0・負の経路の数、Σ|重み| と Σ|和| を歩ごとに記録（判定しない）
  G3 前線は最短路だけ：出発からちょうど t 歩の五角形に着く経路は、すべて戻りがなく、重みはすべて正
  G4 景色（記録）：外周から十分離れた五芒星を出発点にして、t 歩目の（中心からの距離², 経路の和）の並びを景色とし、
     地図の対称（回転と鏡）で移り合わない組の数と、景色の種類の数を 1〜12 歩で記録する（判定しない）。
     外周から離れているかの選別だけ実数の索引（第9章と同じ）
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch15_wave as W
import ch05_views as V
B = W.B

def wgt(w, v, x):
    k = 24 // len(W.nb[w])
    return k - 12 if x == v else k

def enumerate_paths(T):
    """t 歩の全経路を、到着の辺ごとに（和, 絶対値の和, 本数, 0 の本数, 負の本数）で集める。"""
    out = []
    layer = []; agg = {}; stats = Counter(); front_bad = 0              # 1 歩目：五芒星から入って両隣へ 12
    for w in W.START:
        for x in W.nb[w]:
            st = (x, w); layer.append((st, 12, 12, False))
            s = agg.setdefault(st, [0, 0]); s[0] += 12; s[1] += 12; stats["本数"] += 1
            if W.dist[x] != 1: front_bad += 1
    out.append((agg, stats, front_bad))
    for t in range(2, T + 1):
        nxt = []; agg = {}; stats = Counter(); front_bad = 0
        for (w, v), a, b, back in layer:
            for x in W.nb[w]:
                f = wgt(w, v, x); st = (x, w); na = a * f; nb_ = b * abs(f); bk = back or x == v
                nxt.append((st, na, nb_, bk))
                s = agg.setdefault(st, [0, 0]); s[0] += na; s[1] += nb_
                stats["本数"] += 1; stats["0"] += na == 0; stats["負"] += na < 0
                if W.dist[x] == t and (bk or na <= 0): front_bad += 1
        layer = nxt; out.append((agg, stats, front_bad))
    return out

if __name__ == "__main__":
    T = 10
    hist = [W.start_state()]
    for t in range(T): hist.append(W.step(hist[-1]))
    res = enumerate_paths(T)
    print("G1 経路の和（1〜%d 歩）" % T)
    ok1 = True; okn = True
    for t in range(1, T + 1):
        agg = res[t - 1][0]
        wave = {k: v for k, v in hist[t].items()}
        s = {k: v[0] for k, v in agg.items() if v[0]}
        a = {k: v[1] for k, v in agg.items() if v[1]}
        ok1 &= s == wave
        if t >= 2: okn &= a != wave
    print("   すべての歩・すべての辺で、波の値 ＝ 全経路の重みの和:", ok1, "  負の対照（絶対値の積）が 2〜%d 歩のどの歩でも一致しない:" % T, okn)
    print("   →", "OK" if ok1 and okn else "NG")
    print("G2 打ち消し")
    for t in range(1, T + 1):
        agg, st, _ = res[t - 1]
        sa = sum(v[1] for v in agg.values()); ss = sum(abs(v[0]) for v in agg.values())
        print("   %2d 歩：経路 %8d 本  重み 0 %7d  負 %7d   Σ|重み| / Σ|和| = %d / %d" % (t, st["本数"], st["0"], st["負"], sa, ss))
    print("G3 前線は最短路だけ")
    bad = sum(r[2] for r in res)
    fr = sum(1 for t in range(1, T + 1) for (x, w) in res[t - 1][0] if W.dist[x] == t)
    print("   前線の辺（延べ）:", fr, "  そこへ着く経路のうち、戻りを含むか重みが正でないもの:", bad)
    print("   →", "OK" if bad == 0 and fr > 0 else "NG")

    print("G4 景色")
    cells = W.cells
    G = []
    for k in range(10):
        for m in (False, True):
            f = (lambda z, k=k, m=m: B.zrot(B.zconj(z) if m else z, k))
            if all(f(c) in cells for c in cells): G.append(f)
    stars = [g for g in V.stars]
    print("   地図の対称（回転と鏡）:", len(G), "通り   五芒星:", len(stars), "個")
    def view(g, t):
        st = {(w, x): 1 for w in cells if B.norm2(B.zsub(w, g)) == (1, 1) for x in W.nb[w]}
        for _ in range(t): st = W.step(st)
        s = Counter()
        for (w, x), a in st.items(): s[w] += a
        return tuple(sorted((B.norm2(B.zsub(w, g)), v) for w, v in s.items() if v))
    first = None
    for t in range(1, 13):
        inner = [g for g in stars if V.to_outer(g) >= 1.7 * (t + 3)]
        orb = {g: min(f(g) for f in G) for g in inner}
        vw = {g: view(g, t) for g in inner}
        byv = Counter(vw.values()); byo = len(set(orb.values()))
        same_in_orbit = all(vw[g] == vw[h] for g in inner for h in inner if orb[g] == orb[h])
        ok = len(byv) == byo and same_in_orbit
        print("   t = %2d：五芒星 %3d 個  対称で移り合わない組 %3d  景色の種類 %3d  対称で移り合う五芒星は同じ景色 %s" %
              (t, len(inner), byo, len(byv), same_in_orbit))
