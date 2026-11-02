"""終章　円環に乗るもう一つのペンローズ ── 図のデータと検算

使い方:  python3 ch19_ring_rhombus.py <b13_chain_units.py のあるフォルダ>
  黄金のひし形のペンローズ（de Bruijn の指数と窓。窓を整数の組 g / 97 だけずらして、真ん中の決まらない所をなくす）の頂点を、
  (1 − ζ²) 倍（辺の長さが五角形の辺と同じ）して、円環の中心 c0 = (−2, 2, 0, 3) に置く。判定はすべて整数。
  記録：円環の五角形 10 枚の頂点 30 個のうち、ひし形の頂点に乗る数。円環の内側の十角形をふさぐひし形の数。
  あわせて、図のためにひし形の頂点（番地）と残渣を ring_rhombus.json に書き出す。
"""
import sys, math, json
sys.argv = sys.argv[:2]
import ch19_cut as C, ch19_rhombus as Rh
from collections import Counter
V, B = C.V, C.B
N, G = 97, (1, 0, 0, 0)
def judge_g(a):
    r = Rh.idx(a)
    if r not in Rh.WIN: return False
    p = B.zsub(B.zmul((N, 0, 0, 0), Rh.res(a)), G)
    return Rh.inside([B.zmul((N, 0, 0, 0), v) for v in Rh.WIN[r]], p)

if __name__ == "__main__":
    cand = C.candidates(9.0, 2.0)
    DB = {a for a in cand if judge_g(a)}
    c0 = (-2, 2, 0, 3); m = B.zsub((1, 0, 0, 0), B.zt(2))
    PHI2 = B.zmul(B.PHI, B.PHI)
    ring = [B.zadd(c0, B.zmul(PHI2, B.zt(k))) for k in range(10)]
    pents = [[B.zadd(q, B.zt(V.cells[q] + 2 * i)) for i in range(5)] for q in ring]
    PV = {v for p in pents for v in p}
    T = {B.zadd(c0, B.zmul(m, a)): a for a in DB}
    E = [B.zt(k) for k in range(10)]
    S = set(DB); rh = {}
    for a in DB:
        for k1 in range(10):
            for k2 in range(k1 + 1, 10):
                d = (k2 - k1) % 10
                if d in (0, 5): continue
                p1, p2, p3 = B.zadd(a, E[k1]), B.zadd(B.zadd(a, E[k1]), E[k2]), B.zadd(a, E[k2])
                if p1 in S and p2 in S and p3 in S:
                    rh[frozenset([a, p1, p2, p3])] = "thin" if d in (1, 4, 6, 9) else "fat"
    inner = {v for v in PV if B.norm2(B.zsub(v, c0)) == (2, 1)}
    outer = PV - inner
    print("ひし形の頂点（番地の候補から窓で選んだもの）:", len(DB), "  ひし形:", len(rh), dict(Counter(rh.values())))
    print("円環の五角形の頂点", len(PV), "個のうち、ひし形の頂点に乗るもの:", len(PV & set(T)),
          "（外側の 20 個:", len(outer & set(T)), "／内側の十角形の 10 個:", len(inner & set(T)), "）")
    # 図のデータ（表示の選別だけ実数）
    XY = C.XY
    def inpoly(pt, poly):
        s = None
        for i in range(len(poly)):
            x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
            c = (x2 - x1) * (pt[1] - y1) - (y2 - y1) * (pt[0] - x1)
            if abs(c) < 1e-9: continue
            if s is None: s = c > 0
            elif s != (c > 0): return False
        return True
    dec = sorted(inner, key=lambda v: math.atan2(XY(B.zsub(v, c0))[1], XY(B.zsub(v, c0))[0]))
    pp = [[XY(v) for v in p] for p in pents]; dd = [XY(v) for v in dec]
    out = []; indec = Counter()
    for verts, kind in rh.items():
        vs = [B.zadd(c0, B.zmul(m, a)) for a in verts]
        cx = sum(XY(v)[0] for v in vs) / 4; cy = sum(XY(v)[1] for v in vs) / 4
        ind = inpoly((cx, cy), dd)
        if ind or any(inpoly((cx, cy), q) for q in pp):
            if ind: indec[kind] += 1
            src = {B.zadd(c0, B.zmul(m, a)): a for a in verts}
            o = sorted(vs, key=lambda v: math.atan2(XY(v)[1] - cy, XY(v)[0] - cx))
            out.append({"k": kind, "v": [list(v) for v in o], "r": [list(Rh.res(src[v])) for v in o]})
    print("円環と内側の十角形の上に乗るひし形:", len(out), "  そのうち十角形の中:", dict(indec))
    json.dump({"c0": list(c0), "pents": [[list(v) for v in p] for p in pents], "dec": [list(v) for v in dec], "rh": out},
              open("ring_rhombus.json", "w"))
