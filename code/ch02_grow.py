"""第2章　円環をつなぐ ── 検算スクリプト（ch00_ring.py・ch01_join.py と同じフォルダで走る）

使い方:  python3 ch02_grow.py rings_integer.json
  rings_integer.json の先頭25個（手で組んだ三層）から始め、
  (1) 五芒星のまわりに円環5個を置く規則を、止まるまでくり返す
  (2) 円環の中心が作る大きなひし形の、欠けた頂点を埋める
の二段で円環を足し、記録済みの55環・60環と一致するかを確かめる。判定はすべて整数。
"""
import sys, json
from ch00_ring import add, sub, turn, dist2, ZERO, ONE, PHI, PHI2, ring
from ch01_join import overlap, phi_sign

PHI3 = add(PHI2, PHI)                                   # φ³ = φ² + φ
# 第1章の判定で「つながる」となった20通りを、中心の距離²で2つに分ける
from ch01_join import judge
_ok = {sub(turn(PHI2, j), turn(PHI2, k)) for j in range(10) for k in range(10)} - {ZERO}
_ok = [c for c in _ok if judge(c)[0] == "ok"]
CONT = {c for c in _ok if dist2(c, ZERO) == (7, 11)}    # 連続の10方向
SKIP = {c for c in _ok if dist2(c, ZERO) == (3, 4)}     # 2つ飛ばしの10方向

def lt(x, y): return phi_sign(x[0] - y[0], x[1] - y[1]) < 0

class Plane:
    def __init__(self, centers):
        self.rings, self.cells = [], {}
        for c in centers:
            assert self.fits(c), c
            self.put(c)
    def fits(self, c):
        for p, a in ring(c).items():
            if p in self.cells:
                if (self.cells[p] - a) % 2: return False            # 向きの食い違い
                continue
            for q, b in self.cells.items():
                if lt(dist2(p, q), (4, 0)) and overlap(q, b, p, a): return False   # 重なり
        return True
    def put(self, c):
        self.rings.append(c)
        for p, a in ring(c).items(): self.cells.setdefault(p, a)
    def forming(self, need=3):
        """五芒星になりかけの点：その点に五角形がなく、距離 φ の5方向のうち need 枚以上がそろっている点。"""
        out = {}
        for p in self.cells:
            for j in range(10):
                g = sub(p, turn(PHI, j))
                if g in self.cells or g in out: continue
                k = j % 2
                if sum(add(g, turn(PHI, k + 2 * i)) in self.cells for i in range(5)) >= need:
                    out[g] = k
        return out
    def pentagrams(self):
        """五芒星：その点に五角形がなく、距離 φ の5方向（72° おき）に五角形がある点。向き k（0 か 1）も返す。"""
        out = {}
        for p in self.cells:
            for j in range(10):
                g = sub(p, turn(PHI, j))
                if g in self.cells or g in out: continue
                k = j % 2
                if all(add(g, turn(PHI, k + 2 * i)) in self.cells for i in range(5)):
                    out[g] = k
        return out
    def grow(self):
        """五芒星 g のまわり g + φ³ζ^(k+2i) に円環5個。足りない五芒星の分を、5個そろって入るときだけ足す。"""
        added = 0
        while True:
            have = set(self.rings); step = 0
            for g, k in self.pentagrams().items():
                need = [add(g, turn(PHI3, k + 2 * i)) for i in range(5)]
                miss = [c for c in need if c not in have]
                if miss and all(self.fits(c) for c in miss):
                    trial = Plane.__new__(Plane); trial.rings = list(self.rings); trial.cells = dict(self.cells)
                    ok = True
                    for c in miss:
                        if not trial.fits(c): ok = False; break
                        trial.put(c)
                    if ok:
                        for c in miss: self.put(c); have.add(c)
                        step += len(miss)
            added += step
            if step == 0: return added
    def pairs(self, vecs):
        S = set(self.rings)
        return [(a, b) for a in self.rings for b in self.rings if a < b and sub(b, a) in vecs or (a < b and sub(a, b) in vecs)]
    def big_rhombi(self):
        """2つ飛ばしの組 (a, b) を短い対角線とし、a と b の両方に連続でつながる2点を残りの頂点とするひし形。"""
        S = set(self.rings); out = []
        for a, b in self.pairs(SKIP):
            tips = [add(a, e) for e in CONT if sub(add(a, e), b) in CONT]
            assert len(tips) == 2
            out.append((a, b, tips, [t for t in tips if t not in S]))
        return out
    def shells(self):
        d = {}
        for c in self.rings: d[dist2(c, ZERO)] = d.get(dist2(c, ZERO), 0) + 1
        return sorted(d.items(), key=lambda kv: (kv[0][0] + kv[0][1] * 1.618034))   # 並べる順だけ（観察者側）

def report(P, name):
    ps = P.pentagrams()
    sk = P.pairs(SKIP)
    eq = sum(1 for a, b in sk if dist2(a, ZERO) == dist2(b, ZERO))
    print(f"   {name}: 円環 {len(P.rings)}、五角形 {len(P.cells)}、五芒星 {len(ps)}、"
          f"連続 {len(P.pairs(CONT))} 組、2つ飛ばし {len(sk)} 組（中心からの距離が等しい組 {eq}）")

if __name__ == "__main__":
    rec = [tuple(r) for r in json.load(open(sys.argv[1]))["rings"]]
    NEED = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    print("0. つなぎ方の10方向")
    print("   連続:", {dist2(v, ZERO) for v in CONT}, len(CONT), "方向   2つ飛ばし:", {dist2(v, ZERO) for v in SKIP}, len(SKIP), "方向")
    print("   φ³ =", PHI3)

    print("1. 手で組んだ三層（25環）から始める")
    P = Plane(rec[:25]); report(P, "はじめ")
    for g, k in sorted(P.pentagrams().items()):
        pass
    print("2. 五芒星のまわりに円環5個を置く規則をくり返す")
    rounds = []
    while True:
        before = len(P.rings)
        have = set(P.rings); step = []
        for g, k in P.forming(NEED).items():
            need = [add(g, turn(PHI3, k + 2 * i)) for i in range(5)]
            miss = [c for c in need if c not in have]
            if not miss: continue
            trial = Plane.__new__(Plane); trial.rings = list(P.rings); trial.cells = dict(P.cells); ok = True
            for c in miss:
                if not trial.fits(c): ok = False; break
                trial.put(c)
            if ok:
                for c in miss: P.put(c); have.add(c); step.append(c)
        if not step: break
        report(P, f"{len(rounds) + 1} 回目")
        rounds.append(len(step))
    print("   規則が求める新しい円環: 0（止まった）")
    print("   記録済みの55環と一致:", set(P.rings) == set(rec[:55]))
    miss_all = [g for g, k in P.pentagrams().items()
                if any(add(g, turn(PHI3, k + 2 * i)) not in set(P.rings) for i in range(5))]
    print("   円環5個に囲まれていない五芒星:", len(miss_all))

    print("3. 円環の中心が作る大きなひし形")
    R = P.big_rhombi()
    print("   ひし形", len(R), "個。頂点がそろっている", sum(1 for r in R if not r[3]),
          "個、一つ欠けている", sum(1 for r in R if len(r[3]) == 1), "個")
    new = sorted({t for r in R for t in r[3]})
    print("   欠けた頂点:", new)
    print("   欠けた頂点の中心からの距離²:", {dist2(t, ZERO) for t in new})
    for t in new:
        assert P.fits(t); P.put(t)
    report(P, "埋めたあと")
    R = P.big_rhombi()
    print("   ひし形", len(R), "個。頂点がそろっている", sum(1 for r in R if not r[3]), "個")
    print("   記録済みの60環と一致:", set(P.rings) == set(rec))
    print("   その後、規則が求める新しい円環:",
          sum(1 for g, k in P.pentagrams().items() for i in range(5)
              if add(g, turn(PHI3, k + 2 * i)) not in set(P.rings)),
          "  欠けた頂点:", sum(len(r[3]) for r in R))

    print("4. 円環の中心の殻（中心からの距離² ＝ p + qφ と個数）")
    F = [0, 1]
    for _ in range(30): F.append(F[-1] + F[-2])
    for (p, q), n in P.shells():
        tag = ""
        for m in range(1, 6):
            for e in range(1, 20):
                if (p, q) == (m * m * F[e - 1], m * m * F[e]):   # m²·φ^e ＝ m²(F(e−1) + F(e)φ)
                    tag = f"＝ ({m if m > 1 else ''}φ^{e // 2})²" if e % 2 == 0 else f"＝ {m*m}φ^{e}"
        print(f"   {p} + {q}φ : {n} 個 {tag}")
