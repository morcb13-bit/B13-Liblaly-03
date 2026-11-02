"""第5章　10方向に組む（景色と残渣）── 検算スクリプト

使い方:  python3 ch05_views.py <b13_chain_units.py のあるフォルダ> rings_integer.json
  910環の図（b13_chain_units.py の積み上げを 10 方向に組んだもの）の上で、
  五芒星ごとに「そこから見える円環の並び」を数える。

判定は整数（Z[ζ₁₀] の4成分、距離² は p + qφ、比較は phi_lt）。
実数を使うのは、図の外周から離れた五芒星を選ぶための索引（観察者側）だけ。

事前登録
  検定1 五芒星のまわり：外周から離れた五芒星は、距離 φ³ に円環をちょうど5個持つ。5以外が出れば NG
  検定2 72° 回して自分に重なる景色：見る半径 22 で、図の中心の五芒星1個だけ。ほかにも出れば NG
  検定3 景色の類：見る半径を広げると、類が増える。広げても類が増えなければ「景色はどれも同じ」
  検定4 一歩の大きさ：五芒星→円環の一歩のノルムは 1、円環どうしのつなぎ方のノルムは 5
  負の対照：見る半径を φ³ の殻のすぐ外（距離² < 20）に絞ると、景色は 1 類に潰れる
"""
import sys, json, math
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B
from collections import Counter, defaultdict

# ── 910環を組む（star_views.py と同じ手順） ──
T = (-2, 2, 0, 3)     # 距離 φ³
S = (-5, 5, 0, 8)     # 距離 φ⁵
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = B.fits(U)
PHI3 = B.zmul(B.PHI2, B.PHI)

# ── 五芒星（整数）：点に五角形がなく、距離 φ の5方向（72° おき）に五角形がある ──
stars = {}
for p in cells:
    for j in range(10):
        g = B.zsub(p, B.zrot(B.PHI, j))
        if g in cells or g in stars: continue
        if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)):
            stars[g] = j % 2

# ── 外周（観察者側の索引）：一度しか使われない五角形の辺をたどった外側の輪 ──
poly = {q: [B.zadd(q, B.zt(a + 2 * i)) for i in range(5)] for q, a in cells.items()}
use = Counter(); nbr = defaultdict(set)
for q, vs in poly.items():
    for i in range(5): use[frozenset((vs[i], vs[(i + 1) % 5]))] += 1
for e, n in use.items():
    if n == 1:
        u, v = tuple(e); nbr[u].add(v); nbr[v].add(u)
far = max(nbr, key=lambda u: math.hypot(*B.xy(u)))
OB, prev, cur = [], None, far
while True:
    OB.append(B.xy(cur))
    nx = [w for w in nbr[cur] if w != prev]
    if not nx: break
    prev, cur = cur, nx[0]
    if cur == far: break
def to_outer(g):
    p = B.xy(g); m = 1e9
    for i in range(len(OB)):
        a, b = OB[i], OB[(i + 1) % len(OB)]
        dx, dy = b[0] - a[0], b[1] - a[1]
        t = max(0, min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
        m = min(m, math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy))
    return m

# ── 景色（整数）：五芒星から距離² < R2 の円環中心の差を、10回転×鏡映で一番小さい並びにそろえる ──
def view(g, R2):
    rel = [B.zsub(c, g) for c in U if B.phi_lt(B.norm2(B.zsub(c, g)), (R2, 0))]
    best = None
    for k in range(10):
        for m in (False, True):
            sig = tuple(sorted(B.zrot(B.zconj(v) if m else v, k) for v in rel))
            if best is None or sig < best: best = sig
    return best

def sym72(g, R2):
    rel = {B.zsub(c, g) for c in U if B.phi_lt(B.norm2(B.zsub(c, g)), (R2, 0))}
    return {B.zrot(v, 2) for v in rel} == rel

def norm(pq):
    p, q = pq
    return p * p + p * q - q * q

if __name__ == "__main__":
    rec = {tuple(r) for r in json.load(open(sys.argv[2]))["rings"]}
    print("0. 図：円環", len(U), "個、五角形", len(cells), "枚、五芒星", len(stars), "個")
    print("   規則で育てた60環を 36° 回すと、この図にそっくり含まれる:", {B.zrot(r, 1) for r in rec} <= set(U),
          "  回さないと:", len(rec & set(U)), "個")
    print("   積み上げた三角形の段（行ごとの円環の数）:", [len(r) for r in rows])
    d = {g: to_outer(g) for g in stars}

    print("検定1 五芒星から距離 φ³ にある円環の数（外周から 6 以上離れた五芒星）")
    Us = set(U)
    cnt = Counter(sum(B.zadd(g, B.zrot(PHI3, k + 2 * i)) in Us for i in range(5)) for g, k in stars.items() if d[g] >= 6)
    allphi3 = Counter(sum(1 for c in U if B.norm2(B.zsub(c, g)) == B.norm2(PHI3)) for g in stars if d[g] >= 6)
    print("   5方向の席に円環がある数:", dict(cnt), "  距離 φ³ の円環の総数:", dict(allphi3),
          "→", "OK" if set(cnt) == {5} and set(allphi3) == {5} else "NG")

    print("検定2 72° 回して自分に重なる景色（見る半径 22、外周から 22 以上の五芒星）")
    pool = [g for g in stars if d[g] >= 22]
    sym = [g for g in pool if sym72(g, 484)]
    print("   対象", len(pool), "個のうち", len(sym), "個:", sym, "→", "OK" if sym == [B.ZERO] else "NG")

    print("検定3 景色の類（10回転と鏡映でそろえる）")
    for R, R2 in ((4.47, 20), (10, 100), (15, 225), (22, 484)):
        pool = [g for g in stars if d[g] >= R]
        cls = Counter(view(g, R2) for g in pool)
        sizes = sorted(cls.values(), reverse=True)
        print(f"   見る半径 √{R2}: 対象 {len(pool)} 個 → {len(cls)} 類（一番大きい類 {sizes[0]} 個、1個だけの類 {sizes.count(1)}）")

    print("検定4 一歩の大きさ（距離² のノルム p² + pq − q²）")
    print("   五芒星→円環 φ³ :", B.norm2(T), "ノルム", norm(B.norm2(T)))
    print("   五芒星→円環 φ⁵ :", B.norm2(S), "ノルム", norm(B.norm2(S)))
    print("   連続      :", B.N_CONT, "ノルム", norm(B.N_CONT))
    print("   2つ飛ばし :", B.N_SKIP, "ノルム", norm(B.N_SKIP))

    print("5. まったく同じ景色（回さずに重なる）の五芒星は、残渣（もう一つの影）の側で集まる")
    def sig(a):
        """残渣の側へ写す：ζ^k → ζ^(3k)。足し算と回転だけ。"""
        r = B.ZERO
        for k in range(4):
            r = B.zadd(r, B.zmul((a[k], 0, 0, 0), B.zt(3 * k)))
        return r
    PHIF = (1 + 5 ** 0.5) / 2                     # 観察者側の表示だけ
    val = lambda pq: pq[0] + pq[1] * PHIF
    for R2 in (100, 225, 484):
        pool = [g for g in stars if d[g] >= R2 ** 0.5]
        cls = defaultdict(list)
        for g in pool:
            cls[tuple(sorted(B.zsub(c, g) for c in U if B.phi_lt(B.norm2(B.zsub(c, g)), (R2, 0))))].append(g)
        spread = max([B.norm2(sig(B.zsub(x, y))) for gs in cls.values() for x in gs for y in gs],
                     key=val)
        print(f"   見る半径 √{R2}: 対象 {len(pool)} 個 → {len(cls)} 通りの景色。"
              f"同じ景色どうしの残渣の距離² の最大 = {spread}（≈ {val(spread):.4f}）")
    best = None
    for gs in cls.values():
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                re = B.norm2(B.zsub(gs[i], gs[j])); rs = B.norm2(sig(B.zsub(gs[i], gs[j])))
                if best is None or B.phi_lt(best[0], re): best = (re, rs, gs[i], gs[j])
    re, rs, a, b = best
    print("   見る半径 √484 で同じ景色、実の距離がいちばん離れた2個:", a, b)
    print("   実の距離² =", re, f"（≈ {val(re):.1f}）  残渣の距離² =", rs, f"（≈ {val(rs):.4f}）")
    print("   実の距離² のノルム ＝", norm(re), " 残渣の距離² のノルム ＝", norm(rs), "（同じ整数）")
