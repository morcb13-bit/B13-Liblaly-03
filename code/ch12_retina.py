"""第12章　網膜になる ── 検算スクリプト

使い方:  python3 ch12_retina.py <b13_chain_units.py のあるフォルダ>
  網膜は 910環の図の五角形 6030 枚。信号は平面に置いた正五角形（中心・向き・大きさ φ^n）。
  励起：信号の輪郭の点を含む五角形が励起する（点が五角形の内側か辺の上か、を辺ごとの符号で判定）。
  着地：信号の頂点 5 つが、五角形の中心・五芒星・円環の中心のどれかに一致すること。
  判定はすべて整数（4成分の番地。点を 1/m ずらすときは全体を m 倍して整数のまま比べる。符号は p + qφ）。

事前登録（外周から 12 以上内側の五芒星 296 個を中心とする。外周の索引だけは第4章と同じ実数の観察者側）
  検定1 励起：五芒星を中心に、向きを五芒星にそろえた正五角形の輪郭を入れる。励起した五角形の集合が
        中心のまわりの 72° 回転で変わらない（5 回対称に励起する）五芒星の数を、大きさ φ¹〜φ⁵、向き 2 通りで数える。
        すべての五芒星で対称になる大きさと、真ん中の 1 個だけで対称になる大きさが分かれる。分かれなければ NG
  検定2 着地：頂点 5 つが何に着地するかを、大きさ φ⁰〜φ⁷、向き 2 通り（そろえる／36° 回す）で数える。
        φ¹〜φ⁶ では、どちらかの向きで 296 個すべてが 5 つとも着地する。φ⁰ と φ⁷ では、どちらの向きでもそうならない。外れれば NG
  検定3 ずれ：着地していた信号の中心を 1/10 または 1/100 だけずらす。着地は 0 になるが、励起した五角形の集合は変わらない。外れれば NG
"""
import sys
from collections import Counter
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

T = (-2, 2, 0, 3); S0 = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S0, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = set(B.fits(U)); Uset = set(U)
stars = {}
for p in cells:
    for j in range(10):
        g = B.zsub(p, B.zrot(B.PHI, j))
        if g in cells or g in stars: continue
        if all(B.zadd(g, B.zrot(B.PHI, j % 2 + 2 * i)) in cells for i in range(5)): stars[g] = j % 2
STEP = [B.zrot(B.PHI, j) for j in range(10)]
ORI = {c: [j for j in range(10) if B.zadd(c, STEP[j]) in cells][0] % 2 for c in cells}   # 隣の向きの偶奇
R = [B.ONE]
for _ in range(8): R.append(B.zmul(R[-1], B.PHI))                                        # φ^n

def smul(k, a): return tuple(k * x for x in a)
def im_sign(w): return B.phi_sign((w[1], w[2] + w[3]))

# 近くの五角形だけを調べるための升目（升目の割り当てにだけ実数を使う。判定は整数）
GRID = {}
for c in cells:
    x, y = B.xy(c); GRID.setdefault((round(x / 2), round(y / 2)), []).append(c)
def near(p, m):
    x, y = B.xy(p); x /= m; y /= m
    return [c for dx in (-1, 0, 1) for dy in (-1, 0, 1) for c in GRID.get((round(x / 2) + dx, round(y / 2) + dy), [])]

def inside(c, p, m):
    """m 倍した点 p が五角形 c（頂点 c + ζ^(o+2i)、反時計回り）の内側か辺の上か。"""
    o = 1 - ORI[c]
    V = [smul(m, B.zadd(c, B.zt(o + 2 * i))) for i in range(5)]
    return all(im_sign(B.zmul(B.zconj(B.zsub(V[(i + 1) % 5], V[i])), B.zsub(p, V[i]))) >= 0 for i in range(5))

PER = {1: 6, 2: 10, 3: 16, 4: 26, 5: 42, 6: 68}
def excite(G, a, n, m):
    """m 倍した中心 G、向き a、大きさ φ^n の正五角形の輪郭を一辺 PER[n] 点で標本し、励起した五角形を返す。"""
    V = [B.zadd(G, smul(m, B.zrot(R[n], a + 2 * i))) for i in range(5)]
    per, E = PER[n], set()
    for i in range(5):
        A, e = V[i], B.zsub(V[(i + 1) % 5], V[i])
        for s in range(per):
            p = B.zadd(smul(per, A), smul(s, e))
            for c in near(p, m * per):
                if inside(c, p, m * per): E.add(c)
    return E

def symmetric(E, g):
    return all(B.zadd(g, B.zrot(B.zsub(c, g), 2)) in E for c in E)

def kind(p):
    return "五角形" if p in cells else "五芒星" if p in stars else "円環" if p in Uset else "なし"

if __name__ == "__main__":
    sys.path.insert(0, ".")
    import ch04_views as V4
    far = [g for g in stars if V4.to_outer(g) >= 12]
    print("0. 網膜：五角形", len(cells), "枚、五芒星", len(stars), "個（外周から 12 以上内側", len(far), "個）")

    print("検定1 励起（5 回対称に励起する五芒星の数）")
    allsym, onlyc = [], []
    for n in range(1, 6):
        row = []
        for off in (0, 1):
            cnt = sum(symmetric(excite(g, stars[g] + off, n, 1), g) for g in far)
            row.append(cnt)
        print("   φ^%d  そろえる: %3d   36°回す: %3d" % (n, row[0], row[1]))
        for off, v in zip((0, 1), row):
            (allsym if v == len(far) else onlyc if v == 1 else []).append((n, off))
    E0 = excite((0, 0, 0, 0), stars[(0, 0, 0, 0)], 1, 1)
    print("   すべてで対称:", allsym, "  真ん中だけ:", onlyc, "  （φ¹ の励起枚数", len(E0), "）")
    print("   →", "OK" if allsym and onlyc else "NG")

    print("検定2 着地（頂点 5 つが着地した先）")
    ok2 = True
    for n in range(0, 8):
        full = []
        for off, name in ((0, "そろえる"), (1, "36°回す")):
            cnt = Counter()
            for g in far:
                kinds = Counter(kind(B.zadd(g, B.zrot(R[n], stars[g] + off + 2 * i))) for i in range(5))
                cnt[tuple(sorted(kinds.items()))] += 1
            top = cnt.most_common(1)[0]
            all5 = len(cnt) == 1 and dict(top[0]).get("なし", 0) == 0 and len(top[0]) == 1
            full.append(all5)
            print("   φ^%d %s:" % (n, name), {"・".join("%s%d" % kv for kv in k): v for k, v in cnt.most_common(3)})
        if 1 <= n <= 6: ok2 &= any(full)
        else: ok2 &= not any(full)
    print("   →", "OK" if ok2 else "NG")

    print("検定3 ずれ")
    ok3 = True
    for m in (10, 100):
        for n, off in ((1, 0), (2, 1), (4, 0)):
            land, diff = 0, 0
            for g in far:
                a = stars[g] + off
                G = B.zadd(smul(m, g), B.ONE)                                   # 中心を 1/m だけずらす
                land += sum(1 for i in range(5) if all(x % m == 0 for x in B.zadd(G, smul(m, B.zrot(R[n], a + 2 * i))))
                            and kind(tuple(x // m for x in B.zadd(G, smul(m, B.zrot(R[n], a + 2 * i))))) != "なし")
                diff += len(excite(g, a, n, 1) ^ excite(G, a, n, m))
            print("   1/%d ずらす  φ^%d：着地した頂点 %d、励起が変わった五角形 %d" % (m, n, land, diff))
            ok3 &= land == 0 and diff == 0
    print("   →", "OK" if ok3 else "NG")
