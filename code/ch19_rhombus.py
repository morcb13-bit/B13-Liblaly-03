"""終章　ひし形に置き換える（de Bruijn の指数と窓）── 検算スクリプト

使い方:  python3 ch19_rhombus.py <b13_chain_units.py のあるフォルダ>
  借りる道具（de Bruijn 1981）：ひし形版ペンローズの頂点は、指数 r と残渣 w で決まる。γ = 0（10 回対称の位置）。
    指数 r = a0 − a1 + a2 − a3 (mod 5)、残渣 w = Σ a_k ζ^(7k)（ζ = ζ10）
    窓 V1 = 頂点 1, ζ², ζ⁴, ζ⁶, ζ⁸ の五角形、V2 = 頂点 1+ζ², ζ²+ζ⁴, … の五角形、V3 = −V2、V4 = −V1（内部）
  内部かどうかは、辺ごとの外積の符号（整数の組）だけで判定する。

事前登録
  DB1 円環の中心 910 個は、すべて指数 3 で、残渣がすべて V3 の内部にある（ひし形の頂点）
  DB2 五角形の中心 6030 枚の指数は 2 と 4 だけで、指数 2 ⟺ 向きの番号が偶数、指数 4 ⟺ 奇数（例外なし）
  対照 五芒星の中心 416 個では、指数 0（どの窓もない）が出て、窓の外に出るものもある（判定が「外」を返せる）
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch05_views as V
B = V.B
Z7 = [B.zt(7 * k) for k in range(4)]
def res(a):
    r = (0, 0, 0, 0)
    for k in range(4):
        if a[k]: r = B.zadd(r, B.zmul((a[k], 0, 0, 0), Z7[k]))
    return r
def idx(a): return (a[0] - a[1] + a[2] - a[3]) % 5
def imsign(w): return B.phi_sign((w[1], w[2] + w[3]))
def cross(u, v): return imsign(B.zmul(B.zconj(u), v))
V1 = [B.zt(2 * j) for j in range(5)]
V2 = [B.zadd(B.zt(2 * j), B.zt(2 * j + 2)) for j in range(5)]
neg = lambda p: tuple(-x for x in p)
WIN = {1: V1, 2: V2, 3: [neg(p) for p in V2], 4: [neg(p) for p in V1]}
def inside(poly, p):
    s = [cross(B.zsub(poly[(i + 1) % 5], poly[i]), B.zsub(p, poly[i])) for i in range(5)]
    return all(x > 0 for x in s) or all(x < 0 for x in s)
def judge(a):
    r = idx(a)
    return r, (r in WIN and inside(WIN[r], res(a)))

if __name__ == "__main__":
    print("DB1 円環の中心")
    rc = Counter(judge(c) for c in V.U)
    print("   （指数, 窓の内部）の内訳:", dict(rc))
    ok1 = rc == Counter({(3, True): len(V.U)})
    print("   →", "OK" if ok1 else "NG")
    print("DB2 五角形の中心の指数と向き")
    pc = Counter((idx(q), o % 2) for q, o in V.cells.items())
    print("   （指数, 向きの番号の偶奇）の内訳:", dict(pc))
    ok2 = set(pc) == {(2, 0), (4, 1)}
    print("   →", "OK" if ok2 else "NG")
    print("対照 五芒星の中心")
    sc = Counter(judge(s) for s in V.stars)
    print("   （指数, 窓の内部）の内訳:", dict(sc))
    print("   指数 0 が出る:", any(r == 0 for r, _ in sc), "  窓の外が出る:", any((not ins) and r != 0 for r, ins in sc))
