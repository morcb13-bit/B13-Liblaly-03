"""第4章　2と3でつなぐ ── 検算スクリプト（段の形の部分）

使い方:  python3 ch04_shape.py <b13_chain_units.py のあるフォルダ>
  つなぎ方の 2 種類は ch04_join.py、二歩・真横・部品・角は ch04_side.py で調べる。

事前登録
  検定T2 奇数と偶数：奇数の段は「連続」のつなぎ方だけで、折れ角はどこも 108°。
         偶数の段は「2つ飛ばし」と「連続」を交互に並べ、折れ角はどこも 180°（まっすぐ）。外れれば NG
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
    print("検定T2 奇数と偶数")
    ok2 = True
    for n in range(1, 14):
        L = C5.links(B.unit(n))
        kinds = "".join(C5.KIND[v] for v in L)
        angs = sorted({C5.bend(C5.neg(L[i]), L[i + 1]) for i in range(len(L) - 1)})
        want = 108 if n % 2 else 180
        ok2 &= all(a == want for a in angs) and (kinds == "連" * (n - 1) if n % 2 else kinds == ("飛連" * n)[:n - 1])
        print("   段 %2d：円環 %2d 個  つなぎ方 %-12s 折れ角 %s" % (n, n, kinds or "—", angs or "—"))
    print("   →", "OK" if ok2 else "NG")

