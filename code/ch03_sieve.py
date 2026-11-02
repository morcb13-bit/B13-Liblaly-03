"""第3章　2と3で篩う ── 検算スクリプト

使い方:  python3 ch03_sieve.py <b13_chain_units.py のあるフォルダ>
  判定はすべて整数（4成分の番地、長さ² は p + qφ、角は 18° 目盛りの向きの番号）。

事前登録
  検定T3 段の数 n（n = 1..13）は、n ≥ 2 ならすべて 2a + 3b（a, b ≥ 0）と書け、1 だけは書けない。
         段の形で分けると、まっすぐな段（偶数＝2 の倍数）、折れる段のうち 3 の倍数（3, 9）、残り（1, 5, 7, 11, 13）。
         残りから 1 を除いた 5, 7, 11, 13 は素数で、13 までの素数は 2, 3 とこの 4 つですべて。
         2 と 3 で割れない合成数は 25 が最初。外れれば NG
  検定S1 円環で数える：2 で除く段の円環 42 個、3 で除く段（6, 12 は除き済み）の円環 12 個、残る段の円環 37 個
         （段 1 の 1 個と素数の段 36 個）、合計 91。外れれば NG
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
    print("検定T3 2 と 3 で埋まる（エラトステネスの篩）")
    def two_three(n):
        for b in range(n // 3 + 1):
            if (n - 3 * b) % 2 == 0: return (n - 3 * b) // 2, b
        return None
    def is_prime(n): return n > 1 and all(n % p for p in range(2, n))
    straight = {n for n in range(1, 14) if n % 2 == 0 and all(C5.bend(C5.neg(L[i]), L[i + 1]) == 180
                for L in [C5.links(B.unit(n))] for i in range(len(L) - 1))}
    rest = [n for n in range(1, 14) if n not in straight and n % 3 != 0]
    ok4 = True
    for n in range(1, 14):
        ab = two_three(n)
        form = "書けない" if ab is None else " + ".join(["3"] * (ab[1] if n % 2 else 0) + ["2"] * ab[0] + (["3"] * ab[1] if n % 2 == 0 else []))
        kind = "まっすぐ（2 の倍数）" if n in straight else "折れる・3 の倍数" if n % 3 == 0 else ("最初の1" if n == 1 else "折れる・素数")
        print("   段 %2d：%-22s %s" % (n, form, kind))
        ok4 &= (ab is None) == (n == 1)
    ok4 &= straight == {2, 4, 6, 8, 10, 12} and rest == [1, 5, 7, 11, 13] and [n for n in rest if is_prime(n)] == [5, 7, 11, 13]
    ok4 &= [n for n in range(1, 14) if is_prime(n)] == [2, 3, 5, 7, 11, 13]
    first = next(n for n in range(14, 100) if n % 2 and n % 3 and not is_prime(n))
    print("   残った段:", rest, "  2 と 3 で割れない最初の合成数:", first)
    ok4 &= first == 25
    print("   →", "OK" if ok4 else "NG")

    print("検定S1 円環で数える")
    sizes = [len(r) for r in place]
    r2 = sum(sizes[n - 1] for n in range(1, 14) if n % 2 == 0)
    r3 = sum(sizes[n - 1] for n in range(1, 14) if n % 2 and n % 3 == 0)
    rr = [n for n in range(1, 14) if n % 2 and n % 3]
    r1 = sum(sizes[n - 1] for n in rr)
    print("   2 で除く段の円環:", r2, "  3 で除く段の円環:", r3, "  残る段", rr, "の円環:", r1, "（段 1 が", sizes[0], "、素数の段が", r1 - sizes[0], "）  合計:", r2 + r3 + r1)
    oks = (r2, r3, r1, r1 - sizes[0], r2 + r3 + r1) == (42, 12, 37, 36, 91)
    print("   →", "OK" if oks else "NG")
