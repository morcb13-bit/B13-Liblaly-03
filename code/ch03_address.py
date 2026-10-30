"""第3章　番地をふる ── 検算スクリプト（ch00〜ch02 と同じフォルダで走る）

使い方:  python3 ch03_address.py rings_integer.json
平面の位置 (a0, a1, a2, a3) の4つの整数を、それぞれ 13 で割った余りにする。これを「13の番地」と呼ぶ。
判定はすべて整数の足し算・引き算・余り。
"""
import sys, json
from ch00_ring import add, sub, turn, dist2, ZERO, ring

def bal(x, m):
    """m で割った余りを、真ん中を 0 にした平衡の形（13 なら −6〜+6）で返す。"""
    r = x % m
    return r - m if r > m // 2 else r

def addr(p, m=13):
    return tuple(bal(x, m) for x in p)

def distinct(points, m):
    return len({addr(p, m) for p in points})

def phi_mul3(v, m):
    """(x + yφ)·φ³ を m で割った余りで。φ³ = 1 + 2φ。足し算だけ：(x + 2y) + (2x + 3y)φ"""
    x, y = v
    return ((x + y + y) % m, (x + x + y + y + y) % m)

def order_phi3(m):
    v, n = (1, 0), 0
    while True:
        v = phi_mul3(v, m); n += 1
        if v == (1 % m, 0): return n

def norm(pq):
    """Z[φ] の数 p + qφ のノルム p² + pq − q²（整数）。"""
    p, q = pq
    return p * p + p * q - q * q

def factor(n):
    n, out, d = abs(n), [], 2
    while d * d <= n:
        while n % d == 0: out.append(d); n //= d
        d += 1
    if n > 1: out.append(n)
    return out

def irreducible_quartic(p):
    """x⁴ − x³ + x² − x + 1 が p で割った世界で、1次・2次の因数を持たないか（総当たり）。"""
    f = lambda x: (x**4 - x**3 + x*x - x + 1) % p
    if any(f(x) == 0 for x in range(p)): return False
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):        # (x² + ax + b)(x² + cx + d)
                    if ((a + c) % p, (b + d + a * c) % p, (a * d + b * c) % p, (b * d) % p) == \
                       ((-1) % p, 1 % p, (-1) % p, 1 % p):
                        return False
    return True

def lattice_index(vecs):
    """整数ベクトルたちが足し算・引き算で作る格子の、整数4成分全体に対する指数（互除法で消去）。"""
    rows = [list(v) for v in vecs]
    piv = []
    for col in range(4):
        while True:
            nz = [r for r in rows if r[col] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(r[col]))
            a = nz[0]
            for r in nz[1:]:
                q = r[col] // a[col]
                for i in range(4): r[i] -= q * a[i]
        nz = [r for r in rows if r[col] != 0]
        if not nz: return 0
        piv.append(abs(nz[0][col])); rows.remove(nz[0])
    out = 1
    for x in piv: out *= x
    return out

if __name__ == "__main__":
    rings = [tuple(r) for r in json.load(open(sys.argv[1]))["rings"]]
    cells = {}
    for c in rings:
        for p, a in ring(c).items(): cells.setdefault(p, a)
    print("0. 平面：円環", len(rings), "個、五角形", len(cells), "枚")

    print("1. 4つの整数をそれぞれ m で割った余りにする（番地の衝突）")
    for m in (3, 7, 11, 13, 17):
        print(f"   m = {m:2d}: 番地の数 {m**4:6d}  円環 {distinct(rings, m)}/60  五角形 {distinct(cells, m)}/{len(cells)}")
    print(f"   m = 169: 円環 {distinct(rings, 169)}/60  五角形 {distinct(cells, 169)}/{len(cells)}")
    c0 = rings[0]
    print("   最初の円環の中心", c0, "の13の番地:", addr(c0))

    print("   13 ＝ 3² + 2² ＝", 3 * 3 + 2 * 2)
    print("2. 13 で割った世界の性質")
    print("   x² − x − 1 の根が13の世界にあるか:", [x for x in range(13) if (x*x - x - 1) % 13 == 0],
          " 11の世界:", [x for x in range(11) if (x*x - x - 1) % 11 == 0])
    print("   x⁴ − x³ + x² − x + 1 が割れないか:  13:", irreducible_quartic(13), "  11:", irreducible_quartic(11),
          "  3:", irreducible_quartic(3), "  7:", irreducible_quartic(7))
    zero13 = [(p, q) for p in range(13) for q in range(13) if (p*p + p*q - q*q) % 13 == 0]
    zero11 = [(p, q) for p in range(11) for q in range(11) if (p*p + p*q - q*q) % 11 == 0]
    print("   p² + pq − q² が 13 で割り切れる (p, q):", zero13, "  11 なら", len(zero11), "組")

    print("3. 長さの側（距離² のノルム p² + pq − q²）に出る素数")
    rn = sorted({norm(dist2(c, ZERO)) for c in rings})
    print("   円環の中心:", {n: factor(n) for n in rn})
    cn = sorted({norm(dist2(p, ZERO)) for p in cells})
    primes = sorted({f for n in cn for f in factor(n)})
    print("   五角形の中心に出る素数:", primes)
    print("   13 で割り切れるノルム: 円環", sum(n % 13 == 0 for n in rn), " 五角形", sum(n % 13 == 0 for n in cn))

    print("4. φ³ を掛け続けて番地が戻るまで")
    for m in (13, 169, 2197):
        print(f"   {m:5d} で割った世界: {order_phi3(m)} 回")
    v, F = (1, 0), []
    for n in range(1, 29):
        v = (v[0] + v[1] + v[1], v[0] + v[0] + v[1] + v[1] + v[1])     # 割らずに φ³ を掛ける
    print("   28 回掛けた φ⁸⁴ ＝", v[0], "+", v[1], "φ")
    print("   その13の番地:", (v[0] % 13, v[1] % 13), "（(1, 0) なら戻っている）")

    print("5. 足し算で番地の全部に届くか")
    print("   円環の中心60個が足し算・引き算で作る格子の指数:", lattice_index(rings), "（1 なら整数4成分すべてに届く）")
    print("   深さ k（13^k で割る）で番地は 13^(4k) 個。k = 13 なら 13^52 =", 13 ** 52)
