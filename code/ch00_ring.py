"""序章　二つの1を平面に置く ── 検算スクリプト（このファイルだけで走る）

平面の点は整数4つ (a0, a1, a2, a3) = a0 + a1·ζ + a2·ζ² + a3·ζ³ で書く。
ζ は 36° 回った一歩。ζ⁴ = −1 + ζ − ζ² + ζ³ で畳むので、整数は4つで足りる。
判定はすべて整数の一致と足し算・引き算だけで行う。
"""

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))

def turn(a, k=1):
    """ζ を k 回掛ける（36° を k 回）。中身は足し算と引き算だけ。"""
    v = a
    for _ in range(k % 10):
        c0, c1, c2, c3 = v
        v = (-c3, c0 + c3, c1 - c3, c2 + c3)
    return v

def mul(a, b):
    """掛け算 ＝ b を回したものを a の各成分の回数だけ足す。"""
    out = (0, 0, 0, 0)
    for k, n in enumerate(a):
        t = turn(b, k)
        for _ in range(abs(n)):
            out = add(out, t) if n > 0 else sub(out, t)
    return out

def conj(a):
    """鏡に写す：ζ^k → ζ^(10−k)。"""
    out = (0, 0, 0, 0)
    for k, n in enumerate(a):
        t = turn((1, 0, 0, 0), (10 - k) % 10)
        for _ in range(abs(n)):
            out = add(out, t) if n > 0 else sub(out, t)
    return out

def dist2(a, b):
    """距離の2乗を p + qφ の整数の組 (p, q) で返す。"""
    d = sub(a, b)
    z = mul(d, conj(d))
    assert z[1] == 0 and z[3] == -z[2]          # 鏡と掛けると軸の上に落ちる
    return (z[0] - z[2], z[2])                  # z0 + z2(ζ² − ζ³) = (z0 − z2) + z2·φ

ZERO, ONE = (0, 0, 0, 0), (1, 0, 0, 0)
PHI  = (1, 0, 1, -1)                            # φ  = 1 + ζ² − ζ³
PHI2 = add(PHI, ONE)                            # φ² = φ + 1

def ring(c):
    """円環：中心 c のまわり φ²·ζ^k に五角形10枚。位置 → 番地 k。"""
    return {add(c, turn(PHI2, k)): k for k in range(10)}

if __name__ == "__main__":
    print("1. 一歩と回る一歩")
    print("   ζ⁵ =", turn(ONE, 5), "（半周で −1）   ζ¹⁰ =", turn(ONE, 10))
    print("   φ  =", PHI, "  φ·φ =", mul(PHI, PHI), "  φ + 1 =", PHI2,
          "  一致:", mul(PHI, PHI) == PHI2)

    print("2. 円環1個")
    R = ring(ZERO)
    pos = sorted(R, key=R.get)
    d = {dist2(pos[k], pos[(k + 1) % 10]) for k in range(10)}
    print("   五角形", len(R), "枚。隣どうしの距離² =", d, "＝ 1 + φ ＝ φ² → 距離 φ")
    print("   中心からの距離² =", {dist2(p, ZERO) for p in pos}, "＝ 2 + 3φ ＝ φ⁴ → 距離 φ²")

    print("3. 円環5個（真ん中を空けて囲む）")
    # 円環の中心は φ³ の距離に5つ。φ³ = φ² + φ
    PHI3 = add(PHI2, PHI)
    G = (3, 0, 2, -2)                           # 最初の円環の中心（中心からの距離 φ³）
    centers = [turn(G, 2 * m) for m in range(5)]
    print("   円環の中心:", centers)
    print("   中心からの距離² =", {dist2(c, ZERO) for c in centers}, "＝ 5 + 8φ ＝ φ⁶ → 距離 φ³")
    print("   G の距離² と φ³ の距離² が同じ:", dist2(G, ZERO) == dist2(PHI3, ZERO))
    print("   中心の和 =", tuple(sum(c[i] for c in centers) for i in range(4)))
    seats = [p for c in centers for p in ring(c)]
    cells = {}
    for c in centers:
        for p, a in ring(c).items():
            cells.setdefault(p, set()).add(a % 2)
    print("   席", len(seats), "→ 五角形", len(cells), "枚。二つの円環が分け合う五角形",
          len(seats) - len(cells), "枚")
    print("   分け合う五角形で向き（番地の偶奇）が食い違うもの:",
          sum(1 for s in cells.values() if len(s) > 1))
    # 五角形どうしの重なり：中心どうしの距離² が φ²（＝辺で接する距離）より小さい組があるか
    ps = list(cells)
    def less(x, y):                             # p + qφ < r + sφ を整数だけで判定
        p, q = x[0] - y[0], x[1] - y[1]         # p + qφ < 0 か？  2p + q < −q√5
        a, b = 2 * p + q, q
        if a < 0 and b <= 0: return True
        if a >= 0 and b >= 0: return False
        return (a * a > 5 * b * b) if a < 0 else (a * a < 5 * b * b)
    close = sum(1 for i in range(len(ps)) for j in range(i + 1, len(ps))
                if less(dist2(ps[i], ps[j]), (1, 1)))
    print("   φ より近い五角形の組:", close, "（0 なら重ならない）")
    def star_at_center(cells):
        """真ん中が五芒星の穴：原点に五角形がなく、距離 φ の5方向（72° おき）に五角形がある。"""
        five = [turn(PHI, 2 * i) for i in range(5)]
        return ZERO not in cells and all(p in cells for p in five)
    print("   真ん中に五芒星の穴が残る:", star_at_center(cells))

    print("4. 負の対照：円環を中心から φ² の距離に置いた場合")
    G2 = (-2, 0, -1, 1)
    cells2 = {p for m in range(5) for p in ring(turn(G2, 2 * m))}
    print("   中心からの距離² =", dist2(G2, ZERO), "＝ 2 + 3φ ＝ φ⁴ → 距離 φ²")
    print("   五角形", len(cells2), "枚。原点に五角形がある:", ZERO in cells2,
          "  真ん中に五芒星の穴が残る:", star_at_center(cells2))
