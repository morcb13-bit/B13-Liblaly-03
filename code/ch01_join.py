"""第1章　五角形を10個並べる ── 検算スクリプト（ch00_ring.py と同じフォルダで走る）

円環を一つ置き、二つ目の円環を「五角形を少なくとも1枚分け合う」位置にすべて置いてみる。
重なり・向きの食い違いを整数だけで判定し、残るつなぎ方を数える。

重なりの判定：五角形 P1（中心 c1・番地 a1）と P2（c2・a2）の内側が交わるのは、
差 d = c2 − c1 が「P1 と −P2 を足した形」の内側に入るときに限る。
  ・番地の偶奇が違う → −P2 は P1 と同じ向き → 足した形は 2倍の五角形（頂点 2ζ^(a1+2i)）
  ・番地の偶奇が同じ → −P2 は P1 と 36° ずれる → 足した形は十角形（頂点 ζ^k + ζ^(k+1)）
内側かどうかは、辺ごとの「回りの向き」の符号で決める。符号は
  Im(z) の符号 ＝ c1 + (c2 + c3)φ の符号（sin36°>0 を外に括り出した残り）
で、整数だけで判定できる。
"""
from ch00_ring import add, sub, turn, mul, conj, dist2, ZERO, ONE, PHI, PHI2, ring

def phi_sign(p, q):
    """p + qφ の符号（整数だけ）。2(p+qφ) = (2p+q) + q√5"""
    a, b = 2 * p + q, q
    if a >= 0 and b >= 0: return (a > 0 or b > 0) - 0
    if a <= 0 and b <= 0: return -1 if (a < 0 or b < 0) else 0
    s = a * a - 5 * b * b
    if s == 0: return 0
    return (1 if a > 0 else -1) if s > 0 else (1 if b > 0 else -1)

def im_sign(z):
    c0, c1, c2, c3 = z
    return phi_sign(c1, c2 + c3)

def cross_sign(u, v):
    """u から v へ左回りなら +1。Im(ū·v) の符号。"""
    return im_sign(mul(conj(u), v))

def strictly_inside(d, verts):
    """d が凸多角形 verts（左回り）の内側に入るか。辺の上は内側に数えない。"""
    n = len(verts)
    return all(cross_sign(sub(verts[(i + 1) % n], verts[i]), sub(d, verts[i])) > 0 for i in range(n))

def overlap(c1, a1, c2, a2):
    d = sub(c2, c1)
    if d == ZERO:
        return False                                      # 同じ席（分け合い）は別に扱う
    if (a1 - a2) % 2:
        verts = [add(turn(ONE, a1 + 2 * i), turn(ONE, a1 + 2 * i)) for i in range(5)]
    else:
        verts = [add(turn(ONE, k), turn(ONE, k + 1)) for k in range(10)]
    return strictly_inside(d, verts)

def judge(c2, base=ZERO):
    """円環 base と円環 c2 を並べた結果：('ok'|'向き'|'重なり', 分け合う番地の組)"""
    R1, R2 = ring(base), ring(c2)
    shared = [(R1[p], R2[p]) for p in R1 if p in R2]
    if any((x - y) % 2 for x, y in shared):
        return "向き", shared
    for p, a in R1.items():
        for q, b in R2.items():
            if p != q and overlap(p, a, q, b):
                return "重なり", shared
    return "ok", shared

if __name__ == "__main__":
    print("0. 判定器の対照")
    R = ring(ZERO)
    inner = sum(1 for p in R for q in R if p != q and overlap(p, R[p], q, R[q]))
    print("   円環1個の中で重なる組:", inner, "（0 のはず）")
    print("   五角形2枚、偶奇違い・距離1:", overlap(ZERO, 0, ONE, 1), "（重なる）",
          "  偶奇同じ・距離1:", overlap(ZERO, 0, ONE, 0), "（重なる）")
    nb = sub(turn(PHI2, 1), turn(PHI2, 0))                 # 円環で隣り合う二枚の差（辺で接する）
    print("   辺で接する二枚（円環の隣どうし）:", overlap(ZERO, 0, nb, 1), "（重ならない）")
    print("   十角形の頂点ちょうど（偶奇同じ・点で接する）:", overlap(ZERO, 0, add(ONE, turn(ONE, 1)), 0),
          "（重ならない）")
    print("   円環を 36° 回すと自分に重なる:", {turn(p, 1) for p in R} == set(R),
          "  番地は一つ進む:", all(R[turn(p, 1)] == (R[p] + 1) % 10 for p in R))
    print("   鏡に写しても自分に重なる:", {conj(p) for p in R} == set(R))

    print("1. 二つ目の円環を、五角形を分け合う位置にすべて置く")
    cands = {}
    for j in range(10):
        for k in range(10):
            c2 = sub(turn(PHI2, j), turn(PHI2, k))       # 一つ目の番地 j の席に、二つ目の番地 k を重ねる
            if c2 != ZERO:
                cands.setdefault(c2, (j, k))
    print("   置き方:", len(cands), "通り")
    tally = {}
    kinds = {}
    for c2 in cands:
        res, shared = judge(c2)
        tally[res] = tally.get(res, 0) + 1
        if res == "ok":
            gap = (shared[1][0] - shared[0][0]) % 10 if len(shared) == 2 else None
            key = (len(shared), dist2(c2, ZERO), min(gap, 10 - gap) if gap is not None else None)
            kinds[key] = kinds.get(key, 0) + 1
    print("   結果:", tally)
    for (n, d2, gap), cnt in sorted(kinds.items()):
        name = {1: "隣り合う（連続）", 3: "間に2枚（2つ飛ばし）"}.get(gap, f"番地の差 {gap}")
        print(f"   分け合う五角形 {n} 枚・{name}・中心の距離² = {d2[0]} + {d2[1]}φ ： {cnt} 通り")

    print("2. 二つの距離の比")
    print("   距離²の組:", sorted(k[1] for k in kinds))
    (p1, q1), (p2, q2) = sorted(k[1] for k in kinds)
    # (p2 + q2 φ) が (p1 + q1 φ)·φ² に等しいか：φ² = 1 + φ、φ·φ = 1 + φ
    def times_phi(x): return (x[1], x[0] + x[1])         # (x + yφ)φ = y + (x + y)φ
    print("   遠い方の距離² ＝ 近い方の距離² × φ² :", times_phi(times_phi((p1, q1))) == (p2, q2),
          " → 距離の比はちょうど φ")
