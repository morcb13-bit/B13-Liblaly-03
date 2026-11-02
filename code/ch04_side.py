"""第4章　2と3でつなぐ（二歩・真横・部品と桁・角）── 検算スクリプト

使い方:  python3 ch04_side.py <b13_chain_units.py のあるフォルダ>
  第4章の単位（円環を一列につないだ鎖）を、真横から一本の線として見る。
  判定はすべて整数。向きは 18° を一目盛りとした番号（0〜19）で持つ。

事前登録
  検定1 鎖の折れ角：奇数の鎖はどこも 108°、偶数の鎖はどこも 180°。ほかの角が出れば NG
  検定2 二歩の一致：奇数の二歩と偶数の二歩は同じ一歩 w。2k 番目の環は奇偶とも k·w。合わなければ NG
  検定3 真横から：偶数の鎖はすべての環が軸 w の上。奇数の鎖は偶数番目の環が軸の上、奇数番目の環がすべて軸の同じ側。
        鏡に写した奇数の鎖は反対側。外れれば NG
  検定4 部品：奇数は V字（3）と 1、偶数は 2 だけで、3×V字＋桁 ＝ n（桁は −1, 0, +1）。合わなければ NG
  検定5 910環の網の折れ角：連と連・連と飛の角だけが出て、飛と飛は 0。出れば NG
"""
import sys
from collections import Counter
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

# つなぎ方20本の向きの番号（18° 目盛り）。CONT[k] = φ²(ζ⁴−1)ζ^k は 162° + 36k、SKIP[k] = φ²(ζ²−1)ζ^k は 126° + 36k
DIR = {}
for k in range(10):
    DIR[B.CONT[k]] = (9 + 2 * k) % 20
    DIR[B.SKIP[k]] = (7 + 2 * k) % 20
KIND = {**{v: "連" for v in B.CONT}, **{v: "飛" for v in B.SKIP}}

def neg(v): return tuple(-x for x in v)

def im_sign(z):
    """Im(z) の符号。sin36° > 0 を括り出すと c1 + (c2 + c3)φ の符号（整数で判定）。"""
    return B.phi_sign((z[1], z[2] + z[3]))

def side(w, v):
    """v が軸 w の左（+1）・上（0）・右（−1）のどれか。Im(w̄·v) の符号。"""
    return im_sign(B.zmul(B.zconj(w), v))

def bend(a, b):
    """環から見た二本のつなぎの向きの間の角（度）。a・b は環から外へ向かうベクトル。"""
    d = (DIR[a] - DIR[b]) % 20
    return 18 * min(d, 20 - d)

def links(cs): return [B.zsub(cs[i + 1], cs[i]) for i in range(len(cs) - 1)]

def parts(n):
    """部品の並び。奇数は V字（3）と 1、偶数は 2。桁 d は n = 3k + d（d ∈ {−1, 0, +1}）。"""
    if n % 2 == 0:
        return [("2", [i, i + 1]) for i in range(0, n, 2)], None
    k, d = (n + 1) // 3, None
    for d in (0, 1, -1):
        if (n - d) % 3 == 0: k = (n - d) // 3; break
    if d == 0:
        return [("V", [3 * i, 3 * i + 1, 3 * i + 2]) for i in range(k)], 0
    h = k // 2
    left = [("V", [3 * i, 3 * i + 1, 3 * i + 2]) for i in range(h)]
    if d == 1:
        c = 3 * h
        return left + [("1", [c])] + [("V", [c + 1 + 3 * i, c + 2 + 3 * i, c + 3 + 3 * i]) for i in range(h)], 1
    c = 3 * h - 1                                   # 左の最後の環を右の最初と分け合う
    return left + [("V", [c + 3 * i, c + 1 + 3 * i, c + 2 + 3 * i]) for i in range(h)], -1

def is_prime(n): return n > 1 and all(n % p for p in range(2, n))

if __name__ == "__main__":
    print("検定1 鎖の折れ角と、つなぎ方の種類")
    ok1 = True
    for n in range(1, 14):
        L = links(B.unit(n))
        kinds = "".join(KIND[v] for v in L)
        angs = [bend(neg(L[i]), L[i + 1]) for i in range(len(L) - 1)]
        want = 108 if n % 2 else 180
        ok1 &= all(a == want for a in angs)
        print(f"   n={n:2d}  つなぎ方 {kinds or '—':14s} 折れ角 {sorted(set(angs)) or '—'}")
    print("   →", "OK" if ok1 else "NG")

    print("検定2 二歩の一致")
    w = B.zadd(B.C_A, B.C_B)
    print("   奇数の二歩 連+連 =", w, "  偶数の二歩 飛+連 =", B.zadd(B.S_AX, B.C_AX), "  同じ:", w == B.zadd(B.S_AX, B.C_AX))
    ok2 = all(B.unit(n)[2 * k] == tuple(k * x for x in w) for n in range(1, 14) for k in range((n + 1) // 2))
    print("   2k 番目の環が k·w（n = 1..13 のすべて）:", ok2, "→", "OK" if ok2 and w == B.zadd(B.S_AX, B.C_AX) else "NG")
    print("   w の距離² =", B.norm2(w), "  連の距離² =", B.N_CONT, "  飛の距離² =", B.N_SKIP)

    print("検定3 真横から見る（軸 w の左 +1・上 0・右 −1）")
    ok3 = True
    for n in range(1, 14):
        s = [side(w, c) for c in B.unit(n)]
        if n % 2 == 0: ok3 &= all(x == 0 for x in s)
        else: ok3 &= all(s[i] == 0 for i in range(0, n, 2)) and len({s[i] for i in range(1, n, 2)} - {0}) <= 1 and 0 not in s[1::2]
        print(f"   n={n:2d}  {s}")
    m = [side(w, c) for c in [B.ZERO, B.C_B, w]]
    print("   鏡に写した奇数の鎖（連の順を入れ替える）の1番目:", m[1], " もとの鎖の1番目:", side(w, B.C_A), "（逆側）")
    ok3 &= m[1] == -side(w, B.C_A)
    print("   →", "OK" if ok3 else "NG")
    # 偶数の鎖：飛と連は同じ向き（軸の上）で、長さは 1 : φ。足すと w。|w| = |飛|·φ² ＝ |飛| + |連|
    S2, C2, W2 = B.N_SKIP, B.N_CONT, B.norm2(w)
    phi4 = (2, 3)                                   # φ⁴ = 2 + 3φ
    prod = (S2[0] * phi4[0] + S2[1] * phi4[1], S2[0] * phi4[1] + S2[1] * phi4[0] + S2[1] * phi4[1])
    print("   |飛|²·φ⁴ =", prod, " |w|² =", W2, " 一致:", prod == W2, "（|w| = |飛|·φ² = |飛| + |連|）")

    print("検定4 部品に分ける（V＝3、1、2）と桁")
    ok4 = True
    for n in range(1, 30):
        ps, d = parts(n)
        cover = sorted(i for _, ix in ps for i in ix)
        shared = len(cover) - len(set(cover))
        nv = sum(1 for t, _ in ps if t == "V")
        good = set(cover) == set(range(n)) and (n % 2 == 0 or 3 * nv + d == n)
        if n % 2 and n <= 13:
            L = links(B.unit(n))
            for t, ix in ps:                       # V字は 連・連 の 108° の三環か
                if t == "V":
                    good &= bend(neg(L[ix[0]]), L[ix[1]]) == 108
        ok4 &= good
        tag = "素数" if is_prime(n) else ""
        txt = " ".join(t for t, _ in ps)
        print(f"   n={n:2d}  {'偶' if n % 2 == 0 else '奇'}  桁 {'' if d is None else f'{d:+d}':>3s}  重ね {shared}  部品 {txt:28s} {tag}")
    print("   →", "OK" if ok4 else "NG")
    pr = [p for p in range(5, 100) if is_prime(p)]
    print("   5〜97 の素数の、2で割った余り・3の桁:", Counter((p % 2, parts(p)[1]) for p in pr))

    print("検定5 910環の網で、一つの環から出る二本のつなぎの間の角")
    T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
    rows, place, offs = B.build_stack()
    U, seen = [], set()
    for k in range(10):
        t = B.zrot(T if k % 2 == 0 else S, k)
        for c in sum(place, []):
            c2 = B.zadd(B.zrot(c, k), t)
            if c2 not in seen: seen.add(c2); U.append(c2)
    Us = set(U)
    cnt = Counter()
    for c in U:
        out = [v for v in DIR if B.zadd(c, v) in Us]
        for i in range(len(out)):
            for j in range(i + 1, len(out)):
                cnt[("".join(sorted(KIND[out[i]] + KIND[out[j]])), bend(out[i], out[j]))] += 1
    for key in sorted(cnt): print(f"   {key[0]}  {key[1]:3d}°: {cnt[key]}")
    print("   →", "OK" if not any(k[0] == "飛飛" for k in cnt) else "NG")
