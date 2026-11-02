"""第15章　市松の波（原因：和の式と匂いの項）── 検算スクリプト

使い方:  python3 ch15_cause.py <b13_chain_units.py のあるフォルダ>
  ch15_wave.py と同じ地図・同じ規則。s_t(x) は五角形 x に t 歩目に入ってきた値の和。

事前登録
  R1 和の式：0〜71 歩のすべての五角形で  s_{t+1}(x) = Σ_{隣 w} (24/d_w)·s_t(w) − 144·s_{t−1}(x)  が成り立つ。
     負の対照：144 を 143 にすると崩れる
  R2 匂いの項：r_0 = s_0、r_{t+1}(x) = Σ_{隣 w} (12/d_w)·r_t(w)（整数）。
     全部 +1 の出発では届く五角形のどこでも r > 0。
     −1 を混ぜた出発で s = 0 か完全消去になった五角形では、どこでも r = 0。両方成り立てば OK
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch15_wave as W
import ch15_start as P

def sums(st):
    s = Counter()
    for (w, x), a in st.items(): s[w] += a
    return s

def run(signs, T):
    st = {(w, x): signs[P.ORDER.index(w)] for (w, x) in W.start_state()}
    S = [sums(st)]
    for t in range(T): st = W.step(st); S.append(sums(st))
    r = [Counter(S[0])]
    for t in range(T):
        n = Counter()
        for w, a in r[-1].items():
            if a:
                for x in W.nb[w]: n[x] += (12 // len(W.nb[w])) * a
        r.append(n)
    odor = [Counter({c: 1 for c in W.START})]
    for t in range(T):
        n = Counter()
        for c, a in odor[-1].items():
            for x in W.nb[c]: n[x] += a
        odor.append(n)
    return S, r, odor

def check_rec(S, T, k):
    bad = 0
    for t in range(1, T):
        for x in W.cells:
            rhs = sum((24 // len(W.nb[w])) * S[t][w] for w in W.nb[x]) - k * S[t - 1][x]
            if S[t + 1][x] != rhs: bad += 1
    return bad

if __name__ == "__main__":
    T = W.DMAX
    print("R1 和の式")
    S, r, odor = run([1] * 5, T)
    b144, b143 = check_rec(S, T, 144), check_rec(S, T, 143)
    S2, r2, odor2 = run([1, 1, 1, 1, -1], T)
    b144b = check_rec(S2, T, 144)
    print("   全部 +1：式が崩れる (歩, 五角形) の数", b144, "  負の対照（143）:", b143)
    print("   1 枚だけ −1：式が崩れる数", b144b)
    ok1 = b144 == 0 and b144b == 0 and b143 > 0
    print("   →", "OK" if ok1 else "NG")
    print("R2 匂いの項")
    reach = [(t, c) for t in range(T + 1) for c in odor[t] if odor[t][c]]
    pos = all(r[t][c] > 0 for t, c in reach)
    print("   全部 +1：届く五角形（延べ %d）のどこでも r > 0: %s" % (len(reach), pos))
    allok = True
    for bits in [(1, 0, 0, 0, 0), (1, 1, 0, 0, 0), (1, 0, 1, 0, 0)]:
        sg = [(-1 if b else 1) for b in bits]
        S3, r3, od3 = run(sg, T)
        cnt = Counter()
        for t in range(T + 1):
            for c in od3[t]:
                if not od3[t][c]: continue
                if S3[t][c] == 0:
                    cnt["s=0 の五角形"] += 1
                    cnt["そのうち r=0"] += r3[t][c] == 0
                rz = r3[t][c] == 0
                cnt["r=0 の五角形"] += rz
        ok = cnt["s=0 の五角形"] == cnt["そのうち r=0"]
        allok &= ok
        print("   %s：%s" % ("".join("−" if b else "+" for b in bits), dict(cnt)), " s=0 はすべて r=0:", ok)
    ok2 = pos and allok
    print("   →", "OK" if ok2 else "NG")
