"""第15章　市松の波（消えない理由）── 検算スクリプト

使い方:  python3 ch15_why.py <b13_chain_units.py のあるフォルダ>
  ch15_wave.py と同じ地図・同じ規則。届く道のある五角形の、入ってきた値の和 s を調べる。

事前登録
  N1 2 と 3 の冪：各歩で s を割り切る 2 の冪 a と 3 の冪 b を数え、同じ歩ならどの五角形でも同じかを見る
  N2 剰余：s / (2^a·3^b) を 2・3・5・13 で割った余りで分け、余り 0 が一度も出ない法があるかを見る
  N3 最小の |s|：各歩の最小の |s| を 2^a·3^b と比べる
  判定：各歩で 2^a·3^b がどの五角形でも共通で、残りがある法でいつも 0 でないなら OK（71 歩までの観測で s = 0 は起こりえない形）
  対照模型：出発の 5 枚のうち 1 枚の入口だけ −1 にして同じことを数え、s = 0 が出るかを記録する
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch15_wave as W

def v(n, p):
    k = 0
    while n % p == 0: n //= p; k += 1
    return k

def survey(st0, T):
    st = st0; odor = Counter({c: 1 for c in W.START}); rows = []
    for t in range(T + 1):
        if t:
            st = W.step(st); o2 = Counter()
            for c, a in odor.items():
                for x in W.nb[c]: o2[x] += a
            odor = o2
        s = Counter()
        for (w, x), a in st.items(): s[w] += a
        S = [s[c] for c in odor if odor[c]]
        rows.append(S)
    return rows

if __name__ == "__main__":
    T = W.DMAX
    for name, st0 in [("この章の出発（すべて +1）", W.start_state()),
                      ("対照模型（1 枚の入口だけ −1）", {k: (-1 if k[0] == W.START[0] else 1) for k in W.start_state()})]:
        rows = survey(st0, T)
        zero = sum(1 for S in rows for x in S if x == 0)
        common, ab, res, mins = True, [], {m: Counter() for m in (2, 3, 5, 13)}, []
        for t, S in enumerate(rows):
            nz = [x for x in S if x]
            if not nz: ab.append(None); continue
            a = min(v(x, 2) for x in nz); b = min(v(x, 3) for x in nz)
            if any(v(x, 2) != a or v(x, 3) != b for x in nz) or len(nz) != len(S): common = False
            ab.append((a, b)); g = 2 ** a * 3 ** b
            for x in nz:
                for m in res: res[m][(x // g) % m] += 1
            mins.append(min(abs(x) for x in nz) // g)
        print(name)
        print("   届く五角形で s = 0（延べ）:", zero)
        print("   各歩の共通の冪 (a, b)（0〜12 歩）:", ab[:13], "…")
        print("   どの歩でも、すべての五角形の冪が同じ:", common)
        for m in res:
            print("   s/(2^a·3^b) を %2d で割った余り:" % m, dict(sorted(res[m].items())))
        print("   最小の |s| / (2^a·3^b)（0〜12 歩）:", mins[:13])
        nz0 = [m for m in res if res[m][0] == 0]
        print("   余り 0 が一度も出ない法:", nz0)
        if name.startswith("この章"):
            ok = common and bool(nz0) and zero == 0
            print("   → N1〜N3", "OK" if ok else "NG")
