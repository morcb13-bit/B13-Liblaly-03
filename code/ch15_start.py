"""第15章　市松の波（出発の形）── 検算スクリプト

使い方:  python3 ch15_start.py <b13_chain_units.py のあるフォルダ>
  ch15_wave.py と同じ地図・同じ規則。出発の 5 枚（真ん中の五芒星を囲む五角形）に +1 か −1 を置く
  （1 枚の入口はすべて同じ符号）。32 通りを、回転と鏡で 8 組にまとめる。0〜71 歩、届く道のある五角形で数える。

事前登録
  P1 回転：−1 を 1 枚だけ置く 5 通りと、その鏡で、s = 0 の数と完全消去の数がすべて一致すれば OK
  P2 配置：8 組それぞれの s = 0 の数・完全消去の数を記録。全部の符号を反転したものは同じ数になる（確認）
  P3 総和と配置：出発の 5 枚の隣の数が同じなら、−1 の枚数が同じとき出発の総和は同じ。
     −1 が 2 枚（隣り合う／離れた）、3 枚（同様）で数が違えば「配置が効く」、同じなら「総和だけで決まる」と記録
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch15_wave as W
import b13_chain_units as B

ORDER = sorted(W.START, key=lambda c: next(j for j in range(10) if B.zrot(B.PHI, j) == c))   # 五芒星のまわりの順
DIRS = [next(j for j in range(10) if B.zrot(B.PHI, j) == c) for c in ORDER]

def run(signs, T):
    st = {(w, x): signs[ORDER.index(w)] for (w, x) in W.start_state()}
    odor = Counter({c: 1 for c in W.START}); zs = zc = 0
    for t in range(T + 1):
        if t:
            st = W.step(st); o2 = Counter()
            for c, a in odor.items():
                for x in W.nb[c]: o2[x] += a
            odor = o2
        s = Counter(); live = set()
        for (w, x), a in st.items(): s[w] += a; live.add(w)
        for c in odor:
            if not odor[c]: continue
            if c not in live: zc += 1
            elif s[c] == 0: zs += 1
    return zs, zc

def canon(bits):
    vs = []
    for b in (bits, bits[::-1]):
        for r in range(5): vs.append(b[r:] + b[:r])
    return min(vs)

if __name__ == "__main__":
    T = W.DMAX
    print("出発の 5 枚：向き", DIRS, "（36° 単位）  隣の数", [len(W.nb[c]) for c in ORDER],
          "  入口の数の合計", sum(len(W.nb[c]) for c in ORDER))
    res = {}
    for m in range(32):
        bits = tuple((m >> i) & 1 for i in range(5))                 # 1 が −1
        res[bits] = run([(-1 if b else 1) for b in bits], T)
    print("P1 回転（−1 を 1 枚だけ）")
    one = [res[tuple(1 if i == k else 0 for i in range(5))] for k in range(5)]
    print("   −1 の位置 0〜4 の (s = 0, 完全消去):", one)
    ok1 = len(set(one)) == 1
    print("   →", "OK" if ok1 else "NG")
    print("P2 配置（8 組。1 が −1、並びは五芒星のまわりの順）")
    groups = {}
    for bits, r in res.items(): groups.setdefault(canon(bits), []).append(r)
    for g in sorted(groups, key=lambda b: (sum(b), b)):
        rs = groups[g]
        print("   %s  −1 が %d 枚  出発の総和 %+d  通り %2d  (s = 0, 完全消去): %s" %
              ("".join("−" if b else "+" for b in g), sum(g), sum((-1 if b else 1) * len(W.nb[c]) for b, c in zip(g, ORDER)),
               len(rs), sorted(set(rs))))
    inv = all(res[b] == res[tuple(1 - x for x in b)] for b in res)
    print("   全部の符号を反転すると同じ数（確認）:", inv)
    print("P3 総和と配置")
    adj2, sep2 = res[(1, 1, 0, 0, 0)], res[(1, 0, 1, 0, 0)]
    adj3, sep3 = res[(1, 1, 1, 0, 0)], res[(1, 1, 0, 1, 0)]
    same_deg = len({len(W.nb[c]) for c in ORDER}) == 1
    print("   出発の 5 枚の隣の数が同じ:", same_deg)
    print("   −1 が 2 枚：隣り合う", adj2, " 離れた", sep2, " →", "配置が効く" if adj2 != sep2 else "総和だけで決まる")
    print("   −1 が 3 枚：隣り合う", adj3, " 離れた", sep3, " →", "配置が効く" if adj3 != sep3 else "総和だけで決まる")
