"""第18章　風を通す、波紋を広げる（波紋）── 検算スクリプト

使い方:  python3 ch18_ripple.py <b13_chain_units.py のあるフォルダ>
  地図は五角形 6030 枚、規則は第15章。石＝五芒星を囲む 5 枚の入口に 1 を置くこと。
  石 A：真ん中の五芒星、石 B：真ん中と同じ色（市松）で最も近い五芒星 (−10, 2, −5, 8)（距離² 47 + 76φ）。
  （同じ色でないと、二つの波は同じ歩に同じ五角形に来ることがなく、重なりが起きない）

事前登録
  R5 重ね合わせ：A と B に同時に石を落とした波は、0〜150 歩の毎歩・すべての辺で、A だけの波と B だけの波の和にちょうど等しい。
     負の対照：−1・0・+1 に丸める規則では、和に等しくならない歩が出る
     あわせて、二つ石の波の値の和は毎歩 12 倍、二乗和は毎歩 144 倍（0〜150 歩）。二つの波が同じ五角形に同時に来る回数（延べ）を記録
"""
import sys
from collections import Counter
sys.argv = sys.argv[:2]
import ch15_wave as W
import ch17_chaos as K
T = 150
O, G2 = (0, 0, 0, 0), (-10, 2, -5, 8)

def plus(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + v
        if r[k] == 0: del r[k]
    return r

if __name__ == "__main__":
    a, b = K.start_at(O), K.start_at(G2); ab = plus(a, b)
    ok = True; okc = True; meet = 0
    s0, q0 = sum(ab.values()), sum(v * v for v in ab.values())
    for t in range(1, T + 1):
        a, b, ab2 = W.step(a), W.step(b), W.step(ab)
        ok &= ab2 == plus(a, b)
        meet += len({w for (w, v) in a} & {w for (w, v) in b})
        s1, q1 = sum(ab2.values()), sum(v * v for v in ab2.values())
        okc &= s1 == 12 * s0 and q1 == 144 * q0
        s0, q0, ab = s1, q1, ab2
    sa, sb = K.start_at(O), K.start_at(G2); sab = plus(sa, sb); bad = 0
    for t in range(1, T + 1):
        sa, sb, sab = K.step_sign(sa), K.step_sign(sb), K.step_sign(sab)
        if sab != plus(sa, sb): bad += 1
    print("R5 重ね合わせ")
    print("   石 A：", O, "  石 B：", G2)
    print("   0〜%d 歩の毎歩・すべての辺で、二つ石の波 ＝ A だけ ＋ B だけ: %s" % (T, ok))
    print("   二つ石の波の値の和 ×12・二乗和 ×144（毎歩）: %s" % okc)
    print("   二つの波が同じ五角形に同時に来た回数（延べ）:", meet)
    print("   負の対照（−1・0・+1 に丸める規則）で和に等しくならない歩: %d / %d" % (bad, T))
    print("   →", "OK" if ok and okc and bad > 0 else "NG")
