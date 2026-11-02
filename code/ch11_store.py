"""第11章　光を集めて貯める ── 検算スクリプト

使い方:  python3 ch11_store.py <b13_chain_units.py のあるフォルダ>
  地図は第10章と同じ、910環の図の五角形 6030 枚（辺を共有する五角形が隣）。
  光は五角形から隣の五角形へ、一刻に一歩だけ進む（どの光も同じ速さ）。光の量は整数で、分けも増やしもしない。
  光の道は円環（五角形 10 枚の輪）。光は自分の円環の上を決めた向きに回り、別の円環と共有する五角形に来たとき、
  その円環の中心のほうが真ん中に近ければ乗り換える。
  判定はすべて整数（4成分の番地、距離² と向きの比較は p + qφ の符号）。

事前登録
  検定1 貯める輪：真ん中の五芒星を囲む五角形 10 枚は、隣どうしで閉じた輪になる。輪の外から輪へ入る隣は 5 本だけ（入口 5 つ）。
        輪からの段数で層に分けると、どの層の枚数も 5 の倍数。外れれば NG
  検定2 円環の道：910 の円環はすべて、隣どうしで閉じた五角形 10 枚の輪。五角形が属する円環の数は、
        どの五角形でも隣の枚数より 1 少ない。外れれば NG
  検定3 渦：すべての五角形に光を 1 ずつ置き、左回り（反時計回り）に円環を回して乗り換えさせると、6030 すべてが輪に着く。
        どの光も、向きの区分（10 の向き）の変化を足した巻き数が 0 以上（時計回りに巻く光はない）。
        右回りにすると、巻き数の分布はちょうど符号を逆にしたものになる（鏡写し）。外れれば NG
  検定4 集める：輪に着く量を刻ごとに足すと、輪の量は増えるだけで、最後に 6030。5 つの入口から入る量はすべて等しい。
        負の対照：同じ近さの円環が二つあるとき、並びの先の方を選ぶと、入口の量が等しくならない
  検定5 貯める：輪に入った光は、輪の上を一刻に一歩、回り込んできた向きのまま回り続ける。量は増えも減りもしない。
        輪の五角形 0 と 1 の間を一刻に通る量は、左回りなら反時計回りだけ、右回りなら時計回りだけ（回る向きが記憶）。外れれば NG
"""
import sys
from collections import Counter, deque
sys.path.insert(0, sys.argv[1])
import b13_chain_units as B

T = (-2, 2, 0, 3); S0 = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U, seen = [], set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S0, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = set(B.fits(U))
STEP = [B.zrot(B.PHI, j) for j in range(10)]
nb = {c: [B.zadd(c, s) for s in STEP if B.zadd(c, s) in cells] for c in cells}
O = (0, 0, 0, 0)

# 向きの比較：w の虚部 ÷ sin36° を p + qφ で書く（sin(36°k)/sin36° = 0, 1, φ, φ, 1, 0, …）
SIN = [(0, 0), (1, 0), (0, 1), (0, 1)]
def im_sign(w):
    return B.phi_sign((sum(w[i] * SIN[i][0] for i in range(4)), sum(w[i] * SIN[i][1] for i in range(4))))
def turn(c, s):
    """中心から c を見て、一歩 s が反時計回り（+1）か時計回り（−1）か。conj(c)·s の虚部の符号。"""
    return im_sign(B.zmul(B.zconj(c), s))
# 向きの区分：2Re(c ζ^(−k)) を p + qφ で書き、最も大きい k を c の向きとする
COS2 = [(2, 0), (0, 1), (-1, 1), (1, -1), (0, -1), (-2, 0), (0, -1), (1, -1), (-1, 1), (0, 1)]
def proj2(c, k):
    return (sum(c[i] * COS2[(i - k) % 10][0] for i in range(4)), sum(c[i] * COS2[(i - k) % 10][1] for i in range(4)))
def sector(c):
    best = 0
    for k in range(1, 10):
        a, b = proj2(c, k), proj2(c, best)
        if B.phi_sign((a[0] - b[0], a[1] - b[1])) > 0: best = k
    return best

# 貯める輪：真ん中の五芒星に接する五角形 5 枚と、その外側の 5 枚
CEN = [B.zadd(O, B.zrot(B.PHI, j)) for j in range(10) if B.zadd(O, B.zrot(B.PHI, j)) in cells]
dist = {c: 0 for c in CEN}; q = deque(CEN)
while q:
    u = q.popleft()
    for v in nb[u]:
        if v not in dist: dist[v] = dist[u] + 1; q.append(v)
RING = set(CEN) | {c for c in cells if dist[c] == 1}
# 輪を反時計回りに並べる
cyc = [CEN[0]]
while len(cyc) < 10:
    nxt = [v for v in nb[cyc[-1]] if v in RING and v not in cyc]
    if not nxt: break
    cyc.append(nxt[0] if len(cyc) > 1 else max(nxt, key=lambda v: turn(cyc[0], B.zsub(v, cyc[0]))))
POS = {c: i for i, c in enumerate(cyc)}
dR = {c: 0 for c in RING}; q = deque(RING)
while q:
    u = q.popleft()
    for v in nb[u]:
        if v not in dR: dR[v] = dR[u] + 1; q.append(v)

# 円環：U[r] を中心に、五角形 c + φ²ζ^k（k = 0..9、k が増える向きが反時計回り）
from functools import cmp_to_key
N2 = [B.norm2(c) for c in U]
mem = {c: [] for c in cells}
for r, c in enumerate(U):
    for k in range(10): mem[B.zadd(c, B.zrot(B.PHI2, k))].append((r, k))

def ring_order(c, sense, naive=False):
    """c で乗る円環の順位：中心が真ん中に近い方が先。同じ近さなら、回る向きの先にある方が先。"""
    def cmp(x, y):
        a, b = N2[x[0]], N2[y[0]]
        if a != b: return -1 if B.phi_lt(a, b) else 1
        if naive: return 0                                                   # 負の対照：並びの先の方
        sx = im_sign(B.zmul(B.zconj(c), U[x[0]])) * sense
        sy = im_sign(B.zmul(B.zconj(c), U[y[0]])) * sense
        return -1 if sx > sy else (1 if sx < sy else 0)
    return cmp_to_key(cmp)

def ride(c, sense, naive=False):
    """c に置いた光の道（輪に着くまで）。"""
    r, k = min(mem[c], key=ring_order(c, sense, naive))
    path = [c]
    while c not in RING:
        k = (k + sense) % 10; c = B.zadd(U[r], B.zrot(B.PHI2, k)); path.append(c)
        best = min(mem[c], key=ring_order(c, sense, naive))
        if B.phi_lt(N2[best[0]], N2[r]): r, k = best
    return path

def winding(path):
    """向きの区分の変化（−4..+5 に直す）を足した巻き数（10 で一周）。"""
    w = 0
    for a, b in zip(path, path[1:]):
        d = (sector(b) - sector(a)) % 10
        w += d - 10 if d > 5 else d
    return w

if __name__ == "__main__":
    print("0. 地図：五角形", len(cells), "枚、円環", len(U), "個。真ん中の五芒星に接する五角形", len(CEN), "枚")

    print("検定1 貯める輪")
    closed = len(cyc) == 10 and cyc[0] in nb[cyc[-1]]
    gates = [(u, v) for u in cells if u not in RING for v in nb[u] if v in RING]
    L = Counter(dR.values()); layers = [L[i] for i in range(max(L) + 1)]
    print("   輪の五角形", len(RING), "枚、閉じた輪:", closed, "  輪の隣の数:", [len(nb[c]) for c in cyc])
    print("   輪の外から入る隣（入口）:", len(gates), "本  届く五角形:", len(dR), "/", len(cells))
    print("   輪からの段ごとの枚数:", layers)
    blk = [sum(layers[i:i + 5]) for i in range(0, len(layers) - 4, 5)]
    print("   5 段ずつの和:", blk, " 差:", [blk[i + 1] - blk[i] for i in range(len(blk) - 1)])
    print("   →", "OK" if closed and len(gates) == 5 and len(dR) == len(cells) and all(x % 5 == 0 for x in layers) else "NG")

    print("検定2 円環の道")
    closed_r = sum(1 for c in U if all(B.zadd(c, B.zrot(B.PHI2, (k + 1) % 10)) in nb[B.zadd(c, B.zrot(B.PHI2, k))] for k in range(10)))
    memb = Counter((len(mem[c]), len(nb[c])) for c in cells)
    print("   閉じた円環:", closed_r, "/", len(U), "  （属する円環の数, 隣の枚数）:", dict(sorted(memb.items())))
    print("   →", "OK" if closed_r == len(U) and all(m == n - 1 for m, n in memb) else "NG")

    print("検定3 渦")
    P = {s: {c: ride(c, s) for c in cells} for s in (1, -1)}
    WL = Counter(winding(p) for p in P[1].values()); WR = Counter(winding(p) for p in P[-1].values())
    lens = [len(p) - 1 for p in P[1].values()]
    longer = sum(1 for c in cells if len(P[1][c]) - 1 > dR[c])
    print("   輪に着いた光:", sum(1 for p in P[1].values() if p[-1] in RING), "/", len(cells),
          "  歩数の最大", max(lens), "  最短より長く回る光", longer, "個")
    print("   巻き数（10 で一周）の分布  左回り:", dict(sorted(WL.items())))
    print("   右回り:", dict(sorted(WR.items())))
    ok3 = all(p[-1] in RING for p in P[1].values()) and min(WL) >= 0 and all(WR[-w] == n for w, n in WL.items())
    print("   →", "OK" if ok3 else "NG")

    print("検定4 集める")
    arr = Counter(len(p) - 1 for p in P[1].values())
    cum, run = [], 0
    for t in range(max(arr) + 1): run += arr[t]; cum.append(run)
    gate_pos = sorted({POS[v] for u, v in gates})
    ent = Counter(POS[p[-1]] for c, p in P[1].items() if c not in RING)
    entn = Counter(POS[ride(c, 1, naive=True)[-1]] for c in cells if c not in RING)
    print("   輪の量（刻 0, 10, 20, …）:", cum[::10], "… 最後", cum[-1])
    print("   入口ごとに入った量:", [ent[i] for i in gate_pos], "  負の対照（並びの先の方）:", [entn[i] for i in gate_pos])
    ok4 = cum[-1] == len(cells) and all(cum[i] <= cum[i + 1] for i in range(len(cum) - 1)) and \
          len({ent[i] for i in gate_pos}) == 1 and len({entn[i] for i in gate_pos}) > 1
    print("   →", "OK" if ok4 else "NG")

    print("検定5 貯める")
    def run_ring(sense, ticks):
        ring = [0] * 10
        for c in RING: ring[POS[c]] += 1
        arrive = Counter((len(p) - 1, POS[p[-1]]) for c, p in P[sense].items() if c not in RING)
        snaps, flows = [], []
        for t in range(1, ticks):
            flows.append((ring[0] if sense == 1 else 0, ring[1] if sense == -1 else 0))
            r2 = [0] * 10
            for i in range(10): r2[(i + sense) % 10] = ring[i]
            for i in range(10): r2[i] += arrive[(t, i)]
            ring = r2; snaps.append(ring[:])
        return snaps, flows
    T = max(arr) + 20
    snl, fl = run_ring(1, T); snr, fr = run_ring(-1, T)
    full = max(arr) - 1
    keep = all(sum(x) == len(cells) for x in snl[full:]) and all(sum(x) == len(cells) for x in snr[full:])
    print("   すべて輪に入った刻:", full + 1, "  輪の量:", sum(snl[full]), "  その後、量が変わらない:", keep)
    print("   輪の並び  左回り:", snl[full], "  右回り:", snr[full])
    print("   輪の 0 と 1 の間を一刻に通る量（反時計回り, 時計回り）  左回り:", fl[full + 1], "  右回り:", fr[full + 1])
    ok5 = keep and fl[full + 1][1] == 0 and fl[full + 1][0] > 0 and fr[full + 1][0] == 0 and fr[full + 1][1] > 0
    print("   →", "OK" if ok5 else "NG")
