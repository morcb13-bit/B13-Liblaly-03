"""正の対照：記録済みの平面（carrier_1245.json、五角形8520枚）に overlap を当て、重なり 0 を確かめる。
近いものを探すための升目分けだけ実数を使う（観察者側の索引）。判定は整数。"""
import json, sys, math, cmath
from ch01_join import overlap, dist2, sub
cells = {tuple(int(t) for t in k.split(",")): a for k, a in json.load(open(sys.argv[1]))["cells"].items()}
Z = [cmath.exp(1j * math.pi * k / 5) for k in range(4)]
grid = {}
for p in cells:
    z = sum(p[k] * Z[k] for k in range(4)); grid.setdefault((int(z.real // 2), int(z.imag // 2)), []).append(p)
pairs = hits = touch = 0
for (gx, gy), ps in grid.items():
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for p in ps:
                for q in grid.get((gx + dx, gy + dy), ()):
                    if p < q:
                        pairs += 1
                        if overlap(p, cells[p], q, cells[q]): hits += 1
                        if dist2(p, q) == (1, 1): touch += 1
print("五角形", len(cells), "枚。調べた近い組", pairs, "／ 重なり", hits, "／ 辺で接する組（距離 φ）", touch)
