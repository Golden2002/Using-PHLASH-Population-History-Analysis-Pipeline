#!/usr/bin/env python3
import phlash
import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ===============================
# 用户需要修改的部分
# ===============================

# 群体名称 : 对应的 phlash 结果文件
POP_FILES = {
    "A":   "/home/litianxing/100My_Jino/114.PHLASH/results/A/models/phlash_results.pkl",
    "B":   "/home/litianxing/100My_Jino/114.PHLASH/results/B/models/phlash_results.pkl",
    "C":  "/home/litianxing/100My_Jino/114.PHLASH/results/C/models/phlash_results.pkl",
    "D":"/home/litianxing/100My_Jino/114.PHLASH/results/D/models/phlash_results.pkl",
}

# OUTFIG = "phlash_Ne_compare.png"

# 自动生成输出文件名
pop_tag = "_".join(POP_FILES.keys())
# TAG = "mu1.25e-8_gen25"
OUTFIG = f"phlash_Ne_compare__{pop_tag}.png"
Path("plots").mkdir(exist_ok=True)
OUTFIG = Path("plots") / f"phlash_Ne_compare__{pop_tag}.png"

N_T = 1000          # 时间点数
CI_LOW = 5          # credible interval
CI_HIGH = 95

# ===============================
# 读取所有群体的结果
# ===============================

pop_results = {}

for pop, f in POP_FILES.items():
    with open(f, "rb") as fh:
        results = pickle.load(fh)
    pop_results[pop] = results

# ===============================
# 构建统一时间轴 T
# ===============================

all_times = []

for results in pop_results.values():
    for dm in results:
        all_times.append(dm.eta.t[1:])  # 跳过 0

all_times = np.concatenate(all_times)

T = np.geomspace(all_times.min(), all_times.max(), N_T)

# ===============================
# 计算每个群体的 Ne posterior
# ===============================

pop_summary = {}

for pop, results in pop_results.items():
    Nes = np.array([dm.eta(T, Ne=True) for dm in results])

    pop_summary[pop] = {
        "median": np.median(Nes, axis=0),
        "low":    np.percentile(Nes, CI_LOW, axis=0),
        "high":   np.percentile(Nes, CI_HIGH, axis=0),
    }

# ===============================
# 开始绘图
# ===============================

plt.figure(figsize=(10, 7))

colors = plt.cm.tab10.colors  # 自动给不同群体上色

for i, (pop, stat) in enumerate(pop_summary.items()):
    c = colors[i % len(colors)]

    # credible interval
    plt.fill_between(
        T,
        stat["low"],
        stat["high"],
        color=c,
        alpha=0.25
    )

    # median line
    plt.plot(
        T,
        stat["median"],
        color=c,
        lw=2,
        label=pop
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Time (generations)")
plt.ylabel("Effective Population Size (Ne)")
plt.title("Phlash inference: Ne(t) comparison across populations")
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig(OUTFIG, dpi=300)
plt.close()

print(f"[✓] Figure saved to {OUTFIG}")
