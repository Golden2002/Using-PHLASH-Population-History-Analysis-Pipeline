#!/usr/bin/env python3

"""
简化版CCR分析脚本
功能：使用交叉合并率(CCR)方法分析两个群体的分歧时间
特点：
1. 只保留CCR方法，简化代码
2. 支持任意两个群体的分析，通过变量指定
3. 采用文献中的绘图风格
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import argparse
import json
import sys
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any

# ============================
# 设置文献绘图风格
# ============================
try:
    import scienceplots
    plt.style.use('science')
    print("Using scienceplots style")
except ImportError:
    print("scienceplots not installed, using default style")

# 设置字体和数学符号
mpl.rcParams.update({
    'font.size': 12,
    'mathtext.fontset': 'dejavusans',  # 使用系统字体，避免LaTeX问题
    'text.usetex': False,  # 禁用LaTeX
    'figure.dpi': 300,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif'],
})

def divergence_time_ccr(T: np.ndarray, ccr_curve: np.ndarray, threshold: float = 0.5) -> Tuple[Optional[float], float]:
    """
    使用交叉合并率(CCR)计算分歧时间
    CCR通常在群体分离时接近0，在群体合并时接近1
    分离时间通常定义为CCR下降到阈值以下的时间
    """
    # 找到CCR首次低于阈值的时间（从最古老到最近）
    idx = np.where(ccr_curve < threshold)[0]
    if len(idx) == 0:
        # 如果CCR始终高于阈值，尝试其他阈值
        for alt_threshold in [0.6, 0.7, 0.8, 1.0, 1.2]:
            idx = np.where(ccr_curve < alt_threshold)[0]
            if len(idx) > 0:
                print(f"使用替代阈值 {alt_threshold} 检测到分歧时间")
                return T[idx[-1]], alt_threshold
        return None, threshold
    return T[idx[-1]], threshold

def load_phlash_results(result_path: str) -> Any:
    """加载phlash结果文件"""
    with open(result_path, 'rb') as f:
        results = pickle.load(f)
    return results

# 在 divergence_time_ccr 函数之后添加以下函数（大约第 55 行后）

### 2025.12.05第二次修改
def calculate_threshold_range(T: np.ndarray, ccr_curve: np.ndarray,
                            base_threshold: float, delta: float = 0.05) -> Dict[str, Optional[float]]:
    """
    计算选定阈值附近±delta范围内的分歧时间范围

    参数:
    -----------
    T : np.ndarray
        时间点数组
    ccr_curve : np.ndarray
        CCR曲线值数组
    base_threshold : float
        基础阈值（例如0.5, 0.8等）
    delta : float
        阈值波动范围，默认±0.05

    返回:
    -----------
    time_range : dict
        包含以下键的字典:
        - 'lower_time': 阈值-delta对应的分歧时间
        - 'base_time': 基础阈值对应的分歧时间
        - 'upper_time': 阈值+delta对应的分歧时间
        - 'lower_threshold': 阈值-delta
        - 'upper_threshold': 阈值+delta
    """
    thresholds = {
        'lower_threshold': max(0.01, base_threshold - delta),  # 确保阈值不小于0.01
        'base_threshold': base_threshold,
        'upper_threshold': min(1.5, base_threshold + delta)    # 确保阈值不大于1.5
    }

    times = {}
    for key, threshold in thresholds.items():
        idx = np.where(ccr_curve < threshold)[0]
        if len(idx) == 0:
            times[key.replace('threshold', 'time')] = None
        else:
            times[key.replace('threshold', 'time')] = T[idx[-1]]


    return {
        'lower_threshold': thresholds['lower_threshold'],
        'upper_threshold': thresholds['upper_threshold'],
        'base_threshold': thresholds['base_threshold'],
        'lower_time': times['lower_time'],
        'upper_time': times['upper_time'],
        'base_time': times['base_time']
    }

#    return thresholds, times # 2025.12.05第二次修改提示：1.需要返回两个值 2.不能使用字典（可哈希）

def compute_ccr_curve(Ne1: np.ndarray, Ne2: np.ndarray,
                     Ne_combined: np.ndarray, T: np.ndarray) -> np.ndarray:
    """
    计算交叉合并率(CCR)曲线

    参数:
    Ne1, Ne2: 两个独立群体的有效种群大小曲线
    Ne_combined: 合并群体的有效种群大小曲线
    T: 时间点

    返回:
    ccr_curve: CCR曲线
    """
    # CCR公式: 2 * Ne_combined / (Ne1 + Ne2)
    return 2 * Ne_combined / (Ne1 + Ne2)

def analyze_pairwise_ccr(pop1_name: str, pop2_name: str, base_dir: str,
                         output_dir: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    分析两个群体的CCR并绘制图表

    参数:
    -----------
    pop1_name, pop2_name : str
        两个群体的名称
    base_dir : str
        结果文件的基础目录
    output_dir : str
        输出目录，默认为base_dir
    """
    # 检查是否存在合并群体数据
    merged_pop = f"{pop1_name}::{pop2_name}"
    merged_pop_alt = f"{pop2_name}::{pop1_name}"

     # 检查输入
    if not pop1_name or not pop2_name:
        print("错误: 必须指定两个群体名称")
        return None

    print(f"分析群体对: {pop1_name} vs {pop2_name}")

    # 定义可能的文件路径
    pop1_path = Path(base_dir) / pop1_name / "models" / "phlash_results.pkl"
    pop2_path = Path(base_dir) / pop2_name / "models" / "phlash_results.pkl"
    merged_path = Path(base_dir) / merged_pop / "models" / "phlash_results.pkl"
    merged_path_alt = Path(base_dir) / merged_pop_alt / "models" / "phlash_results.pkl"

    # 检查文件是否存在
    if not pop1_path.exists():
        print(f"错误: {pop1_path} 不存在")
        return None
    if not pop2_path.exists():
        print(f"错误: {pop2_path} 不存在")
        return None

    # 尝试查找合并群体文件
    if merged_path.exists():
        merged_pop_used = merged_pop
        merged_path_used = merged_path
        print(f"找到合并群体文件: {merged_path_used}")
    elif merged_path_alt.exists():
        merged_pop_used = merged_pop_alt
        merged_path_used = merged_path_alt
        print(f"找到合并群体文件: {merged_path_used}")
    else:
        print(f"警告: 合并群体文件不存在，无法计算CCR")
        print(f"  请确保已经运行了合并分析: {merged_pop} 或 {merged_pop_alt}")
        print(f"  请检查是否存在: {merged_path} 或 {merged_path_alt}")
        print(f"  当前基础目录: {base_dir}")

        # 列出基础目录下的所有群体目录
        print(f"  基础目录下的群体目录:")
        base_path = Path(base_dir)
        if base_path.exists():
            for item in base_path.iterdir():
                if item.is_dir():
                    print(f"    - {item.name}")

        return None

    # 加载数据
    print(f"加载数据...")
    try:
        pop1_results = load_phlash_results(pop1_path)
        print(f"  群体 {pop1_name}: 加载了 {len(pop1_results)} 个后验样本")
        pop2_results = load_phlash_results(pop2_path)
        print(f"  群体 {pop2_name}: 加载了 {len(pop2_results)} 个后验样本")
        merged_results = load_phlash_results(merged_path_used)
        print(f"  合并群体 {merged_pop_used}: 加载了 {len(merged_results)} 个后验样本")
    except Exception as e:
        print(f"加载数据时出错: {e}")
        return None

#    print(f"加载群体 {pop1_name} 数据...")
#    pop1_results = load_phlash_results(pop1_path)
#    print(f"加载群体 {pop2_name} 数据...")
#    pop2_results = load_phlash_results(pop2_path)
#    print(f"加载合并群体 {merged_pop_used} 数据...")
#    merged_results = load_phlash_results(merged_path_used)

    # 处理数据：生成时间点并计算Ne中位数
    def process_results(results, pop_name=""):
        try:
            times = np.array([dm.eta.t[1:] for dm in results])
            T_min = times.min()
            T_max = times.max()
            T = np.geomspace(T_min, T_max, 1000)
            Nes = np.array([dm.eta(T, Ne=True) for dm in results])
            ne_median = np.median(Nes, axis=0)
            print(f"  群体 {pop_name}: 时间范围 {T_min:.2e} 到 {T_max:.2e}, Ne范围 {ne_median.min():.2e} 到 {ne_median.max():.2e}")
            return T, ne_median, Nes
        except Exception as e:
            print(f"  处理群体 {pop_name} 数据时出错: {e}")
            raise
#        times = np.array([dm.eta.t[1:] for dm in results])
#        T_min = times.min()
#        T_max = times.max()
#        T = np.geomspace(T_min, T_max, 1000)
#        Nes = np.array([dm.eta(T, Ne=True) for dm in results])
#        ne_median = np.median(Nes, axis=0)
#        print(f"  群体 {pop_name}: 时间范围 {T_min:.2e} 到 {T_max:.2e}, Ne范围 {ne_median.min():.2e} 到 {ne_median.max():.2e}")
#        return T, ne_median, Nes

    T1, ne1_median, ne1_all = process_results(pop1_results, pop1_name)
    T2, ne2_median, ne2_all = process_results(pop2_results, pop2_name)
    T_merged, ne_merged_median, ne_merged_all = process_results(merged_results, merged_pop_used)

    # 找到共同的时间范围
    T_min = max(T1.min(), T2.min(), T_merged.min())
    T_max = min(T1.max(), T2.max(), T_merged.max())

    if T_min >= T_max:
        print(f"错误: 时间范围不重叠 (T_min={T_min:.2e}, T_max={T_max:.2e})")
        return None

    T_common = np.geomspace(T_min, T_max, 1000)
    ### 2025.12.05修改，START FILTER: 仅保留 ≤ 1e5 代的数据 ###
    mask = T_common <= 1e4
    if not np.any(mask):
        print("警告: 没有时间点 ≤ 1e5 代 — 检查输入 T_min, T_max")
        return None
    T_common = T_common[mask]
    print(f"共同时间范围 (截断后): {T_common.min():.2e} 到 {T_common.max():.2e}")
    print(f"共同时间范围: {T_min:.2e} 到 {T_max:.2e}")

    # 插值到共同时间点
    ne1_interp = np.interp(np.log(T_common), np.log(T1), np.log(ne1_median))
    ne2_interp = np.interp(np.log(T_common), np.log(T2), np.log(ne2_median))
    ne_merged_interp = np.interp(np.log(T_common), np.log(T_merged), np.log(ne_merged_median))

    # 计算CCR曲线
    ccr_curve = compute_ccr_curve(
        np.exp(ne1_interp),
        np.exp(ne2_interp),
        np.exp(ne_merged_interp),
        T_common
    )

    # 检查CCR曲线的合理性
    ccr_min = ccr_curve.min()
    ccr_max = ccr_curve.max()
    ccr_mean = ccr_curve.mean()

    print(f"\nCCR曲线统计:")
    print(f"  CCR最小值: {ccr_min:.3f}")
    print(f"  CCR最大值: {ccr_max:.3f}")
    print(f"  CCR平均值: {ccr_mean:.3f}")
    print(f"  CCR中位数: {np.median(ccr_curve):.3f}")

    # 计算分歧时间
    t_div_ccr, used_threshold = divergence_time_ccr(T_common, ccr_curve, threshold=0.5)
    ### 2025.12.05更新：增加了可控制时间范围的模块，这里是配套模块
    if t_div_ccr is not None and t_div_ccr > 1e5:
        print("检测到分歧时间超过 1e5 代，不在分析范围内，设为 None")
        t_div_ccr = None

    ### === 2025.12.05第二次更新：计算阈值波动范围内的分歧时间范围 ===
    if t_div_ccr is not None:
        # 计算阈值±0.05范围内的分歧时间
        threshold_range_times = calculate_threshold_range(T_common, ccr_curve,
                                                         used_threshold, delta=0.005)

        # 计算阈值±0.1范围内的分歧时间（可选）
        threshold_range_times_wide = calculate_threshold_range(T_common, ccr_curve,
                                                              used_threshold, delta=0.01)
    else:
        threshold_range_times = None
        threshold_range_times_wide = None
    ### === 新增结束 ===

    # 输出结果
    print(f"\n{'='*50}")
    print(f"CCR分析结果: {pop1_name} vs {pop2_name}")
    print(f"{'='*50}")
    if t_div_ccr:
        print(f"分歧时间 (CCR方法, 阈值={used_threshold}): {t_div_ccr:.2e} 代")
        print(f"分歧时间 (CCR方法, 阈值={used_threshold}): {t_div_ccr:.0f} 代")

        ### === 2025.12.05 第二次更新：输出阈值范围信息 ===
        if threshold_range_times:
            print(f"\n阈值波动分析 (±0.005):")
            print(f"  阈值 {threshold_range_times['lower_threshold']:.3f} 对应的分歧时间: ", end="")
            if threshold_range_times['lower_time']:
                print(f"{threshold_range_times['lower_time']:.3e} 代 ({threshold_range_times['lower_time']:.0f} 代)")
            else:
                print("未检测到")

            print(f"  阈值 {threshold_range_times['upper_threshold']:.3f} 对应的分歧时间: ", end="")
            if threshold_range_times['upper_time']:
                print(f"{threshold_range_times['upper_time']:.3e} 代 ({threshold_range_times['upper_time']:.0f} 代)")
            else:
                print("未检测到")

            # 计算时间范围（如果有上下界）
            valid_times = []
            if threshold_range_times['lower_time']:
                valid_times.append(threshold_range_times['lower_time'])
            if threshold_range_times['upper_time']:
                valid_times.append(threshold_range_times['upper_time'])

            if len(valid_times) >= 2:
                time_min = min(valid_times)
                time_max = max(valid_times)
                print(f"  分歧时间范围: {time_min:.2e} - {time_max:.2e} 代")
                print(f"  分歧时间范围: {time_min:.0f} - {time_max:.0f} 代")
                print(f"  时间跨度: {(time_max/time_min):.2f} 倍")

        # 可选：输出更宽的范围
        if threshold_range_times_wide:
            print(f"\n阈值波动分析 (±0.01):")
            if threshold_range_times_wide['lower_time']:
                print(f"  阈值 {threshold_range_times_wide['lower_threshold']:.2f} 对应: {threshold_range_times_wide['lower_time']:.0f} 代")
            if threshold_range_times_wide['upper_time']:
                print(f"  阈值 {threshold_range_times_wide['upper_threshold']:.2f} 对应: {threshold_range_times_wide['upper_time']:.0f} 代")
        ### === 新增结束 ===

    else:
        print(f"未检测到明显分歧时间 (尝试阈值最高至0.8，CCR始终高于阈值)")
        print("可能原因:")
        print("  1. 两个群体可能没有明显分离")
        print("  2. 分离时间可能非常古老")
        print("  3. 数据可能不足以检测分离事件")

    # ============================
    # 绘制图表
    # ============================
    try:
      fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10),
                                      gridspec_kw={'height_ratios': [2, 1]},
                                      constrained_layout=True)  # 使用constrained_layout替代layout="constrained"

      # 颜色设置
      colors = {'pop1': '#0C5DA5', 'pop2': '#00B945', 'ccr': '#FF9500'}

      # 图1: 两个群体的Ne曲线
      ax1.plot(T1, ne1_median, color=colors['pop1'], linewidth=1.5, label=pop1_name)
      ax1.plot(T2, ne2_median, color=colors['pop2'], linewidth=1.5, label=pop2_name)
      ax1.plot(T_merged, ne_merged_median, color='gray', linewidth=1.5,
               linestyle='--', alpha=0.7, label=f"{merged_pop_used} (together)")

      if t_div_ccr:
          ax1.axvline(x=t_div_ccr, color='red', linestyle='--',
                     linewidth=1.5, alpha=0.8, label=f'Divergence: {t_div_ccr:.2e} gens')

      ax1.set_xscale('log')
      ax1.set_yscale('log')
      ax1.set_xlabel("Time (generations)", fontsize=12)
      ax1.set_ylabel("Effective Population Size ($N_e$)", fontsize=12)
      ax1.set_title(f"{pop1_name} vs {pop2_name}: Population Size Trajectories", fontsize=14, pad=10)
      ax1.legend(loc='best', frameon=True, fancybox=True, framealpha=0.8)
      ax1.grid(True, alpha=0.2, which='both', linestyle=':')

      # 图2: CCR曲线
      ax2.plot(T_common, ccr_curve, color=colors['ccr'], linewidth=2, label='CCR')

      # 添加阈值线
      ax2.axhline(y=0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
      ax2.text(T_common[-1]*0.7, 0.52, 'CCR=0.5', fontsize=10, color='gray')

      if t_div_ccr:
          ax2.axvline(x=t_div_ccr, color='red', linestyle='--',
                     linewidth=1.5, alpha=0.8)
          ax2.text(t_div_ccr*1.1, 0.1, f'{t_div_ccr:.2e} gens',
                  rotation=90, verticalalignment='top',
                  color='red', fontsize=10, alpha=0.8)

      ax2.set_xscale('log')
      ax2.set_xlabel("Time (generations)", fontsize=12)
      ax2.set_ylabel(r"Cross Coalescence Rate (CCR)", fontsize=12)
      ax2.set_title(f"CCR Curve: {pop1_name} vs {pop2_name}", fontsize=14, pad=10)
      ax2.set_ylim(0, max(1.1, ccr_curve.max() * 1.1))
      ax2.legend(loc='upper right', frameon=True, fancybox=True, framealpha=0.8)
      ax2.grid(True, alpha=0.2, which='both', linestyle=':')

      # 设置输出路径
      if output_dir is None:
          output_dir = base_dir

      output_dir = Path(output_dir)
      ### 2025.12.06更新，目录控制功能
      # 创建以群体对命名的子目录
      subdir_name = f"{pop1_name}_vs_{pop2_name}"
      output_subdir = output_dir / subdir_name
      output_subdir.mkdir(parents=True, exist_ok=True)
      output_dir.mkdir(parents=True, exist_ok=True)

      output_path = output_dir / f"ccr_analysis_{pop1_name}_vs_{pop2_name}.png"
#      output_path = output_subdir / f"ccr_analysis_{pop1_name}_vs_{pop2_name}.png" # 如果图片也想保存到子目录的话

      # 保存图片
      fig.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.1)
      plt.close(fig)

      print(f"\n图表已保存至: {output_path}")
    except Exception as e:
        print(f"绘图时出错: {e}")
        print("尝试简化绘图...")

        # 简化绘图
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(T_common, ccr_curve, color='blue', linewidth=2, label='CCR')
        ax.axhline(y=0.5, color='red', linestyle='--', linewidth=1, label='threshold=0.5')

        if t_div_ccr:
            ax.axvline(x=t_div_ccr, color='green', linestyle='--', linewidth=1, label=f'Split Time: {t_div_ccr:.2e}')

        ax.set_xscale('log')
        ax.set_xlabel("Time (generations)")
        ax.set_ylabel("CCR")
        ax.set_title(f"CCR Curve: {pop1_name} vs {pop2_name}")
        ax.legend()
        ax.grid(True, alpha=0.3)

#        output_dir = Path("/home/litianxing/100My_Jino/114.PHLASH")
        output_dir = Path(output_dir)
        ### 2025.12.06 更新的配套部分
        # 创建以群体对命名的子目录
        subdir_name = f"{pop1_name}_vs_{pop2_name}"
        output_subdir = output_dir / subdir_name
        output_subdir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"ccr_analysis_{pop1_name}_vs_{pop2_name}_simple.png"
#        output_path = output_subdir / f"ccr_analysis_{pop1_name}_vs_{pop2_name}_simple.png" # 如果图片也想保存到子目录的话
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        print(f"简化图表已保存至: {output_path}")

    print(f"{'='*50}")

#    # 返回结果
#    return {
#        'pop1': pop1_name,
#        'pop2': pop2_name,
#        'divergence_time': t_div_ccr,
#        'threshold_used': used_threshold,
#        'T_common': T_common,
#        'ccr_curve': ccr_curve,
#        'ne1': np.exp(ne1_interp),
#        'ne2': np.exp(ne2_interp),
#        'ne_merged': np.exp(ne_merged_interp),
#        'ccr_min': ccr_min,  # 缺少这个
#        'ccr_max': ccr_max,  # 缺少这个
#        'output_path': output_path  # 缺少这个
#    }
    ### 2025.12.05第二次修改，配套部分 返回结果
    result = {
        'pop1': pop1_name,
        'pop2': pop2_name,
        'divergence_time': t_div_ccr,
        'threshold_used': used_threshold,
        'T_common': T_common,
        'ccr_curve': ccr_curve,
        'ne1': np.exp(ne1_interp),
        'ne2': np.exp(ne2_interp),
        'ne_merged': np.exp(ne_merged_interp),
        'ccr_min': ccr_min,
        'ccr_max': ccr_max,
        'output_path': str(output_path) if output_path else None,
        'subdir_path': str(output_subdir) # 2025.12.06 更新，如果需要返回子目录路径，可以添加：
    }

    # === 新增：添加阈值范围信息到返回结果 ===
    if threshold_range_times:
        result['threshold_range_±0.05'] = threshold_range_times

    if threshold_range_times_wide:
        result['threshold_range_±0.1'] = threshold_range_times_wide

    # 如果有时间范围，计算平均时间
    if t_div_ccr and threshold_range_times:
        valid_times = []
        if threshold_range_times['lower_time']:
            valid_times.append(threshold_range_times['lower_time'])
        if threshold_range_times['upper_time']:
            valid_times.append(threshold_range_times['upper_time'])

        if valid_times:
            result['time_range_min'] = min(valid_times)
            result['time_range_max'] = max(valid_times)
            result['time_range_avg'] = np.mean(valid_times)

    return result

def batch_analyze_ccr(pop_pairs: List[Tuple[str, str]],
                     base_dir: str, output_dir: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
    """
    批量分析多个群体对的CCR

    参数:
    -----------
    pop_pairs : list of tuples
        群体对列表，例如 [("Jino", "Han_N"), ("Han_N", "Tibetan")]
    base_dir : str
        结果文件的基础目录
    output_dir : str
        输出目录
    """
    results = {}

    for pop1, pop2 in pop_pairs:
        print(f"\n{'='*60}")
        print(f"开始分析: {pop1} vs {pop2}")
        print(f"{'='*60}")

        result = analyze_pairwise_ccr(pop1, pop2, base_dir, output_dir)
        if result is not None:
            key = f"{pop1}_vs_{pop2}"
            results[key] = result

    # 输出汇总信息
    if results:
        print(f"\n{'='*60}")
        print("CCR分析汇总")
        print(f"{'='*60}")

        for key, result in results.items():
            if result['divergence_time']:
                print(f"{key}: {result['divergence_time']:.2e} 代 (阈值={result.get('threshold_used', 0.5)})")
            else:
                print(f"{key}: 未检测到明显分歧时间")

    return results

def save_results_to_json(results: Dict[str, Dict[str, Any]],
                        output_file: str = "ccr_analysis_results.json") -> None:
    """
    将分析结果保存为JSON文件

    参数:
    -----------
    results : dict
        analyze_ccr或batch_analyze_ccr返回的结果字典
    output_file : str
        输出JSON文件路径
    """
    try:
        # 转换numpy数组为列表以便JSON序列化
        def convert_for_json(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.generic):
                return obj.item()
            elif isinstance(obj, Path):  # 2025.12.05第二次新增：处理Path对象
                return str(obj)
            else:
                return obj

        # 深拷贝并转换结果
        json_results = {}
        for key, value in results.items():
            if isinstance(value, dict):
                json_results[key] = {k: convert_for_json(v) for k, v in value.items()}
            else:
                json_results[key] = convert_for_json(value)

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(json_results, f, indent=2, ensure_ascii=False, default=str)

        print(f"\n结果已保存至: {output_file}")
    except Exception as e:
        print(f"保存结果到JSON文件时出错: {e}")

def main():
    # 使用argparse支持命令行参数
    parser = argparse.ArgumentParser(
        description='CCR分析: 使用交叉合并率计算群体分歧时间',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 分析单个群体对
  python ccr_analysis.py --pop1 Jino --pop2 Han_N

  # 分析单个群体对并指定输出目录
  python ccr_analysis.py --pop1 Jino --pop2 Han_N --output_dir ./results

  # 批量分析多个群体对
  python ccr_analysis.py --batch --pop_pairs "Jino:Han_N" "Han_N:Tibetan" "Jino:Tibetan"

  # 显示详细帮助
  python ccr_analysis.py -h
        """
    )

    parser.add_argument('--pop1', type=str, default="Jino", help='第一个群体名称')
    parser.add_argument('--pop2', type=str, default="Han_N", help='第二个群体名称')
    parser.add_argument('--base_dir', type=str,
                       default="/home/litianxing/100My_Jino/114.PHLASH/results",
                       help='结果文件的基础目录')
    parser.add_argument('--output_dir', type=str,
                       help='输出目录，默认为基础目录')
    parser.add_argument('--output_json', type=str,
                       default="ccr_analysis_results.json",
                       help='结果JSON文件路径 (默认: ccr_analysis_results.json)')
    parser.add_argument('--batch', action='store_true',
                       help='批量分析模式，分析多个群体对')
    parser.add_argument('--pop_pairs', type=str, nargs='+',
                       default=["Jino:Han_N", "Han_N:Tibetan", "Jino:Tibetan"],
                       help='批量分析的群体对，格式为"群体1:群体2"，多个用空格分隔')
    parser.add_argument('--config', type=str,
                       help='配置文件路径 (JSON格式)，包含群体对列表和其他参数')

    args = parser.parse_args()

    # 处理配置文件
    if args.config:
        try:
            with open(args.config, 'r', encoding='utf-8') as f:
                config = json.load(f)

            # 从配置文件中读取参数
            if 'pop_pairs' in config:
                args.pop_pairs = config['pop_pairs']
            if 'base_dir' in config:
                args.base_dir = config['base_dir']
            if 'output_dir' in config:
                args.output_dir = config['output_dir']
            if 'output_json' in config:
                args.output_json = config['output_json']

            print(f"已从配置文件加载参数: {args.config}")
        except Exception as e:
            print(f"加载配置文件时出错: {e}")
            return 1

    # 设置输出目录
    if args.output_dir is None:
        args.output_dir = args.base_dir

    if args.batch:
        # 批量分析模式
        pop_pairs = []
        for pair_str in args.pop_pairs:
            if ':' in pair_str:
                pop1, pop2 = pair_str.split(':')
                pop_pairs.append((pop1.strip(), pop2.strip()))
            else:
                print(f"警告: 忽略格式错误的群体对: {pair_str}")

        if pop_pairs:
            print(f"批量分析以下群体对: {pop_pairs}")
            results = batch_analyze_ccr(pop_pairs, args.base_dir, args.output_dir)

            # 保存结果到JSON文件
            if results:
#                save_results_to_json(results, args.output_json)

                ### 2025.12.06 新修改
                # 为每个结果创建对应的子目录保存JSON
                for key, result in results.items():
                    if 'subdir_path' in result:
                        json_path = Path(result['subdir_path']) / args.output_json
                    else:
                        # 如果没有subdir_path，则创建
                        pop1, pop2 = result['pop1'], result['pop2']
                        subdir = Path(args.output_dir) / f"{pop1}_vs_{pop2}"
                        subdir.mkdir(parents=True, exist_ok=True)
                        json_path = subdir / args.output_json
                    save_results_to_json({key: result}, str(json_path))

                # 同时保存一个汇总文件到输出目录
                total_json_path = Path(args.output_dir) / f"summary_{args.output_json}"
                save_results_to_json(results, str(total_json_path))
                print(f"\n批量分析完成!")
                print(f"共分析了 {len(results)} 个群体对")
                print(f"详细结果已保存至: {args.output_json}")
            else:
                print(f"\n批量分析完成，但未获得有效结果")
                return 1
        else:
            print("未提供有效的群体对")
            return 1
    else:
        # 单对分析模式
        print(f"分析群体对: {args.pop1} vs {args.pop2}")
        print("参数设置:")
        print(f"  群体1: {args.pop1}")
        print(f"  群体2: {args.pop2}")
        print(f"  基础目录: {args.base_dir}")
        print(f"  目录是否存在: {os.path.exists(args.base_dir)}")
#        print(f"分析群体对: {args.pop1} vs {args.pop2}")
        result = analyze_pairwise_ccr(args.pop1, args.pop2, args.base_dir, args.output_dir)

        # 示例: 如何访问结果
        if result:
            print(f"\n分析完成!")
            if result['divergence_time']:
                print(f"建议的分歧时间: {result['divergence_time']:.2e} 代 (阈值={result['threshold_used']})")
            else:
                print("未检测到明显分歧时间")

            # 保存单对结果到JSON
#            save_results_to_json({f"{args.pop1}_vs_{args.pop2}": result}, args.output_json)

            ### 2025.12.06 新修改，JSON文件保存到指定目录
            # 使用subdir_path保存JSON

            if 'subdir_path' in result:
                json_path = Path(result['subdir_path']) / args.output_json
            else:
                # 如果没有subdir_path，则创建
                subdir = Path(args.output_dir) / f"{args.pop1}_vs_{args.pop2}"
                subdir.mkdir(parents=True, exist_ok=True)
                json_path = subdir / args.output_json

            save_results_to_json({f"{args.pop1}_vs_{args.pop2}": result}, str(json_path))
            return 0
        else:
            print(f"\n分析失败")
            return 1

    return 0

if __name__ == "__main__":
    exit(main())
