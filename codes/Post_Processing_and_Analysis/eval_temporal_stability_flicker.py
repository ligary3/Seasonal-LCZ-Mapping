# -*- coding: utf-8 -*-
"""
LCZ Evolution (2017–2024) — Windowed processing with persistent cache
- Block-wise statistics (no MemoryError)
- Exports reusable data tables (CSV/NPZ/JSON) for later plotting
- Sankey (aggregate & full-17 with LCZ colors), Transition heatmaps, Area trajectories

Usage:
  1) 首次运行: USE_CACHE=False （默认），生成全部统计数据与图；
  2) 之后:     USE_CACHE=True   直接读取缓存数据，快速重绘。

Author: you
"""

import os, json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from rasterio.vrt import WarpedVRT
# =============== 0) 基础参数（按需修改） ===============
ROOT = r"E:\luwen2\result\年际LCZ\纯净版"      # 2017_ref.tif … 2024_ref.tif 的文件夹
YEARS = list(range(2017, 2025))       # 2017–2024
REP_YEARS = (2017, 2021, 2024)        # 代表年：启动–扩张–稳定
OUTDIR = os.path.join(ROOT, "outputs")
CACHE = os.path.join(OUTDIR, "cache")
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

USE_CACHE = False     # ✅ 想直接复用数据改 True

# Matplotlib 风格（Times / 细线 / 审稿友好）
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.linewidth": 1.0,
    "figure.dpi": 300,
})

HEATMAP_CMAP = "viridis"   # 按你要求：热力图先不变

# LCZ 17类标准配色（按你给的顺序）
LCZ_COLORS = [
    '#8c0000', '#d10000', '#ff0000',
    '#bf4d00', '#ff6600', '#ff9955',
    '#faee05', '#bcbcbc', '#ffccaa',
    '#555555',
    '#006a00', '#00aa00', '#648525',
    '#b9db79', '#000000', '#fbf7ae', '#6a6aff'
]
LCZ_LABELS = [f"LCZ {i}" for i in range(1, 18)]

# 建筑(1–10) vs 自然(11–17)
BUILT_RANGE  = range(1, 11)
NATUR_RANGE  = range(11, 18)

# 若想强制面积（m²），例如 10m → 100*100=100（注意单位是 m²）
PIXEL_AREA_M2_OVERRIDE = None

# 可调：读取块大小（像素）。None 表示优先用数据内置块配置。
READ_BLOCK_SHAPE = None  # 例如 (2048, 2048)


# =============== 1) IO & 基础工具 ===============
def tif_path(year):
    return os.path.join(ROOT, f"{year}_ref.tif")

def get_pixel_area_m2(ds):
    """返回像元面积（m²），若是经纬度坐标系则返回 None。"""
    if PIXEL_AREA_M2_OVERRIDE is not None:
        return float(PIXEL_AREA_M2_OVERRIDE)
    crs = ds.crs
    if crs and crs.is_projected:
        return abs(ds.transform.a * ds.transform.e)  # a=像元宽, e=像元高(负)
    return None

def window_iter(ds):
    """按块遍历窗口。先尝试数据内置块；否则用 READ_BLOCK_SHAPE。"""
    try:
        for _, w in ds.block_windows(1):
            yield w
        return
    except Exception:
        pass
    h, w = ds.height, ds.width
    bh = READ_BLOCK_SHAPE[0] if READ_BLOCK_SHAPE else 2048
    bw = READ_BLOCK_SHAPE[1] if READ_BLOCK_SHAPE else 2048
    for row_off in range(0, h, bh):
        rows = min(bh, h - row_off)
        for col_off in range(0, w, bw):
            cols = min(bw, w - col_off)
            yield Window(col_off=col_off, row_off=row_off, width=cols, height=rows)

def histogram_1to17(arr_block, mask_block=None):
    if mask_block is None:
        valid = (arr_block >= 1) & (arr_block <= 17)
    else:
        valid = mask_block & (arr_block >= 1) & (arr_block <= 17)
    if not np.any(valid):
        return np.zeros(18, dtype=np.int64)
    vals = arr_block[valid].ravel()
    hist = np.bincount(vals, minlength=18).astype(np.int64)
    return hist

def transition_17x17(src_block, dst_block, mask_block=None):
    if mask_block is None:
        valid = (src_block >= 1) & (src_block <= 17) & (dst_block >= 1) & (dst_block <= 17)
    else:
        valid = mask_block & (src_block >= 1) & (src_block <= 17) & (dst_block >= 1) & (dst_block <= 17)
    if not np.any(valid):
        return np.zeros((17, 17), dtype=np.int64)
    s = src_block[valid].ravel().astype(np.int16) - 1
    d = dst_block[valid].ravel().astype(np.int16) - 1
    k = (s * 17 + d).astype(np.int32)
    bc = np.bincount(k, minlength=17*17).astype(np.int64)
    return bc.reshape((17, 17))

def save_matrix_csv(mat, y1, y2, outdir, px_area_m2=None):
    M = mat.astype(np.float64)
    unit = "pixels"
    if px_area_m2 is not None:
        M = M * (px_area_m2 / 1e6)   # km²
        unit = "km2"
    df = pd.DataFrame(M, index=[f"LCZ{i}" for i in range(1,18)],
                         columns=[f"LCZ{i}" for i in range(1,18)])
    fp = os.path.join(outdir, f"transition_{y1}_{y2}_{unit}.csv")
    df.to_csv(fp, encoding="utf-8-sig")
    return fp, unit

def plot_matrix_heatmap(mat, y1, y2, norm_by="row", cmap=HEATMAP_CMAP, outdir=".", px_area_m2=None):
    M = mat.astype(np.float64)
    title_suffix = "Counts"
    value_label = "Pixel count"
    if px_area_m2 is not None:
        M = M * (px_area_m2 / 1e6)  # km²
        title_suffix = "Area (km²)"
        value_label = "Area (km²)"
    if norm_by == "row":
        rs = M.sum(axis=1, keepdims=True); rs[rs==0]=1
        M = M / rs
        title_suffix += " — Row-normalized"; value_label = "Proportion"
    elif norm_by == "global":
        s = M.sum();  M = M/s if s>0 else M
        title_suffix += " — Global-normalized"; value_label = "Proportion"
    fig = plt.figure(figsize=(4.8, 4.2))
    ax = plt.gca()
    im = ax.imshow(M, cmap=cmap, aspect='auto', origin='upper')
    ax.set_title(f"LCZ Transition {y1}→{y2} ({title_suffix})", pad=8)
    ax.set_xticks(np.arange(17)); ax.set_yticks(np.arange(17))
    ax.set_xticklabels([str(i) for i in range(1,18)])
    ax.set_yticklabels([str(i) for i in range(1,18)])
    ax.set_xlabel("To"); ax.set_ylabel("From")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(value_label); cbar.ax.tick_params(length=3)
    plt.tight_layout()
    fp = os.path.join(outdir, f"transition_heatmap_{y1}_{y2}_{('counts' if norm_by is None else norm_by)}.png")
    plt.savefig(fp, dpi=300, bbox_inches="tight"); plt.close(fig)
    return fp


# =============== 2) 桑基：17 类着色（按 LCZ_COLORS） ===============
def hex_to_rgba(hex_color, alpha=0.35):
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"

def sankey_from_mats(mat_AB, mat_BC, yA, yB, yC, 
                     aggregate=True, out_html=True, outdir=".", px_area_m2=None,
                     colorize_full17=True):
    scale = 1.0
    value_suffix = " (pixels)"
    if px_area_m2 is not None:
        scale = px_area_m2 / 1e6  # km²
        value_suffix = " (km²)"

    if aggregate:
        # 2大类聚合
        b2b = mat_AB[0:10, 0:10].sum(); b2n = mat_AB[0:10,10:17].sum()
        n2b = mat_AB[10:17,0:10].sum(); n2n = mat_AB[10:17,10:17].sum()
        B_b = mat_BC[0:10, 0:10].sum(); B_b2n = mat_BC[0:10,10:17].sum()
        B_n2b = mat_BC[10:17,0:10].sum(); B_n = mat_BC[10:17,10:17].sum()

        vals = [b2b,b2n,n2b,n2n,B_b,B_b2n,B_n2b,B_n]
        values = [v*scale for v in vals]
        labels = [f"{yA}:Built", f"{yA}:Natural", f"{yB}:Built", f"{yB}:Natural", f"{yC}:Built", f"{yC}:Natural"]
        A_b_i, A_n_i, B_b_i, B_n_i, C_b_i, C_n_i = range(6)
        sources = [A_b_i, A_b_i, A_n_i, A_n_i, B_b_i, B_b_i, B_n_i, B_n_i]
        targets = [B_b_i, B_n_i, B_b_i, B_n_i, C_b_i, C_n_i, C_b_i, C_n_i]
        node_colors = None
        link_colors = None
    else:
        # 17类完整版
        def expand(mat, ofs_src, ofs_dst):
            s, t, v = [], [], []
            for i in range(17):
                for j in range(17):
                    val = int(mat[i,j])
                    if val>0:
                        s.append(ofs_src+i); t.append(ofs_dst+j); v.append(val*scale)
            return s,t,v
        labels = [f"{yA}:LCZ {i}" for i in range(1,18)] + \
                 [f"{yB}:LCZ {i}" for i in range(1,18)] + \
                 [f"{yC}:LCZ {i}" for i in range(1,18)]
        s1,t1,v1 = expand(mat_AB, 0, 17)
        s2,t2,v2 = expand(mat_BC, 17, 34)
        sources, targets, values = s1+s2, t1+t2, v1+v2

        if colorize_full17:
            node_colors = (LCZ_COLORS * 3)[:len(labels)]
            link_colors = [hex_to_rgba(node_colors[src], 0.35) for src in sources]
        else:
            node_colors = None; link_colors = None

    fig = go.Figure(go.Sankey(
        node=dict(label=labels, pad=12, thickness=16,
                  line=dict(color="rgba(0,0,0,0.35)", width=0.5),
                  color=node_colors),
        link=dict(source=sources, target=targets, value=values, color=link_colors)
    ))
    fig.update_layout(
        title=dict(text=f"LCZ Transitions: {yA}→{yB}→{yC}{value_suffix}", x=0.02, y=0.95),
        font=dict(family="Times New Roman", size=12),
        margin=dict(l=10, r=10, t=40, b=10),
        width=1000, height=500
    )
    tag = 'agg' if aggregate else 'full17'
    unit = 'km2' if (px_area_m2 is not None) else 'pixels'
    os.makedirs(outdir, exist_ok=True)
    if out_html:
        fp = os.path.join(outdir, f"sankey_{yA}_{yB}_{yC}_{tag}_{unit}.html")
        fig.write_html(fp, include_plotlyjs='cdn')
    else:
        try:
            import kaleido  # noqa
            fp = os.path.join(outdir, f"sankey_{yA}_{yB}_{yC}_{tag}_{unit}.png")
            fig.write_image(fp, scale=2)
        except Exception:
            fp = os.path.join(outdir, f"sankey_{yA}_{yB}_{yC}_{tag}_{unit}.html")
            fig.write_html(fp, include_plotlyjs='cdn')
    return fp, (labels, sources, targets, values)


# =============== 3) 统计（带缓存） ===============
def compute_or_load_stats(use_cache=USE_CACHE):
    (yA, yB, yC) = REP_YEARS

    # 缓存路径
    areas_csv = os.path.join(CACHE, "areas_17classes.csv")
    meta_json = os.path.join(CACHE, "meta.json")
    trans_AB_pixels = os.path.join(CACHE, f"transition_{yA}_{yB}_pixels.csv")
    trans_BC_pixels = os.path.join(CACHE, f"transition_{yB}_{yC}_pixels.csv")
    trans_AB_km2    = os.path.join(CACHE, f"transition_{yA}_{yB}_km2.csv")
    trans_BC_km2    = os.path.join(CACHE, f"transition_{yB}_{yC}_km2.csv")

    # -------- 1) 直接用缓存 --------
    if use_cache and all(os.path.exists(p) for p in [areas_csv, meta_json, trans_AB_pixels, trans_BC_pixels]):
        areas_df = pd.read_csv(areas_csv, index_col=0)
        with open(meta_json, "r", encoding="utf-8") as f:
            meta = json.load(f)
        mat_AB = pd.read_csv(trans_AB_pixels, index_col=0).to_numpy()
        mat_BC = pd.read_csv(trans_BC_pixels, index_col=0).to_numpy()
        px_area_m2_use = meta.get("pixel_area_m2", None)
        return areas_df, mat_AB, mat_BC, px_area_m2_use, (areas_csv, trans_AB_pixels, trans_BC_pixels, meta_json)

    # -------- 2) 分块统计年度面积（不要求同网格） --------
    area_hist = {y: np.zeros(18, dtype=np.int64) for y in YEARS}
    px_area_m2_use = None

    for y in YEARS:
        with rasterio.open(tif_path(y)) as ds:
            if px_area_m2_use is None:
                px_area_m2_use = get_pixel_area_m2(ds)  # EPSG:3857@10m -> ~100
            for w in window_iter(ds):
                arr = ds.read(1, window=w, out_dtype='int32', masked=False)
                try:
                    m = ds.read_masks(1, window=w) > 0
                except Exception:
                    m = None
                area_hist[y] += histogram_1to17(arr, m)

    # -------- 3) 两段转移矩阵（把 B、C 对齐到 A 的参考网格） --------
    mat_AB = np.zeros((17,17), dtype=np.int64)
    mat_BC = np.zeros((17,17), dtype=np.int64)

    with rasterio.open(tif_path(yA)) as dA, \
         rasterio.open(tif_path(yB)) as dB_src, \
         rasterio.open(tif_path(yC)) as dC_src:

        # 以 A 为参考网格
        ref_crs = dA.crs
        ref_transform = dA.transform
        ref_width, ref_height = dA.width, dA.height

        # 在内存构造对齐到 A 网格的 VRT（最近邻，不改变类别值）
        vrt_opts = dict(
            crs=ref_crs,
            transform=ref_transform,
            width=ref_width,
            height=ref_height,
            resampling=rasterio.enums.Resampling.nearest,
            dtype='int16',             # 分类整型
            nodata=0,                  # 保持我们代码中的无效=0 约定
            add_alpha=False,
            src_nodata=0
        )
        with WarpedVRT(dB_src, **vrt_opts) as dB, \
             WarpedVRT(dC_src, **vrt_opts) as dC:

            # A→B：按 A 的块窗口遍历，并从对齐后的 B 读取相同窗口
            for w in window_iter(dA):
                a = dA.read(1, window=w, out_dtype='int32', masked=False)
                b = dB.read(1, window=w, out_dtype='int32', masked=False)

                # 掩膜：结合有效值与 mask band（若有）
                try:
                    ma = dA.read_masks(1, window=w) > 0
                except Exception:
                    ma = None
                try:
                    mb = dB.read_masks(1, window=w) > 0
                except Exception:
                    mb = None
                if (ma is not None) and (mb is not None):
                    m = ma & mb
                else:
                    m = None
                mat_AB += transition_17x17(a, b, m)

            # B→C：以对齐后的 B 为“参考”，窗口与 A 一致（因为 B 已对齐到 A）
            for w in window_iter(dB):
                b = dB.read(1, window=w, out_dtype='int32', masked=False)
                c = dC.read(1, window=w, out_dtype='int32', masked=False)
                try:
                    mb = dB.read_masks(1, window=w) > 0
                except Exception:
                    mb = None
                try:
                    mc = dC.read_masks(1, window=w) > 0
                except Exception:
                    mc = None
                if (mb is not None) and (mc is not None):
                    m = mb & mc
                else:
                    m = None
                mat_BC += transition_17x17(b, c, m)

    # -------- 4) 写缓存 --------
    areas_df = pd.DataFrame(
        {f"LCZ{i}":[int(area_hist[y][i]) for y in YEARS] for i in range(1,18)},
        index=YEARS
    )
    areas_df.to_csv(areas_csv, encoding="utf-8-sig")

    meta = {
        "pixel_area_m2": px_area_m2_use,
        "unit_default": "km2" if px_area_m2_use is not None else "pixels",
        "years": YEARS,
        "rep_years": REP_YEARS,
        "align_to": f"{yA} (WarpedVRT, nearest)",
        "notes": "B & C regridded to A via WarpedVRT (nearest) to ensure pixel-wise transitions."
    }
    with open(meta_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    # 保存像元&面积两种单位的转移矩阵
    pd.DataFrame(mat_AB, index=[f"LCZ{i}" for i in range(1,18)],
                        columns=[f"LCZ{i}" for i in range(1,18)]).to_csv(trans_AB_pixels, encoding="utf-8-sig")
    pd.DataFrame(mat_BC, index=[f"LCZ{i}" for i in range(1,18)],
                        columns=[f"LCZ{i}" for i in range(1,18)]).to_csv(trans_BC_pixels, encoding="utf-8-sig")
    if px_area_m2_use is not None:
        pd.DataFrame(mat_AB * (px_area_m2_use/1e6), index=[f"LCZ{i}" for i in range(1,18)],
                                                   columns=[f"LCZ{i}" for i in range(1,18)]).to_csv(trans_AB_km2, encoding="utf-8-sig")
        pd.DataFrame(mat_BC * (px_area_m2_use/1e6), index=[f"LCZ{i}" for i in range(1,18)],
                                                   columns=[f"LCZ{i}" for i in range(1,18)]).to_csv(trans_BC_km2, encoding="utf-8-sig")

    return areas_df, mat_AB, mat_BC, px_area_m2_use, (areas_csv, trans_AB_pixels, trans_BC_pixels, meta_json)


# =============== 4) 桑基/折线导出成数据表（可复用） ===============
def export_sankey_edges(labels, sources, targets, values, out_csv):
    """把桑基的节点与边导出为 CSV，便于后续直接作图。"""
    nodes = pd.DataFrame({"node_id": list(range(len(labels))), "label": labels})
    edges = pd.DataFrame({"source": sources, "target": targets, "value": values})
    # 合并便于阅读
    edges["source_label"] = edges["source"].map(nodes.set_index("node_id")["label"])
    edges["target_label"] = edges["target"].map(nodes.set_index("node_id")["label"])
    edges = edges[["source","target","value","source_label","target_label"]]
    edges.to_csv(out_csv, index=False, encoding="utf-8-sig")
    return out_csv

def plot_area_trajectory_from_df(areas_df, outdir=OUTDIR, px_area_m2=None):
    """areas_df: 行=year, 列=LCZ1..LCZ17（像元数）"""
    years = areas_df.index.astype(int).tolist()
    built = areas_df[[f"LCZ{i}" for i in BUILT_RANGE]].sum(axis=1).values
    natu  = areas_df[[f"LCZ{i}" for i in NATUR_RANGE]].sum(axis=1).values
    ylabel = "Pixel count"
    if px_area_m2 is not None:
        built = built * (px_area_m2/1e6)
        natu  = natu  * (px_area_m2/1e6)
        ylabel = "Area (km²)"
    fig = plt.figure(figsize=(4.8, 3.6))
    ax = plt.gca()
    ax.plot(years, built, marker='o', linewidth=1.6, label='Built (LCZ 1–10)')
    ax.plot(years, natu,  marker='o', linewidth=1.6, label='Natural (LCZ 11–17)')
    ax.set_xlabel("Year"); ax.set_ylabel(ylabel); ax.set_title("Built vs. Natural Area Trajectories", pad=8)
    ax.grid(linestyle='--', linewidth=0.6, alpha=0.6); ax.legend(frameon=False)
    plt.tight_layout()
    fp = os.path.join(outdir, f"area_trajectory_BvN_{'km2' if px_area_m2 else 'pixels'}.png")
    plt.savefig(fp, dpi=300, bbox_inches="tight"); plt.close(fig)
    return fp


# =============== 5) 主流程 ===============
areas_df, mat_AB, mat_BC, px_area_m2_use, cache_files = compute_or_load_stats(USE_CACHE)
yA, yB, yC = REP_YEARS

# —— 转移矩阵热力图（保留三种版本，任选其一入文）——
plot_matrix_heatmap(mat_AB, yA, yB, norm_by=None,      outdir=OUTDIR, px_area_m2=px_area_m2_use)
plot_matrix_heatmap(mat_AB, yA, yB, norm_by="row",     outdir=OUTDIR, px_area_m2=px_area_m2_use)
plot_matrix_heatmap(mat_AB, yA, yB, norm_by="global",  outdir=OUTDIR, px_area_m2=px_area_m2_use)

plot_matrix_heatmap(mat_BC, yB, yC, norm_by=None,      outdir=OUTDIR, px_area_m2=px_area_m2_use)
plot_matrix_heatmap(mat_BC, yB, yC, norm_by="row",     outdir=OUTDIR, px_area_m2=px_area_m2_use)
plot_matrix_heatmap(mat_BC, yB, yC, norm_by="global",  outdir=OUTDIR, px_area_m2=px_area_m2_use)

# —— 桑基：聚合（2大类） + 17类（LCZ配色） —— 及其数据表导出
html_agg, (lab_a, src_a, tgt_a, val_a) = sankey_from_mats(mat_AB, mat_BC, yA, yB, yC,
                                                          aggregate=True,  out_html=True,
                                                          outdir=OUTDIR, px_area_m2=px_area_m2_use)
html_full, (lab_f, src_f, tgt_f, val_f) = sankey_from_mats(mat_AB, mat_BC, yA, yB, yC,
                                                           aggregate=False, out_html=True,
                                                           outdir=OUTDIR, px_area_m2=px_area_m2_use,
                                                           colorize_full17=True)

# 导出桑基边表（复用绘图无需再统计）
sankey_agg_csv  = os.path.join(CACHE, f"sankey_edges_{yA}_{yB}_{yC}_agg_{'km2' if px_area_m2_use else 'pixels'}.csv")
sankey_full_csv = os.path.join(CACHE, f"sankey_edges_{yA}_{yB}_{yC}_full17_{'km2' if px_area_m2_use else 'pixels'}.csv")
export_sankey_edges(lab_a, src_a, tgt_a, val_a, sankey_agg_csv)
export_sankey_edges(lab_f, src_f, tgt_f, val_f, sankey_full_csv)

# —— 面积年际曲线（由缓存的 areas_df 直接绘制）——
traj_png = plot_area_trajectory_from_df(areas_df, outdir=OUTDIR, px_area_m2=px_area_m2_use)

# —— 附：导出 Built/Natural 年度表（便于直接“数据画图”）——
bn_df = pd.DataFrame({
    "Year": areas_df.index.astype(int),
    "Built": areas_df[[f"LCZ{i}" for i in BUILT_RANGE]].sum(axis=1).values,
    "Natural": areas_df[[f"LCZ{i}" for i in NATUR_RANGE]].sum(axis=1).values
})
bn_df.set_index("Year", inplace=True)
bn_df.to_csv(os.path.join(CACHE, f"areas_built_natural_{'km2' if px_area_m2_use else 'pixels'}.csv"), encoding="utf-8-sig")

print("Cache files:")
for p in cache_files:
    print(" ", p)
print("Sankey edges CSV:", sankey_agg_csv, sankey_full_csv)
print("Area BN table:", os.path.join(CACHE, f"areas_built_natural_{'km2' if px_area_m2_use else 'pixels'}.csv"))
print("Done. All figures in:", OUTDIR)
