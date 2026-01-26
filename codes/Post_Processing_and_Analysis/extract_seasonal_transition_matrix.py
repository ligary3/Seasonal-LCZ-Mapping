#导出CSV文件
# -*- coding: utf-8 -*-
"""
LCZ Seasonal Statistics Extractor (2018)
Function: 
  1. Reads 4 seasonal GeoTIFFs (Spring, Summer, Autumn, Winter).
  2. Computes pixel counts per class (Area).
  3. Computes full 17x17 Transition Matrices (Sp->Su, Su->Au, Au->Wi).
  4. Exports all data to CSVs for lightweight plotting later.
"""

import os, json
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window
from rasterio.vrt import WarpedVRT

# =============== 0. 配置 (Configuration) =================
# 输入文件夹
INPUT_DIR = r"E:\luwen2\result\season\HMM\2018season2"

# 文件名映射
FILES = {
    "SPRING": "2018spring.tif",
    "SUMMER": "2018summer.tif",
    "AUTUMN": "2018autumn.tif",
    "WINTER": "2018winter.tif"
}

# 输出文件夹 (会自动创建)
OUTPUT_DIR = os.path.join(INPUT_DIR, "csv_stats_cache")

# 处理块大小 (根据内存调整，2048通常很安全)
BLOCK_SHAPE = (2048, 2048)

# 季节顺序
SEASONS_ORDER = ["SPRING", "SUMMER", "AUTUMN", "WINTER"]

# ================= 1. 工具函数 =================
def get_tif_path(season_key):
    return os.path.join(INPUT_DIR, FILES[season_key])

def window_iter(ds):
    """分块生成器"""
    h, w = ds.height, ds.width
    bh, bw = BLOCK_SHAPE
    for row_off in range(0, h, bh):
        rows = min(bh, h - row_off)
        for col_off in range(0, w, bw):
            cols = min(bw, w - col_off)
            yield Window(col_off=col_off, row_off=row_off, width=cols, height=rows)

def histogram_1to17(arr):
    """统计单张图的直方图 (忽略0值)"""
    valid = (arr >= 1) & (arr <= 17)
    if not np.any(valid):
        return np.zeros(18, dtype=np.int64)
    return np.bincount(arr[valid].ravel(), minlength=18).astype(np.int64)

def transition_17x17(src_arr, dst_arr):
    """统计两张图的转移矩阵"""
    valid = (src_arr >= 1) & (src_arr <= 17) & (dst_arr >= 1) & (dst_arr <= 17)
    if not np.any(valid):
        return np.zeros((17, 17), dtype=np.int64)
    
    # 转换为 0-16 索引
    s = src_arr[valid].ravel().astype(np.int16) - 1
    d = dst_arr[valid].ravel().astype(np.int16) - 1
    
    # 使用 bincount 快速统计
    k = (s * 17 + d).astype(np.int32)
    bc = np.bincount(k, minlength=17*17).astype(np.int64)
    return bc.reshape((17, 17))

# ================= 2. 主逻辑 =================
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # --- Part A: 统计每个季节的面积 (Pixel Counts) ---
    print("Step 1: Computing Area Statistics...")
    area_stats = {}
    
    for season in SEASONS_ORDER:
        path = get_tif_path(season)
        print(f"  Reading {season}...")
        
        hist = np.zeros(18, dtype=np.int64)
        with rasterio.open(path) as ds:
            # 获取像元面积 (如果投影坐标系)
            res = ds.res
            pixel_area_m2 = res[0] * res[1]
            
            for w in window_iter(ds):
                arr = ds.read(1, window=w)
                hist += histogram_1to17(arr)
        
        area_stats[season] = hist[1:] # 取索引 1-17
    
    # 导出面积表
    df_area = pd.DataFrame(area_stats, index=[f"LCZ{i}" for i in range(1, 18)])
    df_area.to_csv(os.path.join(OUTPUT_DIR, "seasonal_areas_pixels.csv"))
    print(f"  Area stats saved. (Pixel area approx: {pixel_area_m2:.2f} m2)")

    # --- Part B: 计算转移矩阵 (Sp->Su, Su->Au, Au->Wi) ---
    print("\nStep 2: Computing Transition Matrices...")
    
    # 定义需要计算的转换对
    pairs = [
        ("SPRING", "SUMMER"),
        ("SUMMER", "AUTUMN"),
        ("AUTUMN", "WINTER")
    ]
    
    for (s1, s2) in pairs:
        print(f"  Processing transition: {s1} -> {s2} ...")
        path1 = get_tif_path(s1)
        path2 = get_tif_path(s2)
        
        matrix = np.zeros((17, 17), dtype=np.int64)
        
        # 使用 WarpedVRT 确保两个影像像元严格对齐 (以 s1 为基准)
        with rasterio.open(path1) as ds1, rasterio.open(path2) as ds2_src:
             vrt_options = {
                'resampling': rasterio.enums.Resampling.nearest,
                'crs': ds1.crs,
                'transform': ds1.transform,
                'height': ds1.height,
                'width': ds1.width,
                'nodata': 0
            }
             
             with WarpedVRT(ds2_src, **vrt_options) as ds2:
                 for w in window_iter(ds1):
                     arr1 = ds1.read(1, window=w)
                     arr2 = ds2.read(1, window=w)
                     matrix += transition_17x17(arr1, arr2)
        
        # 导出矩阵 CSV
        df_mat = pd.DataFrame(matrix, 
                              index=[f"LCZ{i}" for i in range(1, 18)],
                              columns=[f"LCZ{i}" for i in range(1, 18)])
        out_name = f"matrix_{s1}_to_{s2}.csv"
        df_mat.to_csv(os.path.join(OUTPUT_DIR, out_name))
        print(f"    Saved {out_name}")

    # 导出元数据
    with open(os.path.join(OUTPUT_DIR, "meta.json"), "w") as f:
        json.dump({"pixel_area_m2": pixel_area_m2, "crs": str(ds1.crs)}, f)

    print(f"\nAll Done! Data exported to: {OUTPUT_DIR}")