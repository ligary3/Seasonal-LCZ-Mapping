import pandas as pd
import plotly.graph_objects as go
import os
num = 4
# ==========================================
# 0. 配置 (Configuration)
# ==========================================
# CSV 所在的文件夹 (Step 1 生成的)
STATS_DIR = r"E:\luwen2\result\season\HMM\2018season2\csv_stats_cache"
OUTPUT_DIR = os.path.join(STATS_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 过滤阈值：流转像素小于此值的连线不画 (防止图太乱)
# 根据你的图像大小调整，一般 500-2000 比较合适
THRESHOLD = 1000 

# 目标类别: LCZ 9 + LCZ A, B, C, D, F, G (剔除 E)
# Key: CSV列索引(0-16), Value: 显示标签
# LCZ 1-17 对应索引 0-16。例如 LCZ A(11) -> Index 10
TARGET_MAP = {
    8:  "LCZ 9",  # Sparse built
    10: "LCZ A",  # Dense trees
    11: "LCZ B",  # Scattered trees
    12: "LCZ C",  # Bush
    13: "LCZ D",  # Low plants
    15: "LCZ F",  # Bare soil
    16: "LCZ G"   # Water
}

# 颜色配置 (标准 LCZ 颜色或论文常用色)
COLOR_MAP = {
    "LCZ 9": "#ffccaa", # Red
    "LCZ A": "#006a00", # Dark Green
    "LCZ B": "#00aa00", # Green
    "LCZ C": "#648525", # Olive
    "LCZ D": "#b9db79", # Light Green
    "LCZ F": "#fbf7ae", # Khaki/Sand
    "LCZ G": "#6a6aff"  # Blue
}

# ==========================================
# 1. 数据读取与处理
# ==========================================
def get_flows(csv_path):
    """读取矩阵CSV，返回 (源索引, 目标索引, 流量值) 的列表"""
    if not os.path.exists(csv_path):
        print(f"Warning: File not found {csv_path}")
        return []
    
    df = pd.read_csv(csv_path, index_col=0)
    matrix = df.values
    flows = []
    
    # 遍历矩阵，只提取关注类别的流
    valid_indices = sorted(TARGET_MAP.keys())
    
    for src_idx in valid_indices:
        for dst_idx in valid_indices:
            value = matrix[src_idx, dst_idx]
            
            # 过滤掉极小值，保持图表清晰
            if value > THRESHOLD:
                flows.append((src_idx, dst_idx, value))
                
    return flows

# ==========================================
# 2. 构建 Plotly Sankey 数据结构
# ==========================================
print("Preparing Sankey data...")

seasons = ["Spring", "Summer", "Autumn", "Winter"]
node_labels = []
node_colors = []
link_sources = []
link_targets = []
link_values = []
link_colors = []

# --- A. 构建节点 (Nodes) ---
# 节点 ID 映射: (季节索引, LCZ原始索引) -> Plotly全局节点ID
id_map = {}
current_id = 0

# 每一列(季节)都生成一组节点
for season_idx, season_name in enumerate(seasons):
    for lcz_idx, lcz_name in TARGET_MAP.items():
        # 记录映射关系
        id_map[(season_idx, lcz_idx)] = current_id
        
        # 设置节点显示的标签 (例如 "Spring LCZ A")
        # 为了美观，第一列显示全名，后面可以简化
        label = f"{season_name}<br>{lcz_name}" 
        node_labels.append(label)
        
        # 设置节点颜色
        node_colors.append(COLOR_MAP.get(lcz_name, "#888888"))
        
        current_id += 1

# --- B. 构建连线 (Links) ---
# 定义转换过程: (源季节Idx, 目标季节Idx, CSV文件名)
transitions = [
    (0, 1, "matrix_SPRING_to_SUMMER.csv"),
    (1, 2, "matrix_SUMMER_to_AUTUMN.csv"),
    (2, 3, "matrix_AUTUMN_to_WINTER.csv")
]

# 辅助函数: Hex颜色转RGBA (增加透明度)
def hex_to_rgba(hex_code, alpha=0.4):
    hex_code = hex_code.lstrip('#')
    return f"rgba({int(hex_code[:2], 16)}, {int(hex_code[2:4], 16)}, {int(hex_code[4:], 16)}, {alpha})"

for src_season, dst_season, filename in transitions:
    full_path = os.path.join(STATS_DIR, filename)
    flows = get_flows(full_path)
    
    for src_lcz, dst_lcz, value in flows:
        # 获取 Plotly 的节点 ID
        source_id = id_map[(src_season, src_lcz)]
        target_id = id_map[(dst_season, dst_lcz)]
        
        link_sources.append(source_id)
        link_targets.append(target_id)
        link_values.append(value)
        
        # 连线颜色：跟随【源节点】的颜色，但半透明
        lcz_name = TARGET_MAP[src_lcz]
        base_color = COLOR_MAP.get(lcz_name, "#888888")
        link_colors.append(hex_to_rgba(base_color, alpha=0.35))

# ==========================================
# 3. 绘图与导出
# ==========================================
print("Plotting...")

fig = go.Figure(data=[go.Sankey(
    # 节点属性
    node=dict(
        pad=15,             # 节点垂直间距
        thickness=20,       # 节点宽度
        line=dict(color="black", width=0.5), # 节点边框
        label=node_labels,
        color=node_colors,
        hovertemplate='%{label}<br>Pixels: %{value}<extra></extra>' # 悬停提示
    ),
    # 连线属性
    link=dict(
        source=link_sources,
        target=link_targets,
        value=link_values,
        color=link_colors,
        hovertemplate='%{source.label} -> %{target.label}<br>Pixels: %{value}<extra></extra>'
    )
)])

# 全局布局设置
fig.update_layout(
    #title_text="<b>Seasonal Phenological Transitions of Natural LCZ Classes</b>",
    title_x=0.5, # 标题居中
    font_family="Times New Roman",
    font_size=12+num,
    width=1400,   # 图片宽度
    height=700,   # 图片高度
    margin=dict(l=30, r=30, t=60, b=30), # 边距
    plot_bgcolor='white',
    paper_bgcolor='white'
)

# 1. 导出 HTML (交互式，必做)
html_path = os.path.join(OUTPUT_DIR, "Fig_Sankey_Natural.html")
fig.write_html(html_path)
print(f"Interactive HTML saved to: {html_path}")

# 2. 导出 PNG (静态图，需 kaleido)
try:
    png_path = os.path.join(OUTPUT_DIR, "Fig_Sankey_Natural.png")
    # scale=3 相当于 300 DPI
    fig.write_image(png_path, scale=3) 
    print(f"High-res PNG saved to: {png_path}")
except Exception as e:
    print("---------------------------------------------------------")
    print(f"Error saving PNG: {e}")
    print("Tip: Install kaleido via 'pip install -U kaleido' to fix.")
    print("---------------------------------------------------------")

# 如果在 Jupyter 里运行，可以取消下面注释直接显示
# fig.show()