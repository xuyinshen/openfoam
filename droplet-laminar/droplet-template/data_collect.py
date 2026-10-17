import pyvista as pv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ================= 0. 设置 Matplotlib 中文字体 =================
# 防止中文显示为方块 (Ubuntu 推荐 WenQuanYi，兼容 Windows/Mac)
plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei', 'SimHei', 'PingFang SC', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示坐标轴上的负号
# ===============================================================

# ================= 1. 配置区 =================
# 修改为你的 OpenFOAM 案例路径 (如果是当前目录，保持 "./" 即可)
case_path = "./"  
foam_file = os.path.join(case_path, "drop.foam") # OpenFOAM 12 中通常叫 case.foam

# 假设液滴撞击在 Z=0 的平面上。如果你的壁面在其他 Z 坐标，请修改此值
wall_z = 0.0 

# VOF 方法中，通常取 alpha.water = 0.5 作为液滴边界
alpha_threshold = 0.5
# ==============================================

print(f"🚀 正在读取 OpenFOAM 数据: {foam_file}")
if not os.path.exists(foam_file):
    raise FileNotFoundError(f"❌ 找不到文件: {foam_file}。请检查路径和文件名。")

# 初始化 PyVista 的 OpenFOAM Reader
reader = pv.OpenFOAMReader(foam_file)

# 获取所有可用的时间步
time_values = reader.time_values
print(f"✅ 检测到 {len(time_values)} 个时间步。")

# 用于存储结果的列表
results = []

# 【关键】关闭 PyVista 的自动渲染，否则每读一个时间步都会弹窗，导致卡死
pv.OFF_SCREEN = True 

for i, t in enumerate(time_values):
    print(f"处理时间步 {i+1}/{len(time_values)}: t = {t:.4f} s", end="\r")
    
    # 设置当前时间步并读取网格
    reader.set_active_time_value(t)
    raw_mesh = reader.read()
    
    # 【关键修改】OpenFOAM 读取出来是 MultiBlock，包含内部网格和边界 patch。
    # 我们只需要内部流体网格 (internalMesh)，它通常是第一个元素 (索引 0)。
    if isinstance(raw_mesh, pv.MultiBlock):
        mesh = raw_mesh[0]
    else:
        mesh = raw_mesh
    
    # --- 核心后处理逻辑 ---
    
    # 1. 提取液滴区域 (alpha.water > 0.5)
    droplet = mesh.threshold(alpha_threshold, scalars="alpha.water", preference="cell")
    
    # 如果当前时间步没有液滴数据，跳过
    if droplet.n_cells == 0:
        continue
        
    # 2. 【核心技巧】提取最大连通域（主液滴），排除飞溅的卫星液滴！
    droplet = droplet.connectivity(extraction_mode="largest")
    
    # 提取主液滴表面 (用于计算精确的边界框)
    surface = droplet.extract_surface()
    
    # 3. 计算几何参数
    bounds = surface.bounds # 返回 (xmin, xmax, ymin, ymax, zmin, zmax)
    
    # 铺展直径 (Spreading Diameter)：取 X 和 Y 方向的最大跨度
    D_x = bounds[1] - bounds[0]
    D_y = bounds[3] - bounds[2]
    D_spreading = max(D_x, D_y) 
    
    # 液滴高度 (Drop Height)：Z 方向的最大跨度
    H_drop = bounds[5] - bounds[4]
    
    # 4. 计算运动学参数 (质心)
    # 计算每个网格的体积
    droplet_with_vol = droplet.compute_cell_sizes()
    cell_volumes = droplet_with_vol.cell_data["Volume"]
    cell_centers = droplet_with_vol.cell_centers().points
    
    # 体积加权质心
    total_volume = np.sum(cell_volumes)
    centroid = np.sum(cell_centers * cell_volumes[:, np.newaxis], axis=0) / total_volume
    
    # 5. 保存当前时间步的数据
    results.append({
        "Time": t,
        "Spreading_Diameter": D_spreading,
        "Drop_Height": H_drop,
        "Centroid_X": centroid[0],
        "Centroid_Y": centroid[1],
        "Centroid_Z": centroid[2],
        "Total_Volume": total_volume
    })

# ================= 2. 数据保存与可视化 =================
print("\n\n✨ 处理完成！正在保存数据并生成图表...")

# 转换为 Pandas DataFrame
df = pd.DataFrame(results)

# 保存为 CSV 文件 (方便导入 Origin 或 Excel)
csv_filename = "droplet_impact_results.csv"
df.to_csv(csv_filename, index=False)
print(f"📊 数据已保存至: {csv_filename}")

# 绘制结果图表
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("液滴碰撞后处理结果", fontsize=16)

# 1. 铺展直径随时间变化
axs[0, 0].plot(df["Time"], df["Spreading_Diameter"], 'b-', linewidth=2)
axs[0, 0].set_xlabel("时间 (s)")
axs[0, 0].set_ylabel("铺展直径 (m)")
axs[0, 0].set_title("铺展直径随时间变化")
axs[0, 0].grid(True)

# 2. 液滴高度随时间变化
axs[0, 1].plot(df["Time"], df["Drop_Height"], 'r-', linewidth=2)
axs[0, 1].set_xlabel("时间 (s)")
axs[0, 1].set_ylabel("液滴高度 (m)")
axs[0, 1].set_title("液滴高度随时间变化")
axs[0, 1].grid(True)

# 3. 质心 Z 坐标 (反映下落和反弹过程)
axs[1, 0].plot(df["Time"], df["Centroid_Z"], 'g-', linewidth=2)
axs[1, 0].set_xlabel("时间 (s)")
axs[1, 0].set_ylabel("质心高度 Z (m)")
axs[1, 0].set_title("质心高度随时间变化")
axs[1, 0].grid(True)

# 4. 质量守恒检查 (总体积)
axs[1, 1].plot(df["Time"], df["Total_Volume"], 'k-', linewidth=2)
axs[1, 1].set_xlabel("时间 (s)")
axs[1, 1].set_ylabel("总体积 (m³)")
axs[1, 1].set_title("质量守恒 (总体积)")
axs[1, 1].grid(True)

plt.tight_layout()
plt.show()
