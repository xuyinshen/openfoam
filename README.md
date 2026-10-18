# openfoam-engineering

> 基于 OpenFOAM 12 的两个 CFD 工程算例：**液滴撞击壁面的 VOF 两相流**（局部加密网格 + MPI 并行）与**搅拌器滑移网格动网格**（snappyHexMesh + nonConformalCyclic + k-ω SST）。

[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-12-blue.svg)](https://www.openfoam.org/)
[![MPI](https://img.shields.io/badge/parallel-MPI%20(8%20procs)-green.svg)]()
[![Language](https://img.shields.io/badge/language-C%2B%2B%20%2F%20Python-orange.svg)]()

**TL;DR** — Two OpenFOAM-12 cases: (1) `droplet-laminar` — VOF (`interFoam`) simulation of a droplet impacting a flat wall with `topoSet` + `refineMesh` local refinement; (2) `mixer` — sliding-mesh (`solidBody` rotation, 60 rpm) stirred-tank case meshed with `snappyHexMesh`, solved with RAS k-ω SST. Both run in parallel with MPI.

---

## 目录结构

```
openfoam-engineering/
├── droplet-laminar/                 # 算例一：液滴撞击壁面（VOF 两相流）
│   ├── droplet-template/            #   可复算的 OpenFOAM 算例模板
│   │   ├── runcode                  #   一键流水线：画网格 → 求解 → 重组
│   │   ├── data_collect.py          #   PyVista 后处理（定量指标 + 图表）
│   │   ├── U.org / alpha.water.org  #   初始场模板（.org 防止被覆盖）
│   │   ├── 0/  constant/  system/   #   标准 OpenFOAM 三目录
│   │   └── system/topoSetDict, refineMeshDict, setFieldsDict  # 局部加密与初始化
│   ├── visualize.pvsm               #   ParaView 状态文件（一键还原可视化）
│   ├── 2cm.mp4                      #   初始高度 2 cm：撞击后无溅射
│   ├── 3cm.mp4                      #   初始高度 3 cm：发生溅射
│   └── 4cm.mp4                      #   初始高度 4 cm：溅射更剧烈
└── mixer/                           # 算例二：搅拌器动网格
    ├── Allrun                       #   一键执行 = MakeMesh + SolveRun
    ├── MakeMesh                     #   画网格（含滑移界面构造）
    ├── SolveRun                     #   并行求解（须在 MakeMesh 之后）
    ├── Allclean                     #   清理网格与算例产物
    ├── U.org p.org k.org omega.org nut.org   # 初始场模板
    ├── constant/
    │   ├── triSurface/              #   impeller.stl, rotorZone.stl（几何）
    │   ├── dynamicMeshDict          #   solidBody 旋转（60 rpm）
    │   └── momentumTransport        #   RAS kOmegaSST
    ├── system/                      #   blockMeshDict / snappyHexMeshDict /
    │                                #   createBafflesDict / createNonConformalCouplesDict
    └── Example-Oneround.mp4         #   搅拌器附近流线在一转内的变化
```

> **关于仓库体积**：OpenFOAM 的完整求解产物（各时间步场数据、`processor*` 分解目录、加密后的网格）动辄数 GB，不适合入库。
> 因此本仓库只保留**几何 / 算例模板 / 可执行脚本 / 后处理代码 / 结果动画**，网格与时间步由脚本在本地生成，保证仓库轻量且 `git clone` 后即可复算。

---

## 环境依赖

| 依赖 | 版本 / 说明 | 必需 |
|---|---|---|
| OpenFOAM | 12（`interFoam`、`foamRun`、`snappyHexMesh`、`refineMesh`、`createNonConformalCouples` 均来自官方发行版） | ✅ |
| MPI | OpenMPI / MPICH，`mpirun` 需在 `PATH` 中；默认 **8 进程**（可在 `runcode` / `SolveRun` 中改 `-np`） | ✅ |
| ParaView | 用于打开 `visualize.pvsm` 或 `paraFoam` 查看结果 | 推荐 |
| Python | 3.8+，后处理用：`pyvista`、`numpy`、`pandas`、`matplotlib` | 仅 `data_collect.py` 需要 |

```bash
# 后处理脚本依赖
pip install pyvista numpy pandas matplotlib
```

运行前请先 source OpenFOAM 环境：`source /opt/openfoam12/etc/bashrc`（路径按你的安装调整）。

---

## 算例一：`droplet-laminar` — 液滴撞击壁面（VOF 两相流）

### 物理模型

- **求解器**：`interFoam`，VOF（Volume of Fluid）界面捕捉，两相不可压缩；
- **层流**模型（无湍流封闭），驱动力为**重力 + 表面张力**；
- 关注对象：液滴自由下落后**撞击水平桌面**的铺展、回缩与形变过程。

### 网格策略

背景网格由 `blockMesh` 生成，随后对液滴运动与撞击区域做**局部加密**：

```
blockMesh → topoSet → refineMesh -overwrite → setFields
```

`topoSetDict` 框定需要加密的区域，`refineMesh -overwrite` 就地加密，避免全域加密带来的算力浪费——液滴界面与近壁区域保持高分辨率，远场保持粗网格。

### 运行方式

```bash
cd droplet-laminar/droplet-template
sh ./runcode
```

`runcode` 是一条完整流水线，等价于：

```bash
cp -r alpha.water.org 0/alpha.water ; cp -r U.org 0/U      # 恢复初始场
blockMesh                                                   # 背景网格
topoSet && refineMesh -overwrite                            # 局部加密
setFields                                                   # 初始化液滴相分布
decomposePar                                                # 区域分解
mpirun -np 8 interFoam -parallel                            # MPI 并行求解
reconstructPar                                              # 结果重组
```

> 切换撞击高度（2 / 3 / 4 cm）需要修改 `setFieldsDict` 中液滴的初始位置。
> 后续计划提供三套开箱即用的 case 目录（见 Roadmap）。

### 结果与讨论

| 初始高度 | 动画 | 现象 |
|---|---|---|
| 2 cm | [`2cm.mp4`](droplet-laminar/2cm.mp4) | **未发生溅射**。液滴以较低动能撞击壁面，铺展后在**重力与表面张力**共同驱动下完成回缩与形状振荡 |
| 3 cm | [`3cm.mp4`](droplet-laminar/3cm.mp4) | **发生溅射**。撞击动能超过毛细阈值，破碎出卫星液滴；残留在壁面上的主液滴继续在重力与表面张力驱动下形变 |
| 4 cm | [`4cm.mp4`](droplet-laminar/4cm.mp4) | **溅射更剧烈**。卫星液滴数量更多、飞散更远，壁面残留液滴的形状演化过程同样由重力与表面张力主导 |

> GitHub 不支持在 README 中直接内嵌本地视频，点击上表链接即可查看 / 下载动画。

### 定量后处理

```bash
cd droplet-laminar/droplet-template
# 先生成一个 drop.foam / case.foam 空文件供 ParaView/PyVista 读取
python data_collect.py
```

`data_collect.py` 基于 PyVista 逐时间步读取 `alpha.water > 0.5` 的液相区域，并**取最大连通域**以排除飞溅的卫星液滴干扰，输出：

- **铺展直径** `Spreading_Diameter`（m）— 撞击铺展能力
- **液滴高度** `Drop_Height`（m）— 回缩与振荡
- **质心坐标** `Centroid_X/Y/Z` — 体积加权的运动学轨迹
- **总体积** `Total_Volume` — 作为**质量守恒**的检查项

结果写入 `droplet_impact_results.csv`（可直接导入 Origin / Excel），并自动生成 4 张随时间演化的曲线图。
ParaView 用户可直接加载 `visualize.pvsm` 复现论文配图级别的渲染设置。

---

## 算例二：`mixer` — 搅拌器动网格（滑移网格）

### 物理模型

- **求解器**：`foamRun`（OpenFOAM 12 统一求解入口），不可压缩；
- **湍流模型**：RAS **k-ω SST**（`constant/momentumTransport`）；
- **动网格**：`constant/dynamicMeshDict` 采用 `solidBody` 刚体运动，对 `cellZone rotor` 施加绕 **z 轴、60 rpm** 的旋转。

### 网格与滑移界面

```bash
blockMesh                    # 背景网格
surfaceFeatures              # 提取几何特征边
snappyHexMesh -overwrite     # 贴合 impeller 几何的体网格
createBaffles -overwrite     # 生成 rotorInterface 挡板的 slidingRotor / slidingStator
createNonConformalCouples    # 建立绕 z 轴的旋转非 conformal 耦合（滑移界面）
```

几何由 `constant/triSurface/impeller.stl`（叶轮）与 `rotorZone.stl`（旋转区域）驱动 `snappyHexMesh` 生成；`rotorZone` 内网格随叶轮整体刚体旋转，通过**非 conformal cyclic 耦合**在任意转角下与静止区保持通量连续——这是滑移网格（sliding mesh）相比 MRF 更能反映真实瞬态干涉的关键。

### 运行方式

```bash
cd mixer

sh ./Allrun        # 一键执行：MakeMesh + SolveRun
# 或者分步执行
sh ./MakeMesh      # 画网格（含滑移界面构造）
sh ./SolveRun      # 并行求解，必须在 MakeMesh 之后运行
sh ./Allclean      # 清理网格与算例产物
```

`SolveRun` 默认 8 进程：`decomposePar` → `mpirun --oversubscribe -np 8 foamRun -parallel` → `reconstructPar`，每一步日志写入 `log.*`。

### 结果

- [`Example-Oneround.mp4`](mixer/Example-Oneround.mp4)：搅拌器附近**流线在一整转（60 rpm → 1 s）内的演化**，可清晰看到叶轮尾迹与流场结构的周期性响应。

---

## Roadmap

- [ ] `droplet-laminar`：把 `runcode` 拆分为与 `mixer` 一致的 `MakeMesh` / `SolveRun` / `Allrun`，统一两个算例的入口
- [ ] `droplet-laminar`：提供 2 cm / 3 cm / 4 cm 三套开箱即用的 case（当前需手动改 `setFieldsDict`）
- [ ] `droplet-laminar`：补充铺展因子随 Weber 数的定量对比（脚本已具备数据采集能力）
- [ ] `mixer`：补充扭矩 / 功率准数与混合时间的定量后处理

---

## License

[MIT](LICENSE)
