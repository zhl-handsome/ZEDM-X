# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

ZDEM-X 是基于 Route-B 理论的多面体 DEM（离散元方法）接触模型实现。核心算法来自 Feng et al. (2021) 的能量守恒接触理论。

**双实现架构**：
- Python 版本 (`polyline.py`)：几何原型与算法验证
- C++ 版本 (`src/main.cpp`)：高性能 DEM 仿真

## 常用命令

### Python 版本（几何计算）
```bash
# 基本用法
D:/ProgramData/anaconda3/python.exe polyline.py meshA.stl meshB.stl --k 1000.0

# 完整参数
D:/ProgramData/anaconda3/python.exe polyline.py meshA.stl meshB.stl \
  --tA 0.0 0.0 0.0 --qA 1.0 0.0 0.0 0.0 \
  --tB 0.1 0.0 0.0 --qB 1.0 0.0 0.0 0.0 \
  --k 1000.0 --tol 1e-6 --split-contacts \
  --vtk output.vtk
```

### C++ 版本（DEM 仿真）
```bash
# 构建
cmake -S . -B build
cmake --build build --config Debug

# 运行
build/Debug/zdem_cpu --config config/example_sim.txt
```

## 核心架构

### 物理模型

**法向力（能量守恒模型）**：
- 弹性力：`F_n^e = k_n * ΔV^(1/3) * n`
- 阻尼力：`F_n^d = -c_n * (v_r · n) * n`
- 刚度：`k_n = (E1/R1 + E2/R2) * A_n`

**切向力（Cundall-Strack）**：
- 弹性：`F_t^e = -k_t * Δs`
- 阻尼：`F_t^d = -c_t * v_t`
- 库仑约束：`||F_t|| ≤ μ * ||F_n||`

### 几何算法流水线

```
STL 输入 → 三角形-三角形相交 → 线段集合
    → 端点 Snap 合并（Union-Find）→ T-junction 切分
    → 闭合环追踪 → 计算 Sn, Gn, n, xc0 → 接触力计算
```

**关键几何量**：
- `Sn`：向量面积（环路积分）
- `Gn`：几何辅助量
- `nA`：接触法向 = -Sn / ||Sn||
- `xc0`：接触点 = nA × Gn / ||Sn||
- `area`：接触面积 = ||Sn||

### 数据结构（C++）

**粒子状态（SoA 友好）**：
- `pos, vel, omega, quat`：位置/速度/角速度/四元数
- `mass, inv_mass, inertia_body`：质量属性
- `radius, equiv_radius`：外接球/等效半径
- `young, poisson, mu, restitution`：材料参数

**接触对**：
- 切向位移历史 `tangential_disp[pair_key]` 需跨时间步持久化
- `pair_key = (i << 32) | j` 用于标识粒子对

## 代码约定

### 符号映射（论文 → 代码）
| 论文 | Python | C++ | 含义 |
|------|--------|-----|------|
| Γ | `loop_pts` | `loop_pts` | 接触边界闭合环 |
| S_n | `Sn` | `Sn` | 向量面积 |
| G_n | `Gn` | `Gn` | 几何辅助量 |
| n | `nA`, `nB` | `nA` | 接触法向 |
| x_c(0) | `xc0` | `xc0` | 接触点 |
| F_n | `Fn` | `fn` | 法向力大小 |

### 容差处理
- Python: `--tol` 参数或自动估计（网格平均边长 * 0.1）
- C++: `tol = meshA.mean_edge > 0 ? meshA.mean_edge * 0.1 : bbox_diag * 1e-2`

### 输出格式
- VTK 文件（ASCII）：用于 ParaView 可视化
- 包含 CELL_DATA：`id`, `mass`, `radius`, `contact_count`, `velocity`, `omega`, `force`, `torque`

## 配置文件格式（C++）

```text
steps = 1000
dt = 0.001
split_contacts = 1
gravity = 0 -9.81 0
center_mesh = 1
vtk_prefix = particles
output_dir = output
output_interval = 10

particle
stl = geometry/low-poly-banbana.stl
pos = 0 0 0
vel = 1 0 0
quat = 1 0 0 0
omega = 0 0 0
scale = 0.001
density = 2500
young = 1e7
poisson = 0.25
mu = 0.5
restitution = 0.5
end_particle
```

## 参考文档

- 理论基础：`docs/DEM_polyhedron_contact_techdoc.md`
- 中文技术文档：`docs/route_b_多面体dem接触模型技术文档.md`
- 原始论文：`docs/feng2021.pdf`
- 开发指南：`AGENTS.md`
