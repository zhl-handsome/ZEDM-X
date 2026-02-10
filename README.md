# ZDEM-X: Route-B 多面体 DEM 接触模型

基于 Feng et al. (2021) *An energy-conserving contact theory for discrete element modelling of arbitrarily shaped particles* 的 Python 实现。

## 项目简介

本项目实现了用于离散元方法（DEM）的非球形多面体颗粒接触模型。采用 Route-B 理论框架，通过接触边界几何参数化和能量一致性原理计算接触力。

### 核心特性

- ✅ 支持任意凸/凹多面体颗粒
- ✅ 多接触区域自动识别与分离
- ✅ 能量守恒的接触力计算
- ✅ 基于 STL 网格的工程实现
- ✅ VTK 可视化输出

## 理论基础

详见 [技术文档](./route_b_多面体dem接触模型技术文档.md)

**关键公式：**

- 向量面积：$\mathbf{S}_n = \frac{1}{2} \oint_{\Gamma} \mathbf{x} \times d\mathbf{x}$
- 接触法向：$\mathbf{n} = - \frac{\mathbf{S}_n}{\|\mathbf{S}_n\|}$
- 接触力：$\mathbf{f} = - \frac{\partial \Psi(\lambda)}{\partial \lambda} \, \mathbf{n}(\lambda)$
- 最小能量模型：$F_n = \frac{k}{2} \|\mathbf{S}_n\|$

## 环境要求

- Python 3.7+
- NumPy
- （可选）ParaView 用于结果可视化

## 安装

```bash
# 克隆项目
git clone <repository-url>
cd ZDEM-X

# 安装依赖
pip install -r requirements.txt
```

## 使用方法

### 基本用法

```bash
python polyline.py meshA.stl meshB.stl --k 1000.0
```

### 完整参数

```bash
python polyline.py meshA.stl meshB.stl \
  --tA 0.0 0.0 0.0 \              # 颗粒A的平移
  --qA 1.0 0.0 0.0 0.0 \          # 颗粒A的旋转四元数
  --tB 0.1 0.0 0.0 \              # 颗粒B的平移
  --qB 1.0 0.0 0.0 0.0 \          # 颗粒B的旋转四元数
  --k 1000.0 \                    # 接触刚度参数
  --tol 1e-6 \                    # 几何容差
  --split-contacts \              # 启用多接触分离
  --vtk routeB_intersection.vtk \ # 输出VTK文件
  --outA A_transformed.vtk \      # 输出颗粒A网格
  --outB B_transformed.vtk        # 输出颗粒B网格
```

### 输出文件

- `A_transformed.vtk`, `B_transformed.vtk`: 变换后的颗粒网格
- `routeB_intersection.vtk`: 接触环和接触点（带 `contact_id` 标签）

### ParaView 可视化

1. 打开 `routeB_intersection.vtk`
2. 在 Cell Data 中选择 `contact_id` 进行着色
3. 每个闭合环对应一个独立接触

## 算法流程

```
STL 网格输入
    ↓
三角形-三角形相交检测
    ↓
线段集合粗分组（可选）
    ↓
端点 Snap 合并
    ↓
T-junction 切分
    ↓
闭合环追踪
    ↓
计算 Sn, Gn, n, xc(0)
    ↓
接触力计算 (Fn = k/2 * |Sn|)
```

## 项目结构

```
ZDEM-X/
├── README.md                                    # 本文件
├── route_b_多面体dem接触模型技术文档.md          # 理论文档
├── polyline.py                                  # 主程序
├── requirements.txt                             # Python 依赖
├── examples/                                    # 示例文件
│   └── (待添加示例 STL 文件)
├── tests/                                       # 测试用例
│   └── (待添加)
└── docs/                                        # 文档
    ├── Feng - 2021 - An energy-conserving...pdf
    └── feng2021.pdf
```

## 技术特点

### 几何处理

- 使用 UnionFind 进行端点合并
- 空间哈希加速邻域查询
- T-junction 自动检测与切分
- 鲁棒的闭合环提取算法

### 数学模型

- **能量一致性**：力由接触势能导数定义
- **唯一法向**：多点接触时法向由整体几何决定
- **力矩一致**：自动满足角动量守恒

## 开发计划

- [ ] 添加单元测试
- [ ] 提供示例 STL 文件
- [ ] 实现切向摩擦力模型
- [ ] 集成到 DEM 仿真框架
- [ ] 性能优化（Cython/C++）
- [ ] GPU 加速支持

## 引用

如果您在研究中使用本代码，请引用：

```bibtex
@article{feng2021energy,
  title={An energy-conserving contact theory for discrete element modelling of arbitrarily shaped particles: Contact volume based model and computational issues},
  author={Feng, YT and Owen, DRJ and Peri{\'c}, D},
  journal={Computer Methods in Applied Mechanics and Engineering},
  volume={373},
  pages={113493},
  year={2021},
  publisher={Elsevier}
}
```

## 许可证

待定

## 贡献

欢迎提交 Issue 和 Pull Request。

## 联系方式

待补充

---

## 技术规格书与实现步骤（C++ / MPI / CUDA）

### 目标与范围

- 目标：基于 polyline.py 的几何思路与 docs/DEM_polyhedron_contact_techdoc.md 的力学模型，开发高性能非球形多面体 DEM 程序（C++ 实现，支持 MPI + CUDA）。
- 支持对象：凸多面体颗粒（粒-粒）与多面体颗粒–壁面三角网格（粒-壁）。
- 受力汇总：同一颗粒对可能存在多条接触曲线/多接触区域，统一汇总为单一“接触对”的合力与力矩。
- 实施优先级：先完成单线程 CPU 基线；面向千万级颗粒规模进行结构与内存预留，但先用小规模验证正确性。

### 非目标（首版不做）

- 非凸多面体的精确交体（可先以凸分解离线处理）。
- 弹塑性、破碎、粘结等复杂本构。
- 隐式积分或自适应时间步。

### 核心物理模型（来自技术文档）

- 粒-粒法向力：
  - 重叠体积：\(\Delta V\)（对偶变换 + Quickhull + 反极化）。
  - 法向弹性力：\(\mathbf{F}_n^e = k_n\,\Delta V^{1/3}\,\mathbf{n}\)。
  - 阻尼力：\(\mathbf{F}_n^d = -c_n(\mathbf{v}_r\cdot\mathbf{n})\mathbf{n}\)。
- 粒-壁法向力：
  - 通过 Sutherland–Hodgman 裁剪得到重叠面积 \(A_n\)，进入刚度计算。
- 切向力（Cundall–Strack）：
  - \(\mathbf{F}_t^e = -k_t\Delta\mathbf{s}\)，\(\mathbf{F}_t^d=-c_t\mathbf{v}_t\)，并施加库仑约束 \(\|\mathbf{F}_t\|\le \mu\|\mathbf{F}_n\|\)。
- 接触点 \(\mathbf{x}_c\) 用于力矩：\(\mathbf{T} = \mathbf{r}_c\times\mathbf{F}_c\)。

### 数据结构设计（CUDA 友好）

- 颗粒状态（SoA）：
  - `pos[N][3]`, `vel[N][3]`, `omega[N][3]`, `quat[N][4]`
  - `mass[N]`, `inv_mass[N]`, `I_body[N][3]`（主惯量），`I_world[N][9]`
  - `radius[N]`（外接球），`material_id[N]`
- 形状几何（只读常量缓冲）：
  - 预存局部坐标系下的 `verts`, `faces`, `face_normals`，以及面片索引范围（压缩 CSR）。
- 接触对（固定容量数组，避免 GPU 动态分配）：
  - `pair_i[M]`, `pair_j[M]`, `pair_type[M]`（粒-粒/粒-壁）
  - `pair_active[M]`, `pair_n[ M][3 ]`, `pair_xc[ M][3 ]`
  - `pair_tangential[ M ][3 ]`（切向累计位移）
  - `pair_kn[M]`, `pair_kt[M]`, `pair_cn[M]`, `pair_ct[M]`
- 邻居与空间哈希：
  - 规则网格 binning：`cell_start[]`, `cell_count[]`, `particle_ids[]`
  - 固定容量：每个 cell 最大粒子数 `MAX_CELL_OCC`，溢出回退到 CPU 或二次处理。

### 计算流程（每步）

1) 构建邻居列表（空间划分 / cell binning）。
2) 宽相：外接球判定，生成候选对 `pair_list`。
3) 窄相：GJK 判定相交（粒-粒）。
4) 相交后：
   - 粒-粒：对偶变换 + Quickhull 生成交体，得到 \(\Delta V\)、\(\mathbf{n}\)、\(A_n\)、\(\mathbf{x}_c\)。
   - 粒-壁：截交多边形 + Sutherland–Hodgman 得到 \(A_n\)、\(\mathbf{x}_c\)。
5) 力学：根据 \(k_n,k_t,c_n,c_t\) 计算 \(\mathbf{F}_n,\mathbf{F}_t\)，施加到两颗粒并累计力矩。
6) 积分：线动用 Velocity-Verlet/Leapfrog；转动用角动量更新 + 四元数归一化。

### 并行化与内存规划

- CUDA：
  - Kernel 1：构建 cell 列表（prefix-sum 或固定桶）。
  - Kernel 2：宽相对生成（每粒子扫邻域 cell）。
  - Kernel 3：窄相 GJK + 交体/裁剪（按 pair 并行）。
  - Kernel 4：力/力矩累加（原子或分块规约）。
- MPI：
  - 空间域分解（3D 砖块），每个 rank 维护本地粒子 + halo 粒子。
  - 交换边界粒子状态（位置/速度/姿态/半径）与接触对跨域归并。
- 内存：
  - `MAX_PAIR_PER_PARTICLE` 与 `MAX_PAIR_TOTAL` 预估上限；超限时记录并降级处理（例如仅保留最近接触）。
  - 交体构造需要的临时缓冲为固定容量工作区（faces/edges/verts 上限）。
  - 千万级规模预算：优先压缩 SoA（float32/float64 可选），接触对与邻域结构需具备可裁剪策略（稀疏保留与溢出统计）。

### 数值稳定与参数

- 时间步：满足接触刚度与质量的稳定性要求（建议由最大刚度估算）。
- 几何容差：GJK 终止误差、裁剪容差、最小面积阈值。
- 材料参数：\(E, \nu, \mu, \varepsilon\)。

### 验证与基准

- 基础：与 polyline.py 的单对接触几何输出对照（\(\mathbf{n}, \mathbf{x}_c\)）。
- 物理：圆柱颗粒–壁面碰撞解析解对比（见技术文档）。
- 网格无关性：壁面三角网格不同分辨率下力–位移曲线收敛。

### 分阶段实现步骤（里程碑）

1) C++ 核心数学库与几何数据结构（向量/四元数/CSR 网格）。
2) CPU 单线程 DEM 基线（宽相 + GJK + 接触力 + 积分）。
3) 小规模验证集：单对/少量粒子回归测试 + 解析解对比。
4) 加入粒-壁裁剪接触与切向历史量存储。
5) 大规模内存审计：千万级粒子 SoA 预算、接触对上限与溢出策略。
6) CPU 并行（OpenMP）与稳定性回归测试。
7) CUDA 版：邻居列表 + 宽相 + 窄相（分阶段逐核移植）。
8) MPI 域分解 + Halo 交换 + 多 GPU 组合。
9) 性能优化：内存布局、核函数融合、负载均衡、稀疏接触压缩。

### 性能目标（建议）

- 单卡：百万级粒子宽相在亚秒级；窄相与力计算可随接触稀疏线性扩展。
- 多卡：MPI+CUDA 线性弱扩展效率 > 70%。

### 输出与可视化

- 支持输出接触对统计、接触力直方图、VTK 粒子/壁面结果，便于与现有 Python 版本对照。
