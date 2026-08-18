# AGENTS.md - ZDEM-X 开发指南

本文档为 AI 编程代理提供项目规范和最佳实践。

---

## 项目概览

**项目名称**: ZDEM-X  
**核心功能**: 基于 Route-B 理论的多面体 DEM(离散元)接触模型,v8 统一 per-vertex 罚弹簧接触(wall 与粒子-粒子同一模型)

**三实现架构**(同一 v8 物理模型):
- `polyline.py`:Python 几何原型与算法验证(本文下半部分的 Python 风格指南适用于此)
- `zdem_cpu`(C++):`src/main.cpp` + `src/host/` 共享库,高性能仿真与 **golden reference**
- `zdem_gpu`(CUDA):`src/gpu/`,与 CPU 逐条对应;FP64 对照逐位一致,FP32 同判据达标,64 粒子 49× 加速

**技术特点**:
- 任意凸/凹多面体颗粒接触计算(banana STL 有 19/102 面朝内,winding-free 处理)
- STL 网格输入/VTK 可视化输出(CPU/GPU 逐字段一致)
- 能量精确闭合的接触力模型(无跨步接触状态)
- 空间哈希 broadphase(CPU O(N²) 对照;GPU CUB radix sort)

日常开发规范、测试方法论、工作流程规则见 `CLAUDE.md`(权威,持续更新);GPU 设计与验收记录见 `docs/superpowers/specs|plans|notes/`。

---

## 构建与测试命令

### 环境设置

```bash
# Python 环境 (Windows Anaconda)
D:/ProgramData/anaconda3/python.exe -m pip install -r requirements.txt

# 或激活 conda base 环境
D:/ProgramData/anaconda3/Scripts/activate.bat base
```

### C++ / CUDA 构建(CMake ≥ 3.18,MSVC + nvcc)

```bash
# CPU 版
cmake -S . -B build && cmake --build build --config Release
build/Release/zdem_cpu --config config/example_sim.txt

# GPU 版(默认 FP32;FP64 对照构建 -DZDEM_PRECISION=double)
cmake -S . -B build -DZDEM_PRECISION=float && cmake --build build --config Release
build/Release/zdem_gpu --config config/example_sim.txt

# 测试(CTest,5 项烟雾/回归)
cd build && ctest -C Release --output-on-failure
```

config 格式 CPU/GPU 完全兼容;VTK 输出一致,`scratch/` 下分析脚本(analyze_vtk.py / pile_analyze.py / pp_analyze.py / gpu_compare.py)直接复用。

### 运行 Python 原型

```bash
# 基本用法
python polyline.py meshA.stl meshB.stl --k 1000.0

# 完整参数示例
python polyline.py meshA.stl meshB.stl \
  --tA 0.0 0.0 0.0 \
  --qA 1.0 0.0 0.0 0.0 \
  --tB 0.1 0.0 0.0 \
  --qB 1.0 0.0 0.0 0.0 \
  --k 1000.0 \
  --tol 1e-6 \
  --split-contacts \
  --vtk routeB_intersection.vtk \
  --outA A_transformed.vtk \
  --outB B_transformed.vtk
```

### 测试命令

```bash
# C++/CUDA 测试套件(CTest,在 build 目录)
cd build && ctest -C Release --output-on-failure

# 几何单元测试(C++,源 tests/)
ctest -C Release -R geometry

# GPU 对照工具(逐帧 VTK diff,FP64 验收用)
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py <dir_cpu> <dir_gpu> <n_particles>
```

### 代码检查

**当前状态**: 无 linter 配置文件

```bash
# 建议使用 (未配置)
pylint polyline.py
flake8 polyline.py --max-line-length=100
black polyline.py --check  # 代码格式化检查
```

---

## 代码风格指南

### 导入规范

```python
# 标准库
import os
import struct
import argparse

# 第三方库
import numpy as np

# 本地模块 (当前项目为单文件)
# from .module import function
```

**顺序**: 标准库 → 第三方库 → 本地模块，用空行分隔。

### 格式化约定

- **缩进**: 4 空格 (不使用 Tab)
- **行长度**: 无严格限制 (建议 < 120 字符)
- **字符串**: 双引号优先 (`"string"`)
- **编码声明**: `# -*- coding: utf-8 -*-`
- **Shebang**: `#!/usr/bin/env python3`

### 类型注解

**当前状态**: 项目未使用类型注解

建议添加 (保持兼容 Python 3.7+):
```python
def quat_normalize(q: np.ndarray) -> np.ndarray:
    ...

def load_stl(path: str) -> np.ndarray:
    ...
```

### 命名约定

| 类型 | 约定 | 示例 |
|------|------|------|
| 模块 | snake_case | `polyline.py` |
| 函数 | snake_case | `tri_normals()`, `seg_seg_dist2()` |
| 类 | PascalCase | `UnionFind` |
| 常量 | UPPER_SNAKE_CASE | `EPS = 1e-12` |
| 变量 (局部) | snake_case | `tol`, `segments`, `comp_idx` |
| 变量 (数学) | 数学符号 | `nA`, `nB`, `Sn`, `Gn`, `xc0` |
| 私有辅助函数 | _leading_underscore | `_cell_index()` |

**特殊约定** (符合论文符号):
- 数学量保持简洁: `Sn` (向量面积), `Gn` (几何辅助量), `nA`/`nB` (法向)
- 三维数组: `tris`, `pts`
- 几何参数: `tol` (容差), `h` (网格尺寸)

### 注释规范

```python
# 模块级 docstring (文件顶部)
"""
Route-B (Feng 2021) for closed STL contacts, with robust loop reconstruction.

Pipeline (per pair A,B):
1) Triangle-triangle intersections -> raw segment soup (p0,p1)
...
"""

# 函数 docstring
def accumulate_Sn_Gn_from_polyline(loop_pts):
    """
    loop_pts: list/array of shape (M,3), closed (first!=last is ok; we'll wrap)
    Uses same discrete formulas on each edge (p_i, p_{i+1}).
    """
    ...

# 行内注释 (关键逻辑)
tau = tau / lt  # normalize intersection direction
```

**原则**:
- 关键算法步骤必须注释 (如 T-junction split, loop tracing)
- 公式引用论文章节/公式编号
- 行尾注释与代码间至少 2 空格

### 错误处理

```python
# 明确的错误消息
if len(n_tri_bytes) < 4:
    raise ValueError("Invalid STL.")

# 几何退化情况
if ln < 1e-30:  # 使用 EPS 或明确的阈值
    return None, None

# 输入验证
if abs(denom) < 1e-30:
    return False, None
```

**原则**:
- 使用 `ValueError` 处理输入错误
- 几何退化返回 `None` 而非抛出异常
- 数值阈值使用命名常量 (`EPS = 1e-12`)

### 数据结构

**向量/矩阵**: 使用 NumPy 数组
```python
# 正确
pts = np.array(points, dtype=float)
n = np.cross(e1, e2)

# 避免
pts = [list(p) for p in points]  # 避免嵌套列表
```

**字典输出** (DEM 数据结构):
```python
contact = {
    "idA": idA, "idB": idB,
    "cid": int(contact_ids[k]),
    "nVert": int(len(loop)),
    "Sn": Sn.tolist(),  # NumPy → list for JSON serialization
    "nA": nA.tolist(),
    "xc0": None if xc0 is None else xc0.tolist(),
    "Fn": float(Fn),
}
```

### 算法实现约定

**几何容差**:
```python
# 使用一致的容差
tol = 1e-6  # 命令行参数
EPS = 1e-12  # 全局数值阈值

# 距离比较
if np.linalg.norm(p1 - p2) <= tol:
    ...
```

**循环防护**:
```python
# 避免无限循环
steps = 0
while steps < 100000:  # 明确的上限
    steps += 1
    ...
```

**空间哈希**:
```python
# 使用元组键的字典作为网格
grid = {}  # (ix, iy, iz) -> [item_indices]
c = _cell_index(point, h)
grid.setdefault(c, []).append(item)
```

---

## 特定领域约定

### 几何算法

1. **三角形-三角形相交**: 使用平面-线段测试
2. **端点合并**: Union-Find + 空间哈希
3. **T-junction 处理**: 迭代单遍扫描
4. **闭合环提取**: 度优先 (degree=2 节点优先)

### 数学符号映射

| 论文符号 | 代码变量 | 含义 |
|----------|----------|------|
| Γ | `loop_pts` | 接触边界闭合环 |
| **S**_n | `Sn` | 向量面积 |
| **G**_n | `Gn` | 几何辅助量 |
| **n** | `nA`, `nB` | 接触法向 |
| **x**_c(0) | `xc0` | 接触点 |
| λ | `lambda` (概念) | 接触位形参数 |
| Ψ | `Psi` (概念) | 接触势能 |
| F_n | `Fn` | 法向力大小 |

### VTK 输出规范

```python
# 变换后的网格
write_vtk_polydata_tris("A_transformed.vtk", tris, comment="Mesh A")

# 接触环 + 接触点
write_vtk_polydata_polylines(
    "routeB_intersection.vtk",
    loops_pts,       # 闭合环列表
    contact_points,  # xc0 列表
    loop_contact_ids # CELL_DATA contact_id
)
```

---

## 开发流程

### 添加新功能

1. **理论验证**: 参考 `route_b_多面体dem接触模型技术文档.md` 和 `docs/` 中的论文
2. **模块化**: 独立函数 (如当前架构)
3. **测试数据**: 在 `examples/` 目录提供 STL 示例
4. **输出验证**: 使用 ParaView 可视化 VTK 输出

### 性能优化指南

**已落地**:
- C++ 版(`zdem_cpu`)替代 Python 热点
- GPU 版(`zdem_gpu`):64 粒子 49× 加速,交叉点 4–6 粒子(测量方法与规模表见 `docs/superpowers/notes/2026-08-18-gpu-performance.md`)

**后续方向**(GPU B 阶段,10⁵–10⁷ 粒子真实瓶颈出现后):
1. 排序分桶 / shared memory 缓存网格(只动 kernel 内部、不动布局)
2. kernel 分解计时(cudaEvent)定位瓶颈
3. 多 mesh 类型与可迁移粒子

### 常见任务

**修改接触力模型**:
```python
# 当前: 最小能量模型 (Fn = k/2 * |Sn|)
def build_contacts_from_loops(..., k_stiff=1.0):
    ...
    Fn = 0.5 * k_stiff * float(area)  # 修改此处

# 扩展: 添加切向摩擦
# 参考论文 Section X (待补充)
```

**调整容差参数**:
```bash
# 几何鲁棒性问题时调大
python polyline.py ... --tol 1e-5

# 精细网格时调小
python polyline.py ... --tol 1e-7
```

---

## 参考资源

**论文**:
- Feng, Y.T., Owen, D.R.J., & Perić, D. (2021). *An energy-conserving contact theory for discrete element modelling of arbitrarily shaped particles*. Computer Methods in Applied Mechanics and Engineering, 373, 113493.
- 位于 `docs/` 目录

**在线文档**:
- [NumPy API Reference](https://numpy.org/doc/stable/reference/)
- [VTK File Format](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
- [ParaView User Guide](https://docs.paraview.org/)

---

## 许可证与贡献

**许可证**: 待定  
**贡献流程**: 欢迎提交 Issue 和 Pull Request (见 README.md)

---

**文档版本**: 1.1 (2026-08-18,补充 C++/CUDA 实现与真实测试命令;Python 风格指南仍适用于 `polyline.py`)  
**生成工具**: AI Agent (Sisyphus)
