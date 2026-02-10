# AGENTS.md - ZDEM-X 开发指南

本文档为 AI 编程代理提供项目规范和最佳实践。

---

## 项目概览

**项目名称**: ZDEM-X  
**核心功能**: 基于 Route-B 理论的多面体 DEM 接触模型 Python 实现  
**主要文件**: `polyline.py` (839 行)  
**依赖**: NumPy (Python 3.7+)

**技术特点**:
- 任意凸/凹多面体颗粒接触计算
- STL 网格输入/VTK 可视化输出
- 能量守恒的接触力模型
- 多接触区域自动识别

---

## 构建与测试命令

### 环境设置

```bash
# Python 环境 (Windows Anaconda)
D:/ProgramData/anaconda3/python.exe -m pip install -r requirements.txt

# 或激活 conda base 环境
D:/ProgramData/anaconda3/Scripts/activate.bat base
```

### 运行主程序

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

**当前状态**: 测试套件待添加 (见 README.md 开发计划)

```bash
# 单元测试 (待实现)
python -m pytest tests/

# 单个测试文件 (待实现)
python -m pytest tests/test_geometry.py -v

# 特定测试函数 (待实现)
python -m pytest tests/test_geometry.py::test_snap_endpoints -v
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

**当前瓶颈** (README 开发计划):
- 三角形-三角形相交: O(N²) 暴力检测
- T-junction 切分: 迭代扫描

**建议方向**:
1. 空间加速 (BVH / Octree)
2. Cython 重写热点函数
3. C++ 扩展 (`pybind11`)
4. GPU 加速 (CUDA)

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

**文档版本**: 1.0 (2026-02-10)  
**生成工具**: AI Agent (Sisyphus)
