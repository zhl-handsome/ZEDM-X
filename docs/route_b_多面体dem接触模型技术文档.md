# Route‑B 多面体非球形颗粒 DEM 接触模型技术文档

> 基于 Feng et al. (2021) *An energy‑conserving contact theory for discrete element modelling of arbitrarily shaped particles*

---

## 1. 问题背景与目标

在多面体（包括凸/凹）非球形颗粒的 DEM 模拟中，核心困难不在于“是否发生接触”，而在于：

- **接触区域不是点，而是曲线甚至多条曲线**；
- **接触法向、接触点并非几何最近点**；
- **多接触（multi‑contact）必须能量一致、力矩一致**。

Feng (2021) 提出的 **Route‑B 接触理论** 给出了一个：

> **以交叠几何为基础、通过接触边界 Γ 参数化接触位形、并由能量导数定义接触力** 的统一框架。

本文档完整说明：

1. Route‑B 的几何与力学理论；
2. 从 STL 网格到接触环 Γ 的工程实现；
3. 多接触拆分、拓扑重建与闭合环提取；
4. 论文兼容的最小能量模型；
5. DEM 接触数据结构的最终形式。

---

## 2. 几何基础：交叠体 Ω 与接触边界 Γ

### 2.1 交叠体 Ω

设两个闭合多面体颗粒为 \(A, B\)，当发生穿透时，其交叠体定义为：

\[
\Omega = A \cap B
\]

Ω 是一个三维体，其边界由来自 A、B 的部分面片及它们的**交线**组成。

### 2.2 接触边界 Γ

Route‑B 关注的不是 Ω 本身，而是其**接触边界**：

\[
\Gamma = \partial \Omega \cap (\partial A \cap \partial B)
\]

Γ 的性质：

- 是一组 **闭合空间曲线**；
- 可以是 **多个不相交的闭合环**；
- 每一个闭合环在 DEM 中对应一个 **独立的接触（Contact）**。

---

## 3. Route‑B 的核心几何量

### 3.1 接触边界的参数化

对任意一条闭合接触环：

\[
\Gamma = \{ \mathbf{x}(s) \mid s \in [0,L) \}
\]

其中 s 为曲线参数。

---

### 3.2 向量面积 \(\mathbf{S}_n\)

Feng 定义了接触边界的**向量面积**：

\[
\boxed{
\mathbf{S}_n = \frac{1}{2} \oint_{\Gamma} \mathbf{x} \times d\mathbf{x}
}
\]

性质：

- \(|\mathbf{S}_n|\) 等于接触区域的**等效面积**；
- 方向与接触区域的几何法向一致（右手定则）。

---

### 3.3 接触法向

Route‑B 中 **接触法向不是最近点法向**，而是由 Γ 整体几何决定：

\[
\boxed{
\mathbf{n} = - \frac{\mathbf{S}_n}{\|\mathbf{S}_n\|}
}
\]

这一定义保证：

- 多点/多面接触时法向唯一；
- 与接触能量梯度方向一致；
- 自动满足力矩一致性。

---

### 3.4 几何辅助量 \(\mathbf{G}_n\)

论文中定义：

\[
\boxed{
\mathbf{G}_n = -\frac{1}{3} \oint_{\Gamma} \left( \mathbf{x}\cdot\mathbf{x} \, d\mathbf{x} \right)
}
\]

在离散 polyline 情况下，可由边段逐条累加得到。

---

### 3.5 接触点 \(\mathbf{x}_c(0)\)

接触点定义为 **λ=0 时的能量等效作用点**：

\[
\boxed{
\mathbf{x}_c(0) = \frac{\mathbf{n} \times \mathbf{G}_n}{\|\mathbf{S}_n\|}
}
\]

该点用于：

- 施加接触力；
- 计算对颗粒的力矩；
- 保证能量与角动量一致。

---

## 4. 接触位形参数 λ

### 4.1 定义

Route‑B 将接触状态参数化为一个标量 λ：

- λ = 0：当前几何状态；
- λ > 0：沿接触法向 \(\mathbf{n}\) 推进（压入）；
- λ < 0：回退（分离）。

几何量均可写为 λ 的函数：

\[
\Gamma(\lambda),\; \Omega(\lambda),\; \mathbf{n}(\lambda),\; \Psi(\lambda)
\]

---

## 5. 能量一致的接触力定义（论文核心）

### 5.1 接触势能 \(\Psi(\lambda)\)

论文**不固定**具体形式，只要求：

- 连续；
- \(\Psi(0)=0\)（常见约定）；
- 单调递增；
- 与几何接触一致。

---

### 5.2 力的定义（论文给定）

接触力由能量导数定义：

\[
\boxed{
\mathbf{f}(\lambda) = - \frac{\partial \Psi(\lambda)}{\partial \lambda} \, \mathbf{n}(\lambda)
}
\]

这是 Route‑B 与传统 DEM 的本质区别。

---

## 6. 论文兼容的最小能量模型（Minimal Model）

### 6.1 选择（论文允许）

选择最简单的接触势能：

\[
\boxed{
\Psi(\lambda) = \frac{k}{2} \, \Omega(\lambda)
}
\]

其中：

- k：材料刚度参数；
- \(\Omega(\lambda)\)：交叠体体积。

> 该选择：
> - 不引入额外几何假设；
> - 与 Route‑B 完全兼容；
> - 能量守恒。

---

### 6.2 几何关键结果（Route‑B）

在小 λ 情况下：

\[
\boxed{
\frac{d\Omega}{d\lambda} = A_c = \|\mathbf{S}_n\|
}
\]

其中 \(A_c\) 为接触等效面积。

---

### 6.3 最终接触力公式

代入得：

\[
\boxed{
\mathbf{f} = - \frac{k}{2} \, \|\mathbf{S}_n\| \, \mathbf{n}
}
\]

或标量形式：

\[
\boxed{
F_n = \frac{k}{2} \, \|\mathbf{S}_n\|
}
\]

方向沿接触法向。

---

## 7. 多接触与拓扑重建

### 7.1 多接触的来源

- 一个粒子对可产生 **多个不相交接触环**；
- 每个闭合环对应一个独立 Contact；
- 力与接触 **一一对应**。

---

### 7.2 工程实现步骤

1. **Triangle–Triangle Intersection** → 线段集合；
2. **粗分组**（segment‑segment 距离 BFS）；
3. **端点 Snap**（消除浮点误差）；
4. **T‑junction 切分**；
5. **环追踪**（degree=2 优先）；
6. **每个闭合环独立计算 Sn/Gn/n/xc/F**。

---

## 8. DEM 接触数据结构（最终）

```c
struct Contact {
    int    idA, idB;
    int    cid;          // 接触编号

    // 几何（Route‑B）
    Vec3   Sn;
    Vec3   Gn;
    Vec3   nA, nB;
    Vec3   xc0;
    double area;         // |Sn|

    // 力学（最小能量模型）
    double k;
    double Fn;
    Vec3   F_A;
    Vec3   F_B;
};
```

---

## 9. 总结

- Route‑B 提供的是 **几何参数化 + 能量一致框架**；
- 接触力不是经验公式，而是能量导数；
- 多面体、多接触天然支持；
- 本文给出的模型是：
  - **严格论文内**；
  - **最小且闭合**；
  - **工程可实现、可扩展**。

---

> 至此，你已经拥有一套：
> **几何正确、能量一致、可扩展到工业 DEM 的多面体接触模型完整技术路线。**

