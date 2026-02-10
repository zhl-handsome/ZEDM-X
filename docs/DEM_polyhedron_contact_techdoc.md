# DEM 技术文档：多面体颗粒受力计算与颗粒–壁面碰撞（基于所给两篇论文）

> 适用范围：凸多面体颗粒 DEM；壁面以 STL/OBJ 三角面片网格表示；法向力采用能量守恒的“重叠体积”接触理论；切向力采用 Cundall–Strack 弹簧-阻尼+库仑摩擦；壁面接触采用面片投影+多边形裁剪得到“重叠面积”。

---

## 1. 记号与基本量

- 颗粒 1/2 的质心位置：\\(\\mathbf{x}_1,\\mathbf{x}_2\\)；质量：\\(m_1,m_2\\)  
- 接触点（用于相对速度、力矩计算）：\\(\\mathbf{x}_c\\)  
- 质心到接触点向量：\\(\\mathbf{r}_c = \\mathbf{x}_c-\\mathbf{x}\\)（对每个颗粒分别取自身质心）  
- 相对速度（在接触点处）：\\(\\mathbf{v}_r\\)  
- 法向单位向量：\\(\\mathbf{n}\\)（由接触几何确定）  
- 切向相对速度：\\(\\mathbf{v}_t = \\mathbf{v}_r - (\\mathbf{v}_r\\cdot\\mathbf{n})\\mathbf{n}\\)  
- 库仑摩擦系数：\\(\\mu\\)；恢复系数：\\(\\varepsilon\\)  
- 有效质量：\\(m_{\\mathrm{eff}} = (m_1^{-1}+m_2^{-1})^{-1}\\)

---

## 2. 接触几何：多面体–多面体（重叠体积 \\(\\Delta V\\)）

### 2.1 接触检测（宽相+窄相）
1) **宽相**：用外接球（bounding sphere）过滤。若两颗粒质心距离小于外接球半径之和，则进入窄相。  
2) **窄相**：使用 **GJK**（Gilbert–Johnson–Keerthi）算法在 Minkowski difference 上判定是否相交。

### 2.2 计算重叠体积 \\(\\Delta V\\)（精确凸交体）
对已判定相交的凸多面体对：
1) **对偶空间变换（polar transformation）**：每个三角面片 \\(f\\)（外法向 \\(\\mathbf{n}_f\\)）映射为对偶点 \\(\\mathbf{y}\\)。  
2) **增量 Quickhull**：在对偶空间构造凸包。  
3) **反极化（unpolarization）**：将对偶凸包面反变换为交体的顶点，得到交体的三角网格拓扑。  
4) **体积**：以交体顶点均值为内部参考点 \\(\\mathbf{o}\\)，将每个面三角形 \\((\\mathbf{v}_0,\\mathbf{v}_1,\\mathbf{v}_2)\\) 与 \\(\\mathbf{o}\\) 组成四面体并求和：
\\[
V = \\frac{1}{6}\\sum_f (\\mathbf{v}_0-\\mathbf{o})\\cdot\\left((\\mathbf{v}_1-\\mathbf{o})\\times(\\mathbf{v}_2-\\mathbf{o})\\right)
\\]
该 \\(V\\) 即两凸多面体的重叠体积 \\(\\Delta V\\)。

> 备注：工程实现中通常为每个接触对分配预估上限的工作区（faces/edges/verts 缓冲），以避免 GPU 动态分配。

---

## 3. 接触几何：多面体–壁面三角面片（重叠面积 \\(A_n\\)）

### 3.1 为什么需要“重叠面积”与加权
STL/OBJ 壁面通常由大量三角面片组成。直接按“面片是否接触”逐片施加力容易造成 **重复接触（duplicate contact）**，导致不合理的粘附或弹跳。解决思路是显式计算 **颗粒与每个三角面片的重叠面积**，并用该面积参与刚度/力的计算（面积越大贡献越大）。

### 3.2 多边形裁剪：Sutherland–Hodgman（二维）
对单个壁面三角面片：
1) 将三角面片所在平面作为投影平面。  
2) 计算多面体与该平面的截交多边形，并投影到该平面，作为“subject polygon” \\(P\\)。  
3) 将壁面三角形作为“clipping polygon” \\(clip\\)。  
4) 依次用 clip 的 3 条有向边对 P 进行裁剪，得到最终交叠多边形（仍在同一平面）。  
5) 将最终多边形三角剖分求面积，得到该面片的重叠面积 \\(A_n\\)。  

裁剪中，线段 \\(\\overline{\\mathbf{v}_2\\mathbf{v}_3}\\) 与裁剪边直线的交点可写为：
\\[
\\mathbf{I}(\\mathbf{v}_2,\\mathbf{v}_3)=\\mathbf{v}_2+(\\mathbf{v}_3-\\mathbf{v}_2)\\,t
\\]
其中 \\(t\\) 由二维叉积形式给出（见论文公式）。

---

## 4. 法向接触力：能量守恒的重叠体积模型（颗粒–颗粒）

该模型的核心要求：法向力来自 **接触势能函数**（避免能量不一致），并通过重叠体积 \\(\\Delta V\\) 确定。

### 4.1 弹性法向力
\\[
\\mathbf{F}_n^{e} = k_n\\,\\Delta V^{1/3}\\,\\mathbf{n}
\\]
- \\(\\mathbf{n}\\)：单位法向  
- \\(\\Delta V\\)：重叠体积  
- \\(k_n\\)：法向刚度

### 4.2 法向刚度 \\(k_n\\)
\\[
k_n = \\left(\\frac{E_{P1}}{R_{P1}}+\\frac{E_{P2}}{R_{P2}}\\right)A_n
\\]
- \\(E_{Pi}\\)：颗粒材料杨氏模量  
- \\(R_{Pi}\\)：等效半径（与颗粒体积相等的球半径）  
- \\(A_n\\)：投影接触面积（在颗粒–颗粒接触中，由接触几何得到；在颗粒–壁面中可取对应面片重叠面积）

> 直觉：\\(\\Delta V^{1/3}\\) 有“等效穿透深度”的量纲；\\(A_n\\) 则类似接触区域尺度。

### 4.3 粘性阻尼法向力
为描述耗散并抑制振荡，引入：
\\[
\\mathbf{F}_n^{d} = -c_n\\,(\\mathbf{v}_r\\cdot\\mathbf{n})\\,\\mathbf{n}
\\]

阻尼系数与恢复系数 \\(\\varepsilon\\) 的关系（对数衰减形式）：
\\[
c_n = \\frac{2\\ln(\\varepsilon)\\,\\sqrt{k_n m_{\\mathrm{eff}}}}{\\sqrt{\\ln^2(\\varepsilon)+\\pi^2}}
\\]

### 4.4 总法向力
\\[
\\mathbf{F}_n = \\mathbf{F}_n^{e}+\\mathbf{F}_n^{d}
\\]

---

## 5. 切向接触力：Cundall–Strack + 库仑约束

### 5.1 总切向力与摩擦约束
\\[
\\mathbf{F}_t = \\mathbf{F}_t^{e}+\\mathbf{F}_t^{d},\\qquad \\|\\mathbf{F}_t\\|\\le \\mu\\,\\|\\mathbf{F}_n\\|
\\]

### 5.2 弹性切向力（切向“弹簧”）
\\[
\\mathbf{F}_t^{e}=-k_t\\,\\Delta\\mathbf{s}
\\]
- \\(\\Delta\\mathbf{s}\\)：接触期间累计切向位移（需要在每个接触对上持久存储，并在接触结束时清零）

### 5.3 切向刚度 \\(k_t\\)
\\[
k_t = \\left(\\frac{E_{P1}}{2(1+\\nu_{P1})R_{P1}}+\\frac{E_{P2}}{2(1+\\nu_{P2})R_{P2}}\\right)A_n
\\]
- \\(\\nu_{Pi}\\)：泊松比

### 5.4 切向阻尼力
\\[
\\mathbf{F}_t^{d}=-c_t\\,\\mathbf{v}_t,\\qquad
c_t = \\frac{2\\ln(\\varepsilon)\\,\\sqrt{k_t m_{\\mathrm{eff}}}}{\\sqrt{\\ln^2(\\varepsilon)+\\pi^2}}
\\]

---

## 6. 接触合力、力矩与转动动力学

### 6.1 单接触对的合力
\\[
\\mathbf{F}_c = \\mathbf{F}_n + \\mathbf{F}_t
\\]

### 6.2 单接触对力矩
对每个接触对 \\(c\\)：
\\[
\\mathbf{T}_c = \\mathbf{r}_c\\times \\mathbf{F}_c
\\]

总力矩为所有接触对求和：
\\[
\\mathbf{T} = \\sum_c \\mathbf{T}_c
\\]

### 6.3 角动量守恒与惯量张量
\\[
\\dot{\\mathbf{L}}=\\mathbf{T},\\qquad \\mathbf{L}=\\mathbf{I}\\cdot\\boldsymbol{\\Omega}
\\]
多面体的惯量张量随取向变化。体坐标系主惯量为 \\(\\hat{\\mathbf{I}}=\\mathrm{diag}(I_{xx},I_{yy},I_{zz})\\)，全局惯量通过旋转矩阵 \\(\\mathcal{H}\\) 得到：
\\[
\\mathbf{I}=\\mathcal{H}\\,\\hat{\\mathbf{I}}\\,\\mathcal{H}^T
\\]

### 6.4 四元数更新（避免欧拉角奇异）
取四元数 \\(\\mathbf{Q}=(q_0,q_1,q_2,q_3)\\)，由体坐标角速度 \\((\\omega_x,\\omega_y,\\omega_z)\\) 更新 \\(\\dot{\\mathbf{Q}}\\)（见论文给出的分量形式），再进行归一化以保持单位四元数。

---

## 7. 颗粒–壁面碰撞与网格分辨率敏感性（验证用）

### 7.1 圆柱颗粒–壁面碰撞（无摩擦、无重力）
用于验证数值稳定性：设初始角速度为 0；法向恢复系数为 \\(\\varepsilon\\)；碰撞几何由入射角等参数确定。解析解给出碰后角速度与法向速度：
\\[
\\omega_y^{+} = \\frac{mV_z^{-}(1+\\varepsilon)r\\cos(\\alpha+\\theta)}{I_{yy}+mr^2\\cos^2(\\alpha+\\theta)}
\\]
\\[
V_z^{+} = \\omega_y^{+}r\\cos(\\alpha+\\theta) - \\varepsilon V_z^{-}
\\]

### 7.2 壁面网格分辨率对受力的影响
- 壁面用不同分辨率的 STL 三角网格表示，定义无量纲尺度比：
\\[
\\lambda = \\frac{L_c}{L_m}
\\]
其中 \\(L_c\\) 为颗粒特征尺寸（如立方体边长），\\(L_m\\) 为网格面片平均尺寸。
- 颗粒沿平面以恒速滑动，位移无量纲化：
\\[
\\kappa = \\frac{\\Delta x}{L_c}
\\]
观察力–位移曲线在不同 \\(\\lambda\\) 下是否重合，以检验壁面接触力的“网格无关性”。

---

## 8. 实现要点（给 Codex 的可执行提示）

### 8.1 每个时间步的推荐流程（接触对并行）
1) 构建邻居列表（空间 binning）→ 生成 contactPair（一维对列表）。  
2) 宽相过滤（外接球）。  
3) 窄相：GJK 判定相交。  
4) 若颗粒–颗粒相交：  
   - 计算交体（对偶变换 + Quickhull + 反极化）→ 得到 \\(\\Delta V\\)、接触法向 \\(\\mathbf{n}\\)、投影面积 \\(A_n\\)、接触点 \\(\\mathbf{x}_c\\)。  
   - 计算 \\(k_n,k_t,c_n,c_t\\)，再算 \\(\\mathbf{F}_n,\\mathbf{F}_t\\)，施加到两颗粒并累加力矩。  
5) 若颗粒–壁面：  
   - 对每个候选三角面片：投影多边形 → Sutherland–Hodgman 裁剪 → 得到重叠面积 \\(A_n\\)（以及可取多边形重心作为 \\(\\mathbf{x}_c\\)）。  
   - 用 \\(A_n\\) 进入刚度/力模型，计算该面片贡献并累加。  
6) 积分更新：线动用显式积分（如 Velocity-Verlet/Leapfrog）；转动用角动量方程 + 四元数更新。  

### 8.2 状态量存储建议
- 每个“活跃接触对”需要持久存储 \\(\\Delta\\mathbf{s}\\)（切向累计位移）。  
- 若使用 GPU，接触对并行时对粒子力/力矩累加需要原子操作或分块规约。

---

## 9. 参考（来自你提供的两篇论文）
- Xu 等：*Large-Scale Multiple-GPU-based DEM Simulation of Polyhedral Particle Systems*  
- Feng：*An energy-conserving contact theory for discrete element modelling of arbitrarily shaped particles* (2021)

