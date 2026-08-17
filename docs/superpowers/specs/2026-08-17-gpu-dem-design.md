# ZDEM-X GPU 版本设计

日期:2026-08-17
状态:已确认(brainstorming 产出)
前置:v8 统一 per-vertex 罚弹簧接触模型(main.cpp,commit 5c9595e)

## 目标与范围

把 v8 接触模型(DEM 多面体)移植到 GPU。**先中后大**:数据布局与管线按大规模(10⁵–10⁷ 粒子)设计,第一版在中等规模(≤10⁴)落地并验证正确性,大规模优化(B 阶段:排序分桶、shared memory 缓存等)留到真实瓶颈出现后,只动 kernel 内部、不动布局。

硬件/工具:NVIDIA RTX 4090 D(24 GB,CC 8.9),CUDA 12.8,CUB 原语。

非目标(第一版不做):

- 环路遥测管线(三角形相交 → 裁剪 → Sn/Gn)不移植——v8 力不依赖它
- 多 GPU、可迁移粒子、复杂几何(凹槽/黏附)
- 突破 10⁴ 粒子规模的性能调优

## 已确认的关键决策

| 决策点 | 选择 |
|----|----|
| 规模路线 | 先中后大(布局按大规模设计) |
| 精度 | **宏定义可调**:默认 FP32(`using real = float`),`ZDEM_PRECISION=double` 编译 FP64 对照版;能单精度的都用单精度 |
| 架构 | 方案 A:直写 CUDA,SoA + 空间哈希,per-(pair, vertex) 并行,atomicAdd 力累积;排序用 CUB |
| 正确性验收 | 三阶段:FP64 短时逐步对照 + FP32 物理指标 + 中等规模扩展算例 |

## 总体架构

- 新 CMake target `zdem_gpu`(CUDA);`zdem_cpu` 保留为 golden reference 与日常小算例
- config 格式与 CPU 版**完全兼容**(同一批 config 文件直接可用)
- host 侧公共代码(config 解析、STL 加载、VTK 输出)从 main.cpp 抽成 `src/host/` 共享库,两个可执行共链
- 输出 VTK 与现有格式逐字段一致,现有分析脚本(pile_analyze.py / pp_analyze.py / analyze_vtk.py)直接复用
- `real` 贯穿全部物理计算;config 数值、计时等外围用 double 后转换

## 数据布局(为大规模预留,B 阶段不变更)

- **粒子状态 SoA**(device 常驻,`real`):
  - `pos[3N] vel[3N] omega[3N] quat[4N] L[3N]`
  - `mass[N] inv_mass[N] radius[N] equiv_radius[N] inertia_body_inv[9N]`
  - 材料:`young[N] poisson[N] mu[N] restitution[N]`
  - `mesh_index[N]`
- **网格注册表**:全部 STL 的顶点/三角形(体坐标系,center_mesh 之后)拼接平铺 + 偏移表,global memory;接触计算时按粒子姿态旋到世界系(不逐步回写)
- **墙**:世界系三角形静态上传一次;共面分组(quantized n,d)host 启动时一次完成,device 只存组(平面方程 + footprints)
- CUB 排序/扫描 workspace 常驻复用,运行期零逐步分配

## 每步 kernel 序列

与 CPU 循环序完全一致(力 → 积分 → 输出),保证 FP64 逐步对照无相位差:

```
(清力) → build_hash → candidate_pairs → pp_contact → wall_contact → integrate → (可选输出)
```

1. **build_hash**:cell 尺寸 = 2 × max(外接球);粒子 → cell id;CUB radix sort → `cell_start[]/cell_count[]`
2. **candidate_pairs**:每粒子扫 27 邻域 cell,外接球² 重叠 → pair(i<j);两遍法(计数再写)入紧凑数组
3. **pp_contact**(热点,每 pair 一个 CTA):
   - block 内对两 mesh 顶点做对方 AABB 剔除
   - 顶点 × 面并行:内含射线奇偶(+x 射线、y/z AABB 剔除、行列式测试;thread 负责一个面,block 归约交点数)
   - 内含顶点 × 面并行:Ericson 最近点(block 归约取 min,记录最近点)
   - 每内含顶点:罚弹簧 `fn = kh·d^1.5 + Tsuji(vn)`、分布式切向阻尼 + 库仑钳位、推力方向 `-(顶点−最近点)/d`(winding-free)
   - `atomicAdd` 到两粒子力/力矩(公式与 CPU v8 逐条对应)
4. **wall_contact**:每 (粒子, 墙组):平面深度 + `point_in_tri` 足迹判定 → 罚弹簧 → atomicAdd
5. **integrate**:平动半隐式 Euler + 转动中点/指数四元数(per 粒子单线程)
- 诊断:全局 contacts 计数(atomic),日志格式与 CPU 相同字段(去除 loops/segments 列)

## 数值细节

- 数学函数按精度选 `sinf/sqrtf` 或 `sin/sqrt`,不用快速近似 intrinsic;字面量 `real(x)` 包裹
- 容差常量用 `real` 尺度(深度阈值 1e-12 在 FP32 下仍远高于 FLT_MIN,接触深度物理量级 1e-5..1e-3 不受影响)
- **确定性说明**:`atomicAdd` 顺序不定 → 浮点和不可逐位复现;误差 ~ulp(1e-16 相对),验收按容差不按逐位

## 验证计划

| 阶段 | 版本 | 算例 | 通过标准 |
|----|----|----|----|
| 1 算法正确性 | FP64 | 墙单粒子 1k 步;对撞 5k 步;pile12 前 2k 步 | 与 CPU 逐帧对照:位置偏差 < 1e-6 m,能量相对偏差 < 1e-9 |
| 2 物理指标 | FP32(默认) | 墙 e=0.3;对撞;pile6d;pile12 | 与 CPU 同标准:E 单调/完全静止/ALL AT REST;对撞注入 ≤ 0.02 J(单精度噪声余量) |
| 3 中等规模 | FP32 | 新算例 ~64 香蕉随机堆积 | 无爆炸/无穿模;E 单调或恒定;记录性能(步/ms)与 CPU 对比 |

- 工具:`scratch/gpu_compare.py`(逐帧 VTK diff);`tests/` 增 GPU smoke(跑通 + contacts > 0)
- 性能预期:≤12 粒子 GPU 可能不敌 CPU(launch 开销),不作验收项;100+ 粒子应显优,如实记录

## 构建、文件布局与里程碑

- CMake:`enable_language(CUDA)`,`CMAKE_CUDA_ARCHITECTURES` 默认 89;option `ZDEM_PRECISION`(float|double)→ 编译定义 `ZDEM_DOUBLE`;CUDA check 宏(失败带 file:line 退出)
- 文件:
  - `src/host/`:config 解析、STL 加载、VTK 输出(自 main.cpp 抽取)
  - `src/gpu/`:`dem_gpu.cu`(主循环)、kernels(hash / candidates / pp_contact / wall_contact / integrate,`.cuh`)
- 里程碑(每步独立可验证):
  - **M1** host 库抽取,`zdem_cpu` 回归不变(wall smoke + 短 pile)
  - **M2** GPU 墙接触 + 积分:单香蕉落墙 FP64 对照通过
  - **M3** GPU pp 接触:对撞 FP64 对照通过
  - **M4** pile12 FP64 对照 + 全套 FP32 物理指标
  - **M5** 64 粒子堆积算例 + 性能记录

## 错误处理

- 所有 CUDA 调用经 check 宏,非零返回码 → stderr 打印(file:line + cudaGetErrorString)→ 退出码 1
- config 校验与 CPU 版同规则;不支持的参数(未来才有的)明确报错而非忽略
