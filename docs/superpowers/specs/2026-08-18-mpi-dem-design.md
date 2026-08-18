# ZDEM-X MPI 版设计(第一阶段:CPU 通信层)

日期:2026-08-18
状态:已确认(brainstorming 产出)
前置:v8 统一 per-vertex 罚弹簧接触模型(`zdem_cpu`,golden reference)+ GPU 移植(`zdem_gpu`,FP64 parity 验收体系)

## 目标与路线

把 ZDEM-X 并行化到多进程/多节点,最终形态 **MPI × 多 GPU**。本设计只覆盖第一阶段:**zdem_cpu 的 MPI 并行(通信层打通并验证)**;GPU 搬运是后续独立 spec。通信层按设备无关设计(操作 host 缓冲),搬运时复用。

- 拓扑目标:跨节点集群(通信按真实网络设计:聚合打包、通信全封装、非阻塞预留)
- 本地验证:MS-MPI 单机多进程(SDK 环境已确认)
- 参考:open-source-DEM 参考库(`.claude/skills/dem-oss-reference/`)——LIGGGHTS `comm.cpp`(brick 分解、三轴链式迁移)、LAMMPS ghost 体系、PhasicFlow `MPIParallelization/`(多 GPU 形态参考)

**v8 的结构性红利**:接触模型无跨步状态(无切向弹簧历史、无锁定法向)→ ghost 粒子只是只读镜像,LIGGGHTS 类接触历史粒子的反向同步问题不存在。这是选 v8 做并行底座的根本原因。

非目标(第一阶段不做):

- 动态负载均衡/重分区(接口预留)
- 增量 ghost 更新、非阻塞通信 overlap(接口预留)
- CPU 空间哈希 broadphase(GPU 版已有,反向移植属独立工作)
- GPU × MPI、动态墙、粒子插入/删除

## 已确认的关键决策

| 决策点 | 选择 |
|----|----|
| 路线 | 先 CPU 通信层验证,后搬 GPU |
| 拓扑 | 跨节点集群(本地 MS-MPI 验证)|
| 域分解 | 静态 brick(规则长方体)+ 质心越界迁移 |
| 代码组织 | 方案 A:独立 `zdem_mpi` 可执行 + 共享 `zdem_core`;`zdem_cpu`/`zdem_gpu` 一字不动 |
| 接触归属 | Newton-off:每 rank 只算本地粒子为 owner 的对,ghost 只读 |
| 反向通信 | **无**(Newton-off 换来;代价 = 跨域 pair 双算,堆积场景 ~10-20% 计算冗余)|
| 等价验收 | FP64 容差 parity(位置 <1e-6,能量 <1e-9;与 GPU 版同口径)|

## 总体架构

```
zdem_mpi(新):
  src/mpi/decomp.{hpp,cpp}    -- 全局盒 + brick 拓扑(MPI_Cart)+ 粒子归属
  src/mpi/migrate.{hpp,cpp}   -- 质心越界迁移(三轴链式 Sendrecv)
  src/mpi/ghost.{hpp,cpp}     -- ghost 层全量重建(三轴链式,每步)
  src/mpi/gather.{hpp,cpp}    -- 输出步 VTK gather + 日志/能量 reduce
  src/mpi/mpi_main.cpp        -- 主循环(力→积分→迁移→ghost 重建(为下一步)→输出,力/积分序与 CPU 一致)
```

- 复用 `zdem_core`(STL/mesh/config/VTK);网格注册表与墙组**全 rank 复制**(构建一次,banana 102 面量级可忽略)
- `mpiexec -n 1` 走完整 MPI 代码路径 → parity 基线天然存在

## 进程模型与全局盒

- `MPI_Cart_create` 3D 拓扑,px×py×pz 按全局盒长宽比自动分解
- **全局盒**:默认 = 初始粒子包围盒每侧外扩 `mpi_margin`(默认 5×max 外接球;需覆盖粒子整个演化空间——落到墙的堆积场景,z 下界由 margin 给出即可);config 可 `mpi_box = xmin xmax ymin ymax zmin zmax` 显式覆盖
- 粒子越出全局盒 → 明确报错退出(堆积场景天然安全)
- config 新增 `mpi_box` / `mpi_margin` 可选 key,config_io 解析为可选字段,旧可执行(zdem_cpu/gpu)忽略不报错

## 粒子所有权与迁移

- owner = 质心所在 brick;全局 gid = 初始编号,迁移保持
- 迁移 = LAMMPS comm_brick 式三轴链式(x→y→z,每轴 `MPI_Cart_shift` + `MPI_Sendrecv`),角落粒子自动两跳
- 迁移数据 = 粒子全部状态(pos/vel/omega/quat/L + 质量属性 + 材料 + mesh_index + gid)
- 不变量(每步断言,debug 下):全局粒子集合 gid 唯一、无丢失无重复(可用 allreduce 计数校验)

## Ghost 交换

- ghost 深度 = 2×max(外接球)(与 GPU broadphase cell 同口径 → 任意接触对双方至少一方拥有)
- 每步顺序:`力计算 → 积分 → 迁移 → ghost 重建` 下一步力前完成(每步全量重建,第一版)
- ghost 携带(只读):pos/quat/vel/omega/inv_mass/radius/equiv_radius/young/poisson/mu/restitution/mesh_index/gid(≈22 real + 2 int/粒子),三轴链式 Sendrecv,单缓冲打包
- ghost 不写、不积分、不迁移

## 接触计算(Newton-off)

- 每 rank:for 本地粒子 i × (本地 ∪ ghost) 的 j(gid 不同):包围球重叠 → v8 per-vertex 罚力,只累加到 i
- 墙:本地粒子 × 全部墙组(墙全 rank 复制)
- contacts 本地计数,日志 `MPI_Reduce` 求和
- rank 内 broadphase 沿用 CPU O(n²) 双循环(单 rank 10³–10⁴ 可接受;范围外记录)
- 力公式与 CPU v8 **逐条对应**(直接搬 `apply_penalty`/wall 段),不改物理

## 输出与日志

- 输出步:`MPI_Gatherv` 各 rank 本地快照(pos/vel/omega/quat/contact_count + gid)到 rank0 → 按 gid 排序 → 现有 `write_vtk_particles` 写**单一全局帧文件**,文件名与字段与单进程完全一致 → 分析脚本零改动
- 日志:rank0 打印,字段沿用 CPU 格式,contacts 为全局和;每输出步加 `rank_n=` 行(各 rank 粒子数分布,负载诊断)
- 能量:各 rank 本地 KE/PE → `MPI_Reduce` 求和 → 全局 E,判据同单进程

## 验证阶梯

| 阶段 | 算例 | 通过标准 | 结果(2026-08-18 实测)|
|----|----|----|----|
| V0 | `mpiexec -n 1` vs `zdem_cpu`,wall64(40k 步 400 帧)+ pile12(16k 步 160 帧)| 位置 <1e-6,能量 <1e-9(容差;枚举序不同 → ulp 级,不逐位)| **PASS 且字节级一致**(pos 0.000e+00/energy 0.000e+00;`cmp` 全帧一致,强于容差判据;comm_util 重构后复验仍一致)|
| V1 | n=2(x 二分)pile12 16k 全程(跨域接触链)| 同上 | **PASS 且字节级一致**(160 帧;另 4 进程 1x1x4 分解复验一致)|
| V2 | 单粒子跨域飞行 + 对撞跨域 | 轨迹 parity;gid 唯一性断言过 | **PASS 且字节级一致**(fly 20 帧 + pp 80 帧;pile12 n=2 迁移段 rank_n 10,2→12,0 后 158 帧一致;越盒 Abort 与 8 轮迁移上限 Abort 路径实测触发)|
| V3 | pile64(double)× 4 进程 + 更大规模 | 全局 E 单调/恒定,无爆炸穿模,同单进程判据(MPI 版为 CPU 谱系纯 double,无 FP32 变体)| **PASS**:E 92.864→6.895 J,尾段单调缓降(最大帧间回升 0.015 J = 0.22% ≤ ±1%);min_z 0.017 自 20k 步起恒定,无穿模无爆炸;静止比例 60/64 = 93.8% ≥ 90%(4 粒子尾段 |ω|≤3.7 瞬态微振,与 GPU 版同算例表现一致)。**附接触下迁移 parity 窗(24k 步 12 帧)字节级一致**:rank_n 8,24,24,8 → 7,25,25,7(step 20k,contacts=50)→ 6,26,25,7(step 28k,contacts=81),迁移与接触同时活跃且 parity 保持|
| 性能 | 本机 1/2/4 进程 step_ms | 如实记录(单机带宽,不外推集群)| 受控 20k 窗:n=1/2/4 = 5.676/2.678/2.091 ms/步;n=4 vs n=1 = **2.71×**,vs zdem_cpu 11.41 = 5.46×(跨窗口径,zdem_cpu 含环路遥测;披露见 notes);稳态接触段 7.00 ms/步;详见 `../notes/2026-08-18-mpi-performance.md`|

工具:复用 `scratch/gpu_compare.py`(逐帧 VTK diff,输出格式一致故直接可用)。

## 构建与平台

- CMake:`find_package(MPI REQUIRED)` 只挂 `zdem_mpi` target;无 MPI 的机器 zdem_cpu/zdem_gpu 照常构建
- 只用 MPI **C 绑定**(mpi.h;MPICH/OpenMPI/MS-MPI 兼容;不用废弃的 C++ 绑定)
- 本机 MS-MPI(SDK 环境变量已确认);集群模块加载 MPICH/OpenMPI

## 风险与注意

- **brick 尺寸 ≥ ghost 深度**:切分数过多导致薄 brick 会漏 ghost——启动时校验并报错(限制每轴切分数)
- pair 枚举序与单进程不同 → 浮点和 ulp 级差;验收按容差(与 GPU 版 atomicAdd 同口径)
- 跨域 pair 双算的冗余随边界占比上升;大规模下如成瓶颈,后续切 Newton-on + reverse comm(记录为备选)
- 堆积沉底后上稀下密 → 负载不均;静态 brick 第一版接受,动态重分区留接口
- 迁移/ghost 全量重建每步开销:堆积稳态下 ghost 粒子占比小,可接受;增量优化留后
