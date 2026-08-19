# CPU/MPI 性能优化批次落地记录(遥测门控/变换缓存/broadphase 哈希/每 rank 计时)

日期:2026-08-19
代码:分支 `feat/cpu-perf-optimizations`(dffc5ce..),计划见 `docs/superpowers/plans/2026-08-19-cpu-perf-optimizations.md`
来源:扩展性测试报告(`docs/superpowers/notes/2026-08-19-mpi-scaling.md` + `scratch/mpi_scaling_report.md`)识别的优化清单第 1/2/3/4/7 项

## 改动一览

| # | 改动 | 提交 |
|---|------|------|
| 1 | `route_b_telemetry` config 开关,默认关 | 11b31b6 |
| 2 | 每步每粒子世界系三角形缓存 + `pp_contact_pair` 11 参重载 | be7b040 |
| 3 | 末步 ghost 交换 guard + gids/local 等长检查 | 0ea8c8d |
| 4 | 共享 spatial-hash broadphase(O(N²)→O(N)),枚举序保持 | 7c5dbac |
| 5 | MPI 每 rank 计时日志字段 `ms_max`/`comm_ms_max` | 84b1e1c |

**硬约束兑现**:全部改动 VTK 帧输出字节级不变(纯函数外提 + 门控 + 字典序 = 旧双循环序);门控证据:CPU wall64/pile12、MPI n=1 wall64/pile12、MPI n=2 vs n=1、pile64 20k、N=1024 首帧、telemetry-ON pile12,全部 `PASS |dpos|=|dE|=0.000e+00`。

## (a) 遥测成本归因链(引用时注意)

- 2026-08-19 扩展性报告的"CPU 45.0 vs MPI n=1 15.7 ms/步(2.87×)"是**跨实现口径混杂**(CPU 带 Route-B 遥测、MPI 不带),不应引用为加速
- 纯遥测成本:T1 提交实测 4.27×(8k 自由落体窗,368.8→86.3s);HEAD 全批次后复测同窗 **17.2s(OFF)vs 155.2s(ON)= 9.0×**(OFF 路被变换缓存 + broadphase 进一步加速后比值自然放大;同机同放置背靠背)
- commit 11b31b6 message 里的 "2.9x" 引的是最初的口径混杂数——正确引用以本节为准
- 遥测 ON 模式:config 加 `route_b_telemetry = 1`;VTK 输出不受遥测影响(字节级验证);`tri_ms` 字段不再含 transform 时间(外提到计时器外)

## (b) broadphase 诚实数字

- 同二进制配对 A/B 计时:pile64 −11%、ff1024 −10%;**N=256 在噪声带内(−2..4%)**——本仓库规模的堆积布局稀释(256 粒子 / ~1200 cells,27 邻域扫描多为空 cell),哈希开销 ≈ 暴力扫描
- 收益定位 **N≥10³ 起**,面向 Phase B(10⁵–10⁷);`cell = 2·max_radius·(1+2⁻³⁰)`,安全域 **<2e6 cells/轴**(超出则静默漏 pair——注释内有完整 ulp 论证)
- 已知退化:强多粒径(cell 由 max_r 决定,大粒子+海量小粒子 → occupancy 膨胀 → 趋近 O(N²));与 GPU broadphase 同设计
- 序保持不变量:`broadphase_pairs` 输出 = O(N²) i<j 双循环的 pair 集与字典序(距离测试操作数逐点一致);已 ctest 化(`tests/test_broadphase.cpp`:多布局逐元素对照参考双循环 + ghost 语义/退化/相切边界,变异验证断言有牙)

## (c) ms_max / comm_ms_max 口径

- `phys` 段:wtris 缓存 + broadphase + pp + 墙 + 积分;`comm` 段:迁移 + ghost 交换
- **不在任一计时器内**:gather/输出、计时 Allreduce 本身、forces/torques/cc 的 assign
- step 0 恒报 0;区间均值,整数除法(亚毫秒步进位 0——小算例正常);每输出步一次 `MPI_Allreduce`(MAX, LONG_LONG, 2 元)全体 rank 无条件调用,print 仅 rank0
- 动机:混合核平台上区分"核慢"(ms_max 整体抬升)与"算法慢"(phys/comm 结构变化)

## (d) 稀疏场景 transform 取舍

旧代码包围球预检先行,无重叠对时零次 `transform_tris`;新代码每步无条件对全部粒子(含 MPI ghost)各变换一次。自由落体窗实测净收益仍为正(配对 A/B 见 (b)),但 perf 归因时注意:稀疏 + 大 N 时该 O(N·tris)/步是纯增项;Phase B 的 arena 复用可摊销。

## 延后项(Phase B 台账)

2026-08-19 晚间批次 `feat/perf-deferred-cleanup` 已落地 4 项(见 `docs/superpowers/plans/2026-08-19-perf-deferred-cleanup.md`;全部 VTK 字节级不变,MPI V1 n=2 ghost 路径同样字节 PASS):

- ~~MPI 每步 combined 数组整体拷贝 → 索引/指针视图~~ ✅ `broadphase_pairs(local, ghosts*, out)` 双数组接口,组合索引视图,零拷贝
- ~~world_tris 每步 vector<vector> 分配 → arena~~ ✅ `transform_tris_into` + 步循环外持久缓冲,稳态零分配
- ~~末步 `migrate_particles` guard~~ ✅ 与 ghost 交换同一守卫(`step + 1 < cfg.steps`)
- ~~9 参 `pp_contact_pair` 兼容 shim~~ ✅ 已裁,零调用方核实
- broadphase 每步 unordered_map/向量分配 → 复用(**仍开放**;现驱动步循环内仅存的每步堆分配)
- telemetry=1 + split_contacts=1 组合未显式门控过(**仍开放**;门内逐字旧代码,风险近零)

序保持单测已 ctest 化(`tests/test_broadphase.cpp`,ctest 8/8;O(N²) 参考逐元素对照 + 变异验证 282 处断言咬合)。

## 验证基线(全部通过)

ctest 7/7(终审在 HEAD 实跑复核);全部字节 parity 门见"硬约束兑现";GPU 零改动(broadphase.cpp 仅进 zdem_core 宿主库)。
