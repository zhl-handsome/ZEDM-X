# CPU/MPI 性能优化实施计划(2026-08-19 扩展性测试产出)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** 落地扩展性测试报告(scratch/mpi_scaling_report.md)识别的第一、二批优化:CPU 金参考遥测开关(实测 2.87× 口径差)、per-pair 网格变换缓存(~10% 接触相)、两个一行修、共享 broadphase 空间哈希(O(N²)→O(N),#1 算法瓶颈)、MPI 每 rank 计时(测量基础设施)。

**Architecture:** 全部改动落在共享层(src/physics/ + src/host/config_io)与两个驱动(main.cpp / mpi_main.cpp),GPU 不动。**硬约束:所有物理路径改动必须保持 VTK 输出字节级不变**(纯函数外提 + 枚举序保持),用既有 V0 基线门控。

**Tech Stack:** C++17 (MSVC),MPI C 绑定,既有 ctest + scratch/gpu_compare.py 字节对照。

## Global Constraints

- `zdem_cpu`/`zdem_mpi` 的 **VTK 帧输出字节级不变**(除明确声明的 stdout 日志字段);门控:`gpu_compare.py` 对照既有基线 `scratch/out_wall64_cpu`、`scratch/out_pile12_cpu64` PASS 且 |dpos|=0、|dE|=0
- ctest 7/7 不变绿不合并(wall_multiloop 断言 `wall_contact:...n_cross=` 调试行与 `step=0 contacts=1`,**不依赖** Route-B 遥测——已核实 tests/smoke/check_wall_multiloop_regression.py:31-45)
- stdout 日志行的**字段集合与顺序不变**(字段值可以为 0);mpi_smoke 只 grep `step=0 contacts=` 与 `done`
- 新 config key 旧可执行静默忽略;MPI 驱动不消费的新 key 不报错
- 新代码注释 English/ASCII(MSVC C4819)
- GPU target 与代码零改动

---

### Task 1: route_b_telemetry 配置开关(默认关)+ CPU 遥测门控

**Files:**
- Modify: `src/host/config_io.hpp`(SimConfig 加 `int route_b_telemetry = 0;`)
- Modify: `src/host/config_io.cpp`(解析 `route_b_telemetry = 0|1`,模式照抄 mpi_margin 分支)
- Modify: `src/main.cpp` 1091-1248(pp 对循环)
- Test: 既有 ctest + V0 字节对照 + 计时

**改动内容**(遥测块仅喂 stdout 计数器 1288-1298 与计时器,VTK/txt 不消费——已核实):

main.cpp pp 对循环内(包围球重叠的每对):
1. 1102-1186 的块(transform 后的 planes/AABB/tol/tri-tri → `segments`)包进 `if (cfg.route_b_telemetry)`,`segments` 声明留在门外(空向量)
2. 1197 的 `pp_contact_pair` 调用与 n_inc>0 累加块**不动**(1199-1207 的 `continue` 保留)
3. 1209-1247 的 loop-tracing 段:入口加 `if (!cfg.route_b_telemetry) continue;`(在 `segments.empty()` 判断之前——遥测关时 segments 恒空,但显式守卫更清晰)
4. 计数器 `total_segments/total_loops/total_comps/tri_checks/t_tri/t_rebuild` 声明不动,关时输出 0(日志行格式 1288-1298 逐字段保持)

**门控:**
```bash
cmake --build build --config Release 2>&1 | grep -iE "error"   # 无输出
cd build && ctest -C Release 2>&1 | grep -E "passed|failed"    # 7/7
# V0 字节对照(遥测关,VTK 必须与旧基线逐字节一致):
sed 's|output_dir = scratch/out_wall64_cpu|output_dir = scratch/opt_wall64_out|' scratch/wall64_cpu.txt > scratch/opt_wall64.txt
build/Release/zdem_cpu --config scratch/opt_wall64.txt | tail -1
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_wall64_cpu scratch/opt_wall64_out 1 | tail -1   # PASS 0.000e+00
# 遥测开(回归旧行为,stdout 计数非零):
sed -e 's|output_dir = scratch/out_wall64_cpu|output_dir = scratch/opt_wall64_tel_out|' -e 's/$/\nroute_b_telemetry = 1/' scratch/wall64_cpu.txt > scratch/opt_wall64_tel.txt
build/Release/zdem_cpu --config scratch/opt_wall64_tel.txt | grep -m1 "tri_checks="   # tri_checks>0
# 计时(pile256 自由落体 8k 窗,关 vs 开):
sed -e 's|output_dir = scratch/scal_ff_out|output_dir = scratch/opt_ff_off|' scratch/scal_ff.txt > scratch/opt_ff_off.txt
sed -e 's|output_dir = scratch/scal_ff_out|output_dir = scratch/opt_ff_on|' -e 's/$/\nroute_b_telemetry = 1/' scratch/scal_ff.txt > scratch/opt_ff_on.txt
time build/Release/zdem_cpu --config scratch/opt_ff_off.txt   # vs time ... opt_ff_on.txt,记录两个墙钟
```

- [ ] 实现 → [ ] ctest 7/7 → [ ] V0 字节 PASS → [ ] 遥测开 stdout 非零 → [ ] 计时记录(off 预期显著快于 on)→ [ ] commit `perf(cpu): route_b_telemetry config gate, default off (2.9x CPU reference speedup)`

### Task 2: per-pair transform_tris 缓存(共享层重载 + 两驱动接入)

**Files:**
- Modify: `src/physics/pp_contact.hpp` / `pp_contact.cpp`
- Modify: `src/main.cpp`(pp 对循环)
- Modify: `src/mpi/mpi_main.cpp`(本地/ghost 各一份缓存)

**Interfaces:**
```cpp
// pp_contact.hpp 追加(9 参版保留,内部改为委托):
// trisA/trisB = transform_tris(ma, pa.tf) / transform_tris(mb, pb.tf) 的调用方缓存
int pp_contact_pair(const Particle& pa, const Particle& pb,
                    const std::vector<std::array<Vec3, 3>>& trisA,
                    const std::vector<std::array<Vec3, 3>>& trisB,
                    Vec3& f_i, Vec3& t_i, Vec3& f_j, Vec3& t_j,
                    double tangential_damping);
```
pp_contact.cpp:现 73-74 行的两行 transform 删除,9 参版 = `return pp_contact_pair(pa, pb, transform_tris(ma,pa.tf), transform_tris(mb,pb.tf), ...)`(临时物生命周期注意:先存具名局部再传)。**其余代码零改动。**

**调用方:**
- CPU main:每步在 pp 循环前 `std::vector<std::vector<std::array<Vec3,3>>> world_tris(N)`,一次 `transform_tris(meshes[particles[i].mesh_index], particles[i].tf)`;pp 循环内 1105-1106(遥测块,若开)与 pp_contact_pair 都用 `world_tris[i]/[j]`
- MPI mpi_main:每步 `world_tris_local` + `world_tris_ghost` 同样缓存;本地-本地与本地-ghost 调用都用缓存的重载(ghost 粒子包围球预检保持)
- wall_contact_particle 不动(wall tris 每 wall 仅 2 三角形,忽略)

**数值不变性论证(写给审查者):** transform_tris 是纯函数,同输入同输出;缓存只是把每对的重复计算外提,per-pair 数值与累积序完全不变 → VTK 必须字节级一致。

**门控:** build → ctest 7/7 → CPU V0 字节 PASS(Task 1 的 opt_wall64 流程重跑)→ **MPI V0 parity 重跑**(mpi_main 被动):
```bash
mpiexec -n 1 build/Release/zdem_mpi.exe --config scratch/mpi_v0_wall64.txt | tail -1
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_wall64_cpu scratch/mpi_v0_wall_out 1 | tail -1
mpiexec -n 1 build/Release/zdem_mpi.exe --config scratch/mpi_v0_pile12.txt | tail -1
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_pile12_cpu64 scratch/mpi_v0_p12_out 12 | tail -1
```
计时:`time mpiexec -n 1 ... scratch/pile256.txt`(24k 全窗)记录新旧(旧值 377.7s 在案)。

- [ ] 实现 → [ ] ctest → [ ] CPU+MPI V0 字节 PASS → [ ] 计时记录 → [ ] commit `perf(shared): cache per-particle world triangles per step; pp_contact_pair tris overload`

### Task 3: 末步 ghost 交换 guard + gids/local 等长检查

**Files:**
- Modify: `src/mpi/mpi_main.cpp`(~194-196 循环尾):`exchange_ghosts(...)` 包 `if (step + 1 < cfg.steps)`
- Modify: `src/mpi/ghost.cpp`(exchange_ghosts 入口):gids.size() != local.size() 时 fprintf(stderr) + MPI_Abort(风格照抄文件内既有 die/abort 模式)

**门控:** build → ctest 7/7 → MPI V0 parity 字节 PASS(同 Task 2 流程)→ commit `fix(mpi): skip final-step ghost exchange; guard gids/local length mismatch`

---

### Task 4: 共享 broadphase 空间哈希(O(N²) → O(N),枚举序保持)

**Files:**
- Create: `src/physics/broadphase.hpp` / `broadphase.cpp`
- Modify: `src/main.cpp`(pp 对循环 1091-1100 替换)、`src/mpi/mpi_main.cpp`(本地-本地 + 本地-ghost 枚举替换)
- Modify: `CMakeLists.txt`(zdem_core 源列表加 broadphase.cpp)

**Interfaces:**
```cpp
// broadphase.hpp
// Bounding-sphere overlap candidate pairs over ps[0..n) (combined local+ghost
// array in MPI). Cell = 2*max(radius) so any overlapping pair sits in adjacent
// cells; 27-neighborhood scan with j>i dedup, then LEXICOGRAPHIC SORT so the
// output order equals the O(N^2) i<j double loop EXACTLY -- downstream force
// accumulation order and bit-level output are unchanged. Pairs with first
// index >= n_local are dropped (ghost-ghost never computed, matches drivers).
void broadphase_pairs(const std::vector<Particle>& ps, int n_local,
                      std::vector<std::pair<int, int>>& out);
```
实现:unordered_map<array<int64_t,3>(或 hash 组合), vector<int>>;cell 索引 floor((pos-lo)/cell),lo 用扫描时的 min pos(或 0 基,负坐标注意 floor 语义——用 std::floor);遍历每个 i,扫 27 邻 cell 的 occupant j:j>i 且 dist2 ≤ (ri+rj)² 才收;最后 `std::sort(out)`(pair 默认字典序 = 所需枚举序)。
**CPU 消费:** `broadphase_pairs(particles, N, pairs)` 一次,循环改 `for (auto [i,j] : pairs)`,循环体 1093 起原样(包围球重复判断保留无害,或删——删则循环体从 Particle& 引用起)。
**MPI 消费:** 组合数组(local 后接 ghost)一次调用 n_local=local.size();本地-ghost 半循环的包围球预检逻辑保持(或直接信 broadphase——**保守:保留现有预检行,只换枚举来源**)。

**数值不变性:** 同 pair 集(证明:重叠 ⇒ 距离 ≤ 2max_r ≤ cell 长 ⇒ cell 坐标每轴差 ≤1 ⇒ 27 邻域覆盖)+ 同枚举序(sort 后字典序 = i<j 双循环序)→ 字节级一致。

**门控:** build → ctest 7/7 → CPU V0 字节 PASS → MPI V0 + **V1(n=2)** 字节 PASS → pile64 16k 全程 CPU vs 基线 `scratch/out_pile64_cpu64`(若在;否则跑 2k 窗自建基线)字节 PASS → 计时:
```bash
# 自由落体 8k 窗(旧值 CPU 45ms/步级、MPI n=1 4.41ms/步在案):
time build/Release/zdem_cpu --config scratch/opt_ff_off.txt
time mpiexec -n 1 build/Release/zdem_mpi.exe --config scratch/scal_ff.txt
```
- [ ] 实现 → [ ] ctest → [ ] V0/V1 字节 PASS → [ ] 计时 → [ ] commit `perf(shared): spatial-hash broadphase, enumeration order preserved (O(N^2) -> O(N))`

### Task 5: MPI 每 rank 计时(ms_max / comm_ms_max 日志字段)

**Files:**
- Modify: `src/mpi/mpi_main.cpp`

每 rank 两个 steady_clock 累积器:phys_ns(力+积分段)、comm_ns(迁移+ghost 段);输出步:`MPI_Allreduce`(MAX,lua... MPI_LONG_LONG,全体 rank 无条件调用)→ rank0 在既有 `step=... contacts=...` 行尾追加 ` ms_max=<int> comm_ms_max=<int>`(interval 步均值,整除)。**集合卫生:Allreduce 全体调用,print 仅 rank0**(Task 4 教训)。
门控:build → ctest 7/7(mpi_smoke 的 grep 仍命中)→ `mpiexec -n 2 ... mpi_v1_pile12.txt` 日志行含新字段 → V0/V1 parity 字节 PASS(VTK 不受日志影响,跑一遍确认)→ commit `feat(mpi): per-rank step timing in output log (ms_max, comm_ms_max)`

---

## Self-Review

- Spec 覆盖:报告建议 1/2/3 → Task 1/2/3;建议 4 → Task 4;建议 7 → Task 5 ✓
- 占位符:无(全部改动含精确行号与门控命令)
- 类型一致:pp_contact_pair 重载签名 Task 2 内自洽;broadphase_pairs 签名 Task 4 内自洽
- 已知风险内联:Task 2 临时物生命周期;Task 4 负坐标 floor;Task 5 集合卫生
