# 性能延后项清理实施计划(Phase B 前置)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** 落地 2026-08-19 优化批次终审延后清单:末步迁移 guard、9 参 shim 裁剪、world_tris 持久缓冲(零每步分配)、broadphase 双数组接口(消灭 MPI combined 拷贝)、broadphase 序保持单测 ctest 化(把分支最关键不变量从 scratch 手工门变成可重放测试)。

**Architecture:** 全部改动数值中性(分配/接口层,零表达式变化),**VTK 字节级不变**继续作为硬门。测试从 7 项增至 8 项。

**Tech Stack:** C++17 (MSVC),既有 ctest + scratch/gpu_compare.py。

## Global Constraints

- VTK 帧输出字节级不变;门控:`gpu_compare.py` 对照 `scratch/out_wall64_cpu` / `scratch/out_pile12_cpu64` PASS 且 0.000e+00
- ctest 收尾时 8/8(新增 test_broadphase);既有 7 项不破坏
- stdout 日志字段集不变
- 新代码注释 English/ASCII;GPU 代码零改动(dem_state.cu 的 transform_tris 按值调用保持)
- 集合通信卫生:任何新分支条件包住的 MPI 调用,条件必须各 rank 一致

---

### Task 1: 末步迁移 guard + 9 参 shim 裁剪

**Files:**
- Modify: `src/mpi/mpi_main.cpp`(循环尾:现有 `if (step + 1 < cfg.steps) exchange_ghosts(...)` 的 if 扩为同时包住其前的 `migrate_particles(d, local, gids);`——两行同一守卫)
- Modify: `src/physics/pp_contact.hpp`(删 19-22 行 9 参声明,注释块保留在 11 参声明上)
- Modify: `src/physics/pp_contact.cpp`(删 69-78 行 9 参定义)

**事实(已核实):** 9 参版零调用方(main.cpp:1225、mpi_main.cpp:173/196 全走 11 参 tris 重载);`cfg.steps` 各 rank 一致 → guard 集合安全;循环后无 local/ghost 消费者(仅打印 done)。

**门控:** build 无 error → ctest 7/7 → MPI V0 parity 字节 PASS(wall64 + pile12,config scratch/mpi_v0_wall64.txt / mpi_v0_pile12.txt)→ CPU V0 parity 字节 PASS(opt_wall64 流程)→ commit `perf(mpi): skip final-step migration; prune caller-less 9-arg pp_contact_pair shim`

### Task 2: transform_tris_into + world_tris 持久缓冲

**Files:**
- Modify: `src/geometry/mesh_build.hpp` / `.cpp`(加 `void transform_tris_into(const Mesh&, const Transform&, std::vector<Triangle>& out)`——resize 到 tris 数 + 逐元素与现 by-value 版**同一表达式**填充;by-value 版改为委托 `_into`)
- Modify: `src/main.cpp`(~1095-1102:world_tris 提到步循环外持久化;每步 `resize(N)` 后逐粒子 `transform_tris_into(..., world_tris[p])`——vector 容量跨步复用,稳态零分配)
- Modify: `src/mpi/mpi_main.cpp`(~128-140:wtris_local/wtris_ghost 同样持久化,每步在迁移/ghost 重建后 resize+refill)

**数值不变性:** 每元素表达式逐字相同(委托重构),仅消除堆分配;容量复用不改值。

**门控:** build → ctest 7/7 → CPU+MPI V0 字节 PASS → commit `perf(shared): persistent world-tris buffers via transform_tris_into (zero steady-state alloc)`

### Task 3: broadphase 双数组接口(消灭 MPI combined 拷贝)

**Files:**
- Modify: `src/physics/broadphase.hpp` / `.cpp`
- Modify: `src/main.cpp` / `src/mpi/mpi_main.cpp`(调用点;删 combined 数组构建)

**Interfaces:**
```cpp
// 新签名(旧 combined 签名删除,调用方仅两个驱动):
// local 后可选 ghosts;输出 pair 用组合索引(ghost j -> local.size()+gj),
// first 恒 < local.size(),字典序 = 旧 [locals;ghosts] 组合数组序 —— 不变量不变。
// cell 插入遍历 local 后 ghosts;外层 pair 发射只跑 local i → ghost-ghost 结构性不生成。
void broadphase_pairs(const std::vector<Particle>& local,
                      const std::vector<Particle>* ghosts,
                      std::vector<std::pair<int, int>>& out);
```
内部距离检查处按组合索引取 Particle(j < n_local ? local[j] : (*ghosts)[j - n_local])。cell 尺寸 = 2·max(全 local+ghost radius)(与现实现同,每次调用重算)。CPU 传 ghosts=nullptr。

**门控:** build → ctest 7/7 → CPU+MPI V0 字节 PASS → **MPI V1(n=2)字节 PASS**(ghost 路径接口改动)→ commit `perf(shared): two-array broadphase interface, no combined copy`

### Task 4: broadphase 序保持单测(ctest 化)

**Files:**
- Create: `tests/test_broadphase.cpp`(框架仿 tests/test_geometry.cpp:独立 main,断言失败打印并返回非零)
- Modify: `CMakeLists.txt`:
```cmake
add_executable(test_broadphase
    tests/test_broadphase.cpp
)
target_link_libraries(test_broadphase PRIVATE zdem_core)
add_test(NAME test_broadphase COMMAND test_broadphase)
```

**测试内容**(对 Task 3 的最终 API;确定性 seed,不依赖随机设备):
1. **序保持核心断言**:多组布局(N ∈ {0,1,2,5,40},单粒径与混合粒径、均匀与聚集、含负坐标、含重合位置),参考实现 = O(N²) i<j 双循环 + **逐字相同**的距离谓词 `if (dot(dpos,dpos) > rsum*rsum) continue;`;断言 `broadphase_pairs` 输出与参考 pair 列表**逐元素相等**(集合与字典序双重锁定)
2. **ghost 语义**:local=N_l、ghosts=N_g;断言所有 pair first < N_l;local-ghost 对组合索引正确(与参考 combined 循环对照);ghost-ghost 零出现
3. **退化**:N<2 → 空;全零半径 → 空(现实现文档化行为);ghosts=nullptr 路径(CPU 形态)
4. **相切边界**:恰好 dist == rsum 的对(构造对称两点)必须保留(非严格 `>` 语义)

**门控:** build → `ctest -C Release -R broadphase` PASS → 全量 ctest 8/8 → commit `test(shared): broadphase order-preservation unit test (pair set + lexicographic order vs O(N^2) reference)`

---

## Self-Review

- 覆盖:终审延后清单 5 项全落(guard=Task 1、shim=Task 1、arena=Task 2、索引视图=Task 3、单测=Task 4)
- 顺序依赖:Task 4 的 API 依赖 Task 3 完成;Task 1/2 独立
- 数值风险点内联:Task 2 委托重构的表达式同一性;Task 3 组合索引取 Particle 的分支
