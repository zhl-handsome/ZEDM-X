# ZDEM-X MPI CPU 通信层 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 zdem_cpu 的 v8 接触模型并行化到 MPI(静态 brick 分解 + ghost 交换 + 质心迁移,Newton-off),FP64 容差 parity 验证,为 MPI×多 GPU 铺通信层。

**Architecture:** 新 `zdem_mpi` 可执行(共享 `zdem_core` + 新 `src/mpi/` 通信层 + 从 main.cpp 抽出的共享物理/构建层 `src/host/sim_build`、`src/physics/`);每 rank 持有 brick 内粒子,ghost 只读镜像(v8 无接触状态),跨域 pair 双算免反向通信;输出 Gatherv 到 rank0 写全局 VTK(格式与单进程逐字段一致)。

**Tech Stack:** MPI-3 **C 绑定**(MS-MPI 本地验证 / MPICH/OpenMPI 集群)、CMake `find_package(MPI)`、C++17、Python 对照工具(复用 `scratch/gpu_compare.py`)。

**设计文档:** `docs/superpowers/specs/2026-08-18-mpi-dem-design.md`

## Global Constraints

- `zdem_cpu` / `zdem_gpu` 行为不变:每任务后 `cd build && ctest -C Release` 必须 5/5(或 6/6 含 gpu_smoke)通过
- 物理公式与 CPU v8 **逐条对应**——通过共享抽取(Task 2)结构性保证,不复制粘贴公式
- 只用 MPI C 绑定(`mpi.h`),不用已废弃的 C++ 绑定;打包用 `MPI_BYTE` + 整结构 memcpy(**同构集群假设**:同一可执行文件、同一 ABI,异构集群范围外)
- Newton-off:每 rank 只累加本地粒子的力;ghost 深度 = `2 × max(外接球)`
- config 兼容:新增 `mpi_box` / `mpi_margin` 可选 key;现有 config_io 对未知 key 静默忽略(已确认现状),旧可执行不受影响
- MPI 版为 CPU 谱系**纯 double**,无 FP32 变体
- 主循环序:清力 → pp(本地×(本地∪ghost))→ 墙(本地×全墙)→ 积分 → **输出** → 迁移 → ghost 重建(为下一步)。输出在迁移前:force/cc 数组不随迁移重排,物理状态等价
- 测试命令 `D:/ProgramData/anaconda3/python.exe`;MPI 结果目录 `scratch/out_mpi_*`;`mpiexec` 在 `D:\Program Files\Microsoft MPI\Bin\mpiexec.exe`(已在 PATH)
- 分支:`feat/mpi-cpu-comm`,每任务一 commit;完成后 ff merge 回 main、删分支、烟雾(固定流程,不呈现选项)
- Parity 容差:位置 < 1e-6 m,能量相对 < 1e-9(不逐位;pair 枚举序不同 → ulp 级差)

---

### Task 1: config 增加 mpi_box / mpi_margin 可选字段

**Files:**
- Modify: `src/host/config_io.hpp`(SimConfig 增字段)
- Modify: `src/host/config_io.cpp`(解析两个 key)
- Test: `scratch/mpi_keys_check.txt`(临时 config,手工验证)

**Interfaces:**
- Produces(后续所有任务依赖):
```cpp
// config_io.hpp —— SimConfig 内追加:
std::array<double, 6> mpi_box{0, 0, 0, 0, 0, 0};  // xmin xmax ymin ymax zmin zmax
bool has_mpi_box = false;
double mpi_margin = -1.0;   // <0 = 未设置(用默认 5*max_radius)
```

- [ ] **Step 1: 在 config_io.hpp 的 SimConfig 末尾(`std::vector<WallInit> wall_inits;` 之后)加上述三个字段**(文件顶部补 `#include <array>`)

- [ ] **Step 2: config_io.cpp 解析链(`tangential_damping` 分支后)追加**

```cpp
        } else if (key == "mpi_box") {
            iss >> cfg.mpi_box[0] >> cfg.mpi_box[1] >> cfg.mpi_box[2]
                >> cfg.mpi_box[3] >> cfg.mpi_box[4] >> cfg.mpi_box[5];
            cfg.has_mpi_box = true;
        } else if (key == "mpi_margin") {
            iss >> cfg.mpi_margin;
        }
```

- [ ] **Step 3: 构建 + 新旧 config 双验证**

```bash
cmake --build build --config Release 2>&1 | grep -iE "error"        # 无输出
# 旧可执行容忍新 key(未知 key 静默忽略是现状):
printf 'steps = 1\ndt = 1e-5\nmpi_box = -1 1 -1 1 0 2\nmpi_margin = 0.5\n\nparticle\nstl = geometry/low-poly-banbana.stl\npos = 0 0 0.1\nscale = 0.001\ndensity = 2500\nend_particle\n' > scratch/mpi_keys_check.txt
build/Release/zdem_cpu.exe --config scratch/mpi_keys_check.txt | tail -1   # done(不报错)
cd build && ctest -C Release 2>&1 | grep -E "tests passed|failed"
```
Expected: `done`;测试全过(5/5 或 6/6)

- [ ] **Step 4: Commit**

```bash
git add src/host/config_io.hpp src/host/config_io.cpp && git commit -m "feat(mpi): optional mpi_box/mpi_margin config keys"
```

---

### Task 2: 抽取共享构建层与物理层(zdem_core)

**Files:**
- Create: `src/host/sim_build.hpp`, `src/host/sim_build.cpp`(粒子/墙/网格注册表构建,自 main.cpp 1129-1257 抽取;`struct Wall` 自 main.cpp:30 移入)
- Create: `src/physics/pp_contact.hpp`, `src/physics/pp_contact.cpp`(单 pair 接触力,自 main.cpp 1370-1490 抽取)
- Create: `src/physics/wall_contact.hpp`, `src/physics/wall_contact.cpp`(单粒子墙接触,自 main.cpp 1540-1690 抽取)
- Create: `src/physics/integrate.hpp`, `src/physics/integrate.cpp`(自 main.cpp 1034-1102 原样)
- Modify: `src/main.cpp`(删除被抽走的段,改调用;主循环逻辑不变)
- Modify: `CMakeLists.txt`(zdem_core 增 4 个 .cpp)

**Interfaces:**
- Produces(Task 3-7 依赖):
```cpp
// sim_build.hpp
struct Wall {                       // 自 main.cpp:30 原样移动
    Mesh mesh;
    Transform tf;
    double mu = 0.5, restitution = 0.5;
    std::vector<Vec3> tri_normals;
};
struct SimBuild {
    std::vector<Mesh> meshes;          // 去重网格注册表(build_sim 内按 stl 路径去重,同 main.cpp 现逻辑)
    std::vector<Particle> particles;   // config 序;gid == 下标
    std::vector<Wall> walls;
    double max_radius = 0.0;           // max(p.radius),ghost 深度/margin 用
};
bool build_sim(const SimConfig& cfg, SimBuild& out, std::string& err);

// physics/pp_contact.hpp
// 单 pair:v8 内含检测 + per-vertex 罚力。返回内含顶点总数 n_inc(0 = 无接触)。
// f_i/t_i/f_j/t_j 为输出(先置零再累加):owner=i 侧与 other=j 侧两半。
// 调用方取舍:CPU 双侧都加;MPI Newton-off 只取 own 半。
int pp_contact_pair(const Particle& pa, const Particle& pb,
                    const Mesh& ma, const Mesh& mb,
                    Vec3& f_i, Vec3& t_i, Vec3& f_j, Vec3& t_j,
                    double tangential_damping);

// physics/wall_contact.hpp
// 单粒子 × 全部墙:分组/穿透/罚力自 main.cpp 1540-1690 原样。
// 返回本粒子产生接触的墙组数(contact_counts 语义,同 CPU)。
int wall_contact_particle(const Particle& p, const Mesh& pmesh,
                          const std::vector<Wall>& walls,
                          double tangential_damping,
                          Vec3& f, Vec3& t);

// physics/integrate.hpp
void integrate_particle(Particle& p, const Vec3& force, const Vec3& torque,
                        const Vec3& gravity, double dt);
```

- [ ] **Step 1: 抽 integrate(最简单,先建立流程)**

`src/physics/integrate.{hpp,cpp}`:把 main.cpp 1034-1102 的 `integrate_particle` **逐字**移入(去 `static`),hpp 只放声明。main.cpp 该处改 `#include "physics/integrate.hpp"` 并删除原函数体。

- [ ] **Step 2: 抽 pp_contact**

`pp_contact.cpp`:把 main.cpp 1370-1490 段的双层循环**体内**逻辑(内含扫描 A→B/B→A、`incA/incA_cp/incA_d` 收集、`apply_penalty` lambda 1452-1485)移入 `pp_contact_pair`,改造点:
- lambda `apply_penalty` 变 .cpp 内 `static` 函数,签名 `(const Vec3& wv, const Vec3& cp, double d, const Particle& own, const Particle& oth, const Mesh& /*mesh 常数已在调用前算好*/, int n_inc, double tangential_damping, Vec3& f_own, Vec3& t_own, Vec3& f_oth, Vec3& t_oth)` —— 原 lambda 捕获的 `k_hr/kt_r/m_eff/m_eff_v/e_r/ct_v/mu_r` 常数块在 pair 函数内计算(公式逐条照搬 1420-1450)
- 函数开头 `f_i=t_i=f_j=t_j=Vec3{}`;A→B 命中调 apply_penalty(own=a, oth=b),B→A 命中调 apply_penalty(own=b, oth=a)(**保持 CPU 的调用序:先全部 incA 再全部 incB,1430-1490 的循环序原样**)
- 返回 `n_inc = incA.size() + incB.size()`
main.cpp 主循环改为:
```cpp
Vec3 f_i, t_i, f_j, t_j;
if (pp_contact_pair(pa, pb, meshes[pa.mesh_index], meshes[pb.mesh_index],
                    f_i, t_i, f_j, t_j, cfg.tangential_damping) > 0) {
    forces[i] += f_i; torques[i] += t_i;      // 加法序与原 apply_penalty 内联累积一致:
    forces[j] += f_j; torques[j] += t_j;      // CPU 原序是逐顶点交錯加两侧,此处改为整块加 —— 见 Step 5 逐位验证,
    contact_counts[i]++; contact_counts[j]++; // 若逐位验证失败则改回逐顶点接口(见 Step 5 note)
    total_contacts++;
}
```
**注意**:CPU 原实现里每个顶点的力是交错累加到 forces[i]/forces[j] 的(i 侧加一点、j 侧加一点、再 i 侧…);整块累加改变了浮点加序 → **可能产生 ulp 级差**。若 Step 5 的逐位对照失败,备选方案:pp_contact_pair 改为回调/输出逐顶点增量列表 `std::vector<std::array<Vec3,2>>`(保持交错序),CPU 调用方按原序回放。默认先试整块(代码最简),逐位验证不过再切换。

- [ ] **Step 3: 抽 wall_contact 与 sim_build**

`wall_contact.cpp`:main.cpp 1540-1690 的 `for i` 循环体 + `for w` 循环体(WallPlaneGroup 分组、穿透检测、罚力)移入 `wall_contact_particle`,原样保留内部结构(含 `WallPlaneGroup` 局部 struct 与 quantized 分组);捕获变量 `cfg.tangential_damping` 参数化。返回值 = 有 any_force 的墙组数(同 CPU contact_counts 语义)。
`sim_build.{hpp,cpp}`:main.cpp 1129-1257 的粒子构建(load_stl→scale→build_mesh→compute_mass_properties→radius/equiv_radius→mesh 注册表按 stl 路径去重)与墙构建段移入 `build_sim`;`struct Wall`(main.cpp:30-37)移入 sim_build.hpp;`max_radius` = max(p.radius)。main.cpp include 并调用,删除原段。
CMakeLists `add_library(zdem_core ...)` 追加 `src/host/sim_build.cpp src/physics/pp_contact.cpp src/physics/wall_contact.cpp src/physics/integrate.cpp`。

- [ ] **Step 4: 构建 + ctest**

```bash
cmake --build build --config Release 2>&1 | grep -iE "error"    # 无输出
cd build && ctest -C Release 2>&1 | grep -E "tests passed|failed"
```
Expected: 全过

- [ ] **Step 5: CPU 逐位回归(抽取不改行为的硬证据)**

```bash
cd .. && build/Release/zdem_cpu.exe --config scratch/wall64_cpu.txt >/dev/null 2>&1
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_wall64_cpu scratch/out_wall64_cpu 1 | tail -2
# 自对照恒过;真正的基线:与 git stash 前的旧版输出对照 —— 用下法:
build/Release/zdem_cpu.exe --config scratch/mpi_v0_base.txt   # 见 Task 4 Step 1 生成;先跳到那里建好再回来
```
正确的逐位验证流程:
```bash
# 1) 用抽取前代码跑基线(现在还没改?已改 —— 用 git):
git stash && cmake --build build --config Release >/dev/null 2>&1 && \
  build/Release/zdem_cpu.exe --config scratch/mpi_v0_base.txt >/dev/null 2>&1 && git stash pop && \
  cmake --build build --config Release >/dev/null 2>&1 && \
  build/Release/zdem_cpu.exe --config scratch/mpi_v0_new.txt >/dev/null 2>&1 && \
  D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/mpi_v0_base_out scratch/mpi_v0_new_out 1 | tail -2
```
(mpi_v0_base.txt 的 output_dir=scratch/mpi_v0_base_out,mpi_v0_new.txt 同 config 改 output_dir=scratch/mpi_v0_new_out;config 内容 = wall64_cpu.txt 改 steps=40000 已是、只改目录名)
Expected: `RESULT: PASS`,且 max|dpos| = 0(**逐位**;若非 0 但 <1e-6,记录数值并切换 Step 2 的备选接口重验)

- [ ] **Step 6: Commit**

```bash
git add -A && git commit -m "refactor: extract shared sim build + contact physics into zdem_core"
```

---

### Task 3: CMake MPI target + decomp + 初始分布(0 步骨架)

**Files:**
- Create: `src/mpi/decomp.hpp`, `src/mpi/decomp.cpp`
- Create: `src/mpi/mpi_main.cpp`(骨架:Init → config → build_sim → 过滤本 rank → 打印 → Finalize)
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: Task 1 `SimConfig.mpi_box/mpi_margin`、Task 2 `SimBuild/build_sim`
- Produces:
```cpp
// decomp.hpp
#include <mpi.h>
struct Decomp {
    double box_lo[3], box_hi[3];     // 全局盒
    int dims[3];                     // px py pz
    int coords[3];                   // 本 rank 笛卡尔坐标
    MPI_Comm cart = MPI_COMM_NULL;   // 3D 拓扑(reorder=false,rank 与 COMM_WORLD 一致)
    int rank = 0, nprocs = 1;
    double sub_lo[3], sub_hi[3];     // 本 brick
    double ghost_depth = 0.0;        // 2*max_radius
};
// cfg.mpi_box 优先;否则初始粒子包围盒每侧外扩 margin(默认 5*max_radius,可 mpi_margin 覆盖)
// dims 自动分解:贪心——从 {1,1,1} 起,重复把剩余进程数乘到"盒最长轴"对应的维度,
//   使 dims[a]*box_len[a] 尽量均衡;然后校验每轴 brick 厚 >= ghost_depth,
//   违反则该轴 dims 减半重试(仍违反且 dims==1 → 报错退出:盒子太小/进程太多)
Decomp make_decomp(const SimConfig& cfg, const SimBuild& sim);
int owner_rank(const Decomp& d, const Vec3& pos);   // pos→全局 rank(floor 到 brick,越界 clamp 后报错由调用方)
bool in_sub(const Decomp& d, const Vec3& pos);      // 本 rank 归属(边界=左闭右开,与 owner_rank 同 floor 语义)
```

- [ ] **Step 1: CMake 挂 MPI target**

CMakeLists 末尾(gpu_smoke test 之后):
```cmake
find_package(MPI REQUIRED)   # 只影响下面这个 target;无 MPI 环境时上面全部照常
add_executable(zdem_mpi
    src/mpi/mpi_main.cpp
    src/mpi/decomp.cpp
)
target_link_libraries(zdem_mpi PRIVATE zdem_core MPI::MPI_C)
```

- [ ] **Step 2: 实现 decomp.cpp**

`make_decomp`:`MPI_Dims_create` 不合需求(不看盒长比),按接口注释手写贪心分解;`MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods={0,0,0}, 0, &cart)`(非周期);`MPI_Cart_coords` 取 coords;brick 边界按盒长均分。`owner_rank`/`in_sub` 用同一表达式 `floor((pos[a]-lo[a])/box_len[a]*dims[a])` 截断到 [0,dims-1](越全局盒:报错——mpi_main 检测后 exit 1)。
**边界一致性关键**:owner 判定与 ghost/migrate 的邻居判定必须用同一 floor 表达式,写成 `decomp.cpp` 内 static helper `cell_index(pos, axis)`,两个公开函数共用。

- [ ] **Step 3: mpi_main 骨架**

```cpp
// src/mpi/mpi_main.cpp -- MPI DEM driver(v8,Newton-off,静态 brick)
#include <cstdio>
#include <string>
#include <vector>
#include "mpi.h"
#include "host/config_io.hpp"
#include "host/sim_build.hpp"
#include "mpi/decomp.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    // --config 解析(与 zdem_cpu 相同的 argv 约定)
    ... parse config_path ...
    SimConfig cfg; parse_config_file(config_path, cfg);
    SimBuild sim; std::string err;
    if (!build_sim(cfg, sim, err)) { if (rank==0) fprintf(stderr,...); MPI_Finalize(); return 1; }
    Decomp d = make_decomp(cfg, sim);
    // 初始分布:全 rank build 后本地过滤(gid = 原 config 下标,随粒子携带)
    std::vector<Particle> local; std::vector<int> gids;
    for (int i = 0; i < (int)sim.particles.size(); ++i)
        if (in_sub(d, sim.particles[i].tf.pos)) { local.push_back(sim.particles[i]); gids.push_back(i); }
    // 守恒断言:allreduce sum(local) == 全局 N
    long long loc = local.size(), tot = 0;
    MPI_Allreduce(&loc, &tot, 1, MPI_LONG_LONG, MPI_SUM, d.cart);
    if (tot != (long long)sim.particles.size()) { /* rank0 报错,Abort */ }
    if (d.rank == 0)
        std::printf("zdem_mpi: nprocs=%d dims=%dx%dx%d box=[%.3f %.3f]-[%.3f %.3f]-[%.3f %.3f] ghost=%.3f N=%lld\n",
                    d.nprocs, d.dims[0], d.dims[1], d.dims[2], ...);
    std::printf("  rank %d: nlocal=%zu gids=[%d..%d]\n", d.rank, local.size(),
                local.empty()?-1:gids.front(), local.empty()?-1:gids.back());
    MPI_Finalize(); return 0;
}
```
(steps 字段暂不消费;主循环 Task 4 加)

- [ ] **Step 4: 构建 + 1/2/3 进程冒烟**

```bash
cmake --build build --config Release 2>&1 | grep -iE "error"     # 无输出
build/Release/zdem_mpi.exe --config scratch/pile12_v8.txt | head -2          # 单进程直接跑(MPI_Init 本地)
mpiexec -n 2 build/Release/zdem_mpi.exe --config scratch/pile12_v8.txt
mpiexec -n 3 build/Release/zdem_mpi.exe --config scratch/pile12_v8.txt
```
Expected: 各 rank `nlocal` 之和 = 12;dims 合理(n=3 → 3×1×1);exit 0

- [ ] **Step 5: ctest(确认未破坏其它 target)+ Commit**

```bash
cd build && ctest -C Release 2>&1 | grep -E "tests passed|failed"
cd .. && git add -A && git commit -m "feat(mpi): cart decomp, global box, initial distribution skeleton"
```

---

### Task 4: 主循环(Newton-off)+ gather 输出 + V0 parity

**Files:**
- Create: `src/mpi/gather.hpp`, `src/mpi/gather.cpp`
- Create: `src/mpi/ghost.hpp`, `src/mpi/ghost.cpp`(本任务为**空实现占位**:容器与接口,exchange 为 no-op——n=1 不需要 ghost;Task 5 填实现)
- Modify: `src/mpi/mpi_main.cpp`(完整主循环)
- Test: `scratch/mpi_v0_wall64.txt`、`scratch/mpi_v0_pile12.txt`

**Interfaces:**
- Consumes: Task 2 全部物理接口、Task 3 Decomp
- Produces:
```cpp
// gather.hpp
struct FrameSnapshot {            // MPI_BYTE 打包单元(POD,同构集群假设)
    Particle p; int gid; int cc;
    Vec3 force, torque;
};
void gather_write_frame(const Decomp& d, const SimConfig& cfg, const SimBuild& sim,
                        const std::vector<Particle>& local, const std::vector<int>& gids,
                        const std::vector<int>& cc, const std::vector<Vec3>& forces,
                        const std::vector<Vec3>& torques, int step);
int  reduce_add(const Decomp& d, int local_val);                 // MPI_Reduce SUM → rank0 返回,他 rank 0
void log_rank_n(const Decomp& d, const std::vector<int>& local_sizes_hint);  // rank0 打印 rank_n=a,b,c(用 Allgatherv)

// ghost.hpp(本任务占位)
struct GhostLayer {
    std::vector<Particle> particles;
    std::vector<int> gids;
};
void exchange_ghosts(const Decomp& d, const std::vector<Particle>& local,
                     const std::vector<int>& gids, GhostLayer& out);   // Task 4: no-op(清空 out)
```
- `gather_write_frame` 实现:本地打包 `FrameSnapshot` 数组 → `MPI_Gatherv`(MPI_BYTE,counts/displs 由各 rank `local.size()` Allgather 得)→ rank0 按 gid `std::sort` → 组装 `particles/forces/torques/cc` 调 `write_vtk_particles`(文件名构造与 main.cpp 1695-1701 相同:`vtk_prefix_%06d.vtk`)。

- [ ] **Step 1: 建 V0 config**

```bash
sed -e 's|output_dir = scratch/out_wall64_cpu|output_dir = scratch/mpi_v0_wall_out|' scratch/wall64_cpu.txt > scratch/mpi_v0_wall64.txt
sed -e 's/^steps = 250000/steps = 2000/' -e 's|output_dir = scratch/out_pile12_v8|output_dir = scratch/mpi_v0_p12_out|' scratch/pile12_v8.txt > scratch/mpi_v0_pile12.txt
```

- [ ] **Step 2: 实现主循环(mpi_main.cpp)**

```cpp
// 循环体(在 Task 3 骨架 initial 过滤后):
GhostLayer ghost; exchange_ghosts(d, local, gids, ghost);   // n=1: 空
std::vector<Vec3> forces, torques; std::vector<int> cc;
int total_contacts = 0;
for (int step = 0; step < cfg.steps; ++step) {
    forces.assign(local.size(), Vec3{}); torques.assign(local.size(), Vec3{}); cc.assign(local.size(), 0);
    // ---- pp:Newton-off ----
    for (int li = 0; li < (int)local.size(); ++li) {
        // 本地 j(gid 大于本地 i 的 gid,每本地对只算一次,双侧累加,同 CPU 序)
        for (int lj = li + 1; lj < (int)local.size(); ++lj) {
            Vec3 f_i,t_i,f_j,t_j;
            if (pp_contact_pair(local[li], local[lj], sim.meshes[local[li].mesh_index],
                                sim.meshes[local[lj].mesh_index], f_i,t_i,f_j,t_j,
                                cfg.tangential_damping) > 0) {
                forces[li]+=f_i; torques[li]+=t_i; forces[lj]+=f_j; torques[lj]+=t_j;
                cc[li]++; cc[lj]++; total_contacts++;
            }
        }
        // ghost j(只取 own 半)
        for (int gj = 0; gj < (int)ghost.particles.size(); ++gj) {
            if (ghost.gids[gj] == gids[li]) continue;    // 不可能(迁移保证唯一),防御
            const Particle& g = ghost.particles[gj];
            Vec3 dv = g.tf.pos - local[li].tf.pos;
            double rs = local[li].radius + g.radius;
            if (dot(dv,dv) > rs*rs) continue;
            Vec3 f_i,t_i,f_j,t_j;
            if (pp_contact_pair(local[li], g, sim.meshes[local[li].mesh_index],
                                sim.meshes[g.mesh_index], f_i,t_i,f_j,t_j,
                                cfg.tangential_damping) > 0) {
                forces[li]+=f_i; torques[li]+=t_i; cc[li]++; total_contacts++;
            }
        }
    }
    // ---- 墙(本地粒子 × 全墙)----
    for (int li = 0; li < (int)local.size(); ++li) {
        Vec3 f, t;
        cc[li] += wall_contact_particle(local[li], sim.meshes[local[li].mesh_index],
                                        sim.walls, cfg.tangential_damping, f, t) > 0 ? 1 : 0;
        forces[li] += f; torques[li] += t;
        // 注意保持 CPU 语义:wall_contact_particle 返回组数,cc += (组数>0 ? 1 : 0)——CPU 每墙一个计数;
        // 对照 main.cpp 1676-1690 的 contact_counts 语义,若 CPU 是 +=wall组数 则照抄
    }
    // ---- 积分(本地)----
    for (int li = 0; li < (int)local.size(); ++li)
        integrate_particle(local[li], forces[li], torques[li], cfg.gravity, cfg.dt);
    // ---- 输出 ----
    if (step % cfg.output_interval == 0) {
        gather_write_frame(d, cfg, sim, local, gids, cc, forces, torques, step);
        if (d.rank == 0) std::printf("step=%d contacts=%d\n", step, reduce_add(d, total_contacts));
        // total_contacts 为累计?——CPU 语义是每步值:对照 main.cpp 日志行,若为每步计数则每步重置
    }
    total_contacts = 0;   // 依 CPU 语义调整位置
    // ---- 迁移(Task 6)与 ghost 重建(Task 5)在此插入 ----
}
```
**执行者注意**:两处 `对照 CPU` 注释必须先读 main.cpp 对应行确认语义(wall cc 累加方式、contacts 每步/累计),照抄,不要猜。

- [ ] **Step 3: 实现 gather.cpp(如 Interfaces 所述)+ ghost.cpp 占位**

- [ ] **Step 4: V0 parity 门槛(n=1)**

```bash
cmake --build build --config Release 2>&1 | grep -iE "error"
mpiexec -n 1 build/Release/zdem_mpi.exe --config scratch/mpi_v0_wall64.txt | tail -2
mpiexec -n 1 build/Release/zdem_mpi.exe --config scratch/mpi_v0_pile12.txt | tail -2
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_wall64_cpu scratch/mpi_v0_wall_out 1 | tail -2
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/out_pile12_cpu64 scratch/mpi_v0_p12_out 12 | tail -2
```
Expected: 两条 `RESULT: PASS`(pos ≤1e-6, energy ≤1e-9)。**先跑 wall64 再跑 pile12**;FAIL 则按 GPU 版调试方法(能量序列定位注入源)
注:pile12 对照基线 out_pile12_cpu64 是 16k 步版,mpi_v0_pile12 是 2k 步、interval=2000 → 帧集合不一致会让 gpu_compare 报 frame-set mismatch。修正:pile12 基线用 `mpiexec -n 1` 之外的 zdem_cpu 跑同一 2k config 建 `scratch/mpi_v0_p12_cpu_out`,或把 mpi_v0_pile12.txt 改成 16k 步复用现基线(推荐后者,改 steps=16000, interval 不变)。

- [ ] **Step 5: ctest + Commit**

```bash
cd build && ctest -C Release 2>&1 | grep -E "tests passed|failed"
cd .. && git add -A && git commit -m "feat(mpi): Newton-off main loop, gather output, V0 single-rank parity"
```

---

### Task 5: ghost 三轴链式交换 + mpi_smoke ctest + V1 parity

**Files:**
- Modify: `src/mpi/ghost.cpp`(真实现)
- Create: `tests/smoke/check_mpi_smoke.py`
- Modify: `CMakeLists.txt`(注册 mpi_smoke)
- Test: `scratch/mpi_v1_pile12.txt`

**Interfaces:**
- Consumes: Task 3 `Decomp`(cart/coords/ghost_depth)、Task 4 GhostLayer
- Produces: `exchange_ghosts` 真实现(签名不变)

- [ ] **Step 1: 实现 exchange_ghosts(三轴链式,LAMMPS comm_brick 式)**

```cpp
// 伪代码级完整逻辑(执行者按此写):
// 工作集 W = local(带 gid)——注意链式:轴 k 交换收到的 ghost 加入 W 参与轴 k+1 的选择
// for axis in {0,1,2}:
//   MPI_Cart_shift(cart, axis, 1, &src, &dst)     // dst = +axis 邻居,src = -axis 邻居
//   to_hi  = [w in W | w.pos[axis] >  sub_hi[axis] - ghost_depth]   // 发给 +邻居
//   to_lo  = [w in W | w.pos[axis] <  sub_lo[axis] + ghost_depth]   // 发给 -邻居
//   MPI_Sendrecv(send to dst: pack(to_hi), recv from src: probe->recv -> got_lo,
//                send to src: pack(to_lo), recv from dst: probe->recv -> got_hi)
//   // MS-MPI 无 MPI_Sendrecv 变长问题:用 MPI_Probe + MPI_Get_count 两段式,或
//   // 先 Alltoall 交换 count 再定长 Sendrecv(推荐后者,简单可移植):
//   //   先与两个邻居各交换 count(MPI_Sendrecv 1 int),再交换数据
//   W += got_lo + got_hi(标记为 ghost)
// 返回 out.particles/gids = W 中非本地的全部(去重:同 gid 可能经不同路径到达?
//   ——链式选择条件是坐标邻近,同粒子经两条路径到达同一 rank 当且仅当它同时在两个
//   发送窗口内;brick ≥ ghost_depth 校验保证窗口不重叠 → 无重复。加 assert(开发期):
//   out.gids 无重复)
```
打包:`struct PackedPart { Particle p; int gid; };`(POD,MPI_BYTE)。

- [ ] **Step 2: mpi_smoke ctest**

`tests/smoke/check_mpi_smoke.py`(仿 check_gpu_smoke.py):
```python
# --exe zdem_mpi 路径 --root 仓库根 --mpiexec "D:/Program Files/Microsoft MPI/Bin/mpiexec.exe"(where mpiexec 探测,无则 SKIP)
# 改写 config/wall_multiloop_regression.txt 的 output_dir → build/test_output/mpi_smoke
# 跑 [mpiexec, "-n", "2", exe, "--config", cfg]
# 断言:returncode 0、"done" in stdout、r"step=0 contacts=(\d+)" >= 1
```
CMakeLists(python 分支内追加):
```cmake
    add_test(
        NAME mpi_smoke
        COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/tests/smoke/check_mpi_smoke.py
                --exe $<TARGET_FILE:zdem_mpi>
                --root ${CMAKE_SOURCE_DIR}
    )
```

- [ ] **Step 3: V1 parity(n=2,跨域接触)**

```bash
# x 二分:n=2 时 dims 由盒长比自动定(pile12 盒 x 最长 → 2×1×1,x 切,粒子跨域接触)
sed -e 's/^steps = 16000/steps = 2000/' -e 's|output_dir = scratch/mpi_v0_p12_out|output_dir = scratch/mpi_v1_p12_out|' scratch/mpi_v0_pile12.txt > scratch/mpi_v1_pile12.txt
mpiexec -n 2 build/Release/zdem_mpi.exe --config scratch/mpi_v1_pile12.txt | tail -2
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/mpi_v0_p12_out scratch/mpi_v1_p12_out 12 | tail -2
```
Expected: `RESULT: PASS`(≤1e-6/1e-9;跨域 pair ulp 级差,容差内)
注:V1 基线用 V0 的 n=1 输出(已过 parity),等价于对照 CPU。

- [ ] **Step 4: ctest(7 项)+ Commit**

```bash
cd build && ctest -C Release --output-on-failure 2>&1 | tail -4
cd .. && git add -A && git commit -m "feat(mpi): chained ghost exchange, mpi smoke test, V1 two-rank parity"
```

---

### Task 6: 粒子迁移 + V2(跨域飞行/对撞)

**Files:**
- Create: `src/mpi/migrate.hpp`, `src/mpi/migrate.cpp`
- Modify: `src/mpi/mpi_main.cpp`(积分输出后插 `migrate_particles` + `exchange_ghosts`)
- Test: `scratch/mpi_v2_fly.txt`、`scratch/mpi_v2_pp.txt`

**Interfaces:**
- Produces:
```cpp
// migrate.hpp
// 质心越界迁移(三轴链式,同 ghost 的交换模式但只发本地粒子中越出本 brick 该轴范围的;
// 收到的进入本地集继续参与下一轴)。单轮结束后若本地仍有越界粒子(一步跨多格),
// 重复整轮,上限 8 轮;仍越界或越出全局盒 → rank0 报错 Abort。
// 每步结束(仅 debug 或每 K 步)调 assert_global_count 可选校验。
void migrate_particles(const Decomp& d, std::vector<Particle>& local, std::vector<int>& gids);
void assert_global_count(const Decomp& d, long long expected_total);  // Allreduce SUM,不符则 Abort
```

- [ ] **Step 1: V2 config**

```bash
# (a) 单粒子跨域飞行:全局盒 x [-0.5,0.5],2 rank x 二分,vx=1 → 0.4s 内飞过边界(x=0)
cat > scratch/mpi_v2_fly.txt <<'EOF'
steps = 40000
dt = 0.00001
gravity = 0 0 0
center_mesh = 1
vtk_prefix = particle
output_dir = scratch/mpi_v2_fly_out
output_interval = 2000
mpi_box = -0.5 0.5 -0.5 0.5 -0.5 0.5
particle
stl = geometry/low-poly-banbana.stl
pos = -0.2 0 0
vel = 1 0 0
quat = 1 0 0 0
scale = 0.001
density = 2500
young = 1e7
poisson = 0.25
mu = 0.5
restitution = 0.5
end_particle
EOF
# (b) 对撞跨域:pp_v8full.txt 改 steps=5000 + mpi_box + 输出目录 scratch/mpi_v2_pp_out
```

- [ ] **Step 2: 实现 migrate.cpp(如 Interfaces 所述;交换模式复用 ghost.cpp 的 count-then-data Sendrecv,抽成两文件共用的 static helper 或直接重复——建议抽 `src/mpi/comm_util.hpp`: `sendrecv_particles(...)`)**

- [ ] **Step 3: 主循环接线 + V2 验证**

```bash
# 输出之后:migrate_particles(d, local, gids); exchange_ghosts(d, local, gids, ghost);
mpiexec -n 2 build/Release/zdem_mpi.exe --config scratch/mpi_v2_fly.txt | grep -E "step=(0|10000|20000|30000|39000)|done"
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/mpi_v2_fly_cpu_out scratch/mpi_v2_fly_out 1 | tail -2
mpiexec -n 2 build/Release/zdem_mpi.exe --config scratch/mpi_v2_pp.txt | tail -1
D:/ProgramData/anaconda3/python.exe scratch/gpu_compare.py scratch/mpi_v2_pp_cpu_out scratch/mpi_v2_pp_out 2 | tail -2
```
(单进程 CPU 基线:zdem_cpu 跑同 config 改 output_dir,飞行算例无墙无接触纯积分——CPU 可跑)
Expected: 飞行:x 从 -0.2 线性到 +0.2(gpu_compare PASS);对撞:PASS(粒子相遇在 x=0 边界,跨域 pair 双算 + 迁移同时生效);`rank_n` 日志显示迁移后粒子数变化正确

- [ ] **Step 4: ctest + Commit**

```bash
cd build && ctest -C Release 2>&1 | grep -E "tests passed|failed"
cd .. && git add -A && git commit -m "feat(mpi): centroid migration, V2 cross-domain flight and collision parity"
```

---

### Task 7: V3 验收(pile64 × 4 进程)+ 性能记录

**Files:**
- Create: `scratch/mpi_v3_pile64.txt`(pile64.txt 改 output_dir + steps=250000 不变)
- Create: `docs/superpowers/notes/2026-08-18-mpi-performance.md`
- Modify: `docs/superpowers/specs/2026-08-18-mpi-dem-design.md`(验证表勾选结果)

- [ ] **Step 1: V3 算例**

```bash
sed -e 's|output_dir = scratch/out_pile64|output_dir = scratch/mpi_v3_p64_out|' scratch/pile64.txt > scratch/mpi_v3_pile64.txt
mpiexec -n 4 build/Release/zdem_mpi.exe --config scratch/mpi_v3_pile64.txt > scratch/mpi_v3_run.log 2>&1; echo "EXIT=$?"
D:/ProgramData/anaconda3/python.exe scratch/pile_analyze.py scratch/mpi_v3_p64_out 64 | tail -8
```
Expected: E 单调降后恒定(±1%)、min_z 稳定无穿模、ALL AT REST 或静止比例 ≥ 90%(同单进程判据;4 进程 dims 自动 2×2×1)

- [ ] **Step 2: 性能记录(1/2/4 进程,20k 步无输出)**

```bash
for n in 1 2 4; do
  sed -e 's/^steps = 250000/steps = 20000/' -e 's/^output_interval = 2000/output_interval = 25000/' \
      -e 's|output_dir = scratch/mpi_v3_p64_out|output_dir = scratch/mpi_t'"$n"'_out|' scratch/mpi_v3_pile64.txt > scratch/mpi_t$n.txt
  echo "=== n=$n ==="; time mpiexec -n $n build/Release/zdem_mpi.exe --config scratch/mpi_t$n.txt | tail -1
done
```
ms/step = (real - 启动)/20000(启动≈0.5s,mpiexec 额外开销记录);对照 zdem_cpu 单进程 11.41 ms/步(既有 GPU notes 数据)写入 `docs/superpowers/notes/2026-08-18-mpi-performance.md`:表(n, ms/step, 加速比, rank_n 分布)+ 结论(单机共享内存上限,通信占比)+ 已知瓶颈(O(n²) broadphase、Newton-off 冗余)

- [ ] **Step 3: Commit + 收尾(fixed 流程)**

```bash
git add -A && git commit -m "test(mpi): pile64 four-rank acceptance and performance notes (V3)"
# 然后 finishing-a-development-branch 固定路径:ctest 全绿 → checkout main → merge --ff-only → 删分支 → 烟雾(ctest + mpiexec -n 2 smoke)
```

---

## Self-Review 记录

- **Spec 覆盖**:§总体架构(Task 2 抽取 + Task 3 骨架)、§进程模型/全局盒(Task 3)、§所有权与迁移(Task 6)、§Ghost(Task 5)、§Newton-off 接触(Task 4)、§输出/日志(Task 4 gather + Task 5 reduce/mpi_smoke)、§验证阶梯 V0-V3(Task 4/5/6/7)、§构建(Task 3 CMake)——全覆盖;spec"每步顺序"修正为"输出在迁移前"(force/cc 数组不随迁移重排),已在 Global Constraints 说明理由
- **Placeholder 扫描**:Task 4 主循环中两处"对照 CPU 语义"标注了精确行号(main.cpp 1676-1690)并要求执行者先读后抄——这是引用而非占位;Task 2 Step 2 的备选接口(逐顶点增量)完整给出
- **类型一致**:`SimBuild/Decomp/GhostLayer/FrameSnapshot/pp_contact_pair/wall_contact_particle/migrate_particles` 各任务间签名一致;`Vec3/Mesh/Particle/Wall` 来自现有 core/host 头
- **已知风险内联**:Task 2 整块累加的浮点序风险(备选方案给出)、Task 4 frame-set 基线匹配(pile12 用 16k 版)、Task 5 链式去重论证(brick≥ghost_depth 保证)
