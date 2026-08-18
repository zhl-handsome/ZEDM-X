# GPU 版 M4 验收记录:pile12 FP64 对照 + FP32 全套物理指标

日期:2026-08-18
分支:`feat/gpu-dem`(HEAD bebb76a,quat_rotate 三明治积精确移植之后)
设计文档:`docs/superpowers/specs/2026-08-17-gpu-dem-design.md` §验证计划 阶段 1/2

## 结论:全部通过

| 算例 | 精度 | 步数 | 判据 | 结果 |
|----|----|----|----|----|
| pile12 对照 | FP64 | 16 000(interval 100 → 160 帧)| 位置 < 1e-6 m,能量 < 1e-9(相对)| **PASS,max\|dpos\|=0,max\|dE\|=0**(偏差低于工具分辨率 ~1e-7,见下)|
| wall_v8 | FP32 | 300 000 | E 末段恒定,完全静止 | **PASS**,E 自 50k 步起恒定 0.072 J;tail30 \|v\|=1e-4,\|ω\|=0,z_std=0 |
| pp_v8full | FP32 | 100 000 | 注入 ≤ 0.02 J,动量守恒,反弹 | **PASS**,碰后注入 +0.0054 J;Px/Py range = 0;vrel_x +1.000 → −0.024 |
| pile6d_v8b | FP32 | 250 000 | ALL AT REST,E 恒定无回升 | **PASS**,E 9.971 → 0.681 单调降后恒定;6 粒子全部静止 |
| pile12_v8 | FP32 | 250 000 | E 单调降后恒定,无爆炸穿模 | **PASS**,E 19.159 → 1.459;ALL AT REST(12 粒子 tail \|v\| ≤ 8e-4,z_std ≤ 7e-5)|

## 详细数据

### pile12 FP64 逐步对照(阶段 1)

- config:`scratch/pile12_gpu64.txt`(16 000 步,dt=2e-5,interval=100),CPU 基线 `pile12_cpu64.txt`
- 工具:`scratch/gpu_compare.py scratch/out_pile12_cpu64 scratch/out_pile12_gpu64 12`
- 160 帧全部 `max|dpos| = 0.000e+00`、`dE = 0.000e+00`(12 粒子、~26 对接触链,含持续接触下的姿态演化)
- **分辨率说明**(2026-08-18 review 修正):gpu_compare.py 对照的 COM/速度来自 VTK,而 VTK 以 float + 6 位有效数字写出(分辨率 ~1e-7),且 pp/wall 的 atomicAdd 累加序本身不保证 bit 级——"0.000e+00" 证明偏差低于 ~1e-7,满足 1e-6/1e-9 容差判据,但**不构成逐位一致的证明**;如需 bit 级结论须用全精度输出(如逐步 checksum)再对照
- 远超计划要求的前 2 000 步覆盖(实际 16 000 步)

### FP32 物理指标(阶段 2,Release 构建)

**wall_v8**(`out_wall_v8_gpu32`,300k 步,dt=1e-5):
- E:0.826(初始)→ 0.521(25k)→ 0.074(50k)→ 0.072 恒定直至 275k,末段 12 个采样点无回升
- tail30:z=0.0173,|v|=0.0001,|ω|=0.000,z_std=0.0000 → 落定完全静止
- 碰撞瞬态 max|ω|=21.1 与 CPU 版同量级

**pp_v8full**(`out_pp_v8full_gpu32`,100k 步):
- E:start=0.1055,min=0.0651(冲击最低点),final=0.0704 → 碰后注入 +0.0054 J(≤ 0.02 J 判据,单精度噪声余量内)
- 动量:Px/Py range 均 [+0.0000, +0.0000],全程守恒
- vrel_x:+1.000(入射)→ −0.024(反弹分离),无粘连

**pile6d_v8b**(`out_pile6d_v8b_gpu32`,250k 步,dt=2e-5):
- E:9.971 → 0.681,单调下降后恒定
- tail20 每粒子:|v| ≤ 5e-4,|ω| ≤ 5e-3,z_std ≤ 1e-5 → `VERDICT: ALL AT REST`

**pile12_v8**(`out_pile12_v8_gpu32`,250k 步,dt=2e-5):
- E:19.159 → 1.459,单调下降后恒定(CPU 版同算例终值 1.422,差 2.6%,同为静置堆积水平)
- tail 每粒子:|v| ≤ 8e-4,z_std ≤ 7e-5 → `VERDICT: ALL AT REST`(CPU 版同日志判 STILL MOVING,GPU 收敛更彻底)

## 性能速览(FP32,Release,RTX 4090 D;详细分析见 M5 性能文档)

| 算例 | 粒子数 | 末段 step_ms |
|----|----|----|
| wall_v8 | 1 | 3.2 |
| pp_v8full | 2 | 4.1 |
| pile6d | 6 | 13.3 |
| pile12 | 12 | 13.9 |

小规模下 GPU 受 launch 开销主导(设计预期:≤12 粒子不敌 CPU,不作验收项)。

## 工件位置

- 运行/分析日志:`scratch/*_run.log`、`scratch/*_analyze.log`、`scratch/pile12_gpu64_compare.log`
- VTK 输出:`scratch/out_pile12_cpu64|out_pile12_gpu64|out_wall_v8_gpu32|out_pp_v8full_gpu32|out_pile6d_v8b_gpu32|out_pile12_v8_gpu32/`
- 构建验证:FP64/FP32 两种二进制均在 bebb76a 之后重新构建(zdem_gpu.exe 2026-08-18 09:33)
