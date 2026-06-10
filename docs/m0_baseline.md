# M0 Baseline

Date: 2026-06-10

## Purpose

This document records the reproducible baseline for the current C++ DEM prototype before the M1 math and geometry extraction.

## Commands

```powershell
cmake -S . -B build
cmake --build build
.\build\Debug\zdem_cpu.exe --config config/example_sim.txt
```

## Expected Behavior

- The executable builds as `zdem_cpu`.
- The example simulation prints progress lines containing `step=`.
- The run exits with code 0.
- The run prints `done` before exit.

## Known Current Limitations

- `src/main.cpp` contains math, geometry, contact detection, force calculation, integration, IO, and CLI code in one file.
- The baseline uses a brute-force particle-pair loop and triangle checks inside candidate pairs.
- Particle-wall contact and tangential history exist in the prototype but are not isolated behind stable interfaces.
- There is no test target in CMake before M1.
