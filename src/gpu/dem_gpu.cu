// src/gpu/dem_gpu.cu -- GPU DEM driver skeleton (v8 penalty contact)
// Real driver logic starts in Task 3; this file only verifies the CUDA
// build skeleton: real type switch, CUDA_CHECK macro, zdem_core linkage.
#include <cstdio>
#include <cuda_runtime.h>

#include "gpu/cuda_check.hpp"
#include "gpu/real.hpp"
#include "host/config_io.hpp"

__host__ inline const char* precision_name() {
    return sizeof(real) == 8 ? "double" : "float";
}

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;
    std::printf("zdem_gpu (CUDA, precision=%s)\n", precision_name());
    return 0;
}
