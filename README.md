# RabbitBAM-X

An ultra-fast I/O framework for BAM files, supporting both x86 and Sunway platforms.

## Features

- High-performance BAM file reading and writing
- Multi-threaded processing with OpenMP (x86) or slave cores (Sunway)
- Support for both x86 and Sunway many-core architectures
- Optimized compression/decompression using libdeflate

## Platform Support

- **x86**: Uses OpenMP for parallel processing
- **Sunway**: Uses slave cores (athread) for parallel processing

## Dependencies

- CMake 3.10 or newer
- GCC 8.5.0 or newer (for x86)
- SWGCC/SWG++ (for Sunway platform)
- htslib
- libdeflate
- zlib

## Installation

### Step 1: Build External Dependencies

Before building RabbitBAM-X, you need to build the external dependencies (htslib and libdeflate).

**For x86 Platform:**
```bash
bash build_ext_x86.sh
```

**For Sunway Platform:**
```bash
bash build_ext_sw.sh
```

### Step 2: Build RabbitBAM-X

**For x86 Platform:**
```bash
mkdir build_x86 && cd build_x86
cmake .. -DPLATFORM=x86
make
```

**For Sunway Platform:**
```bash
mkdir build_sunway && cd build_sunway
cmake .. -DPLATFORM=sunway
make
```

**Optional: Enable SWLU Profiling (Sunway only):**
```bash
mkdir build_sunway && cd build_sunway
cmake .. -DPLATFORM=sunway -DUSE_SWLU=ON
make
```

## Usage

After building, the executable `RabbitBAM-X` will be in the build directory.

**For x86 Platform:**
```bash
./RabbitBAM-X -h
```

**For Sunway Platform:**

When running on Sunway platform using `bsub`, you need to set the `LD_LIBRARY_PATH` environment variable:

```bash
bsub -q q_sw_expr -I -b -J test -N 1 -np 1 -cgsp 64 -share_size 12000 -cache_size 32 \
     -o output.log \
     -x LD_LIBRARY_PATH=./ext/htslib-1.20:./ext/libdeflate-1.20/build \
     ./RabbitBAM-X api_test -i input.bam -o output.bam --nr 1 --nw 1

note: sunway platform only support compress level = 1 with cache_size = 32KB && sunway_api_test para with --tgs
```

Note: The `-x LD_LIBRARY_PATH=./ext/htslib-1.20:./ext/libdeflate-1.20/build` parameter is required for Sunway platform to locate the shared libraries.

## Citation

If you use RabbitBAM-X in your research, please cite:

```
Yan, L., Zhao, Z., Yin, Z., Zhang, T., Yang, Y., Zhu, F., ... & Liu, W. (2025). 
RabbitBAM: Accelerating BAM File Manipulation on Multi-Core Platforms. 
IEEE Transactions on Computational Biology and Bioinformatics.
```
