# 原生 Fortran 后端（HASH-HP）

CNCHASH 基于 Modern Fortran + OpenMP 后端构建。所有数值计算（网格搜索、
S/P 振幅、不确定性分析、速度模型表）都在 `libhashhp` 中运行；Python
层负责 API、文件处理和批量调度。如果没有编译原生库，包会给出明确的
构建错误。

## 架构

```
HASH_complete/            # 原始 Fortran HASH v1.2（不可变上游，golden reference）
cnchash/                  # Python 包
  ├── backend/
  │   ├── base.py         # HashBackend 接口
  │   └── fortran_backend.py  # libhashhp 的 ctypes 绑定
  ├── driver.py           # run_hash / run_hash_batch
  ├── io.py               # 震相/台站文件处理
  └── utils.py            # 纯 Python 工具
src/hash_hp/              # Modern Fortran 核心
  ├── hash_kinds.f90      # real64/int32 类型定义
  ├── hash_geometry.f90   # TO_CAR/FPCOOR/cross
  ├── hash_rotation.f90   # 旋转网格缓存
  ├── hash_runtime.f90    # 网格缓存 + OpenMP 线程控制
  ├── hash_focalmc.f90    # FOCALMC，旋转方向 OpenMP 并行
  ├── hash_misfit.f90     # GET_GAP/GET_MISF
  ├── hash_uncertainty.f90# MECH_ROT/MECH_AVG/MECH_PROB
  ├── hash_amplitude.f90  # FOCALAMP_MC/GET_MISF_AMP
  ├── hash_velocity.f90   # LAYERTRACE/MK_TABLE/GET_TTS
  ├── hash_batch.f90      # 完整事件流水线 + 事件级并行批处理
  └── hash_c_api.f90      # 稳定的 ISO_C_BINDING ABI（ctypes，无需 f2py）
```

## 设计规则

- **内核无 I/O。** `src/hash_hp` 只接收数组并返回数组；文件、格式和
  异常都由 Python 处理。
- **热循环零分配。** 每个事件的 workspace 复用；旋转网格每个 `dang`
  只构建一次，之后共享只读。
- **禁止嵌套 OpenMP。** 批量模式要么按事件并行、要么按旋转并行，绝不
  同时嵌套。
- **默认确定性。** 可接受机制按网格顺序取前 `maxout` 个
  （`selection=0`）。传入 `selection=1` 可复现原始 HASH 的随机选择。
- **默认不开 `-ffast-math`。** 科学计算正确性优先；如需开启显式设置
  `CNCHASH_FAST_MATH=ON`。
- **`HASH_complete/` 保持不动。** 它是 golden reference。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

默认为可移植构建。设置 `-DCNCHASH_NATIVE_MARCH=ON` 可生成仅在当前
机器上运行、约快 2.5 倍的本机优化版本。

## 后端选择

```python
from cnchash import run_hash, available_backends, get_backend_info

run_hash(az, the, pol, qual, backend="auto")     # 默认（fortran）
run_hash(az, the, pol, qual, backend="fortran")  # 显式指定
print(available_backends())
print(get_backend_info())
```

环境变量：`CNCHASH_BACKEND`、`CNCHASH_NUM_THREADS`、
`CNCHASH_HASHHP_LIB`（显式指定库路径）。

## 批量模式

`run_hash_batch` 一次调用将大量事件送入 Fortran 核心（CSR 风格扁平
数组），大型目录使用事件级 OpenMP 并行：

```python
from cnchash import run_hash_batch

results = run_hash_batch(events, nmc=30, backend="fortran", num_threads=16)
```

## 与 HASH v1.2 的有意差异

| 项目 | HASH v1.2 | HASH-HP |
|------|-----------|---------|
| GET_GAP 空隙 | INTEGER（隐式类型截断） | REAL（空隙为连续值） |
| GET_MISF 断层面向量 | 弧度被误当作度数（单位 bug） | 已修复：向 FPCOOR 传度数 |
| GET_MISF_AMP stdr | 振幅台站累计权重 | 相同 |
| nf > maxout 时的选择 | `rand(0)` 随机 | 默认确定性 |
| velocity GET_TTS 深度检查 | `d(nd0)`（未初始化） | `d(ndep)` |
| velocity 射线缓冲区 | 固定 10001（可能溢出） | 20001 堆分配 |

所有差异都是 bug 修复。测试保证可接受机制集合和首选解在相同输入下
与原始 HASH v1.2 一致（见 tests/test_accuracy.py），仅剩聚类边界上
单精度与双精度带来的微小差异。
