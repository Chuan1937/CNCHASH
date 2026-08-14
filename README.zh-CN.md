# CNCHASH

基于 **Modern Fortran/OpenMP 原生后端** + 简洁 Python 前端的 HASH 地震
震源机制解算高性能实现。

**English** | [中文](./README.zh-CN.md)

![Python](https://img.shields.io/badge/python-3.10+-orange.svg)
![License](https://img.shields.io/badge/license-BSD%203--blue.svg)
![Fortran](https://img.shields.io/badge/fortran-modern%20fortran-orange.svg)
![Numpy](https://img.shields.io/badge/numpy-1.20+-yellow.svg)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)

所有数值计算（网格搜索、S/P 振幅、不确定性分析、速度模型表）都在
Modern Fortran/OpenMP 核心（`libhashhp`）中运行；Python 层负责数据、
文件格式和 API。OpenMP 线程与事件级批量并行可跨核扩展。

## 性能

30 台站、30 次 MC 试验、dang=5 度（详见
[docs/benchmarks.md](docs/benchmarks.md)）：

| 后端构建               | ms/事件 |
|------------------------|---------|
| 可移植，1 线程         | ~49     |
| 可移植，4 线程         | ~19     |
| 原生（AVX2），1 线程   | ~23     |
| 原生（AVX2），16 线程  | ~7      |
| 原版 HASH v1.2         | ~145    |

原生构建（`-DCNCHASH_NATIVE_MARCH=ON`）在 1 线程下约快 2.5 倍；它是
机器相关的，只应在编译它的机器上使用。

## 准确性验证

![准确性验证](docs/images/accuracy_verification.png)

**主要结果：**（见 tests/test_accuracy.py）
- 无极性错误时合成恢复：10/10 个机制，Kagan 旋转角中位数约 8°
  （dang=5 网格离散化下限）
- 约 10% 极性翻转时：9/10 保持在 25° 接受角内
- 相同输入下可接受机制集合与原版 HASH v1.2 完全一致（golden-reference 测试）

**注意：** 走向差异（40-80°）是正常的——震源机制有两个正交节面，都
能满足极性数据。

## 快速开始（源码构建）

CNCHASH 2.x 是源码安装包：原生后端需要用 CMake 编译一次。需要支持
OpenMP 的 Fortran 编译器（gfortran >= 9）。

```bash
git clone https://github.com/Chuan1937/CNCHASH.git
cd CNCHASH

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# 可选：本机性能版（约快 2.5 倍，仅限当前机器）
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCNCHASH_NATIVE_MARCH=ON

export CNCHASH_HASHHP_LIB=$PWD/build/libhashhp.so
pip install -e .
```

没有 Fortran 编译器时安装仍会成功，但后端不可用；`get_backend_info()`
会报告库路径，`available_backends()` 保持为空。

### 基本用法（仅 P 波极性）

```python
from cnchash import run_hash
import numpy as np

# 方位角、离源角、极性、质量
p_azi = np.array([45.0, 135.0, 225.0, 315.0])
p_the = np.array([30.0, 45.0, 60.0, 75.0])
p_pol = np.array([1, -1, 1, -1])  # 1=向上, -1=向下
p_qual = np.array([0, 0, 0, 0])

result = run_hash(p_azi, p_the, p_pol, p_qual)

print(f"走向: {result['strike_avg']:.1f}")
print(f"倾角: {result['dip_avg']:.1f}")
print(f"滑动角: {result['rake_avg']:.1f}")
print(f"质量等级: {result['quality']}")
```

### 加入 S/P 振幅比

```python
from cnchash import run_hash_with_amp

# S/P 比值（log10 尺度），0.0 = 无数据
sp_amp = np.array([0.3, -0.2, 0.5, 0.0])

result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp)

print(f"极性失配率: {result['mfrac']*100:.1f}%")
print(f"振幅失配: {result['mavg']:.2f}")
```

## 后端选择

```python
from cnchash import run_hash, run_hash_batch, available_backends, get_backend_info

# backend="auto"（默认）使用原生后端
result = run_hash(p_azi, p_the, p_pol, p_qual, backend="auto", num_threads=16)

# 批量模式：一次原生调用处理大量事件（事件级 OpenMP）
results = run_hash_batch(events, backend="fortran", num_threads=16)

print(available_backends())   # ['fortran']
print(get_backend_info())
```

后端架构、构建说明和设计规则见
[docs/native_backend.md](docs/native_backend.md)。

## 特性

- 震源机制网格搜索
- Monte Carlo 不确定性分析
- S/P 振幅比约束
- 质量等级评定（A-D, E, F）
- 多种震相文件格式
- Modern Fortran/OpenMP 原生后端
- 事件级批量并行
- 核心算法与原版 HASH v1.2 一致（golden-reference 测试）

## 文档

完整文档见 https://chuan1937.github.io/CNCHASH/，包括：
- 安装与源码构建
- 使用与后端选择
- 原生后端架构
- 性能基准
- API 参考

## 运行测试

```bash
pytest tests/          # parity, determinism, batch, velocity suites
jupyter notebook HASH_Tests.ipynb
```

## 项目结构

```
cnchash/
├── backend/       # 原生后端绑定（ctypes）
├── driver.py      # 主驱动（run_hash, run_hash_batch）
├── io.py          # 文件读写
└── utils.py       # 工具

src/hash_hp/       # Modern Fortran/OpenMP 核心（libhashhp）
HASH_complete/     # 原版 Fortran HASH v1.2（不可变 golden reference）
benchmarks/        # 性能基准套件
tests/             # parity, determinism, batch, velocity 测试
```

## 作者

He XingChen

## 许可证

BSD 3-Clause

## 参考文献

Hardebeck, Jeanne L. and Peter M. Shearer, A new method for determining first-motion
focal mechanisms, Bulletin of the Seismological Society of America, 92,
2264-2276, 2002.

Hardebeck, Jeanne L. and Peter M. Shearer, Using S/P Amplitude Ratios to
Constrain the Focal Mechanisms of Small Earthquakes, Bulletin of the
Seismological Society of America, 93, 2434-2444, 2003.
