# 安装

## 环境要求

- Python >= 3.10
- NumPy >= 1.20.0
- 支持 OpenMP 的 Fortran 编译器（gfortran >= 9），用于构建原生后端

## 安装（源码构建）

CNCHASH 2.x 是源码安装包：原生后端需要用 CMake 编译一次。

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

包会自动检测 `libhashhp`（也可通过 `CNCHASH_HASHHP_LIB` 显式指定）。

## 开发环境

```bash
pip install -e ".[dev]"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## 验证安装

```python
from cnchash import run_hash, get_backend_info
print(get_backend_info())
```

`get_backend_info()` 报告后端类型、OpenMP 线程数和库路径。如果没有编译
原生库，会给出明确的构建错误提示。
