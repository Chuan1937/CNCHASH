# API 参考

## 主要函数

### run_hash()

根据 P 波初动极性确定震源机制。

```python
result = run_hash(p_azi, p_the, p_pol, p_qual, **kwargs)
```

**参数：**

| 名称 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `p_azi` | ndarray | 必填 | 方位角（度），形状 (nsta,) 或 (nsta, nmc) |
| `p_the` | ndarray | 必填 | 离源角（度） |
| `p_pol` | ndarray | 必填 | 极性：1=向上, -1=向下 |
| `p_qual` | ndarray | 必填 | 质量：0=清晰, 1=模糊 |
| `dang` | float | 5.0 | 网格角度步长 |
| `nmc` | int | 30 | Monte Carlo 试验次数 |
| `maxout` | int | 500 | 最大输出机制数 |
| `badfrac` | float | 0.1 | 允许的坏极性比例 |
| `npolmin` | int | 8 | 最少极性数 |
| `max_agap` | float | 90.0 | 最大方位角空隙 |
| `max_pgap` | float | 60.0 | 最大离源角空隙 |

**返回值：**

| 键 | 类型 | 说明 |
|----|------|------|
| `success` | bool | 是否求解成功 |
| `strike_avg` | float | 走向（度） |
| `dip_avg` | float | 倾角（度） |
| `rake_avg` | float | 滑动角（度） |
| `quality` | str | A、B、C、D、E 或 F |
| `mfrac` | float | 失配比例（0-1） |
| `prob` | float | 解的概率 |
| `stdr` | float | 台站分布比 |
| `nout2` | int | 解的个数 |

---

### run_hash_with_amp()

使用极性 + S/P 振幅比确定机制。

```python
result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp, **kwargs)
```

**额外参数：**

| 名称 | 类型 | 说明 |
|------|------|------|
| `sp_amp` | ndarray | S/P 比值（log10），0.0 = 无数据 |

**额外返回值：**

| 键 | 类型 | 说明 |
|----|------|------|
| `mavg` | float | 振幅失配（log10） |
| `npol` | int | 极性数量 |
| `nspr` | int | S/P 比值数量 |

---

### run_hash_batch() / run_hash_batch_with_amp()

一次原生调用批量处理大量事件（事件级 OpenMP 并行）。

```python
results = run_hash_batch(events, nmc=30, backend="fortran", num_threads=16)
```

### run_hash_from_file()

处理 HASH 输入文件中的事件。

```python
results = run_hash_from_file("example.inp")
```

---

## 质量等级

| 等级 | 判据 |
|------|------|
| **A** | prob > 0.8, var ≤ 25°, misfit ≤ 15%, stdr ≥ 0.5 |
| **B** | prob > 0.6, var ≤ 35°, misfit ≤ 20%, stdr ≥ 0.4 |
| **C** | prob > 0.5, var ≤ 45°, misfit ≤ 30%, stdr ≥ 0.3 |
| **D** | 有解但不满足 C 判据 |
| **E** | 方位角或离源角空隙过大 |
| **F** | 未找到可接受机制 |

---

## 模块

| 模块 | 函数 | 说明 |
|------|------|------|
| `driver.py` | `run_hash`, `run_hash_with_amp`, `run_hash_batch`, `run_hash_batch_with_amp`, `run_hash_from_file` | 主要入口 |
| `backend/fortran_backend.py` | `run_event`, `run_event_amp`, `run_batch`, `build_velocity_table`, `get_tts` | 原生后端绑定 |
| `io.py` | `read_phase_file`, `read_station_file`, `read_polarity_reversal_file`, `read_velocity_model`, `write_mechanism_output` | 文件读写 |
| `utils.py` | `fp_coord_angles_to_vectors`, `fp_coord_vectors_to_angles` | 坐标转换 |

---

## 底层后端方法

这些方法位于 Fortran 后端实例上（通过 `backend.get_backend("fortran")`
获取），测试套件用它们将结果钉在原版 HASH 上：

### run_event()

单个事件的完整极性流水线（等价于不经过输入归一化的 `run_hash`）。

### run_event_amp()

单个事件的完整极性 + S/P 振幅流水线。

### run_batch() / run_batch_amp()

批量流水线（CSR 风格扁平数组），供 `run_hash_batch` 和
`run_hash_batch_with_amp` 使用。

### mech_prob()

从一组可接受机制中求首选机制和概率。

```python
result = backend.mech_prob(nf, faults, slips, cangle=45.0, prob_max=0.1)
```

### get_misfit() / get_misfit_amp() / get_gap()

给定机制下的失配比例、台站分布比和方位角/离源角空隙。

### build_velocity_table() / get_tts()

一维速度模型的离源角表构建与插值。

### get_rotation_grid()

给定 `dang` 的旋转网格缓存（b1/b2/b3 方向余弦）。

---

## 离源角 API

原生速度表之上的公开封装（`cnchash.takeoff`）：

| 接口 | 说明 |
|------|------|
| `TakeoffTable(depth, velocity, params=None, backend="auto")` | 由一维速度模型构建离源角表（MK_TABLE）；按模型内容缓存复用 |
| `TakeoffTable.takeoff(distance, source_depth)` | 单点查询；超出可用范围抛 `TakeoffRangeError` |
| `TakeoffTable.takeoff_batch(distances, source_depths)` | 向量化查询（支持广播）；超出范围返回 NaN |
| `compute_takeoff_angles(depth, velocity, distances, source_depths, ...)` | 便捷函数：建表 + 批量查询一步完成 |
| `DEFAULT_TABLE_PARAMS` | 默认表参数（距离 0-120 km 步长 1，深度 0-35 km 步长 1，nump 1000） |

离源角约定：自竖直方向量起，0° = 垂直向上，90° = 水平，180° = 垂直向下。

---

## 完整原版 HASH 覆盖

Python 接口完整覆盖原版 Fortran 功能：

| 原版 Fortran | Python 接口 |
|--------------|-------------|
| driver1（FPFIT 文件内离源角） | `read_phase_file` 自动解析方位角/离源角/不确定度；`process_event` 按台站不确定度 Monte Carlo |
| driver2（速度模型表） | `velocity_models` 参数 + `TakeoffTable`；深度扰动 + 模型轮换 |
| driver3（S/P 振幅比） | `read_amp_file` + `read_statcor_file` + `ratmin` SNR 过滤；自动走 `run_hash_with_amp` |
| driver5（SIMULPS 离源角） | `read_simul_takeoff_file` + `simul` 参数；按事件 ID/台站名匹配极性 |
| CHECK_POL / GETSTAT | `read_polarity_reversal_file` / `read_station_file`（含 TRI 格式） |
| MK_TABLE / GET_TTS | `TakeoffTable` / `compute_takeoff_angles` |
| MECH_ROT | `kagan_angle(n1, s1, n2, s2)`（原生绑定） |
| 核心版本 | `backend.get_backend("fortran").native_version()` |
