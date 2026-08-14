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
