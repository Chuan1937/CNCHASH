# 使用指南

## 基本用法

根据 P 波初动极性确定震源机制：

```python
from cnchash import run_hash
import numpy as np

# 台站数据
p_azi = np.array([45, 135, 225, 315, 0, 90, 180, 270])    # 方位角（度）
p_the = np.array([30, 45, 60, 75, 40, 50, 55, 65])         # 离源角（度）
p_pol = np.array([1, -1, 1, -1, 1, -1, 1, -1])             # 极性：1=向上, -1=向下
p_qual = np.zeros(8)                                        # 质量：0=清晰, 1=模糊

# 运行 HASH
result = run_hash(p_azi, p_the, p_pol, p_qual)

# 输出
print(f"走向: {result['strike_avg']:.1f}°")
print(f"倾角: {result['dip_avg']:.1f}°")
print(f"滑动角: {result['rake_avg']:.1f}°")
print(f"质量等级: {result['quality']}")
```

## 加入 S/P 振幅比

加入 S/P 振幅比可以获得更好的约束：

```python
from cnchash import run_hash_with_amp

# S/P 比值（log10 尺度）
# 0.0 = 无数据，典型范围：-1.0 到 2.0
sp_amp = np.array([0.3, -0.2, 0.5, 0.0, 0.4, -0.1, 0.6, 0.2])

result = run_hash_with_amp(p_azi, p_the, p_pol, sp_amp)

print(f"极性失配率: {result['mfrac']*100:.1f}%")
print(f"振幅失配: {result['mavg']:.2f}")
```

## 从输入文件处理

批量处理 HASH 输入文件中的事件：

```python
from cnchash import run_hash_from_file

results = run_hash_from_file("example.inp")

for result in results:
    if result['success']:
        print(f"S={result['strike_avg']:.0f}° D={result['dip_avg']:.0f}° "
              f"R={result['rake_avg']:.0f}° Q={result['quality']}")
```

### HASH 输入文件（原始 Fortran 格式）

原始 Fortran HASH 程序使用固定格式输入文件（`example.inp`），每行一个
数值，不支持注释。

**格式 1 结构：**

| 行 | 内容 | 说明 |
|----|------|------|
| 1 | polfile | 极性反转文件 |
| 2 | phasefile | 震相数据文件 |
| 3 | outfile1 | 输出文件（机制解） |
| 4 | outfile2 | 输出文件（可接受节面） |
| 5-14 | parameters | 算法参数（见下） |

**算法参数（第 5-14 行）：**

| 行 | 参数 | 说明 |
|----|------|------|
| 5 | npolmin | 所需最少极性数 |
| 6 | max_agap | 最大方位角空隙（度） |
| 7 | max_pgap | 最大离源角空隙（度） |
| 8 | dang | 网格角度步长（度） |
| 9 | nmc | Monte Carlo 试验次数 |
| 10 | maxout | 最大输出机制数 |
| 11 | badfrac | 允许的坏极性比例 |
| 12 | delmax | 最大台站距离（km） |
| 13 | cangle | 聚类角度（度） |
| 14 | prob_max | 概率阈值 |

**其他格式变体：**

| 格式 | 文件结构 |
|------|----------|
| 2/4 | stationfile, polfile, phasefile, out1, out2, params |
| 3 | stationfile, polfile, statcor, amp, phasefile, out1, params |
| 5 | polfile, simulfile, phasefile, out1, out2, params |

### 支持的文件格式

#### 震相文件（`read_phase_file`）

支持多种格式，自动检测：

| 格式 | 示例 | 说明 |
|------|------|------|
| 1 | north1.phase | 2 位年份，压缩格式 |
| 2 | north2.phase | 4 位年份，独立台站行 |
| 3 | north4.phase | 8 位日期格式 |
| 4 | north5.simul | SIMUL2000 格式 |

**格式 1**（2 位年份）：
```
YY MDD hhmmss lat lon depth mag ... event_id
STATION polcode ...
...
event_id
```

**格式 2**（4 位年份）：
```
YYYY MDDHHMM lat lon depth mag ... event_id
STATION NETWORK COMPONENT ONSET POLARITY
```

#### 台站文件（`read_station_file`）

```
STATION COMPONENT lat lon elevation
```

示例：
```
PAS   HHZ   34.1512  -118.1567   0.405
CAL   HHZ   34.1296  -117.9258   0.757
```

#### 极性反转文件（`read_polarity_reversal_file`）

记录台站极性随时间反转的信息：

```
STATION start_date end_date
```

日期为 YYYYMMDD 格式。`end_date = 0` 表示持续至今。

示例：
```
PAS  20100101  20121231
CAL  20150601  0
```

#### 速度模型文件（`read_velocity_model`）

```
depth velocity
```

示例：
```
0.0   5.5
5.0   5.8
10.0  6.2
20.0  6.8
```

### 直接读取文件

```python
from cnchash.io import (
    read_phase_file,
    read_station_file,
    read_polarity_reversal_file,
    read_velocity_model,
    read_hash_input_file
)

# 读取震相数据
events = read_phase_file("north1.phase")

# 读取台站坐标
stations = read_station_file("scsn.stations")

# 读取极性反转
reversals = read_polarity_reversal_file("scsn.pol")

# 读取速度模型
depths, velocities = read_velocity_model("vz.layer")

# 解析输入文件
params = read_hash_input_file("example.inp")
```
